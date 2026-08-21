//! Pure-Rust model run (plan.md, Phase 12: "Retire Fortran From The Rust
//! Build").
//!
//! This module replaces the last Fortran runtime authority: it wires the
//! Phase 9 geometry (surrounding points, upwind neighbours, mask/Neumann
//! refinement, boundary support-point mapping, observation weights) into the
//! Phase 11 solver port, and ports the boundary-condition update
//! (`update_boundary_points`, `update_wind_field`, `make_theta_grid`,
//! `build_boundary_support_points_spectra`, `update_boundaries`), the
//! observation-point update (`update_obs_points`) and the directionally
//! integrated output (`directional_spreading`) that were still Fortran.
//!
//! # Scope and deviations (recorded in the Phase 12 parity report)
//!
//! * Mesh input is NetCDF only (`nc_read_net` port in `crate::mesh`). The
//!   structured index/mask and ASCII mesh readers, and the file-backed
//!   `fw`/`fwig`/`u10`/`u10dir` interpolation (`read_interpolate_map_input`
//!   → `triintfast`/Triangle/`kdtree2`), have no checked-in testcase and
//!   stay Fortran; the Rust run rejects them with a clear error.
//! * Vegetation (`ja_vegetation == 1`) is likewise rejected (its NetCDF
//!   reader is not ported). No checked-in case enables it.
//! * Wind growth is ported for the *uniform* value path; a wind *list*
//!   (`windlistfile`) is rejected for the same reason.
//!
//! The output writer (`crate::output`) is reused unchanged: the model builds
//! the same [`crate::capture`] structs the Fortran capture stream used to
//! produce, so the NetCDF schema and fill logic stay byte-identical.

use std::path::Path;

use anyhow::{bail, Context, Result};

use crate::capture::{HisRecord, MapRecord, StaticHis, StaticMap};
use crate::date::date_to_iso8601;
use crate::geometry::{self, DomainGeometry};
use crate::input::SnapWaveInput;
use crate::interp;
use crate::mesh::Mesh;
use crate::solver;
use crate::state::ModelState;
use crate::text_input::{deg2rad_f32, pi_f32, BoundaryInput, ParsedTextInputs};

/// `FILL_VALUE` of `snapwave_data.f90` (`-999999.0`, real*4).
const FILL_VALUE: f32 = -999_999.0;
/// Water density [kg/m³] — `snapwave_data::rho`.
const RHO: f32 = 1025.0;
/// Gravitational acceleration [m/s²] — `snapwave_data::g`.
const G: f32 = 9.813;

/// `rad2deg = 180d0 / pi` of `snapwave_data.f90`: a double-precision
/// division (`180d0`) truncated back to `real*4` on the parameter's kind.
pub fn rad2deg_f32() -> f32 {
    (180.0f64 / (pi_f32() as f64)) as f32
}

/// Fortran `mod(a, m)`: sign follows `a` (truncated division).
fn fortran_mod(a: f32, m: f32) -> f32 {
    a % m
}

/// Fortran `modulo(a, m)`: result in `[0, m)`.
fn fortran_modulo(a: f32, m: f32) -> f32 {
    a.rem_euclid(m)
}

/// Fortran `nint(x)`: round half away from zero.
fn nint(x: f32) -> i32 {
    x.round() as i32
}

/// `mod2` of `snapwave_boundaries.f90`: `mod(a,b)` with 0 → b and
/// negatives wrapped into `1..=b`.
fn mod2(a: i32, b: i32) -> i32 {
    let mut c = a % b;
    if c == 0 {
        c = b;
    }
    if c < 0 {
        c += b;
    }
    c
}

/// `weighted_average` of `snapwave_boundaries.f90` (`iopt = 2`: angles in
/// radians). All arithmetic is `real*4`.
fn weighted_average_angle(val1: f32, val2: f32, fac: f32) -> f32 {
    let (u1, v1) = (val1.cos(), val1.sin());
    let (u2, v2) = (val2.cos(), val2.sin());
    let u = u1 * fac + u2 * (1.0 - fac);
    let v = v1 * fac + v2 * (1.0 - fac);
    v.atan2(u)
}

/// The full solver-state of one model run, owned by Rust (the explicit
/// replacement of the `snapwave_data` globals for the solver path).
pub struct Model {
    config: SnapWaveInput,
    mesh: Mesh,

    // ---- boundary forcing series (Rust-owned, converted to radians) ----
    nwbnd: usize,
    ntwbnd: usize,
    t_bwv: Vec<f32>,
    hs_bwv: Vec<f32>,
    tp_bwv: Vec<f32>,
    wd_bwv: Vec<f32>,
    ds_bwv: Vec<f32>,
    zs_bwv: Vec<f32>,

    // ---- geometry ----
    ntheta: usize,
    ntheta360: usize,
    dhdx: Vec<f32>,
    dhdy: Vec<f32>,
    w360: Vec<f32>,
    prev360: Vec<i32>,
    ds360: Vec<f32>,
    theta360: Vec<f32>,
    msk: Vec<i32>,
    neumannconnected: Vec<i32>,
    nmindbnd: Vec<i32>,
    ind1: Vec<i32>,
    ind2: Vec<i32>,
    fac: Vec<f32>,

    // ---- observation interpolation ----
    wobs: Vec<f64>,
    irefobs: Vec<i32>,
    /// Observation points as parsed (x, y, name), for `xobs`/`yobs`/`nameobs`.
    obs_points: Vec<(f64, f64, String)>,

    // ---- friction / wind ----
    fw: Vec<f32>,
    fw_ig: Vec<f32>,
    u10: Vec<f32>,
    u10dir: Vec<f32>,
    windspread360: Vec<f32>,

    // ---- boundary per-time state (update_boundary_points) ----
    hst_bwv: Vec<f32>,
    tpt_bwv: Vec<f32>,
    wdt_bwv: Vec<f32>,
    dst_bwv: Vec<f32>,
    zst_bwv: Vec<f32>,
    eet_bwv: Vec<f32>,
    tpmean_bwv: f32,
    zsmean_bwv: f32,
    wdmean_bwv: f32,
    /// `itwbndlast` (1-based; starts at 2 like `read_boundary_data`).
    itwbndlast: usize,
    /// `Tpini`, mutable because `update_boundaries` overwrites it with
    /// `tpmean_bwv` when wind is off.
    tpini: f32,

    // ---- solver state (per timestep) ----
    dtheta_rad: f32,
    depth: Vec<f32>,
    theta: Vec<f32>,
    w: Vec<f32>,
    prev: Vec<i32>,
    ds: Vec<f32>,
    windspreadfac: Vec<f32>,
    ee: Vec<f32>,
    ee_ig: Vec<f32>,
    aa: Vec<f32>,
    wsor_e: Vec<f32>,
    wsor_a: Vec<f32>,
    swe: Vec<f32>,
    swa: Vec<f32>,
    kwav: Vec<f32>,
    cg: Vec<f32>,
    ctheta: Vec<f32>,
    ctheta_ig: Vec<f32>,
    sinhkh: Vec<f32>,
    hmx: Vec<f32>,
    hmx_ig: Vec<f32>,
    kwav_ig: Vec<f32>,
    cg_ig: Vec<f32>,
    sig: Vec<f32>,
    tp: Vec<f32>,
    h: Vec<f32>,
    h_ig: Vec<f32>,
    dw: Vec<f32>,
    df: Vec<f32>,
    f: Vec<f32>,
    thetam: Vec<f32>,
    dveg: Vec<f32>,
    fx: Vec<f32>,
    fy: Vec<f32>,

    // ---- observation outputs ----
    hm0obs: Vec<f32>,
    zsobs: Vec<f32>,
    tpobs: Vec<f32>,
    hm0igobs: Vec<f32>,
    dwobs: Vec<f32>,
    dfobs: Vec<f32>,
    stobs: Vec<f32>,
    swobs: Vec<f32>,
    wdobs: Vec<f32>,
    dirsprobs: Vec<f32>,
}

impl Model {
    /// Build the Rust-owned model state from the parsed configuration, the
    /// Rust-read mesh and the Rust-parsed auxiliary text inputs.
    pub fn new(config: &SnapWaveInput, mesh: &Mesh, text: &ParsedTextInputs) -> Result<Model> {
        // plan.md Phase 12 deviations: the paths that stay Fortran are
        // rejected up front with a clear error (no checked-in testcase).
        if !crate::state::rust_owns_gridfile(&config.grid.gridfile) {
            bail!(
                "gridfile '{}' is not a NetCDF mesh; the Rust build only supports NetCDF \
                 meshes (structured/ASCII mesh readers stay Fortran)",
                config.grid.gridfile
            );
        }
        if config.vegetation.ja_vegetation == 1 {
            bail!("ja_vegetation = 1 is not supported by the Rust build (vegetation NetCDF reader stays Fortran)");
        }
        if !config.wind.windlistfile.is_empty() {
            bail!("windlistfile is not supported by the Rust build (wind list reader stays Fortran)");
        }

        let no_nodes = mesh.no_nodes;
        let dtheta_deg = config.grid.dtheta;
        let ntheta360 = (360.0f32 / dtheta_deg).round() as usize;
        let ntheta = (config.grid.sector / dtheta_deg).round() as usize;
        let dtheta_rad = dtheta_deg * deg2rad_f32();
        let deg2rad = deg2rad_f32();

        // ---- boundary series (converted to radians, time-major flat) ----
        let (nwbnd, ntwbnd, t_bwv, hs_bwv, tp_bwv, wd_bwv, ds_bwv, zs_bwv, x_bwv, y_bwv) =
            match &text.boundary {
                BoundaryInput::None => (0usize, 0usize, Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new()),
                BoundaryInput::Single(j) => (
                    1usize,
                    j.len(),
                    j.t.clone(),
                    j.hs.clone(),
                    j.tp.clone(),
                    j.wd_rad(),
                    j.ds_rad(),
                    j.zs.clone(),
                    Vec::new(),
                    Vec::new(),
                ),
                BoundaryInput::Timeseries(s) => (
                    s.nwbnd,
                    s.ntwbnd,
                    s.t.clone(),
                    s.hs.clone(),
                    s.tp.clone(),
                    s.wd_rad(),
                    s.ds_rad(),
                    s.zs.clone(),
                    s.x.clone(),
                    s.y.clone(),
                ),
            };

        // ---- derived geometry (Phase 9 ports) ----
        let enclosure = text.enclosure.as_ref().map(|p| (p.x.as_slice(), p.y.as_slice()));
        let neumann = text.neumann.as_ref().map(|p| (p.x.as_slice(), p.y.as_slice()));
        let DomainGeometry { dhdx, dhdy, w360, prev360, ds360, msk, neumannconnected, nmindbnd, .. } =
            geometry::compute_domain_geometry(
                &mesh.x,
                &mesh.y,
                &mesh.zb,
                mesh.sferic,
                &mesh.face_nodes,
                mesh.no_faces,
                dtheta_deg,
                config.boundary.tol,
                enclosure,
                neumann,
                &mesh.msk,
            );

        // theta360 (real*4 radians): `theta360 = (0.5*dtheta + (i-1)*dtheta)`,
        // the degree value evaluated in real*4, then `* deg2rad`.
        let theta360: Vec<f32> = (0..ntheta360)
            .map(|i| (0.5f32 * dtheta_deg + (i as f32) * dtheta_deg) * deg2rad)
            .collect();

        // Boundary support-point mapping (find_boundary_indices).
        let (ind1, ind2, fac) = if nwbnd > 0 {
            let b = geometry::find_boundary_indices(&mesh.x, &mesh.y, &msk, &x_bwv, &y_bwv, nwbnd);
            (b.ind1, b.ind2, b.fac)
        } else {
            (Vec::new(), Vec::new(), Vec::new())
        };

        // Observation interpolation weights (make_map_fm) and metadata.
        let (wobs, irefobs, obs_points) = match &text.obs {
            Some(o) => {
                let xobs: Vec<f64> = o.points.iter().map(|p| p.x).collect();
                let yobs: Vec<f64> = o.points.iter().map(|p| p.y).collect();
                let map = interp::make_map_fm(&mesh.x, &mesh.y, &mesh.face_nodes, mesh.no_faces, &xobs, &yobs);
                let pts = o.points.iter().map(|p| (p.x, p.y, p.name.clone())).collect();
                (map.w, map.iref, pts)
            }
            None => (Vec::new(), Vec::new(), Vec::new()),
        };
        let nobs = obs_points.len();

        // ---- friction (uniform; `where(zsini - zb > fwcutoff) fw = 0`) ----
        let fw0: f32 = config
            .physics
            .fw
            .parse()
            .with_context(|| format!("fw '{}' is not a uniform value (file-backed friction stays Fortran)", config.physics.fw))?;
        let fw0_ig: f32 = config.physics.fw_ig.parse().with_context(|| {
            format!("fwig '{}' is not a uniform value (file-backed friction stays Fortran)", config.physics.fw_ig)
        })?;
        let fwcutoff = config.physics.fwcutoff;
        let zsini = config.physics.zsini;
        let fw: Vec<f32> = mesh
            .zb
            .iter()
            .map(|&zb| if zsini - zb > fwcutoff { 0.0 } else { fw0 })
            .collect();
        let fw_ig: Vec<f32> = mesh
            .zb
            .iter()
            .map(|&zb| if zsini - zb > fwcutoff { 0.0 } else { fw0_ig })
            .collect();

        // ---- wind (uniform; read_wind_data else-branch) ----
        let wind = if config.wind.enabled { 1 } else { 0 };
        let (u10, u10dir) = if wind == 1 {
            let u10_0: f32 = config
                .wind
                .u10
                .parse()
                .with_context(|| format!("u10 '{}' is not a uniform value (file-backed wind stays Fortran)", config.wind.u10))?;
            let u10dir_0: f32 = config.wind.u10dir.parse().with_context(|| {
                format!("u10dir '{}' is not a uniform value (file-backed wind stays Fortran)", config.wind.u10dir)
            })?;
            // `u10dir = mod(270 - u10dir, 360) * deg2rad` (nautical coming-from
            // degrees → cartesian going-to radians).
            let dir = fortran_mod(270.0 - u10dir_0, 360.0) * deg2rad;
            (vec![u10_0; no_nodes], vec![dir; no_nodes])
        } else {
            (vec![0.0; no_nodes], vec![0.0; no_nodes])
        };

        // windspread360 initialized to zero (matches the wind==0 path; the
        // wind==1 distribution is computed per-timestep in update_wind_field).
        let windspread360 = vec![0.0f32; ntheta360 * no_nodes];

        // ---- solver state arrays ----
        let ntheta_no = ntheta * no_nodes;
        let depth = vec![0.0f32; no_nodes];
        let theta = vec![0.0f32; ntheta];
        let w = vec![0.0f32; 2 * ntheta_no];
        let prev = vec![0i32; 2 * ntheta_no];
        let ds = vec![0.0f32; ntheta_no];
        let windspreadfac = vec![0.0f32; ntheta_no];
        let ee = vec![0.0f32; ntheta_no];
        let ee_ig = vec![0.0f32; ntheta_no];
        let aa = vec![0.0f32; ntheta_no];
        let wsor_e = vec![0.0f32; ntheta_no];
        let wsor_a = vec![0.0f32; ntheta_no];
        let swe = vec![0.0f32; no_nodes];
        let swa = vec![0.0f32; no_nodes];
        let kwav = vec![0.0f32; no_nodes];
        let cg = vec![0.0f32; no_nodes];
        let ctheta = vec![0.0f32; ntheta_no];
        let ctheta_ig = vec![0.0f32; ntheta_no];
        let sinhkh = vec![0.0f32; no_nodes];
        let hmx = vec![0.0f32; no_nodes];
        let hmx_ig = vec![0.0f32; no_nodes];
        let kwav_ig = vec![0.0f32; no_nodes];
        let cg_ig = vec![0.0f32; no_nodes];
        let sig = vec![0.0f32; no_nodes];
        let tp = vec![config.physics.Tpini; no_nodes];
        let h = vec![0.0f32; no_nodes];
        let h_ig = vec![0.0f32; no_nodes];
        let dw = vec![0.0f32; no_nodes];
        let df = vec![0.0f32; no_nodes];
        let f = vec![0.0f32; no_nodes];
        let thetam = vec![0.0f32; no_nodes];
        let dveg = vec![0.0f32; no_nodes];
        let fx = vec![0.0f32; no_nodes];
        let fy = vec![0.0f32; no_nodes];

        // Boundary per-time work arrays.
        let hst_bwv = vec![0.0f32; nwbnd];
        let tpt_bwv = vec![0.0f32; nwbnd];
        let wdt_bwv = vec![0.0f32; nwbnd];
        let dst_bwv = vec![0.0f32; nwbnd];
        let zst_bwv = vec![0.0f32; nwbnd];
        let eet_bwv = vec![0.0f32; ntheta * nwbnd];

        // Observation outputs.
        let hm0obs = vec![FILL_VALUE; nobs];
        let zsobs = vec![FILL_VALUE; nobs];
        let tpobs = vec![FILL_VALUE; nobs];
        let hm0igobs = vec![FILL_VALUE; nobs];
        let dwobs = vec![FILL_VALUE; nobs];
        let dfobs = vec![FILL_VALUE; nobs];
        let stobs = vec![FILL_VALUE; nobs];
        let swobs = vec![FILL_VALUE; nobs];
        let wdobs = vec![FILL_VALUE; nobs];
        let dirsprobs = vec![FILL_VALUE; nobs];

        Ok(Model {
            config: config.clone(),
            mesh: mesh.clone(),
            nwbnd,
            ntwbnd,
            t_bwv,
            hs_bwv,
            tp_bwv,
            wd_bwv,
            ds_bwv,
            zs_bwv,
            ntheta,
            ntheta360,
            dhdx,
            dhdy,
            w360,
            prev360,
            ds360,
            theta360,
            msk,
            neumannconnected,
            nmindbnd,
            ind1,
            ind2,
            fac,
            wobs,
            irefobs,
            obs_points,
            fw,
            fw_ig,
            u10,
            u10dir,
            windspread360,
            hst_bwv,
            tpt_bwv,
            wdt_bwv,
            dst_bwv,
            zst_bwv,
            eet_bwv,
            tpmean_bwv: config.physics.Tpini,
            zsmean_bwv: 0.0,
            wdmean_bwv: 0.0,
            itwbndlast: 2,
            tpini: config.physics.Tpini,
            dtheta_rad,
            depth,
            theta,
            w,
            prev,
            ds,
            windspreadfac,
            ee,
            ee_ig,
            aa,
            wsor_e,
            wsor_a,
            swe,
            swa,
            kwav,
            cg,
            ctheta,
            ctheta_ig,
            sinhkh,
            hmx,
            hmx_ig,
            kwav_ig,
            cg_ig,
            sig,
            tp,
            h,
            h_ig,
            dw,
            df,
            f,
            thetam,
            dveg,
            fx,
            fy,
            hm0obs,
            zsobs,
            tpobs,
            hm0igobs,
            dwobs,
            dfobs,
            stobs,
            swobs,
            wdobs,
            dirsprobs,
        })
    }

    /// Run the model end to end: write the map/history NetCDF output files.
    /// `map_path`/`his_path` are `None` when that output family is disabled.
    pub fn run(&mut self, map_path: Option<&Path>, his_path: Option<&Path>) -> Result<()> {
        let map_file_nonempty = map_path.is_some();
        let his_file_nonempty = his_path.is_some();
        let nobs = self.wobs.len() / 4;
        let his_interval = self.config.output.his_interval;
        let map_interval = self.config.output.map_interval;
        let ja_save_each_iter = self.config.output.ja_save_each_iter;

        // Static map / history state (built once).
        let static_map = if map_file_nonempty {
            Some(self.build_static_map())
        } else {
            None
        };
        let static_his = if his_file_nonempty && nobs > 0 {
            Some(self.build_static_his())
        } else {
            None
        };

        let mut map_records: Vec<MapRecord> = Vec::new();
        let mut his_records: Vec<HisRecord> = Vec::new();

        let mut schedule = ModelState::new(&self.config);
        while schedule.is_running() {
            schedule.advance_iteration();

            self.update_boundary_conditions(schedule.t);
            let iter_outputs = self.compute_wave_field(schedule.t);

            // Per-iteration map output (ja_save_each_iter == 1): the solver
            // emits one record per sweep iteration at time+iter, matching
            // `ncoutput_update_map(time + iter, iter)` in the Fortran.
            if map_file_nonempty {
                for out in &iter_outputs {
                    map_records.push(self.build_iter_map_record(out));
                }
            }

            if schedule.should_output_his(his_file_nonempty, nobs as i32) {
                schedule.record_his_output(his_interval);
                self.update_obs_points();
                his_records.push(self.build_his_record(schedule.t, schedule.his_output_count));
            }

            if schedule.should_output_map(ja_save_each_iter, map_file_nonempty) {
                schedule.record_map_output(map_interval);
                map_records.push(self.build_map_record(schedule.t, schedule.map_output_count));
            }

            schedule.advance_time();
        }

        if let Some(sm) = &static_map {
            if let Some(path) = map_path {
                crate::output::write_map(path, &self.config, sm, &map_records)?;
            }
        }
        if let Some(sh) = &static_his {
            if let Some(path) = his_path {
                crate::output::write_his(path, &self.config, sh, &his_records)?;
            }
        }
        Ok(())
    }

    // ------------------------------------------------------------------
    // Static output state
    // ------------------------------------------------------------------

    fn build_static_map(&self) -> StaticMap {
        let mesh = &self.mesh;
        // face_nodes(1:max_nodes,:) section — node-major within a face.
        let mut face_nodes = Vec::with_capacity(mesh.max_nodes * mesh.no_faces);
        for f in 0..mesh.no_faces {
            for j in 0..mesh.max_nodes {
                face_nodes.push(mesh.face_nodes[f * 4 + j]);
            }
        }
        StaticMap {
            no_nodes: mesh.no_nodes,
            no_faces: mesh.no_faces,
            max_nodes: mesh.max_nodes,
            ntheta: self.ntheta,
            sferic: mesh.sferic,
            tref_iso8601: date_to_iso8601(&self.config.time.tref)
                .unwrap_or_else(|_| self.config.time.tref.clone()),
            libvers: crate::output::LIBVERS.to_string(),
            x: mesh.x.clone(),
            y: mesh.y.clone(),
            zb: mesh.zb.clone(),
            face_nodes,
            fw: self.fw.clone(),
            fw_ig: self.fw_ig.clone(),
            veg: None,
        }
    }

    fn build_static_his(&self) -> StaticHis {
        let names: Vec<String> = self.obs_points.iter().map(|(_, _, n)| n.clone()).collect();
        StaticHis {
            tref_iso8601: date_to_iso8601(&self.config.time.tref)
                .unwrap_or_else(|_| self.config.time.tref.clone()),
            libvers: crate::output::LIBVERS.to_string(),
            nobs: names.len(),
            xobs: self.obs_points.iter().map(|(x, _, _)| *x).collect(),
            yobs: self.obs_points.iter().map(|(_, y, _)| *y).collect(),
            names,
        }
    }
}

impl Model {
    // ------------------------------------------------------------------
    // update_boundary_conditions (plan.md Phase 12 port)
    // ------------------------------------------------------------------

    fn update_boundary_conditions(&mut self, t: f64) {
        self.update_boundary_points(t);
        if self.config.wind.enabled {
            self.update_wind_field(t);
        }

        let thetamean: f32;
        if self.ntwbnd > 0 && !self.config.wind.enabled {
            thetamean = self.wdmean_bwv;
            self.make_theta_grid(self.wdmean_bwv);
        } else {
            // wind == 1 (u10dmean) — set by update_wind_field.
            thetamean = self.u10dmean();
            self.make_theta_grid(thetamean);
        }

        self.build_boundary_support_points_spectra(thetamean);
        self.update_boundaries();
    }

    /// `update_boundary_points(t)`: interpolate the boundary time series to
    /// time `t`, then compute `tpmean_bwv`/`zsmean_bwv`/`wdmean_bwv`/`depth`.
    fn update_boundary_points(&mut self, t: f64) {
        if self.ntwbnd == 0 {
            self.depth = self
                .mesh
                .zb
                .iter()
                .map(|&zb| (self.config.physics.zsini - zb).max(self.config.physics.hmin))
                .collect();
            return;
        }
        let nwbnd = self.nwbnd;
        let t_bwv = &self.t_bwv;
        let (hs_bwv, tp_bwv, wd_bwv, ds_bwv, zs_bwv) =
            (&self.hs_bwv, &self.tp_bwv, &self.wd_bwv, &self.ds_bwv, &self.zs_bwv);

        let mut itb = self.itwbndlast;
        while itb <= self.ntwbnd {
            let t_cur = t_bwv[itb - 1] as f64;
            if t_cur > t || itb == self.ntwbnd {
                let denom = t_bwv[itb - 1] - t_bwv[itb - 2];
                let tbfac = ((t - t_bwv[itb - 2] as f64) / (denom as f64)) as f32;
                let f = 1.0 - tbfac;
                for ib in 0..nwbnd {
                    let j1 = (itb - 2) * nwbnd + ib;
                    let j2 = (itb - 1) * nwbnd + ib;
                    let hs = hs_bwv[j1] + (hs_bwv[j2] - hs_bwv[j1]) * tbfac;
                    let tps = tp_bwv[j1] + (tp_bwv[j2] - tp_bwv[j1]) * tbfac;
                    let dsp = ds_bwv[j1] + (ds_bwv[j2] - ds_bwv[j1]) * tbfac;
                    let zst = zs_bwv[j1] + (zs_bwv[j2] - zs_bwv[j1]) * tbfac;
                    let wd = weighted_average_angle(wd_bwv[j1], wd_bwv[j2], f);
                    self.hst_bwv[ib] = hs;
                    self.tpt_bwv[ib] = tps;
                    self.wdt_bwv[ib] = wd;
                    self.dst_bwv[ib] = dsp;
                    self.zst_bwv[ib] = zst;
                }
                self.itwbndlast = itb;
                break;
            }
            itb += 1;
        }

        let tpmean = self.tpt_bwv.iter().sum::<f32>() / nwbnd as f32;
        let zsmean = self.zst_bwv.iter().sum::<f32>() / nwbnd as f32;
        self.tpmean_bwv = tpmean;
        self.zsmean_bwv = zsmean;

        let hmin = self.config.physics.hmin;
        self.depth = self
            .mesh
            .zb
            .iter()
            .map(|&zb| (zsmean - zb).max(hmin))
            .collect();

        let hsum = self.hst_bwv.iter().sum::<f32>();
        let wdmean = (self.wdt_bwv.iter().zip(self.hst_bwv.iter()).map(|(w, h)| w.sin() * h).sum::<f32>() / hsum)
            .atan2(
                self.wdt_bwv.iter().zip(self.hst_bwv.iter()).map(|(w, h)| w.cos() * h).sum::<f32>() / hsum,
            );
        self.wdmean_bwv = wdmean;
    }

    fn u10dmean(&self) -> f32 {
        let num = self.u10.iter().zip(self.u10dir.iter()).map(|(u, d)| d.sin() * u).sum::<f32>();
        let den = self.u10.iter().zip(self.u10dir.iter()).map(|(u, d)| d.cos() * u).sum::<f32>();
        num.atan2(den)
    }

    /// `update_wind_field(t)`: uniform wind (`ntu10bnd == 1`), and rebuild
    /// the 360° wind-spread distribution.
    fn update_wind_field(&mut self, _t: f64) {
        let dtheta = self.dtheta_rad;
        let ntheta360 = self.ntheta360;
        let no_nodes = self.mesh.no_nodes;
        // windspread360(k) for each node: (cos(theta360 - u10dir))^2,
        // clipped at zero, normalized and /dtheta.
        for k in 0..no_nodes {
            let u10dirk = self.u10dir[k];
            let mut sum = 0.0f32;
            let mut w = Vec::with_capacity(ntheta360);
            for itheta in 0..ntheta360 {
                let c = (self.theta360[itheta] - u10dirk).cos();
                let mut v = c * c;
                if c < 0.0 {
                    v = 0.0;
                }
                w.push(v);
                sum += v;
            }
            for itheta in 0..ntheta360 {
                self.windspread360[k * ntheta360 + itheta] = if sum > 0.0 { w[itheta] / sum / dtheta } else { 0.0 };
            }
        }
    }

    /// `make_theta_grid(central_theta)`: build `theta`, and select
    /// `w`/`prev`/`ds`/`windspreadfac` from the 360° tables.
    fn make_theta_grid(&mut self, central_theta: f32) {
        let ntheta = self.ntheta;
        let ntheta360 = self.ntheta360;
        let no_nodes = self.mesh.no_nodes;
        let ind = nint(central_theta / self.dtheta_rad) - (ntheta as i32) / 2;

        let mut i360 = vec![0usize; ntheta];
        for itheta in 0..ntheta {
            i360[itheta] = (mod2(itheta as i32 + 1 + ind, ntheta360 as i32) - 1) as usize;
        }

        for itheta in 0..ntheta {
            let i360k = i360[itheta];
            self.theta[itheta] = self.theta360[i360k];
            for k in 0..no_nodes {
                let wbase = k * (2 * ntheta360) + i360k * 2;
                let dbase = k * ntheta360 + i360k;
                let tbase = k * (2 * ntheta) + itheta * 2;
                self.w[tbase] = self.w360[wbase];
                self.w[tbase + 1] = self.w360[wbase + 1];
                self.prev[tbase] = self.prev360[wbase];
                self.prev[tbase + 1] = self.prev360[wbase + 1];
                self.ds[k * ntheta + itheta] = self.ds360[dbase];
                self.windspreadfac[k * ntheta + itheta] = self.windspread360[dbase];
            }
        }

        if self.config.wind.enabled {
            let mwind = self.config.wind.mwind;
            for k in 0..no_nodes {
                let u10dirk = self.u10dir[k];
                let mut sum = 0.0f32;
                let mut w = Vec::with_capacity(ntheta);
                for itheta in 0..ntheta {
                    let c = (self.theta[itheta] - u10dirk).cos();
                    let mut v = c.powi(mwind);
                    if c < 0.0 {
                        v = 0.0;
                    }
                    w.push(v);
                    sum += v;
                }
                for itheta in 0..ntheta {
                    self.windspreadfac[k * ntheta + itheta] = if sum > 0.0 {
                        w[itheta] / sum / self.dtheta_rad
                    } else {
                        0.0
                    };
                }
            }
        }
    }

    /// `build_boundary_support_points_spectra`: JONSWAP-like directional
    /// spectrum on the boundary support points.
    fn build_boundary_support_points_spectra(&mut self, thetamean: f32) {
        let ntheta = self.ntheta;
        let nwbnd = self.nwbnd;
        let pi = pi_f32();
        for ib in 0..nwbnd {
            let e0 = 0.0625 * RHO * G * self.hst_bwv[ib] * self.hst_bwv[ib];
            let ms = 1.0 / (self.dst_bwv[ib] * self.dst_bwv[ib]) - 1.0;
            let mut dist = vec![0.0f32; ntheta];
            let mut sum = 0.0f32;
            for itheta in 0..ntheta {
                let c = (self.theta[itheta] - thetamean).cos();
                // `sign(1.0, c) * abs(c)**ms` (Fortran `sign` treats -0.0 as +).
                let v = if c >= 0.0 { c.abs().powf(ms) } else { -c.abs().powf(ms) };
                let arg = (fortran_mod(pi + self.theta[itheta] - thetamean, 2.0 * pi) - pi).abs();
                if arg > 0.999 * pi / 2.0 {
                    dist[itheta] = 0.0;
                } else {
                    dist[itheta] = v;
                }
                sum += dist[itheta];
            }
            for itheta in 0..ntheta {
                self.eet_bwv[ib * ntheta + itheta] = dist[itheta] / sum * e0 / self.dtheta_rad;
            }
        }
    }

    /// `update_boundaries`: set the boundary-node spectra (`ee`) and peak
    /// period (`Tp`) from the support-point spectra, and update `Tpini`.
    fn update_boundaries(&mut self) {
        if self.nwbnd == 0 {
            return;
        }
        let ntheta = self.ntheta;
        let nb = self.nmindbnd.len();
        for i in 0..nb {
            let k = self.nmindbnd[i] as usize - 1;
            let i1 = (self.ind1[i] - 1) as usize;
            let i2 = (self.ind2[i] - 1) as usize;
            let fac = self.fac[i];
            for itheta in 0..ntheta {
                self.ee[itheta + k * ntheta] = self.eet_bwv[i1 * ntheta + itheta] * fac
                    + self.eet_bwv[i2 * ntheta + itheta] * (1.0 - fac);
            }
            if !self.config.wind.enabled {
                self.tp[k] = self.tpmean_bwv;
            } else {
                self.tp[k] = self.tpt_bwv[i1] * fac + self.tpt_bwv[i2] * (1.0 - fac);
            }
        }
        if !self.config.wind.enabled {
            self.tpini = self.tpmean_bwv;
        }
    }

    // ------------------------------------------------------------------
    // compute_wave_field (Phase 11 solver port, wired with real geometry)
    // ------------------------------------------------------------------

    /// Run one solver step and return any per-iteration map outputs
    /// (`ja_save_each_iter == 1`).
    fn compute_wave_field(&mut self, time: f64) -> Vec<crate::solver::IterOutput> {
        let ig = self.config.physics.ig;
        let wind = if self.config.wind.enabled { 1 } else { 0 };
        let msk: Vec<i8> = self.msk.iter().map(|&m| m as i8).collect();

        let mut iter_outputs: Vec<crate::solver::IterOutput> = Vec::new();

        // IG arrays stay zero when ig == 0 (the Fortran sets them to zero
        // in that branch); the ig == 1 path is not exercised by any
        // checked-in testcase.
        solver::compute_wave_field(
            time,
            self.config.time.restart,
            ig,
            wind,
            0, // ja_vegetation (rejected earlier)
            self.config.output.ja_save_each_iter,
            self.ntheta,
            self.mesh.no_nodes,
            0, // no_secveg
            &self.mesh.x,
            &self.mesh.y,
            &self.dhdx,
            &self.dhdy,
            &msk,
            &self.neumannconnected,
            &self.theta,
            self.thetamean(),
            self.tpmean_bwv,
            &self.depth,
            &self.fw,
            &self.fw_ig,
            &self.w,
            &self.ds,
            &self.prev,
            self.config.time.dt,
            RHO,
            self.config.physics.alpha,
            self.config.physics.gamma,
            self.config.physics.gammax,
            &self.u10,
            self.config.time.niter,
            self.config.time.crit,
            self.config.physics.upwindref,
            self.tpini,
            &self.windspreadfac,
            self.config.physics.jadcgdx,
            self.config.physics.sigmin,
            self.config.physics.sigmax,
            self.config.physics.c_dispT,
            &self.kwav_ig,
            &self.cg_ig,
            &self.ctheta_ig,
            &self.hmx_ig,
            &[],
            &[],
            &[],
            &[],
            &mut self.kwav,
            &mut self.cg,
            &mut self.ctheta,
            &mut self.ee,
            &mut self.ee_ig,
            &mut self.sinhkh,
            &mut self.hmx,
            &mut self.tp,
            &mut self.sig,
            &mut self.aa,
            &mut self.wsor_e,
            &mut self.wsor_a,
            &mut self.swe,
            &mut self.swa,
            &mut self.h,
            &mut self.h_ig,
            &mut self.dw,
            &mut self.df,
            &mut self.f,
            &mut self.thetam,
            &mut self.dveg,
            &mut self.fx,
            &mut self.fy,
            &mut iter_outputs,
        );
        iter_outputs
    }

    fn thetamean(&self) -> f32 {
        if self.ntwbnd > 0 && !self.config.wind.enabled {
            self.wdmean_bwv
        } else {
            self.u10dmean()
        }
    }

    // ------------------------------------------------------------------
    // update_obs_points + output records
    // ------------------------------------------------------------------

    /// `update_obs_points` of `snapwave_obspoints.f90`.
    fn update_obs_points(&mut self) {
        let nobs = self.wobs.len() / 4;
        if nobs == 0 {
            return;
        }
        let ntheta = self.ntheta;
        let dtheta_dp = self.dtheta_rad as f64;
        let rad2deg_dp = rad2deg_f32() as f64;
        let sqrt2 = 2.0f32.sqrt();
        let ig = self.config.physics.ig == 1;
        let wind = self.config.wind.enabled;

        let cos_theta: Vec<f64> = self.theta.iter().map(|&t| (t as f64).cos()).collect();
        let sin_theta: Vec<f64> = self.theta.iter().map(|&t| (t as f64).sin()).collect();

        for iobs in 0..nobs {
            self.hm0obs[iobs] = FILL_VALUE;
            self.zsobs[iobs] = FILL_VALUE;
            self.tpobs[iobs] = FILL_VALUE;
            self.hm0igobs[iobs] = FILL_VALUE;
            self.dwobs[iobs] = FILL_VALUE;
            self.dfobs[iobs] = FILL_VALUE;
            self.stobs[iobs] = FILL_VALUE;
            self.swobs[iobs] = FILL_VALUE;
            self.wdobs[iobs] = FILL_VALUE;
            self.dirsprobs[iobs] = FILL_VALUE;

            if self.irefobs[4 * iobs] > 0 {
                self.hm0obs[iobs] = 0.0;
                self.zsobs[iobs] = 0.0;
                self.tpobs[iobs] = 0.0;
                self.dwobs[iobs] = 0.0;
                self.dfobs[iobs] = 0.0;
                self.dirsprobs[iobs] = 0.0;
                let mut hm0x_sum = 0.0f64;
                let mut hm0y_sum = 0.0f64;
                let mut m0_obs = 0.0f64;
                let mut a1_obs = 0.0f64;
                let mut b1_obs = 0.0f64;
                if ig {
                    self.hm0igobs[iobs] = 0.0;
                }
                if wind {
                    self.swobs[iobs] = 0.0;
                    self.stobs[iobs] = 0.0;
                }

                for ip in 0..4 {
                    let k = self.irefobs[ip + 4 * iobs].max(1) as usize - 1;
                    let weight = self.wobs[ip + 4 * iobs];
                    let weight_sp = weight as f32;
                    self.hm0obs[iobs] += weight_sp * self.h[k];
                    self.zsobs[iobs] += weight_sp * (self.depth[k] + self.mesh.zb[k]);
                    self.tpobs[iobs] += weight_sp * self.tp[k];
                    self.dwobs[iobs] += weight_sp * self.dw[k];
                    self.dfobs[iobs] += weight_sp * self.df[k];
                    hm0x_sum += weight * (self.h[k] as f64) * (self.thetam[k] as f64).cos();
                    hm0y_sum += weight * (self.h[k] as f64) * (self.thetam[k] as f64).sin();
                    for itheta in 0..ntheta {
                        let energy_bin = (self.ee[itheta + k * ntheta] as f64).max(0.0);
                        let energy_weight = weight * energy_bin * dtheta_dp;
                        m0_obs += energy_weight;
                        a1_obs += energy_weight * cos_theta[itheta];
                        b1_obs += energy_weight * sin_theta[itheta];
                    }
                    if ig {
                        self.hm0igobs[iobs] += weight_sp * self.h_ig[k];
                    }
                    if wind {
                        self.swobs[iobs] += weight_sp * self.swe[k];
                        self.stobs[iobs] += weight_sp * self.swa[k];
                    }
                }

                self.hm0obs[iobs] *= sqrt2;
                if ig {
                    self.hm0igobs[iobs] *= sqrt2;
                }
                // `wdobs = real(mod(270 - atan2(...)*rad2deg + 360, 360), sp)`:
                // computed in double precision, truncated to real*4 once.
                let ang_deg = 270.0_f64 - hm0y_sum.atan2(hm0x_sum) * rad2deg_dp + 360.0_f64;
                self.wdobs[iobs] = (ang_deg % 360.0) as f32;
                if m0_obs > 0.0 {
                    let a1 = a1_obs / m0_obs;
                    let b1 = b1_obs / m0_obs;
                    let mut r1 = (a1 * a1 + b1 * b1).sqrt();
                    if r1 > 1.0 {
                        r1 = 1.0;
                    }
                    if r1 < 0.0 {
                        r1 = 0.0;
                    }
                    self.dirsprobs[iobs] = ((2.0 * (1.0 - r1)).sqrt() * rad2deg_dp) as f32;
                } else {
                    self.dirsprobs[iobs] = FILL_VALUE;
                }
            }
        }
    }

    /// Directional spreading for one node's directional spectrum
    /// (`directional_spreading` in `snapwave_ncoutput.F90`), used by
    /// `map_dirspr`. `ee` is the `ntheta`-bin spectrum of one node.
    fn directional_spreading_of(&self, ee: &[f32]) -> f32 {
        let ntheta = self.ntheta;
        let dtheta_dp = self.dtheta_rad as f64;
        let rad2deg_dp = rad2deg_f32() as f64;
        if ntheta == 0 || self.dtheta_rad <= 0.0 {
            return FILL_VALUE;
        }
        let mut m0 = 0.0f64;
        let mut a1 = 0.0f64;
        let mut b1 = 0.0f64;
        for itheta in 0..ntheta {
            let mut energy = ee[itheta] as f64;
            if energy < 0.0 {
                energy = 0.0;
            }
            let weight = energy * dtheta_dp;
            let c = (self.theta[itheta] as f64).cos();
            let s = (self.theta[itheta] as f64).sin();
            m0 += weight;
            a1 += weight * c;
            b1 += weight * s;
        }
        if m0 > 0.0 {
            let a1 = a1 / m0;
            let b1 = b1 / m0;
            let mut r1 = (a1 * a1 + b1 * b1).sqrt();
            if r1 > 1.0 {
                r1 = 1.0;
            }
            if r1 < 0.0 {
                r1 = 0.0;
            }
            ((2.0 * (1.0 - r1)).sqrt() * rad2deg_dp) as f32
        } else {
            FILL_VALUE
        }
    }

    fn build_map_record(&self, t: f64, ntmapout: i32) -> MapRecord {
        let o = &self.config.output;
        let wind = self.config.wind.enabled;
        let ig = self.config.physics.ig == 1;
        let hmin = self.config.physics.hmin;
        let no_nodes = self.mesh.no_nodes;
        let rad2deg = rad2deg_f32();

        // `where (depth < hmin)` fill mask applied exactly like
        // `ncoutput_capture_update_map`.
        let fill_masked = |v: &[f32]| -> Vec<f32> {
            v.iter()
                .zip(self.depth.iter())
                .map(|(&x, &d)| if d < hmin { FILL_VALUE } else { x })
                .collect()
        };

        let depth = if o.map_depth == 1 { Some(self.depth.clone()) } else { None };
        let hm0 = if o.map_Hm0 == 1 {
            Some(self.h.iter().map(|&h| h * 2.0f32.sqrt()).collect())
        } else {
            None
        };
        let hm0_ig = if ig && o.map_Hig == 1 {
            Some(fill_masked(&self.h_ig.iter().map(|&h| h * 2.0f32.sqrt()).collect::<Vec<f32>>()))
        } else {
            None
        };
        let tp = if o.map_Tp == 1 { Some(fill_masked(&self.tp)) } else { None };
        let wd = if o.map_dir == 1 {
            Some(fill_masked(&self.thetam.iter().map(|&th| fortran_modulo(270.0 - th * rad2deg + 360.0, 360.0)).collect::<Vec<f32>>()))
        } else {
            None
        };
        let wdspr = if o.map_dirspr == 1 {
            let v: Vec<f32> = (0..no_nodes)
                .map(|k| self.directional_spreading_of(&self.ee[k * self.ntheta..(k + 1) * self.ntheta]))
                .collect();
            Some(fill_masked(&v))
        } else {
            None
        };
        let cg = if o.map_cg == 1 { Some(fill_masked(&self.cg)) } else { None };
        let dw = if o.map_Dw == 1 { Some(fill_masked(&self.dw)) } else { None };
        let df = if o.map_Df == 1 { Some(fill_masked(&self.df)) } else { None };
        let sw = if wind && o.map_SwE == 1 { Some(fill_masked(&self.swe)) } else { None };
        let st = if wind && o.map_SwA == 1 { Some(fill_masked(&self.swa)) } else { None };
        let sig = if wind && o.map_sig == 1 { Some(fill_masked(&self.sig)) } else { None };
        let u10 = if wind && o.map_u10 == 1 { Some(fill_masked(&self.u10)) } else { None };
        let u10dir = if wind && o.map_u10 == 1 {
            Some(fill_masked(&self.u10dir.iter().map(|&d| fortran_modulo(270.0 - d * rad2deg + 360.0, 360.0)).collect::<Vec<f32>>()))
        } else {
            None
        };
        let dveg = None; // ja_vegetation == 1 rejected earlier
        let ee = if o.map_ee == 1 { Some(self.ee.clone()) } else { None };
        let ctheta = if o.map_ctheta == 1 { Some(self.ctheta.clone()) } else { None };
        let theta_deg = if o.map_ee == 1 || o.map_ctheta == 1 {
            Some(self.theta.iter().map(|&th| fortran_modulo(270.0 - th * rad2deg, 360.0)).collect())
        } else {
            None
        };

        MapRecord {
            t,
            ntmapout,
            depth,
            hm0,
            hm0_ig,
            tp,
            wd,
            wdspr,
            cg,
            dw,
            df,
            sw,
            st,
            sig,
            u10,
            u10dir,
            dveg,
            ee,
            ctheta,
            theta_deg,
        }
    }

    fn build_his_record(&self, t: f64, nthisout: i32) -> HisRecord {
        let ig = self.config.physics.ig == 1;
        let wind = self.config.wind.enabled;
        HisRecord {
            t,
            nthisout,
            zs: self.zsobs.clone(),
            hm0: self.hm0obs.clone(),
            tp: self.tpobs.clone(),
            wavdir: self.wdobs.clone(),
            dirspr: self.dirsprobs.clone(),
            hm0ig: if ig { Some(self.hm0igobs.clone()) } else { None },
            dw: self.dwobs.clone(),
            df: self.dfobs.clone(),
            sw: if wind { Some(self.swobs.clone()) } else { None },
            st: if wind { Some(self.stobs.clone()) } else { None },
        }
    }

    /// Build a per-iteration map record from an [`IterOutput`] snapshot
    /// (the `ja_save_each_iter == 1` path). The fill masking and direction
    /// conversions are identical to [`Model::build_map_record`]; the static
    /// fields (depth, `theta`) come from the model, the solver state from
    /// the snapshot.
    fn build_iter_map_record(&self, out: &crate::solver::IterOutput) -> MapRecord {
        let o = &self.config.output;
        let wind = self.config.wind.enabled;
        let ig = self.config.physics.ig == 1;
        let hmin = self.config.physics.hmin;
        let rad2deg = rad2deg_f32();
        let s = &out.snapshot;

        let fill_masked = |v: &[f32]| -> Vec<f32> {
            v.iter()
                .zip(self.depth.iter())
                .map(|(&x, &d)| if d < hmin { FILL_VALUE } else { x })
                .collect()
        };

        let depth = if o.map_depth == 1 { Some(self.depth.clone()) } else { None };
        let hm0 = if o.map_Hm0 == 1 {
            Some(s.h.iter().map(|&h| h * 2.0f32.sqrt()).collect())
        } else {
            None
        };
        let hm0_ig = if ig && o.map_Hig == 1 {
            Some(fill_masked(&s.h_ig.iter().map(|&h| h * 2.0f32.sqrt()).collect::<Vec<f32>>()))
        } else {
            None
        };
        let tp = if o.map_Tp == 1 { Some(fill_masked(&s.tp)) } else { None };
        let wd = if o.map_dir == 1 {
            Some(fill_masked(&s.thetam.iter().map(|&th| fortran_modulo(270.0 - th * rad2deg + 360.0, 360.0)).collect::<Vec<f32>>()))
        } else {
            None
        };
        let wdspr = if o.map_dirspr == 1 {
            // Per-iteration directional spreading from the snapshot's ee.
            let v: Vec<f32> = (0..s.h.len())
                .map(|k| self.directional_spreading_of(&s.ee[k * self.ntheta..(k + 1) * self.ntheta]))
                .collect();
            Some(fill_masked(&v))
        } else {
            None
        };
        let cg = if o.map_cg == 1 { Some(fill_masked(&s.cg)) } else { None };
        let dw = if o.map_Dw == 1 { Some(fill_masked(&s.dw)) } else { None };
        let df = if o.map_Df == 1 { Some(fill_masked(&s.df)) } else { None };
        let sw = if wind && o.map_SwE == 1 { Some(fill_masked(&s.swe)) } else { None };
        let st = if wind && o.map_SwA == 1 { Some(fill_masked(&s.swa)) } else { None };
        let sig = if wind && o.map_sig == 1 { Some(fill_masked(&s.sig)) } else { None };
        let u10 = if wind && o.map_u10 == 1 { Some(fill_masked(&self.u10)) } else { None };
        let u10dir = if wind && o.map_u10 == 1 {
            Some(fill_masked(&self.u10dir.iter().map(|&d| fortran_modulo(270.0 - d * rad2deg + 360.0, 360.0)).collect::<Vec<f32>>()))
        } else {
            None
        };
        let dveg = None;
        let ee = if o.map_ee == 1 { Some(s.ee.clone()) } else { None };
        let ctheta = if o.map_ctheta == 1 { Some(s.ctheta.clone()) } else { None };
        let theta_deg = if o.map_ee == 1 || o.map_ctheta == 1 {
            Some(self.theta.iter().map(|&th| fortran_modulo(270.0 - th * rad2deg, 360.0)).collect())
        } else {
            None
        };

        MapRecord {
            t: out.time,
            ntmapout: out.ntmapout,
            depth,
            hm0,
            hm0_ig,
            tp,
            wd,
            wdspr,
            cg,
            dw,
            df,
            sw,
            st,
            sig,
            u10,
            u10dir,
            dveg,
            ee,
            ctheta,
            theta_deg,
        }
    }
}
