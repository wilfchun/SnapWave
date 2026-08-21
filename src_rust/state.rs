//! Explicit Rust state for the non-solver model data (plan.md Phase 8)
//! and the Rust-owned time loop / output scheduling (plan.md Phase 10).
//!
//! # What this module is
//!
//! Phases 3–7 made Rust parse and validate everything the model reads:
//! the configuration (`input`), the auxiliary text inputs (`text_input`)
//! and the UGRID mesh (`mesh`). Until now the run path still had the
//! Fortran readers re-read those files into the `snapwave_data` module
//! globals. This module closes that gap:
//!
//! * [`DomainState`] composes the already-Rust-owned pieces — config,
//!   mesh, boundary forcing, wind forcing, observation points, enclosure
//!   and Neumann polylines — plus the [`RuntimeState`] scalars the
//!   Fortran timestep loop initialises (the Phase 10 seam where output
//!   scheduling becomes Rust-owned). It is the explicit replacement of
//!   the `snapwave_data` globals for the migrated subsystems
//!   (Phase 8, step 1).
//! * [`FfiState`] owns the buffers that cross the FFI boundary with
//!   Fortran-compatible widths and layouts — `real*8`/`real*4`/
//!   `integer*1`/`integer`, column-major (Phase 8, step 2; the layout
//!   facts are pinned in `ffi_layout`) — and hands them to the coarse
//!   Fortran entry point `snapwave_run_capture_state_c` as one
//!   `#[repr(C)]` struct (Phase 8, steps 4–5: allocation ownership for
//!   the non-solver arrays moves to Rust; Fortran *associates* its
//!   module globals with this memory via `c_f_pointer` instead of
//!   reading files).
//! * [`ModelState`] (Phase 10) owns the time-loop scheduling state —
//!   current time, iteration counter, next-output times, output counters
//!   — and provides the scheduling predicates that mirror the Fortran
//!   `run_time_loop` logic exactly. The Rust `execute()` function now
//!   drives the loop, calling the coarse Fortran entry points
//!   `snapwave_timestep_c` / `snapwave_capture_{map,his}_c` /
//!   `snapwave_finalize_capture_c` (plan.md Phase 10, steps 1, 4).
//!
//! # What still comes from Fortran
//!
//! The Fortran side remains the runtime authority for everything that is
//! *derived* from this state rather than read from disk: surrounding
//! points / upwind-neighbour weights, the enclosure→mask refinement
//! (which **writes into the `msk` buffer handed over here** — the buffer
//! is Fortran-owned content from the handoff on, while `Mesh.msk` keeps
//! the as-read values), `make_map_fm` / `find_boundary_indices`
//! interpolation weights, the wind and `fw`/`fwig` value-or-file inputs
//! (their file-backed branch needs `triintfast` mesh interpolation,
//! Phase 9) and vegetation input. Those stay callable exactly as before;
//! the legacy `make` oracle keeps reading every file itself.
//!
//! # FFI struct contract
//!
//! [`SnapWaveStateC`] must stay field-for-field identical to the
//! `snapwave_state_t` `bind(C)` type in `src/snapwave_c_api.f90`
//! (AGENTS.md: FFI signatures match exactly). Absent data is
//! represented by zero extents plus empty (non-null, zero-length)
//! buffers, never by null pointers.

use std::ffi::c_int;

use anyhow::{bail, Result};

use crate::input::SnapWaveInput;
use crate::mesh::Mesh;
use crate::text_input::{BoundaryInput, ParsedTextInputs, Polyline};

/// Width of the Fortran `nameobs` global (`character*32`); every
/// observation-point name crosses as one 32-byte, blank-padded record.
const WIDTH_NAMEOBS: usize = 32;

// ----------------------------------------------------------------------
// Runtime state (plan.md Phase 8, step 1; the Phase 10 seam)
// ----------------------------------------------------------------------

/// The scheduling scalars `run_time_loop` in `src/snapwave_c_api.f90`
/// initialises from the configuration. Rust owns the *values* from here
/// on (they feed diagnostics and the Phase 8 handoff summary); the
/// authoritative loop that steps them stays Fortran until plan.md
/// Phase 10 makes output scheduling Rust-owned.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RuntimeState {
    pub tstart: f64,
    pub tstop: f64,
    pub timestep: f32,
    pub map_interval: f32,
    pub his_interval: f32,
    /// Both output schedules start at `tstart` (loop init).
    pub next_map_output: f64,
    pub next_his_output: f64,
    pub map_output_count: i32,
    pub his_output_count: i32,
    /// `max(1.0d-6, abs(dble(timestep))*1.0d-6)` of the Fortran loop.
    pub output_tol: f64,
}

impl RuntimeState {
    /// Derive the scheduling scalars exactly the way `run_time_loop`
    /// does, so the Rust-side summary cannot drift from the run.
    pub fn new(cfg: &SnapWaveInput) -> Self {
        RuntimeState {
            tstart: cfg.time.tstart,
            tstop: cfg.time.tstop,
            timestep: cfg.time.timestep,
            map_interval: cfg.output.map_interval,
            his_interval: cfg.output.his_interval,
            next_map_output: cfg.time.tstart,
            next_his_output: cfg.time.tstart,
            map_output_count: 0,
            his_output_count: 0,
            output_tol: 1.0e-6_f64.max((cfg.time.timestep as f64).abs() * 1.0e-6),
        }
    }
}

// ----------------------------------------------------------------------
// Model state (plan.md Phase 10: Rust-owned time loop and output scheduling)
// ----------------------------------------------------------------------

/// The time-loop scheduling state, moved from Fortran's `run_time_loop`
/// to Rust (plan.md Phase 10, steps 1 and 4). Rust now owns the loop
/// orchestration and output scheduling; Fortran is called as a coarse
/// numerical kernel: one timestep at a time, plus capture calls when
/// Rust decides output is due.
///
/// # Scheduling invariants (mirror `run_time_loop` in `snapwave_c_api.f90`)
///
/// 1. Both output schedules start at `tstart` (first iteration outputs).
/// 2. Output fires when `t >= next_*_output - output_tol`.
/// 3. After output, `next_*_output` advances by the interval in a
///    `do while` loop (handles timestep > interval).
/// 4. History output additionally requires `nobs > 0` and a non-empty
///    `his_filename`.
/// 5. Map output is suppressed when `ja_save_each_iter != 0`.
/// 6. Time advances by `timestep` after output checks.
/// 7. The loop runs while `t <= tstop`.
#[derive(Debug, Clone)]
pub struct ModelState {
    /// Current model time (seconds since reference).
    pub t: f64,
    /// Iteration counter (1-based, incremented at the top of each iteration).
    pub it: i32,
    /// End time; the loop stops when `t > tstop`.
    pub tstop: f64,
    /// Time step (seconds).
    pub timestep: f32,
    /// Next scheduled map output time.
    pub next_map_output: f64,
    /// Next scheduled history output time.
    pub next_his_output: f64,
    /// Map output record counter (1-based, passed to `ncoutput_update_map`).
    pub map_output_count: i32,
    /// History output record counter (1-based, passed to `ncoutput_update_his`).
    pub his_output_count: i32,
    /// `max(1.0d-6, abs(dble(timestep))*1.0d-6)` — the Fortran loop's
    /// floating-point tolerance for output-time comparisons.
    pub output_tol: f64,
}

impl ModelState {
    /// Initialise the scheduling state exactly the way `run_time_loop` does.
    pub fn new(cfg: &SnapWaveInput) -> Self {
        ModelState {
            t: cfg.time.tstart,
            it: 0,
            tstop: cfg.time.tstop,
            timestep: cfg.time.timestep,
            next_map_output: cfg.time.tstart,
            next_his_output: cfg.time.tstart,
            map_output_count: 0,
            his_output_count: 0,
            output_tol: 1.0e-6_f64.max((cfg.time.timestep as f64).abs() * 1.0e-6),
        }
    }

    /// Whether the time loop should continue (`t <= tstop`).
    pub fn is_running(&self) -> bool {
        self.t <= self.tstop
    }

    /// Advance the iteration counter (called at the top of each iteration,
    /// before the solver step).
    pub fn advance_iteration(&mut self) {
        self.it += 1;
    }

    /// Advance model time by one timestep (called at the bottom of each
    /// iteration, after output checks).
    pub fn advance_time(&mut self) {
        self.t += self.timestep as f64;
    }

    /// Whether map output should fire at the current time. Mirrors the
    /// Fortran condition:
    ///   `ja_save_each_iter == 0 .and. map_filename /= '' .and.
    ///    t >= next_map_output - output_tol`
    pub fn should_output_map(&self, ja_save_each_iter: i32, map_file_nonempty: bool) -> bool {
        ja_save_each_iter == 0 && map_file_nonempty && self.t >= self.next_map_output - self.output_tol
    }

    /// Whether history output should fire at the current time. Mirrors the
    /// Fortran condition:
    ///   `his_filename /= '' .and. nobs > 0 .and.
    ///    t >= next_his_output - output_tol`
    pub fn should_output_his(&self, his_file_nonempty: bool, nobs: i32) -> bool {
        his_file_nonempty && nobs > 0 && self.t >= self.next_his_output - self.output_tol
    }

    /// Record a map output event: increment the counter and advance
    /// `next_map_output` by `map_interval` in a `do while` loop (handles
    /// the case where the timestep is larger than the output interval).
    pub fn record_map_output(&mut self, map_interval: f32) {
        self.map_output_count += 1;
        let interval = map_interval as f64;
        while self.next_map_output <= self.t + self.output_tol {
            self.next_map_output += interval;
        }
    }

    /// Record a history output event: increment the counter and advance
    /// `next_his_output` by `his_interval` in a `do while` loop.
    pub fn record_his_output(&mut self, his_interval: f32) {
        self.his_output_count += 1;
        let interval = his_interval as f64;
        while self.next_his_output <= self.t + self.output_tol {
            self.next_his_output += interval;
        }
    }
}

// ----------------------------------------------------------------------
// Domain state
// ----------------------------------------------------------------------

/// The non-solver model data, composed from the Rust-owned parses of
/// Phases 3–7. Everything the Fortran core needs *as data* (rather than
/// as derived geometry or file-backed interpolation) is reachable from
/// here, so new Rust code has no reason to touch `snapwave_data`
/// globals for these subsystems (Phase 8 acceptance).
#[derive(Debug)]
pub struct DomainState<'a> {
    /// The resolved configuration. Read by the Phase 10 output
    /// scheduling / solver-state work; kept here so the state struct is
    /// the single composition point (tests and diagnostics use it).
    #[allow(dead_code)]
    pub config: &'a SnapWaveInput,
    pub mesh: &'a Mesh,
    pub text: &'a ParsedTextInputs,
    pub runtime: RuntimeState,
}

impl<'a> DomainState<'a> {
    pub fn new(config: &'a SnapWaveInput, mesh: &'a Mesh, text: &'a ParsedTextInputs) -> Self {
        DomainState { config, mesh, text, runtime: RuntimeState::new(config) }
    }
}

/// Does the Rust mesh reader own this gridfile? Mirrors the extension
/// dispatch of `initialize_snapwave_domain` *verbatim*, including the
/// quirks: the two characters after the **first** `.` decide (`a.ncd`
/// takes the NetCDF branch), and a name without any `.` compares its
/// first two characters (Fortran `index` returns 0, so `gridfile(1:2)`
/// is inspected).
pub fn rust_owns_gridfile(gridfile: &str) -> bool {
    let after_dot = match gridfile.find('.') {
        Some(j) => &gridfile[j + 1..],
        None => gridfile,
    };
    after_dot.get(..2) == Some("nc")
}

// ----------------------------------------------------------------------
// FFI handoff (plan.md Phase 8, steps 2, 4 and 5)
// ----------------------------------------------------------------------

/// `bind(C)` mirror of `snapwave_state_t` in `src/snapwave_c_api.f90`.
/// Field order and widths must stay identical on both sides; absent
/// data is zero extents + empty buffers (never null pointers).
#[repr(C)]
pub struct SnapWaveStateC {
    // ---- mesh (nc_read_net + the two post-processing steps, Rust-owned)
    pub no_nodes: c_int,
    pub no_faces: c_int,
    pub max_nodes: c_int,
    pub sferic: c_int,
    /// `x` (`real*8`), length `no_nodes`.
    pub x: *const f64,
    /// `y` (`real*8`), length `no_nodes`.
    pub y: *const f64,
    /// `zb` (`real*4`) after `zb = -posdwn*zb`, length `no_nodes`.
    pub zb: *const f32,
    /// `msk` (`integer*1`), length `no_nodes`. **Fortran writes into
    /// this buffer** (enclosure/Neumann mask refinement); see the module
    /// docs.
    pub msk: *const i8,
    /// `face_nodes` (`integer`), shape `(4, no_faces)` column-major —
    /// byte-identical to `Mesh::face_nodes` (pinned in `ffi_layout`).
    pub face_nodes: *const i32,
    // ---- boundary enclosure polyline (`x_bndenc`, `y_bndenc`)
    pub n_bndenc: c_int,
    pub x_bndenc: *const f64,
    pub y_bndenc: *const f64,
    // ---- Neumann polyline (`x_neu`, `y_neu`)
    pub n_neu: c_int,
    pub x_neu: *const f64,
    pub y_neu: *const f64,
    // ---- observation points (`xobs`, `yobs`, `nameobs`)
    pub nobs: c_int,
    pub xobs: *const f64,
    pub yobs: *const f64,
    /// `nameobs` (`character*32`): `nobs` records of 32 blank-padded
    /// bytes, in point order.
    pub names: *const u8,
    // ---- boundary forcing series (`*_bwv`)
    /// 0 = none, 1 = single-point JONSWAP, 2 = space/time-varying.
    pub boundary_mode: c_int,
    pub nwbnd: c_int,
    pub ntwbnd: c_int,
    /// `x_bwv`/`y_bwv` (`real*8`), length `nwbnd` (timeseries mode).
    pub x_bwv: *const f64,
    pub y_bwv: *const f64,
    /// `t_bwv` (`real*4`), length `ntwbnd`.
    pub t_bwv: *const f32,
    /// `hs_bwv`/`tp_bwv`/`zs_bwv` (`real*4`), shape `(nwbnd, ntwbnd)`
    /// column-major — byte-identical to the time-major Rust layout.
    pub hs_bwv: *const f32,
    pub tp_bwv: *const f32,
    /// `wd_bwv` (`real*4`), already converted: `(270 - wd) * deg2rad`.
    pub wd_bwv: *const f32,
    /// `ds_bwv` (`real*4`), already converted: `ds * deg2rad`.
    pub ds_bwv: *const f32,
    pub zs_bwv: *const f32,
}

/// Rust-owned buffers backing [`SnapWaveStateC`]. Build one from a
/// [`DomainState`], keep it alive for the duration of the Fortran call,
/// and hand over `c_state()`.
///
/// The conversions here are the *only* place where the logical Rust
/// state meets Fortran widths and memory order (Phase 8, steps 2–4):
/// `msk` narrows `i32 → integer*1`, the boundary series flatten into
/// `(nwbnd, ntwbnd)` column-major buffers, directions/spreadings are
/// pre-converted with the exact `snapwave_data` recipes, and names pack
/// into `character*32` records.
#[derive(Debug)]
pub struct FfiState {
    // mesh
    x: Vec<f64>,
    y: Vec<f64>,
    zb: Vec<f32>,
    msk: Vec<i8>,
    face_nodes: Vec<i32>,
    no_nodes: c_int,
    no_faces: c_int,
    max_nodes: c_int,
    sferic: c_int,
    // polylines
    n_bndenc: c_int,
    x_bndenc: Vec<f64>,
    y_bndenc: Vec<f64>,
    n_neu: c_int,
    x_neu: Vec<f64>,
    y_neu: Vec<f64>,
    // observation points
    nobs: c_int,
    xobs: Vec<f64>,
    yobs: Vec<f64>,
    names: Vec<u8>,
    // boundary series
    boundary_mode: c_int,
    nwbnd: c_int,
    ntwbnd: c_int,
    x_bwv: Vec<f64>,
    y_bwv: Vec<f64>,
    t_bwv: Vec<f32>,
    hs_bwv: Vec<f32>,
    tp_bwv: Vec<f32>,
    wd_bwv: Vec<f32>,
    ds_bwv: Vec<f32>,
    zs_bwv: Vec<f32>,
}

impl FfiState {
    /// Convert the logical Rust state into Fortran-compatible buffers.
    /// Fails (rather than truncating) on data that cannot cross the
    /// boundary, e.g. mask values outside `integer*1`.
    pub fn build(state: &DomainState) -> Result<Self> {
        let mesh = state.mesh;

        // msk narrows to integer*1 (its Fortran declaration).
        let mut msk = Vec::with_capacity(mesh.no_nodes);
        for &m in &mesh.msk {
            let Ok(narrow) = i8::try_from(m) else {
                bail!("mesh mask value {m} does not fit the Fortran integer*1 global 'msk'");
            };
            msk.push(narrow);
        }

        // face_nodes: (4, no_faces) column-major == node-major within
        // face (ffi_layout pins the equivalence); handed as stored.
        if mesh.face_nodes.len() != 4 * mesh.no_faces {
            bail!(
                "mesh face_nodes has {} elements, expected 4 * no_faces = {}",
                mesh.face_nodes.len(),
                4 * mesh.no_faces
            );
        }

        let (n_bndenc, x_bndenc, y_bndenc) = polyline_buffers(state.text.enclosure.as_ref());
        let (n_neu, x_neu, y_neu) = polyline_buffers(state.text.neumann.as_ref());

        let nobs = state.text.obs.as_ref().map_or(0, |o| o.points.len()) as c_int;
        let (xobs, yobs, names) = match &state.text.obs {
            Some(o) => {
                let mut names = Vec::with_capacity(o.points.len() * WIDTH_NAMEOBS);
                for p in &o.points {
                    // The parser already truncated to character*32; pad
                    // the record with blanks exactly like Fortran
                    // character assignment would.
                    let bytes = p.name.as_bytes();
                    debug_assert!(bytes.len() <= WIDTH_NAMEOBS);
                    names.extend_from_slice(&bytes[..WIDTH_NAMEOBS.min(bytes.len())]);
                    names.resize(names.len() + WIDTH_NAMEOBS - WIDTH_NAMEOBS.min(bytes.len()), b' ');
                }
                (o.points.iter().map(|p| p.x).collect(), o.points.iter().map(|p| p.y).collect(), names)
            }
            None => (Vec::new(), Vec::new(), Vec::new()),
        };

        let (boundary_mode, nwbnd, ntwbnd, x_bwv, y_bwv, t_bwv, hs_bwv, tp_bwv, wd_bwv, ds_bwv, zs_bwv) =
            match &state.text.boundary {
                BoundaryInput::None => (0, 0, 0, Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new(), Vec::new()),
                BoundaryInput::Single(j) => {
                    // hs_bwv(1, nrec): the flat JONSWAP columns already
                    // are the (1, nrec) column-major layout.
                    (
                        1,
                        1,
                        as_c_int(j.len())?,
                        Vec::new(),
                        Vec::new(),
                        j.t.clone(),
                        j.hs.clone(),
                        j.tp.clone(),
                        j.wd_rad(),
                        j.ds_rad(),
                        j.zs.clone(),
                    )
                }
                BoundaryInput::Timeseries(s) => {
                    // (nwbnd, ntwbnd) column-major == time-major rows
                    // (ffi_layout pins the equivalence); handed as stored.
                    (
                        2,
                        as_c_int(s.nwbnd)?,
                        as_c_int(s.ntwbnd)?,
                        s.x.clone(),
                        s.y.clone(),
                        s.t.clone(),
                        s.hs.clone(),
                        s.tp.clone(),
                        s.wd_rad(),
                        s.ds_rad(),
                        s.zs.clone(),
                    )
                }
            };

        Ok(FfiState {
            x: mesh.x.clone(),
            y: mesh.y.clone(),
            zb: mesh.zb.clone(),
            msk,
            face_nodes: mesh.face_nodes.clone(),
            no_nodes: as_c_int(mesh.no_nodes)?,
            no_faces: as_c_int(mesh.no_faces)?,
            max_nodes: as_c_int(mesh.max_nodes)?,
            sferic: mesh.sferic,
            n_bndenc,
            x_bndenc,
            y_bndenc,
            n_neu,
            x_neu,
            y_neu,
            nobs,
            xobs,
            yobs,
            names,
            boundary_mode,
            nwbnd,
            ntwbnd,
            x_bwv,
            y_bwv,
            t_bwv,
            hs_bwv,
            tp_bwv,
            wd_bwv,
            ds_bwv,
            zs_bwv,
        })
    }

    /// The `#[repr(C)]` view handed to `snapwave_run_capture_state_c`.
    /// All pointers stay valid for as long as `self` is alive; empty
    /// buffers intentionally yield dangling-but-aligned non-null
    /// pointers so the Fortran side can associate zero-extent arrays
    /// uniformly.
    pub fn c_state(&self) -> SnapWaveStateC {
        SnapWaveStateC {
            no_nodes: self.no_nodes,
            no_faces: self.no_faces,
            max_nodes: self.max_nodes,
            sferic: self.sferic,
            x: self.x.as_ptr(),
            y: self.y.as_ptr(),
            zb: self.zb.as_ptr(),
            msk: self.msk.as_ptr(),
            face_nodes: self.face_nodes.as_ptr(),
            n_bndenc: self.n_bndenc,
            x_bndenc: self.x_bndenc.as_ptr(),
            y_bndenc: self.y_bndenc.as_ptr(),
            n_neu: self.n_neu,
            x_neu: self.x_neu.as_ptr(),
            y_neu: self.y_neu.as_ptr(),
            nobs: self.nobs,
            xobs: self.xobs.as_ptr(),
            yobs: self.yobs.as_ptr(),
            names: self.names.as_ptr(),
            boundary_mode: self.boundary_mode,
            nwbnd: self.nwbnd,
            ntwbnd: self.ntwbnd,
            x_bwv: self.x_bwv.as_ptr(),
            y_bwv: self.y_bwv.as_ptr(),
            t_bwv: self.t_bwv.as_ptr(),
            hs_bwv: self.hs_bwv.as_ptr(),
            tp_bwv: self.tp_bwv.as_ptr(),
            wd_bwv: self.wd_bwv.as_ptr(),
            ds_bwv: self.ds_bwv.as_ptr(),
            zs_bwv: self.zs_bwv.as_ptr(),
        }
    }
}

/// Polyline buffers: absent polylines cross as zero extent + empty
/// buffers (`n_bndenc = 0` / `n_neu = 0` is exactly what the Fortran
/// readers leave behind for a disabled `encfile`/`neumannfile`).
fn polyline_buffers(poly: Option<&Polyline>) -> (c_int, Vec<f64>, Vec<f64>) {
    match poly {
        Some(p) => (p.len() as c_int, p.x.clone(), p.y.clone()),
        None => (0, Vec::new(), Vec::new()),
    }
}

fn as_c_int(n: usize) -> Result<c_int> {
    let Ok(v) = c_int::try_from(n) else {
        bail!("array length {n} exceeds the 32-bit Fortran global dimension");
    };
    Ok(v)
}

// ----------------------------------------------------------------------
// Tests: focused coverage of the shape/indexing conversions that cross
// the FFI boundary (plan.md Phase 8 acceptance).
// ----------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::input::parse_str;
    use crate::text_input::{
        parse_boundary_timeseries, parse_jonswap, parse_obs_points, parse_polyline, JonswapSeries,
        WindInput,
    };

    fn toy_mesh() -> Mesh {
        Mesh {
            no_nodes: 4,
            no_faces: 2,
            max_nodes: 3,
            sferic: 0,
            x: vec![0.0, 10.0, 0.0, 10.0],
            y: vec![0.0, 0.0, 10.0, 10.0],
            zb: vec![5.0, 4.0, 3.0, 2.0],
            msk: vec![1, 1, 2, 3],
            face_nodes: vec![1, 2, 3, -999, 1, 3, 4, -999],
        }
    }

    #[test]
    fn runtime_state_mirrors_the_fortran_loop_init() {
        let cfg = parse_str("timestep = 12.5\nmap_interval = 3600\n").unwrap();
        let rt = RuntimeState::new(&cfg);
        assert_eq!(rt.timestep, 12.5);
        assert_eq!(rt.map_interval, 3600.0);
        assert_eq!(rt.his_interval, 12.5, "his_interval defaults to the timestep");
        assert_eq!(rt.next_map_output, rt.tstart);
        assert_eq!(rt.next_his_output, rt.tstart);
        assert_eq!(rt.map_output_count, 0);
        assert_eq!(rt.his_output_count, 0);
        // max(1e-6, |timestep|*1e-6): for a small timestep the absolute
        // floor wins; for a large one the relative term.
        let small = RuntimeState::new(&parse_str("timestep = 0.001\n").unwrap());
        assert_eq!(small.output_tol, 1.0e-6);
        let large = RuntimeState::new(&parse_str("timestep = 3600\n").unwrap());
        assert_eq!(large.output_tol, 3600.0 * 1.0e-6);
    }

    #[test]
    fn gridfile_dispatch_mirrors_the_fortran_extension_check() {
        assert!(rust_owns_gridfile("mesh.nc"), "'nc' extension is Rust-owned");
        assert!(rust_owns_gridfile("run/mesh.nc"));
        assert!(rust_owns_gridfile(".nc"), "bare '.nc' works");
        assert!(rust_owns_gridfile("a.ncd"), "quirk: only two chars after the first dot are compared");
        assert!(rust_owns_gridfile("ncmesh"), "quirk: no dot inspects the first two characters");
        assert!(!rust_owns_gridfile("mesh.txt"));
        assert!(!rust_owns_gridfile(".txt"));
        assert!(!rust_owns_gridfile("mesh.NET"));
        assert!(!rust_owns_gridfile("n"));
        assert!(!rust_owns_gridfile("snapshot.vnc"));
        assert!(!rust_owns_gridfile(""));
    }

    #[test]
    fn ffi_state_converts_widths_and_layouts() {
        let cfg = parse_str("timestep = 60\n").unwrap();
        let mesh = toy_mesh();
        let text = ParsedTextInputs {
            obs: Some(parse_obs_points("1.0 2.0 'north'\n3.0 4.0\n").unwrap()),
            boundary: BoundaryInput::None,
            wind: WindInput::Uniform(crate::text_input::parse_uniform_wind("0.0", "270.0")),
            enclosure: Some(parse_polyline("0 0\n10 0\n").unwrap()),
            neumann: None,
        };
        let state = DomainState::new(&cfg, &mesh, &text);
        let ffi = FfiState::build(&state).unwrap();

        assert_eq!(ffi.no_nodes, 4);
        assert_eq!(ffi.no_faces, 2);
        assert_eq!(ffi.msk, vec![1i8, 1, 2, 3], "msk narrows to integer*1");
        assert_eq!(ffi.face_nodes, mesh.face_nodes, "face_nodes layout handed as stored");
        assert_eq!(ffi.n_bndenc, 2);
        assert_eq!(ffi.x_bndenc, vec![0.0, 10.0]);
        assert_eq!(ffi.n_neu, 0);
        assert!(ffi.x_neu.is_empty(), "absent neumann crosses as empty buffer");
        assert_eq!(ffi.nobs, 2);
        assert_eq!(ffi.xobs, vec![1.0, 3.0]);
        // nameobs: two 32-byte blank-padded records.
        assert_eq!(ffi.names.len(), 64);
        assert_eq!(&ffi.names[..5], b"north");
        assert!(ffi.names[5..32].iter().all(|&b| b == b' '), "padded with blanks, not NULs");
        assert_eq!(&ffi.names[32..39], b"station");
        assert_eq!(ffi.boundary_mode, 0);
        assert_eq!(ffi.nwbnd, 0);
        assert_eq!(ffi.ntwbnd, 0);

        let c = ffi.c_state();
        assert!(!c.x.is_null());
        assert!(!c.x_neu.is_null(), "empty buffers still yield non-null pointers");
        // The pointers actually point at the owned buffers.
        assert_eq!(unsafe { *c.msk.add(2) }, 2i8);
    }

    #[test]
    fn ffi_state_rejects_masks_that_do_not_fit_integer1() {
        let cfg = parse_str("").unwrap();
        let mut mesh = toy_mesh();
        mesh.msk[1] = 300;
        let text = ParsedTextInputs {
            obs: None,
            boundary: BoundaryInput::None,
            wind: WindInput::Uniform(crate::text_input::parse_uniform_wind("0.0", "270.0")),
            enclosure: None,
            neumann: None,
        };
        let err = FfiState::build(&DomainState::new(&cfg, &mesh, &text)).expect_err("must reject");
        assert!(format!("{err}").contains("integer*1"), "error was: {err}");
    }

    #[test]
    fn ffi_state_single_point_jonswap() {
        let cfg = parse_str("").unwrap();
        let mesh = toy_mesh();
        let j = parse_jonswap("0 1.0 10.0 270. 30. 0.5\n3600 1.2 10.0 90. 30. 0.5\n").unwrap();
        let text = ParsedTextInputs {
            obs: None,
            boundary: BoundaryInput::Single(j),
            wind: WindInput::Uniform(crate::text_input::parse_uniform_wind("0.0", "270.0")),
            enclosure: None,
            neumann: None,
        };
        let ffi = FfiState::build(&DomainState::new(&cfg, &mesh, &text)).unwrap();
        assert_eq!(ffi.boundary_mode, 1);
        assert_eq!(ffi.nwbnd, 1);
        assert_eq!(ffi.ntwbnd, 2);
        assert!(ffi.x_bwv.is_empty(), "single-point mode has no support-point coordinates");
        assert_eq!(ffi.t_bwv, vec![0.0f32, 3600.0]);
        assert_eq!(ffi.hs_bwv.len(), 2, "hs_bwv(1, nrec) flat");
        // wd already converted: (270-270)*deg2rad = 0 and (270-90)*deg2rad.
        let deg2rad = crate::text_input::deg2rad_f32();
        assert_eq!(ffi.wd_bwv[0], 0.0f32);
        assert_eq!(ffi.wd_bwv[1], (270.0f32 - 90.0f32) * deg2rad);
        assert_eq!(ffi.ds_bwv, vec![30.0f32 * deg2rad, 30.0f32 * deg2rad]);
    }

    #[test]
    fn ffi_state_timeseries_layout_and_conversion() {
        let cfg = parse_str("").unwrap();
        let mesh = toy_mesh();
        let bnd = "10 20\n30 40\n";
        let bhs = "0.0 1.0 2.0\n100.0 3.0 4.0\n";
        let btp = "0.0 5.0 6.0\n100.0 7.0 8.0\n";
        let bwd = "0.0 270.0 250.0\n100.0 270.0 250.0\n";
        let bds = "0.0 10.0 20.0\n100.0 10.0 20.0\n";
        let bzs = "0.0 0.0 0.0\n100.0 1.0 1.0\n";
        let s = parse_boundary_timeseries(bnd, bhs, btp, bwd, bds, bzs).unwrap();
        let text = ParsedTextInputs {
            obs: None,
            boundary: BoundaryInput::Timeseries(s),
            wind: WindInput::Uniform(crate::text_input::parse_uniform_wind("0.0", "270.0")),
            enclosure: None,
            neumann: None,
        };
        let ffi = FfiState::build(&DomainState::new(&cfg, &mesh, &text)).unwrap();
        assert_eq!(ffi.boundary_mode, 2);
        assert_eq!(ffi.nwbnd, 2);
        assert_eq!(ffi.ntwbnd, 2);
        assert_eq!(ffi.x_bwv, vec![10.0, 30.0]);
        assert_eq!(ffi.t_bwv, vec![0.0f32, 100.0]);

        // The time-major Rust buffer must equal the (nwbnd, ntwbnd)
        // column-major Fortran memory: hs_bwv(ib, itb) at (ib-1) +
        // nwbnd*(itb-1) — the equivalence ffi_layout pins.
        let cm = crate::ffi_layout::ColMajor::new(&[2, 2]);
        let expect_hs = |itb: i64, ib: i64, v: f32| {
            let off = cm.offset_fortran(&[ib, itb]).unwrap();
            assert_eq!(ffi.hs_bwv[off], v, "hs_bwv(ib={ib}, itb={itb})");
        };
        expect_hs(1, 1, 1.0);
        expect_hs(1, 2, 2.0);
        expect_hs(2, 1, 3.0);
        expect_hs(2, 2, 4.0);

        let deg2rad = crate::text_input::deg2rad_f32();
        assert_eq!(ffi.wd_bwv[cm.offset_fortran(&[1, 1]).unwrap()], 0.0f32);
        assert_eq!(ffi.wd_bwv[cm.offset_fortran(&[2, 1]).unwrap()], (270.0f32 - 250.0f32) * deg2rad);
        assert_eq!(ffi.ds_bwv[cm.offset_fortran(&[2, 1]).unwrap()], 20.0f32 * deg2rad);
    }

    #[test]
    fn jonswap_conversion_matches_the_fortran_recipe() {
        // read_boundary_data_singlepoint: wd_bwv = (270.0 - wd_bwv) *
        // deg2rad with the single-precision snapwave parameter.
        let j: JonswapSeries = parse_jonswap("0 1 1 135. 45. 0\n").unwrap();
        let deg2rad = crate::text_input::deg2rad_f32();
        assert_eq!(j.wd_rad()[0], (270.0f32 - 135.0f32) * deg2rad);
        assert_eq!(j.ds_rad()[0], 45.0f32 * deg2rad);
    }

    // ------------------------------------------------------------------
    // Phase 10: ModelState scheduling tests
    // ------------------------------------------------------------------

    #[test]
    fn model_state_initialises_like_the_fortran_loop() {
        let cfg = parse_str("timestep = 60\nmap_interval = 3600\nhis_interval = 1800\n").unwrap();
        let m = ModelState::new(&cfg);
        assert_eq!(m.t, cfg.time.tstart);
        assert_eq!(m.it, 0);
        assert_eq!(m.tstop, cfg.time.tstop);
        assert_eq!(m.timestep, 60.0);
        assert_eq!(m.next_map_output, cfg.time.tstart);
        assert_eq!(m.next_his_output, cfg.time.tstart);
        assert_eq!(m.map_output_count, 0);
        assert_eq!(m.his_output_count, 0);
        assert_eq!(m.output_tol, 60.0f64 * 1.0e-6);
        assert!(m.is_running());
    }

    #[test]
    fn model_state_output_tol_floor() {
        // For very small timesteps the absolute floor (1e-6) wins.
        let cfg = parse_str("timestep = 0.001\n").unwrap();
        let m = ModelState::new(&cfg);
        assert_eq!(m.output_tol, 1.0e-6);
    }

    #[test]
    fn model_state_iteration_and_time_advance() {
        let cfg = parse_str("timestep = 10\n").unwrap();
        let mut m = ModelState::new(&cfg);
        assert_eq!(m.it, 0);
        m.advance_iteration();
        assert_eq!(m.it, 1);
        assert_eq!(m.t, cfg.time.tstart);
        m.advance_time();
        assert_eq!(m.t, cfg.time.tstart + 10.0);
    }

    #[test]
    fn model_state_output_predicates() {
        let cfg = parse_str("timestep = 60\nmap_interval = 3600\nhis_interval = 1800\n").unwrap();
        let mut m = ModelState::new(&cfg);

        // At t=tstart, both should fire (next_*_output == tstart).
        assert!(m.should_output_map(0, true), "map output at t=tstart");
        assert!(m.should_output_his(true, 5), "his output at t=tstart with nobs>0");

        // Suppressed conditions.
        assert!(!m.should_output_map(1, true), "ja_save_each_iter suppresses map");
        assert!(!m.should_output_map(0, false), "empty map_file suppresses map");
        assert!(!m.should_output_his(false, 5), "empty his_file suppresses his");
        assert!(!m.should_output_his(true, 0), "nobs==0 suppresses his");

        // After recording, next_output advances.
        m.record_map_output(3600.0);
        assert_eq!(m.map_output_count, 1);
        assert!(m.next_map_output > m.t, "next_map_output advanced past current t");
        assert!(!m.should_output_map(0, true), "map output not due again yet");

        m.record_his_output(1800.0);
        assert_eq!(m.his_output_count, 1);
        assert!(m.next_his_output > m.t);
        assert!(!m.should_output_his(true, 5));
    }

    #[test]
    fn model_state_handles_timestep_larger_than_interval() {
        // timestep=7200, map_interval=3600: the do-while loop should
        // advance next_map_output past t+output_tol so only one record
        // is emitted per iteration.
        let cfg = parse_str("timestep = 7200\nmap_interval = 3600\n").unwrap();
        let mut m = ModelState::new(&cfg);

        // First iteration at t=tstart: output fires.
        assert!(m.should_output_map(0, true));
        m.record_map_output(3600.0);
        assert_eq!(m.map_output_count, 1);
        // next_map_output should now be > tstart + output_tol, so no
        // further output this iteration.
        assert!(!m.should_output_map(0, true));

        // Advance time by one timestep.
        m.advance_time();
        // t = tstart + 7200. next_map_output should be tstart + 3600
        // (or tstart + 7200 if the do-while ran twice). Either way,
        // output should fire again.
        assert!(m.should_output_map(0, true));
    }

    #[test]
    fn model_state_stops_at_tstop() {
        let cfg = parse_str("timestep = 10\n").unwrap();
        let mut m = ModelState::new(&cfg);
        // Manually advance past tstop.
        m.t = m.tstop + 1.0;
        assert!(!m.is_running());
    }
}
