//! Rust side of the Phase 9 geometry comparison hook (plan.md Phase 9).
//!
//! [`check`] compares the Rust-computed derived geometry — surrounding
//! points and upwind neighbours (`crate::geometry`), observation
//! interpolation weights (`crate::interp::make_map_fm`) and boundary
//! support-point mapping (`crate::geometry::find_boundary_indices`) —
//! against the canonical dump produced by the Fortran facade hook
//! `snapwave_geometry_dump_c` (see `src/snapwave_c_api.f90`).
//!
//! # Dump format
//!
//! Line-oriented `section <name>` blocks of `key value` lines, exactly like
//! the Phase 6 text dump (`crate::text_compare`): array keys are followed by
//! a count and then one value per line, reals as IEEE-754 bit patterns
//! (`real*8` = 16 hex digits, `real*4` = 8), integers decimal, scalars one
//! value. The array order is the Fortran column-major do-loop order and the
//! Rust layouts are chosen to match it byte-for-byte, so no reshuffling is
//! needed here (see `crate::ffi_layout` for the underlying layout facts).
//!
//! Reals are compared bit-exact first and within a small relative tolerance
//! otherwise: the geometry involves `hypot`/`tan`/`atan2`/`cos`/`sin`, whose
//! libm results can differ in the last ulp between the Fortran and Rust
//! runtimes while being numerically identical (the same rationale as the
//! Phase 3/4 and Phase 6 comparisons). Integers are compared exactly — an
//! index or mask that differs by even one is a real divergence.

use std::collections::BTreeMap;

use anyhow::{bail, Context, Result};

use crate::geometry::{self, BoundaryIndices, DomainGeometry};
use crate::input::SnapWaveInput;
use crate::interp::{self, MapFmResult};
use crate::mesh::Mesh;
use crate::text_input::{BoundaryInput, ParsedTextInputs};

/// Everything the comparison needs from the Rust side, computed by
/// [`compute_geometry`].
#[derive(Debug)]
pub struct GeometryState {
    pub domain: DomainGeometry,
    pub obs: Option<MapFmResult>,
    pub boundary: BoundaryIndices,
    /// `nwbnd` of the boundary section (0 = none, 1 = single point,
    /// >1 = space/time-varying support points).
    pub nwbnd: usize,
}

/// Compute the Rust-side geometry from the Rust-owned mesh and text inputs,
/// mirroring `initialize_snapwave_domain` + `read_obs_points` +
/// `read_boundary_data` (geometry only).
pub fn compute_geometry(mesh: &Mesh, config: &SnapWaveInput, text: &ParsedTextInputs) -> GeometryState {
    let enclosure = text.enclosure.as_ref().map(|p| (p.x.as_slice(), p.y.as_slice()));
    let neumann = text.neumann.as_ref().map(|p| (p.x.as_slice(), p.y.as_slice()));

    let domain = geometry::compute_domain_geometry(
        &mesh.x,
        &mesh.y,
        &mesh.zb,
        mesh.sferic,
        &mesh.face_nodes,
        mesh.no_faces,
        config.grid.dtheta,
        config.boundary.tol,
        enclosure,
        neumann,
        &mesh.msk,
    );

    let obs = text.obs.as_ref().map(|o| {
        let xobs: Vec<f64> = o.points.iter().map(|p| p.x).collect();
        let yobs: Vec<f64> = o.points.iter().map(|p| p.y).collect();
        interp::make_map_fm(&mesh.x, &mesh.y, &mesh.face_nodes, mesh.no_faces, &xobs, &yobs)
    });

    let (nwbnd, boundary) = match &text.boundary {
        BoundaryInput::Timeseries(s) => (
            s.nwbnd,
            geometry::find_boundary_indices(&mesh.x, &mesh.y, &domain.msk, &s.x, &s.y, s.nwbnd),
        ),
        BoundaryInput::Single(_) => (
            1,
            geometry::find_boundary_indices(&mesh.x, &mesh.y, &domain.msk, &[], &[], 1),
        ),
        BoundaryInput::None => {
            (0, geometry::find_boundary_indices(&mesh.x, &mesh.y, &domain.msk, &[], &[], 0))
        }
    };

    GeometryState { domain, obs, boundary, nwbnd }
}

// ----------------------------------------------------------------------
// Dump parsing and comparison
// ----------------------------------------------------------------------

type Section = BTreeMap<String, Vec<String>>;
type Dump = BTreeMap<String, Section>;

/// Array keys in the geometry dump (followed by a count and one value per
/// line). Must match `dump_geometry_globals` in `src/snapwave_c_api.f90`.
const ARRAY_KEYS: [&str; 16] = [
    "kp",
    "dhdx",
    "dhdy",
    "w360",
    "prev360",
    "ds360",
    "msk",
    "neumannconnected",
    "nmindbnd",
    "neubnd",
    "wobs",
    "irefobs",
    "nrefobs",
    "ind1_bwv_cst",
    "ind2_bwv_cst",
    "fac_bwv_cst",
];

/// Compare `state` against the Fortran dump; returns the number of scalar /
/// element values compared. Errors name every mismatching field.
pub fn check(state: &GeometryState, dump_text: &str) -> Result<usize> {
    let dump = parse_dump(dump_text)?;
    let mut mismatches: Vec<String> = Vec::new();
    let mut compared = 0usize;

    // ---- domain --------------------------------------------------------------
    let d = section(&dump, "domain");
    compared += cmp_scalar_usize(&mut mismatches, "domain.no_nodes", state.domain.no_nodes, expect_scalar(d, "no_nodes")?);
    compared += cmp_scalar_usize(&mut mismatches, "domain.np", state.domain.np, expect_scalar(d, "np")?);
    compared +=
        cmp_scalar_usize(&mut mismatches, "domain.ntheta360", state.domain.ntheta360, expect_scalar(d, "ntheta360")?);
    compared += cmp_scalar_usize(&mut mismatches, "domain.nb", state.domain.nb, expect_scalar(d, "nb")?);
    compared += cmp_scalar_usize(&mut mismatches, "domain.nnmb", state.domain.nnmb, expect_scalar(d, "nnmb")?);

    compared += cmp_i32(&mut mismatches, "domain.kp", &state.domain.kp, dump_i32(d, "kp")?);
    compared += cmp_f32(&mut mismatches, "domain.dhdx", &state.domain.dhdx, dump_f32(d, "dhdx")?);
    compared += cmp_f32(&mut mismatches, "domain.dhdy", &state.domain.dhdy, dump_f32(d, "dhdy")?);
    compared += cmp_f32(&mut mismatches, "domain.w360", &state.domain.w360, dump_f32(d, "w360")?);
    compared += cmp_i32(&mut mismatches, "domain.prev360", &state.domain.prev360, dump_i32(d, "prev360")?);
    compared += cmp_f32(&mut mismatches, "domain.ds360", &state.domain.ds360, dump_f32(d, "ds360")?);
    compared += cmp_i32(&mut mismatches, "domain.msk", &state.domain.msk, dump_i32(d, "msk")?);
    compared += cmp_i32(
        &mut mismatches,
        "domain.neumannconnected",
        &state.domain.neumannconnected,
        dump_i32(d, "neumannconnected")?,
    );
    compared += cmp_i32(&mut mismatches, "domain.nmindbnd", &state.domain.nmindbnd, dump_i32(d, "nmindbnd")?);
    if state.domain.nnmb > 0 {
        compared += cmp_i32(&mut mismatches, "domain.neubnd", &state.domain.neubnd, dump_i32(d, "neubnd")?);
    }

    // ---- obs -----------------------------------------------------------------
    let o = section(&dump, "obs");
    let n_obs = expect_scalar(o, "nobs")?;
    compared += 1;
    match &state.obs {
        None => {
            if n_obs != 0 {
                mismatches.push(format!("obs.nobs: rust has no observation points but Fortran computed {n_obs}"));
            }
        }
        Some(map) => {
            let n2 = map.iref.len() / 4;
            if n_obs != n2 {
                mismatches.push(format!("obs.nobs: rust {n2} vs fortran {n_obs}"));
            }
            if n_obs == n2 && n2 > 0 {
                compared += cmp_f64(&mut mismatches, "obs.wobs", &map.w, dump_f64(o, "wobs")?);
                compared += cmp_i32(&mut mismatches, "obs.irefobs", &map.iref, dump_i32(o, "irefobs")?);
                compared += cmp_i32(&mut mismatches, "obs.nrefobs", &map.nref, dump_i32(o, "nrefobs")?);
            }
        }
    }

    // ---- boundary -------------------------------------------------------------
    let b = section(&dump, "boundary");
    compared += cmp_scalar_usize(&mut mismatches, "boundary.nwbnd", state.nwbnd, expect_scalar(b, "nwbnd")?);
    compared += cmp_scalar_usize(&mut mismatches, "boundary.nb", state.domain.nb, expect_scalar(b, "nb")?);
    if state.nwbnd > 0 && state.domain.nb > 0 {
        compared += cmp_i32(&mut mismatches, "boundary.ind1", &state.boundary.ind1, dump_i32(b, "ind1_bwv_cst")?);
        compared += cmp_i32(&mut mismatches, "boundary.ind2", &state.boundary.ind2, dump_i32(b, "ind2_bwv_cst")?);
        compared += cmp_f32(&mut mismatches, "boundary.fac", &state.boundary.fac, dump_f32(b, "fac_bwv_cst")?);
    }

    if !mismatches.is_empty() {
        bail!(
            "Rust and Fortran geometry disagree ({} of {} values):\n  {}",
            mismatches.len(),
            compared,
            mismatches.join("\n  ")
        );
    }
    Ok(compared)
}

fn section<'a>(dump: &'a Dump, name: &str) -> &'a Section {
    dump.get(name).unwrap_or(&EMPTY)
}

static EMPTY: Section = BTreeMap::new();

fn cmp_scalar_usize(mismatches: &mut Vec<String>, label: &str, rust: usize, fortran: usize) -> usize {
    if rust != fortran {
        mismatches.push(format!("{label}: rust {rust} vs fortran {fortran}"));
    }
    1
}

fn cmp_i32(mismatches: &mut Vec<String>, label: &str, rust: &[i32], fortran: Vec<i32>) -> usize {
    let n = rust.len().min(fortran.len());
    let mut compared = 0;
    for i in 0..n {
        compared += 1;
        if rust[i] != fortran[i] {
            mismatches.push(format!("{label}[{i}]: rust {} vs fortran {}", rust[i], fortran[i]));
        }
    }
    if rust.len() != fortran.len() {
        mismatches.push(format!("{label}: length {} vs {}", rust.len(), fortran.len()));
    }
    compared
}

fn cmp_f64(mismatches: &mut Vec<String>, label: &str, rust: &[f64], fortran: Vec<f64>) -> usize {
    let n = rust.len().min(fortran.len());
    let mut compared = 0;
    for i in 0..n {
        compared += 1;
        if !f64_close(rust[i], fortran[i]) {
            mismatches.push(format!(
                "{label}[{i}]: rust {} (0x{:016x}) vs fortran {} (0x{:016x})",
                rust[i],
                rust[i].to_bits(),
                fortran[i],
                fortran[i].to_bits()
            ));
        }
    }
    if rust.len() != fortran.len() {
        mismatches.push(format!("{label}: length {} vs {}", rust.len(), fortran.len()));
    }
    compared
}

fn cmp_f32(mismatches: &mut Vec<String>, label: &str, rust: &[f32], fortran: Vec<f32>) -> usize {
    let n = rust.len().min(fortran.len());
    let mut compared = 0;
    for i in 0..n {
        compared += 1;
        if !f32_close(rust[i], fortran[i]) {
            mismatches.push(format!(
                "{label}[{i}]: rust {} (0x{:08x}) vs fortran {} (0x{:08x})",
                rust[i],
                rust[i].to_bits(),
                fortran[i],
                fortran[i].to_bits()
            ));
        }
    }
    if rust.len() != fortran.len() {
        mismatches.push(format!("{label}: length {} vs {}", rust.len(), fortran.len()));
    }
    compared
}

fn expect_scalar(sec: &Section, key: &str) -> Result<usize> {
    let Some(vals) = sec.get(key) else {
        bail!("dump is missing scalar '{key}'");
    };
    if vals.is_empty() {
        bail!("dump scalar '{key}' has no value");
    }
    vals[0].parse::<usize>().with_context(|| format!("dump scalar '{key}' = '{}' is not an integer", vals[0]))
}

fn dump_f64(sec: &Section, key: &str) -> Result<Vec<f64>> {
    let Some(vals) = sec.get(key) else {
        bail!("dump is missing array '{key}'");
    };
    vals.iter()
        .map(|v| {
            u64::from_str_radix(v.trim(), 16)
                .map(f64::from_bits)
                .with_context(|| format!("dump array '{key}' value '{v}' is not a real*8 bit pattern"))
        })
        .collect()
}

fn dump_f32(sec: &Section, key: &str) -> Result<Vec<f32>> {
    let Some(vals) = sec.get(key) else {
        bail!("dump is missing array '{key}'");
    };
    vals.iter()
        .map(|v| {
            u32::from_str_radix(v.trim(), 16)
                .map(f32::from_bits)
                .with_context(|| format!("dump array '{key}' value '{v}' is not a real*4 bit pattern"))
        })
        .collect()
}

fn dump_i32(sec: &Section, key: &str) -> Result<Vec<i32>> {
    let Some(vals) = sec.get(key) else {
        bail!("dump is missing array '{key}'");
    };
    vals.iter()
        .map(|v| {
            v.trim()
                .parse::<i32>()
                .with_context(|| format!("dump array '{key}' value '{v}' is not an integer"))
        })
        .collect()
}

/// Parse the dump text into sections (same grammar as `text_compare`).
fn parse_dump(dump_text: &str) -> Result<Dump> {
    let mut dump: Dump = BTreeMap::new();
    let mut current: Option<String> = None;
    let mut lines = dump_text.lines().peekable();

    while let Some(line) = lines.next() {
        let line = line.trim_end_matches('\r');
        if line.trim().is_empty() {
            continue;
        }
        if let Some(name) = line.strip_prefix("section ") {
            current = Some(name.trim().to_string());
            dump.entry(name.trim().to_string()).or_default();
            continue;
        }
        let Some(sec_name) = &current else {
            continue;
        };
        let Some((key, rest)) = line.split_once(' ') else {
            continue;
        };
        let key = key.trim();
        let rest = rest.trim();
        let section = dump.get_mut(sec_name).expect("section registered above");
        if ARRAY_KEYS.contains(&key) {
            let count: usize = rest
                .parse()
                .with_context(|| format!("dump array '{key}' count '{rest}' is not an integer"))?;
            let mut vals = Vec::with_capacity(count);
            for _ in 0..count {
                let Some(v) = lines.next() else {
                    bail!("dump truncated in section '{sec_name}' array '{key}'");
                };
                vals.push(v.trim_end_matches('\r').to_string());
            }
            section.insert(key.to_string(), vals);
        } else {
            section.insert(key.to_string(), vec![rest.to_string()]);
        }
    }
    Ok(dump)
}

fn f32_close(a: f32, b: f32) -> bool {
    if a.to_bits() == b.to_bits() {
        return true;
    }
    let scale = a.abs().max(b.abs()).max(1e-30);
    (a - b).abs() <= 1e-6 * scale
}

fn f64_close(a: f64, b: f64) -> bool {
    if a.to_bits() == b.to_bits() {
        return true;
    }
    let scale = a.abs().max(b.abs()).max(1e-300);
    (a - b).abs() <= 1e-9 * scale
}
