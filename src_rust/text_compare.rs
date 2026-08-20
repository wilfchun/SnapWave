//! Rust side of the Phase 6 comparison hook (plan.md Phase 6, step 2).
//!
//! [`check`] compares the Rust-parsed auxiliary text inputs
//! ([`crate::text_input::ParsedTextInputs`]) against the canonical dump
//! produced by the Fortran facade hook `snapwave_text_dump_c` (see
//! `src/snapwave_c_api.f90`).
//!
//! # Dump format
//!
//! The dump is line-oriented: `section <name>` blocks of `key value` lines.
//! Array keys (`x`, `y`, `name`, `t`, `hs`, `tp`, `wd`, `ds`, `zs`, `u10`,
//! `u10dir`) are followed by a count and then one value per line; reals are
//! IEEE-754 bit patterns in zero-padded hex (`real*8` = 16 hex digits,
//! `real*4` = 8), names are trimmed text. Scalar keys (`n`, `mode`,
//! `nwbnd`, `ntwbnd`, `ntu10bnd`) carry a single value.
//!
//! Reals are compared with a small relative tolerance (bit-exactness is not
//! expected: the `wd`/`ds` conversions involve `atan`/`deg2rad` arithmetic
//! whose last-ulp result can differ between the Fortran and Rust runtimes
//! while being numerically identical for parsing purposes). Integers and
//! names are compared exactly.

use std::collections::BTreeMap;

use anyhow::{bail, Context, Result};

use crate::text_input::{
    BoundaryInput, ParsedTextInputs, Polyline, WindInput,
};

/// One dump section: variable name -> raw value lines (scalars have one
/// line, arrays have `count` lines).
type Section = BTreeMap<String, Vec<String>>;
type Dump = BTreeMap<String, Section>;

/// Array keys in the dump (followed by a count and one value per line).
const ARRAY_KEYS: [&str; 11] = ["x", "y", "name", "t", "hs", "tp", "wd", "ds", "zs", "u10", "u10dir"];

/// Compare `rust` against the Fortran dump; returns the number of scalar /
/// element values compared. Errors name every mismatching field.
pub fn check(rust: &ParsedTextInputs, dump_text: &str) -> Result<usize> {
    let dump = parse_dump(dump_text)?;
    let mut mismatches: Vec<String> = Vec::new();
    let mut compared = 0usize;

    // ---- observation points -------------------------------------------------
    let obs = section(&dump, "obs");
    let n_obs = expect_scalar(obs, "n")?;
    match &rust.obs {
        None => {
            if n_obs != 0 {
                mismatches.push(format!("obs: rust has no observation file but Fortran read {n_obs} points"));
            }
            compared += 1;
        }
        Some(points) => {
            compare_len(&mut mismatches, "obs", points.len(), n_obs);
            compared += 1;
            if n_obs == points.len() {
                compared += cmp_f64(&mut mismatches, "obs.x", &points.x(), dump_f64(obs, "x")?);
                compared += cmp_f64(&mut mismatches, "obs.y", &points.y(), dump_f64(obs, "y")?);
                compared += cmp_str(&mut mismatches, "obs.name", &points.names(), dump_str(obs, "name")?);
            }
        }
    }

    // ---- boundary conditions ------------------------------------------------
    let bnd = section(&dump, "boundary");
    let mode = expect_scalar_str(bnd, "mode")?;
    match (&rust.boundary, mode) {
        (BoundaryInput::None, "none") => {
            compared += 1;
        }
        (BoundaryInput::None, other) => {
            mismatches.push(format!("boundary.mode: rust has no boundary but Fortran reported '{other}'"));
        }
        (BoundaryInput::Single(j), "single") => {
            let nwbnd = expect_scalar(bnd, "nwbnd")?;
            let ntwbnd = expect_scalar(bnd, "ntwbnd")?;
            compared += 2;
            if nwbnd != 1 {
                mismatches.push(format!("boundary.nwbnd: rust single-point expects 1, Fortran reported {nwbnd}"));
            }
            compare_len(&mut mismatches, "boundary.ntwbnd", j.len(), ntwbnd);
            if ntwbnd == j.len() {
                compared += cmp_f32(&mut mismatches, "boundary.t", &j.t, dump_f32(bnd, "t")?);
                compared += cmp_f32(&mut mismatches, "boundary.hs", &j.hs, dump_f32(bnd, "hs")?);
                compared += cmp_f32(&mut mismatches, "boundary.tp", &j.tp, dump_f32(bnd, "tp")?);
                compared += cmp_f32(&mut mismatches, "boundary.wd", &j.wd_rad(), dump_f32(bnd, "wd")?);
                compared += cmp_f32(&mut mismatches, "boundary.ds", &j.ds_rad(), dump_f32(bnd, "ds")?);
                compared += cmp_f32(&mut mismatches, "boundary.zs", &j.zs, dump_f32(bnd, "zs")?);
            }
        }
        (BoundaryInput::Timeseries(s), "timeseries") => {
            let nwbnd = expect_scalar(bnd, "nwbnd")?;
            let ntwbnd = expect_scalar(bnd, "ntwbnd")?;
            compared += 2;
            compare_len(&mut mismatches, "boundary.nwbnd", s.nwbnd, nwbnd);
            compare_len(&mut mismatches, "boundary.ntwbnd", s.ntwbnd, ntwbnd);
            if nwbnd == s.nwbnd && ntwbnd == s.ntwbnd {
                compared += cmp_f64(&mut mismatches, "boundary.x", &s.x, dump_f64(bnd, "x")?);
                compared += cmp_f64(&mut mismatches, "boundary.y", &s.y, dump_f64(bnd, "y")?);
                compared += cmp_f32(&mut mismatches, "boundary.t", &s.t, dump_f32(bnd, "t")?);
                compared += cmp_f32(&mut mismatches, "boundary.hs", &s.hs, dump_f32(bnd, "hs")?);
                compared += cmp_f32(&mut mismatches, "boundary.tp", &s.tp, dump_f32(bnd, "tp")?);
                compared += cmp_f32(&mut mismatches, "boundary.wd", &s.wd_rad(), dump_f32(bnd, "wd")?);
                compared += cmp_f32(&mut mismatches, "boundary.ds", &s.ds_rad(), dump_f32(bnd, "ds")?);
                compared += cmp_f32(&mut mismatches, "boundary.zs", &s.zs, dump_f32(bnd, "zs")?);
            }
        }
        (BoundaryInput::Single(_), other) => {
            mismatches.push(format!("boundary.mode: rust parsed a single-point JONSWAP file but Fortran reported '{other}'"));
        }
        (BoundaryInput::Timeseries(_), other) => {
            mismatches.push(format!("boundary.mode: rust parsed time-series files but Fortran reported '{other}'"));
        }
    }

    // ---- wind ----------------------------------------------------------------
    let wind = section(&dump, "wind");
    let ntu10bnd = expect_scalar(wind, "ntu10bnd")?;
    let wmode = expect_scalar_str(wind, "mode")?;
    compared += 1;
    match (&rust.wind, wmode) {
        (WindInput::Uniform(u), "uniform") => {
            if ntu10bnd != 1 {
                mismatches.push(format!("wind.ntu10bnd: rust uniform expects 1, Fortran reported {ntu10bnd}"));
            }
            let d_u10 = dump_f32(wind, "u10")?;
            let d_u10dir = dump_f32(wind, "u10dir")?;
            // File-backed wind (None) cannot be compared without Phase 9
            // interpolation; skip those values.
            if let Some(u10) = u.u10 {
                compared += cmp_f32(&mut mismatches, "wind.u10", &[u10], d_u10);
            }
            if let Some(dir_deg) = u.u10dir_deg {
                let dir_rad = (270.0f32 - dir_deg) * crate::text_input::deg2rad_f32();
                compared += cmp_f32(&mut mismatches, "wind.u10dir", &[dir_rad], d_u10dir);
            }
        }
        (WindInput::List(list), "list") => {
            compare_len(&mut mismatches, "wind.ntu10bnd", list.len(), ntu10bnd);
            if ntu10bnd == list.len() {
                compared += cmp_f32(&mut mismatches, "wind.t", &list.t(), dump_f32(wind, "t")?);
            }
        }
        (WindInput::Uniform(_), other) => {
            mismatches.push(format!("wind.mode: rust parsed uniform wind but Fortran reported '{other}'"));
        }
        (WindInput::List(_), other) => {
            mismatches.push(format!("wind.mode: rust parsed a wind list but Fortran reported '{other}'"));
        }
    }

    // ---- enclosure + neumann polylines ---------------------------------------
    compared += cmp_polyline(&mut mismatches, "enc", rust.enclosure.as_ref(), section(&dump, "enc"))?;
    compared += cmp_polyline(&mut mismatches, "neu", rust.neumann.as_ref(), section(&dump, "neu"))?;

    if !mismatches.is_empty() {
        bail!(
            "Rust and Fortran text-input parses disagree ({} of {} values):\n  {}",
            mismatches.len(),
            compared,
            mismatches.join("\n  ")
        );
    }
    Ok(compared)
}

impl crate::text_input::ObsPoints {
    fn x(&self) -> Vec<f64> {
        self.points.iter().map(|p| p.x).collect()
    }
    fn y(&self) -> Vec<f64> {
        self.points.iter().map(|p| p.y).collect()
    }
    fn names(&self) -> Vec<String> {
        self.points.iter().map(|p| p.name.clone()).collect()
    }
}

fn cmp_polyline(
    mismatches: &mut Vec<String>,
    label: &str,
    rust: Option<&Polyline>,
    sec: &Section,
) -> Result<usize> {
    let n = expect_scalar(sec, "n")?;
    let mut compared = 1;
    match rust {
        None => {
            if n != 0 {
                mismatches.push(format!("{label}: rust has no polyline but Fortran read {n} points"));
            }
        }
        Some(p) => {
            compare_len(mismatches, &format!("{label}.n"), p.len(), n);
            if n == p.len() {
                compared += cmp_f64(mismatches, &format!("{label}.x"), &p.x, dump_f64(sec, "x")?);
                compared += cmp_f64(mismatches, &format!("{label}.y"), &p.y, dump_f64(sec, "y")?);
            }
        }
    }
    Ok(compared)
}

fn compare_len(mismatches: &mut Vec<String>, label: &str, rust: usize, fortran: usize) {
    if rust != fortran {
        mismatches.push(format!("{label}: rust {rust} vs fortran {fortran}"));
    }
}

fn cmp_str(mismatches: &mut Vec<String>, label: &str, rust: &[String], fortran: Vec<&str>) -> usize {
    let n = rust.len().min(fortran.len());
    let mut compared = 0;
    for i in 0..n {
        compared += 1;
        if rust[i] != fortran[i] {
            mismatches.push(format!("{label}[{i}]: rust '{}' vs fortran '{}'", rust[i], fortran[i]));
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

fn section<'a>(dump: &'a Dump, name: &str) -> &'a Section {
    dump.get(name).unwrap_or(&EMPTY)
}

// A shared empty section so missing sections fall through to "missing key"
// errors rather than panicking.
static EMPTY: Section = BTreeMap::new();

fn expect_scalar(sec: &Section, key: &str) -> Result<usize> {
    let Some(vals) = sec.get(key) else {
        bail!("dump is missing scalar '{key}'");
    };
    if vals.is_empty() {
        bail!("dump scalar '{key}' has no value");
    }
    vals[0].parse::<usize>().with_context(|| format!("dump scalar '{key}' = '{}' is not an integer", vals[0]))
}

fn expect_scalar_str<'a>(sec: &'a Section, key: &str) -> Result<&'a str> {
    let Some(vals) = sec.get(key) else {
        bail!("dump is missing scalar '{key}'");
    };
    if vals.is_empty() {
        bail!("dump scalar '{key}' has no value");
    }
    Ok(vals[0].as_str())
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

fn dump_str<'a>(sec: &'a Section, key: &str) -> Result<Vec<&'a str>> {
    let Some(vals) = sec.get(key) else {
        bail!("dump is missing array '{key}'");
    };
    Ok(vals.iter().map(|s| s.as_str()).collect())
}

/// Parse the dump text into sections (see the module docs).
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
            continue; // before any section: ignore
        };
        let Some((key, rest)) = line.split_once(' ') else {
            continue; // line without "key value" form: ignore
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
