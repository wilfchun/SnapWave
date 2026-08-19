//! Rust side of the temporary Phase 3 Fortran comparison hook
//! (plan.md Phase 3, step 5).
//!
//! [`check`] compares a [`SnapWaveInput`] parse result against the
//! canonical `key=value` dump produced by `snapwave_read_input_dump_c`
//! (`src/snapwave_c_api.f90`), which runs the legacy Fortran reader and
//! dumps every `snapwave_data` global it sets. Both value tables must
//! stay in sync; the comparison reports missing/extra keys so a mismatch
//! is caught by `tests/input_parse.rs`.
//!
//! Dump conventions (mirrored from the Fortran side):
//! integers decimal, real*4/real*8 as zero-padded IEEE-754 bit patterns in
//! hex, logicals 1/0, characters trimmed raw. Reals are compared
//! bit-exact first and within a 1e-6 relative tolerance otherwise: the
//! `sigmin`/`sigmax` defaults involve `atan(1.0)`, whose libm results may
//! differ in the last ulp between the Fortran and Rust runtimes while
//! being numerically identical for configuration purposes.
//!
//! This module — and the facade hook it mirrors — is temporary scaffolding
//! and gets removed once plan.md Phase 4 retires the Fortran input reader
//! from the Rust path.

use std::collections::BTreeMap;

use anyhow::{bail, Result};

use crate::input::SnapWaveInput;

/// Expected value of one dump key, rendered the way the Fortran hook
/// renders it.
enum Expected {
    Int(i32),
    Real4(f32),
    Real8(f64),
    Flag(bool),
    Str(String),
}

/// Compare `cfg` against a Fortran dump; returns the number of values
/// compared. Errors name every mismatching key.
pub fn check(cfg: &SnapWaveInput, dump_text: &str) -> Result<usize> {
    let dump: BTreeMap<&str, &str> = dump_text
        .lines()
        .filter_map(|line| line.split_once('='))
        .map(|(k, v)| (k, v))
        .collect();

    let expected = expected_values(cfg);
    let mut mismatches: Vec<String> = Vec::new();

    for (key, exp) in &expected {
        let Some(raw) = dump.get(key) else {
            mismatches.push(format!("{key}: missing from Fortran dump"));
            continue;
        };
        match exp {
            Expected::Int(v) => {
                if raw.parse::<i32>().as_ref() != Ok(v) {
                    mismatches.push(format!("{key}: rust {v} vs fortran '{raw}'"));
                }
            }
            Expected::Real4(v) => match u32::from_str_radix(raw, 16) {
                Ok(bits) => {
                    let f = f32::from_bits(bits);
                    if !f32_close(*v, f) {
                        mismatches.push(format!(
                            "{key}: rust {} (0x{:08x}) vs fortran {} (0x{raw})",
                            v, v.to_bits(), f
                        ));
                    }
                }
                Err(_) => mismatches.push(format!("{key}: fortran dump value '{raw}' is not a real*4 bit pattern")),
            },
            Expected::Real8(v) => match u64::from_str_radix(raw, 16) {
                Ok(bits) => {
                    let f = f64::from_bits(bits);
                    if !f64_close(*v, f) {
                        mismatches.push(format!(
                            "{key}: rust {} (0x{:016x}) vs fortran {} (0x{raw})",
                            v, v.to_bits(), f
                        ));
                    }
                }
                Err(_) => mismatches.push(format!("{key}: fortran dump value '{raw}' is not a real*8 bit pattern")),
            },
            Expected::Flag(v) => {
                let want = if *v { "1" } else { "0" };
                if *raw != want {
                    mismatches.push(format!("{key}: rust {want} vs fortran '{raw}'"));
                }
            }
            Expected::Str(v) => {
                if *raw != v {
                    mismatches.push(format!("{key}: rust '{v}' vs fortran '{raw}'"));
                }
            }
        }
    }

    for key in dump.keys() {
        if !expected.iter().any(|(k, _)| k == key) {
            mismatches.push(format!("{key}: unexpected in Fortran dump (not compared on the Rust side)"));
        }
    }

    if !mismatches.is_empty() {
        bail!(
            "Rust and Fortran input parses disagree ({} of {} values):\n  {}",
            mismatches.len(),
            expected.len(),
            mismatches.join("\n  ")
        );
    }
    Ok(expected.len())
}

/// The Rust parse result rendered exactly like the Fortran dump (see
/// `snapwave_read_input_dump_c`); keys are the `snapwave_data` global
/// names. Keep in sync with the Fortran dump list.
fn expected_values(cfg: &SnapWaveInput) -> Vec<(&'static str, Expected)> {
    vec![
        // ---- time / control
        ("trefstr", Expected::Str(cfg.time.tref.clone())),
        ("tstartstr", Expected::Str(cfg.time.tstart_str.clone())),
        ("tstopstr", Expected::Str(cfg.time.tstop_str.clone())),
        ("tstart", Expected::Real8(cfg.time.tstart)),
        ("tstop", Expected::Real8(cfg.time.tstop)),
        ("timestep", Expected::Real4(cfg.time.timestep)),
        ("dt", Expected::Real4(cfg.time.dt)),
        ("niter", Expected::Int(cfg.time.niter)),
        ("crit", Expected::Real4(cfg.time.crit)),
        ("restart", Expected::Flag(cfg.time.restart)),
        // ---- grid / domain
        ("mmax", Expected::Int(cfg.grid.mmax)),
        ("nmax", Expected::Int(cfg.grid.nmax)),
        ("dx", Expected::Real4(cfg.grid.dx)),
        ("dy", Expected::Real4(cfg.grid.dy)),
        ("x0", Expected::Real4(cfg.grid.x0)),
        ("y0", Expected::Real4(cfg.grid.y0)),
        ("rotation", Expected::Real4(cfg.grid.rotation)),
        ("posdwn", Expected::Real4(cfg.grid.posdwn)),
        ("sferic", Expected::Int(cfg.grid.sferic)),
        ("dtheta", Expected::Real4(cfg.grid.dtheta)),
        ("sector", Expected::Real4(cfg.grid.sector)),
        ("gridfile", Expected::Str(cfg.grid.gridfile.clone())),
        ("depfile", Expected::Str(cfg.grid.depfile.clone())),
        ("mskfile", Expected::Str(cfg.grid.mskfile.clone())),
        ("indfile", Expected::Str(cfg.grid.indfile.clone())),
        ("upwfile", Expected::Str(cfg.grid.upwfile.clone())),
        // ---- boundary forcing
        ("jonswapfile", Expected::Str(cfg.boundary.jonswapfile.clone())),
        ("bndfile", Expected::Str(cfg.boundary.bndfile.clone())),
        ("encfile", Expected::Str(cfg.boundary.encfile.clone())),
        ("neumannfile", Expected::Str(cfg.boundary.neumannfile.clone())),
        ("bhsfile", Expected::Str(cfg.boundary.bhsfile.clone())),
        ("btpfile", Expected::Str(cfg.boundary.btpfile.clone())),
        ("bwdfile", Expected::Str(cfg.boundary.bwdfile.clone())),
        ("bdsfile", Expected::Str(cfg.boundary.bdsfile.clone())),
        ("bzsfile", Expected::Str(cfg.boundary.bzsfile.clone())),
        ("obsfile", Expected::Str(cfg.boundary.obsfile.clone())),
        ("tol", Expected::Real4(cfg.boundary.tol)),
        // ---- wind
        ("u10str", Expected::Str(cfg.wind.u10.clone())),
        ("u10dirstr", Expected::Str(cfg.wind.u10dir.clone())),
        ("windlistfile", Expected::Str(cfg.wind.windlistfile.clone())),
        ("mwind", Expected::Int(cfg.wind.mwind)),
        ("wind", Expected::Flag(cfg.wind.enabled)),
        // ---- output
        ("map_filename", Expected::Str(cfg.output.map_file.clone())),
        ("his_filename", Expected::Str(cfg.output.his_file.clone())),
        ("map_interval", Expected::Real4(cfg.output.map_interval)),
        ("his_interval", Expected::Real4(cfg.output.his_interval)),
        ("map_dep", Expected::Int(cfg.output.map_depth)),
        ("map_Hm0", Expected::Int(cfg.output.map_Hm0)),
        ("map_Hig", Expected::Int(cfg.output.map_Hig)),
        ("map_Tp", Expected::Int(cfg.output.map_Tp)),
        ("map_dir", Expected::Int(cfg.output.map_dir)),
        ("map_dirspr", Expected::Int(cfg.output.map_dirspr)),
        ("map_cg", Expected::Int(cfg.output.map_cg)),
        ("map_Dw", Expected::Int(cfg.output.map_Dw)),
        ("map_Df", Expected::Int(cfg.output.map_Df)),
        ("map_SwE", Expected::Int(cfg.output.map_SwE)),
        ("map_SwA", Expected::Int(cfg.output.map_SwA)),
        ("map_sig", Expected::Int(cfg.output.map_sig)),
        ("map_u10", Expected::Int(cfg.output.map_u10)),
        ("map_Dveg", Expected::Int(cfg.output.map_Dveg)),
        ("map_ee", Expected::Int(cfg.output.map_ee)),
        ("map_ctheta", Expected::Int(cfg.output.map_ctheta)),
        ("ja_save_each_iter", Expected::Int(cfg.output.ja_save_each_iter)),
        // ---- diagnostics
        ("writetestfiles", Expected::Flag(cfg.diagnostics.writetestfiles)),
        // ---- vegetation
        ("ja_vegetation", Expected::Int(cfg.vegetation.ja_vegetation)),
        ("vegmapfile", Expected::Str(cfg.vegetation.vegmapfile.clone())),
        // ---- solver physics knobs
        ("gamma", Expected::Real4(cfg.physics.gamma)),
        ("alpha", Expected::Real4(cfg.physics.alpha)),
        ("gammax", Expected::Real4(cfg.physics.gammax)),
        ("hmin", Expected::Real4(cfg.physics.hmin)),
        ("fwcutoff", Expected::Real4(cfg.physics.fwcutoff)),
        ("fwstr", Expected::Str(cfg.physics.fw.clone())),
        ("fw_igstr", Expected::Str(cfg.physics.fw_ig.clone())),
        ("Tpini", Expected::Real4(cfg.physics.Tpini)),
        ("zsini", Expected::Real4(cfg.physics.zsini)),
        ("sigmin", Expected::Real4(cfg.physics.sigmin)),
        ("sigmax", Expected::Real4(cfg.physics.sigmax)),
        ("jadcgdx", Expected::Int(cfg.physics.jadcgdx)),
        ("c_dispT", Expected::Real4(cfg.physics.c_dispT)),
        ("ig", Expected::Int(cfg.physics.ig)),
        ("upwindref", Expected::Int(cfg.physics.upwindref)),
    ]
}

/// Bit-exact or within a tiny relative tolerance (libm `atan` ulp drift on
/// the sigmin/sigmax defaults; see module docs).
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
