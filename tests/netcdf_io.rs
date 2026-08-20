//! Phase 7 integration tests: NetCDF mesh input (plan.md Phase 7, step 2).
//!
//! Runs the wrapper's `--compare-mesh` mode on curated testcases: the mode
//! reads the mesh NetCDF in Rust (`src_rust/mesh.rs`) and compares it
//! against the unchanged Fortran `nc_read_net` reader through the temporary
//! `snapwave_mesh_dump_c` hook, so exit 0 pins that the Rust port agrees
//! with the numerical oracle for real, checked-in meshes (node coordinates,
//! bathymetry, mask, face connectivity and the `sferic` fix).
//!
//! The map/history writers (step 3-4) are pinned end-to-end by
//! `tests/regression.rs`: the wrapper now writes those files itself
//! (capture + `src_rust/output.rs`), and the regression suite compares them
//! against the committed baselines and the live Fortran oracle.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

fn stdout_of(out: &Output) -> String {
    String::from_utf8_lossy(&out.stdout).into_owned()
}

fn stderr_of(out: &Output) -> String {
    String::from_utf8_lossy(&out.stderr).into_owned()
}

fn wrapper() -> Command {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_snapwave"));
    cmd.env("OMP_NUM_THREADS", "1");
    cmd
}

/// Recursively copy a directory tree (used to prepare a temp run dir so
/// checked-in inputs are never modified).
fn copy_tree(src: &Path, dst: &Path) {
    fs::create_dir_all(dst).expect("create copy dir");
    for entry in fs::read_dir(src).expect("read source dir") {
        let entry = entry.expect("dir entry");
        let from = entry.path();
        let to = dst.join(entry.file_name());
        if from.is_dir() {
            copy_tree(&from, &to);
        } else {
            fs::copy(&from, &to).expect("copy file");
        }
    }
}

/// Replace `\` with `/` in a copied `*.inp` (testcases are authored on
/// Windows; the Linux readers only accept `/`).
fn normalize_inp(path: &Path) {
    let text = fs::read_to_string(path).expect("read inp");
    let normalized = text.replace('\\', "/");
    if normalized != text {
        fs::write(path, normalized).expect("write normalized inp");
    }
}

/// Run `--compare-mesh` on a copied run directory; returns the temp dir so
/// failures can be inspected (removed on success).
fn run_compare_mesh(case: &str, run_src: &str, inp_rel: &str) -> PathBuf {
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let src = manifest.join(run_src);
    assert!(src.is_dir(), "run directory not found: {}", src.display());

    let work = std::env::temp_dir().join(format!(
        "snapwave_mesh_{}_{}_{}",
        case,
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_micros())
            .unwrap_or(0)
    ));
    if work.exists() {
        fs::remove_dir_all(&work).expect("clean stale work dir");
    }
    copy_tree(&src, &work);

    let inp = work.join(inp_rel);
    assert!(inp.is_file(), "input not found in copy: {}", inp.display());
    normalize_inp(&inp);

    let out = wrapper()
        .current_dir(inp.parent().expect("inp parent"))
        .arg("--compare-mesh")
        .arg(&inp)
        .output()
        .expect("failed to run wrapper binary");

    assert_eq!(
        out.status.code(),
        Some(0),
        "--compare-mesh failed for {case} ({})\n--- stdout ---\n{}\n--- stderr ---\n{}",
        inp.display(),
        stdout_of(&out),
        stderr_of(&out)
    );

    let _ = fs::remove_dir_all(&work);
    work
}

#[test]
fn rust_mesh_reader_matches_fortran_for_checked_in_cases() {
    run_compare_mesh("case31", "testcases/31_linear_shoaling_refraction/run", "coarse/SnapWave.inp");
    run_compare_mesh("case32", "testcases/32_curvi_island/run", "SnapWave.inp");
    run_compare_mesh("case33", "testcases/33_circle_reef/run", "SnapWave.inp");
    run_compare_mesh("case45", "testcases/45_haringvliet/run/snapwave", "snapwave.inp");
}

#[test]
fn compare_mesh_requires_its_input_file() {
    // Missing input path must be a wrapper usage error (exit 2), not a panic.
    let out = wrapper().arg("--compare-mesh").output().expect("run wrapper");
    assert_eq!(out.status.code(), Some(2), "missing input must exit 2");
    assert!(stderr_of(&out).contains("missing input"), "stderr was:\n{}", stderr_of(&out));
}
