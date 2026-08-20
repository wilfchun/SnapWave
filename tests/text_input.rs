//! Phase 6 integration tests: auxiliary text-input readers in Rust
//! (plan.md Phase 6).
//!
//! For each curated testcase, the test copies the run directory to a temp
//! location (normalizing Windows `\` separators in the copied `*.inp`), runs
//! the wrapper's `--compare-text` mode there, and asserts it exits 0. The
//! mode parses the auxiliary text files in Rust and compares them against
//! the unchanged Fortran readers through the temporary `snapwave_text_dump_c`
//! hook, so exit 0 pins that the Rust parsers agree with the numerical
//! oracle for real, checked-in files.
//!
//! Coverage per case:
//!   * 31 (coarse): obs points without names, single-point boundary
//!     time series, boundary enclosure + Neumann polyline;
//!   * 32: single-point JONSWAP boundary + enclosure;
//!   * 33: JONSWAP + enclosure whose file leads with a blank line;
//!   * 45 (haringvliet): obs points with quoted names + enclosure.
//!
//! The trailing-quote obs-name quirk of `43_st_croix` and the multi-point
//! time series are covered by unit tests in `src_rust/text_input.rs`; they
//! are not re-run here to keep the mesh-dependent oracle hook fast.

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

/// Run `--compare-text` on a copied run directory; returns the temp dir so
/// failures can be inspected (removed on success).
fn run_compare_text(case: &str, run_src: &str, inp_rel: &str) -> PathBuf {
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let src = manifest.join(run_src);
    assert!(src.is_dir(), "run directory not found: {}", src.display());

    let work = std::env::temp_dir().join(format!(
        "snapwave_text_{}_{}_{}",
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
        .arg("--compare-text")
        .arg(&inp)
        .output()
        .expect("failed to run wrapper binary");

    assert_eq!(
        out.status.code(),
        Some(0),
        "--compare-text failed for {case} ({})\n--- stdout ---\n{}\n--- stderr ---\n{}",
        inp.display(),
        stdout_of(&out),
        stderr_of(&out)
    );

    let _ = fs::remove_dir_all(&work);
    work
}

#[test]
fn rust_text_parsers_match_fortran_for_checked_in_cases() {
    run_compare_text("case31", "testcases/31_linear_shoaling_refraction/run", "coarse/SnapWave.inp");
    run_compare_text("case32", "testcases/32_curvi_island/run", "SnapWave.inp");
    run_compare_text("case33", "testcases/33_circle_reef/run", "SnapWave.inp");
    run_compare_text("case45", "testcases/45_haringvliet/run/snapwave", "snapwave.inp");
}

#[test]
fn compare_text_requires_its_input_file() {
    // Missing input path must be a wrapper usage error (exit 2), not a panic.
    let out = wrapper().arg("--compare-text").output().expect("run wrapper");
    assert_eq!(out.status.code(), Some(2), "missing input must exit 2");
    assert!(stderr_of(&out).contains("missing input"), "stderr was:\n{}", stderr_of(&out));
}
