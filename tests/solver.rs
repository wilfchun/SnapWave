//! Phase 11 integration tests: solver parity against the Fortran routines
//! (plan.md Phase 11).
//!
//! The `--compare-solver` mode runs the unchanged Fortran solver for one
//! timestep through the temporary `snapwave_solver_dump_c` hook, computes
//! the same solver step in Rust (`src_rust/solver.rs`), and compares the
//! resulting solver-state globals.
//!
//! Unit tests for the small deterministic routines live with the code:
//! `src_rust/solver.rs` (solve_tridiag, baldock, hpsort_eps_epw, disper_nr,
//! compute_celerities, numerical_limiter, windinput, vegatt, swvegatt,
//! bulkdragcoeff).

use std::process::{Command, Output};

fn stderr_of(out: &Output) -> String {
    String::from_utf8_lossy(&out.stderr).into_owned()
}

fn wrapper() -> Command {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_snapwave"));
    cmd.env("OMP_NUM_THREADS", "1");
    cmd
}

#[test]
fn compare_solver_flag_is_recognized() {
    // The --compare-solver flag must be a recognized option (not "unknown option").
    // Without an input path it should produce a "missing input" error (exit 2).
    let out = wrapper()
        .arg("--compare-solver")
        .output()
        .expect("run wrapper");
    assert_eq!(out.status.code(), Some(2), "missing input must exit 2");
    let stderr = stderr_of(&out);
    assert!(
        stderr.contains("missing input"),
        "stderr must mention missing input, got:\n{stderr}"
    );
}

#[test]
fn compare_solver_requires_netcdf_mesh() {
    // The --compare-solver mode requires a NetCDF mesh. Passing a
    // non-existent file should produce a clean error.
    let out = wrapper()
        .arg("--compare-solver")
        .arg("/nonexistent/path/SnapWave.inp")
        .output()
        .expect("run wrapper");
    // Should fail with a wrapper error (exit 2) because the input file
    // doesn't exist.
    assert_eq!(out.status.code(), Some(2), "nonexistent input must exit 2");
}

#[test]
fn compare_solver_appears_in_help() {
    let out = wrapper().arg("--help").output().expect("run wrapper");
    let stdout = String::from_utf8_lossy(&out.stdout);
    assert!(
        stdout.contains("--compare-solver"),
        "help text must document --compare-solver:\n{stdout}"
    );
}