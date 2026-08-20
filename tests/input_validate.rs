//! Phase 4 validation tests (plan.md Phase 4, step 4).
//!
//! Covers:
//!   - bad output intervals (negative/zero with output enabled → exit 2;
//!     accepted without output file);
//!   - optional output settings (disabled outputs with bad intervals
//!     accepted; intervals default to timestep);
//!   - missing input file (wrapper exit 2, already covered by mwe/cli
//!     tests — re-asserted here for the "validation before domain init"
//!     acceptance);
//!   - the resolved-config handoff produces a working run (implicitly
//!     covered by the regression tests, which now exercise the Phase 4
//!     config path).

mod support;

use std::fs;
use std::path::PathBuf;
use std::process::{Command, Output};

fn wrapper() -> Command {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_snapwave"));
    cmd.env("OMP_NUM_THREADS", "1");
    cmd
}

fn stderr_of(out: &Output) -> String {
    String::from_utf8_lossy(&out.stderr).into_owned()
}

fn stdout_of(out: &Output) -> String {
    String::from_utf8_lossy(&out.stdout).into_owned()
}

/// Every bad-interval combination that Fortran rejects (and Rust must
/// reject before the Fortran core runs).
#[test]
fn bad_output_intervals_fail_as_wrapper_errors() {
    let cases: &[(&str, &str, &str)] = &[
        ("negative map_interval with map_file", "map_file = out.nc\nmap_interval = -1\n", "map_interval"),
        ("zero map_interval with map_file", "map_file = out.nc\nmap_interval = 0\n", "map_interval"),
        ("negative his_interval with his_file", "his_file = out.nc\nhis_interval = -0.5\n", "his_interval"),
        ("zero his_interval with his_file", "his_file = out.nc\nhis_interval = 0\n", "his_interval"),
        ("both intervals bad", "map_file = m.nc\nhis_file = h.nc\nmap_interval = -1\nhis_interval = 0\n", "map_interval"),
    ];

    let dir = std::env::temp_dir().join(format!("snapwave_phase4_intervals_{}", std::process::id()));
    fs::create_dir_all(&dir).expect("create temp dir");

    for (label, content, needle) in cases {
        let inp = dir.join("SnapWave.inp");
        fs::write(&inp, content).expect("write input");
        let out = wrapper().arg(&inp).output().expect("failed to run wrapper");

        let stderr = stderr_of(&out);
        assert_eq!(
            out.status.code(),
            Some(2),
            "{label}: expected wrapper exit code 2\n--- stdout ---\n{}\n--- stderr ---\n{stderr}",
            stdout_of(&out)
        );
        assert!(stderr.contains("error"), "{label}: stderr was:\n{stderr}");
        assert!(stderr.contains(needle), "{label}: error must name '{needle}':\n{stderr}");
        assert!(
            !stdout_of(&out).contains("Reading input file"),
            "{label}: the Fortran reader must not have been invoked"
        );
    }

    let _ = fs::remove_dir_all(&dir);
}

/// Non-positive intervals are accepted when the corresponding output file
/// is empty (Fortran behaviour preserved by the Rust parser).
#[test]
fn nonpositive_intervals_accepted_without_output_file() {
    let dir = std::env::temp_dir().join(format!("snapwave_phase4_optional_{}", std::process::id()));
    fs::create_dir_all(&dir).expect("create temp dir");

    let cases: &[&str] = &[
        "map_interval = -1\n",
        "his_interval = 0\n",
        "map_interval = -1\nhis_interval = 0\n",
    ];

    for (i, content) in cases.iter().enumerate() {
        let inp = dir.join(format!("SnapWave_{i}.inp"));
        fs::write(&inp, content).expect("write input");
        let out = wrapper().arg("--compare-input").arg(&inp).output().expect("failed to run wrapper");
        assert_eq!(
            out.status.code(),
            Some(0),
            "case {i}: non-positive interval without output file must be accepted\n\
             --- stdout ---\n{}\n--- stderr ---\n{}",
            stdout_of(&out),
            stderr_of(&out)
        );
    }

    let _ = fs::remove_dir_all(&dir);
}

/// The resolved-config handoff (Phase 4) is verified by `--compare-input`
/// for every checked-in testcase input. This test re-runs the comparison
/// on the MWE input to confirm the new hook works end-to-end.
#[test]
fn resolved_config_handoff_matches_for_mwe() {
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let inp = manifest.join("testcases/31_linear_shoaling_refraction/run/coarse/SnapWave.inp");
    assert!(inp.is_file(), "MWE input not found at {}", inp.display());

    let out = wrapper().arg("--compare-input").arg(&inp).output().expect("failed to run wrapper");
    assert_eq!(
        out.status.code(),
        Some(0),
        "--compare-input (Phase 3 + Phase 4) failed for MWE\n\
         --- stdout ---\n{}\n--- stderr ---\n{}",
        stdout_of(&out),
        stderr_of(&out)
    );
}

/// The wrapper rejects a missing input file before any Fortran code runs
/// (RunContext validation, exit 2). This is the "missing required files"
/// case from step 4.
#[test]
fn missing_input_file_fails_before_fortran_runs() {
    let out = wrapper()
        .arg("/nonexistent/path/to/SnapWave.inp")
        .output()
        .expect("failed to run wrapper");
    assert_eq!(out.status.code(), Some(2), "missing input must exit 2");
    let stderr = stderr_of(&out);
    assert!(stderr.contains("SnapWave.inp"), "error must mention the input path: {stderr}");
    assert!(
        !stdout_of(&out).contains("Reading input file"),
        "Fortran must not have been invoked"
    );
}