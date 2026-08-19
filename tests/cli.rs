//! Phase-2 CLI integration tests (plan.md, Phase 2, step 4).
//!
//! Run the built wrapper binary as a subprocess and pin the wrapper-owned
//! process behaviour: `--help`/`--version` flags and exit codes, usage
//! errors, graceful handling of unusable input paths (missing file,
//! directory, invalid UTF-8) and relative input paths on a real testcase.
//!
//! Embedded NUL bytes cannot be passed through `Command` (execve rejects
//! them), so that path is covered by the unit tests on the FFI name
//! conversion in `src_rust/run_context.rs` instead.

mod support;

use std::fs;
use std::process::{Command, Output};

use support::harness::CaseSpec;

/// Same coarse testcase as the MWE regression target (canonical spec in
/// `regression.rs`; duplicated here so the Phase-1 harness stays untouched).
const CASE_31_COARSE: CaseSpec = CaseSpec {
    name: "31_coarse_cli",
    testcase_dir: "testcases/31_linear_shoaling_refraction",
    run_subdir: "run/coarse",
    map_file: "shoalref_coarse_neu_map.nc",
    his_file: "shoalref_coarse_neu_his.nc",
    expected_map_frames: None,
    expected_his_frames: None,
    use_map_baseline: false,
    use_his_baseline: false,
};

fn run_wrapper(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_snapwave"))
        .args(args)
        .output()
        .unwrap_or_else(|e| panic!("failed to run wrapper binary: {e}"))
}

fn stderr_of(out: &Output) -> String {
    String::from_utf8_lossy(&out.stderr).into_owned()
}

fn stdout_of(out: &Output) -> String {
    String::from_utf8_lossy(&out.stdout).into_owned()
}

fn assert_exit(code: i32, out: &Output, what: &str) {
    assert_eq!(
        out.status.code(),
        Some(code),
        "{what}: expected exit code {code}, got {:?}\n--- stdout ---\n{}\n--- stderr ---\n{}",
        out.status.code(),
        stdout_of(out),
        stderr_of(out)
    );
}

#[test]
fn help_flag_prints_usage_and_succeeds() {
    for flag in ["-h", "--help"] {
        let out = run_wrapper(&[flag]);
        assert_exit(0, &out, flag);
        let stdout = stdout_of(&out);
        assert!(stdout.contains("Usage:"), "{flag} output must show usage:\n{stdout}");
        assert!(stdout.contains("--verbose"), "{flag} output must document --verbose:\n{stdout}");
        assert!(stdout.contains("--version"), "{flag} output must document --version:\n{stdout}");
        assert!(stderr_of(&out).is_empty(), "{flag} must not write to stderr");
    }
}

#[test]
fn version_flag_prints_version_and_succeeds() {
    for flag in ["-V", "--version"] {
        let out = run_wrapper(&[flag]);
        assert_exit(0, &out, flag);
        let stdout = stdout_of(&out);
        assert!(
            stdout.contains(env!("CARGO_PKG_VERSION")),
            "{flag} output must contain the crate version:\n{stdout}"
        );
        assert!(stderr_of(&out).is_empty(), "{flag} must not write to stderr");
    }
}

#[test]
fn no_arguments_is_a_usage_error() {
    let out = run_wrapper(&[]);
    assert_exit(2, &out, "missing argument");
    let stderr = stderr_of(&out);
    assert!(stderr.contains("missing input"), "stderr was:\n{stderr}");
    assert!(stderr.contains("--help"), "stderr should point at --help:\n{stderr}");
}

#[test]
fn unknown_flag_is_a_usage_error() {
    let out = run_wrapper(&["--frobnicate"]);
    assert_exit(2, &out, "unknown flag");
    assert!(stderr_of(&out).contains("unknown option"), "stderr was:\n{}", stderr_of(&out));
}

#[test]
fn extra_positional_argument_is_a_usage_error() {
    let out = run_wrapper(&["a.inp", "b.inp"]);
    assert_exit(2, &out, "two positional arguments");
    assert!(
        stderr_of(&out).contains("unexpected extra argument"),
        "stderr was:\n{}",
        stderr_of(&out)
    );
}

#[test]
fn missing_input_file_is_a_wrapper_error() {
    let out = run_wrapper(&["/nonexistent/path/to/SnapWave.inp"]);
    assert_exit(2, &out, "missing input file");
    let stderr = stderr_of(&out);
    assert!(
        stderr.contains("/nonexistent/path/to/SnapWave.inp"),
        "error must name the missing path, stderr was:\n{stderr}"
    );
    assert!(!stderr.contains("panicked"), "wrapper must fail cleanly, stderr was:\n{stderr}");
}

#[test]
fn directory_input_is_a_wrapper_error() {
    let dir = std::env::temp_dir().join(format!("snapwave_cli_dir_{}", std::process::id()));
    fs::create_dir_all(&dir).expect("create temp dir");
    let out = run_wrapper(&[dir.to_str().expect("temp dir path is UTF-8")]);
    assert_exit(2, &out, "directory as input");
    let stderr = stderr_of(&out);
    assert!(stderr.contains("not a file"), "stderr was:\n{stderr}");
    assert!(!stderr.contains("panicked"), "wrapper must fail cleanly, stderr was:\n{stderr}");
    fs::remove_dir_all(&dir).ok();
}

#[test]
#[cfg(unix)]
fn invalid_utf8_input_argument_fails_gracefully() {
    use std::ffi::OsString;
    use std::os::unix::ffi::OsStringExt;
    // Non-UTF-8 arguments must be reported, not panic on (std::env::args
    // would have panicked; args_os + lossy display must not).
    let bad = OsString::from_vec(b"/nonexistent/\xffSnapWave.inp".to_vec());
    let out = Command::new(env!("CARGO_BIN_EXE_snapwave"))
        .arg(&bad)
        .output()
        .unwrap_or_else(|e| panic!("failed to run wrapper binary: {e}"));
    assert_exit(2, &out, "invalid UTF-8 argument");
    let stderr = stderr_of(&out);
    assert!(!stderr.contains("panicked"), "wrapper must not panic, stderr was:\n{stderr}");
    assert!(
        stderr.contains("SnapWave.inp"),
        "error should name the (lossily displayed) path, stderr was:\n{stderr}"
    );
}

/// plan.md Phase 2, step 4: a *relative* input path must resolve against
/// the process working directory (not the binary's location), and the full
/// coarse testcase must run through it. Also pins `--verbose` output.
#[test]
fn relative_input_path_runs_coarse_testcase() {
    let prepared = CASE_31_COARSE
        .prepare("cli")
        .unwrap_or_else(|e| panic!("preparing temp copy of the coarse testcase: {e:#}"));

    let rel = prepared
        .inp
        .strip_prefix(&prepared.work)
        .expect("input file must live under the temp work dir");
    assert!(rel.is_relative(), "stripped path must be relative: {}", rel.display());

    // Run from the temp work dir (not the run dir) so the relative path
    // resolution is actually exercised. Single OpenMP thread, matching the
    // regression harness convention.
    let out = Command::new(env!("CARGO_BIN_EXE_snapwave"))
        .current_dir(&prepared.work)
        .env("OMP_NUM_THREADS", "1")
        .arg("--verbose")
        .arg(rel)
        .output()
        .unwrap_or_else(|e| panic!("failed to run wrapper binary: {e}"));

    if !out.status.success() {
        panic!(
            "wrapper failed with status {:?}\n--- stdout ---\n{}\n--- stderr ---\n{}\n(run kept at {})",
            out.status.code(),
            stdout_of(&out),
            stderr_of(&out),
            prepared.work.display()
        );
    }

    // --verbose must have printed the run context (to stderr).
    let stderr = stderr_of(&out);
    assert!(stderr.contains("run directory"), "--verbose stderr was:\n{stderr}");
    assert!(stderr.contains("input (resolved)"), "--verbose stderr was:\n{stderr}");

    // Outputs land in <testcase root>/output relative to the run dir.
    for file in [prepared.output_dir.join(CASE_31_COARSE.map_file), prepared.output_dir.join(CASE_31_COARSE.his_file)] {
        let len = fs::metadata(&file)
            .unwrap_or_else(|e| panic!("expected output {} missing: {e} (run kept at {})", file.display(), prepared.work.display()))
            .len();
        assert!(len > 0, "output file {} is empty (run kept at {})", file.display(), prepared.work.display());
    }

    prepared.cleanup();
}
