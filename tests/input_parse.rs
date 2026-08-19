//! Phase 3 integration tests: `SnapWave.inp` parsing in Rust (plan.md
//! Phase 3).
//!
//! 1. `--compare-input` runs *both* parsers — the Rust parser and the
//!    legacy Fortran reader via the temporary facade hook
//!    (`snapwave_read_input_dump_c`) — on every checked-in testcase input
//!    and asserts the results match. This pins the Phase 3 acceptance
//!    criteria "Rust can parse every checked-in SnapWave.inp" and
//!    "Rust parse results match Fortran globals".
//! 2. Quirky-but-legal inputs (duplicate keys, missing `=`, indented keys,
//!    unknown keywords) are cross-checked the same way, proving the Rust
//!    parser preserves the legacy quirks rather than merely the common
//!    case.
//! 3. Invalid inputs fail as *wrapper* errors (exit 2) with the Fortran
//!    core never being invoked (Phase 3 acceptance: "Wrapper failures for
//!    invalid input are Rust errors, not Fortran stop").
//!
//! The checked-in testcase files are used read-only and directly (no temp
//! copies): the parse comparison resolves no paths but the input file
//! itself, so Windows-authored `\` separators in values are irrelevant
//! here.

mod support;

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use support::harness::INPUT_NAMES;

fn wrapper() -> Command {
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_snapwave"));
    // Fortran module init may spin an OpenMP team; keep it deterministic.
    cmd.env("OMP_NUM_THREADS", "1");
    cmd
}

fn stdout_of(out: &Output) -> String {
    String::from_utf8_lossy(&out.stdout).into_owned()
}

fn stderr_of(out: &Output) -> String {
    String::from_utf8_lossy(&out.stderr).into_owned()
}

/// Every checked-in testcase input (the six file names the legacy reader
/// probes), anywhere under `testcases/`.
fn find_checked_in_inputs() -> Vec<PathBuf> {
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("testcases");
    assert!(root.is_dir(), "testcases directory not found at {}", root.display());
    let mut found = Vec::new();
    walk(&root, &mut found);
    assert!(!found.is_empty(), "no testcase inputs found under {}", root.display());
    found.sort();
    found
}

fn walk(dir: &Path, out: &mut Vec<PathBuf>) {
    for entry in fs::read_dir(dir).unwrap_or_else(|e| panic!("reading {}: {e}", dir.display())) {
        let entry = entry.expect("dir entry");
        let path = entry.path();
        if path.is_dir() {
            walk(&path, out);
        } else if INPUT_NAMES.contains(&entry.file_name().to_string_lossy().as_ref()) {
            out.push(path);
        }
    }
}

fn run_compare(inp: &Path) -> Output {
    wrapper().arg("--compare-input").arg(inp).output().expect("failed to run wrapper binary")
}

#[test]
fn rust_parse_matches_fortran_globals_for_all_checked_in_inputs() {
    let files = find_checked_in_inputs();
    // 12 SnapWave.inp + 1 snapwave.inp at the time of writing; the exact
    // count is not the point, coverage of every checked-in file is.
    assert!(
        files.len() >= 13,
        "expected at least 13 checked-in inputs (12 SnapWave.inp + snapwave.inp), found {}: {files:#?}",
        files.len()
    );
    for inp in &files {
        let out = run_compare(inp);
        assert_eq!(
            out.status.code(),
            Some(0),
            "--compare-input failed for {}\n--- stdout ---\n{}\n--- stderr ---\n{}",
            inp.display(),
            stdout_of(&out),
            stderr_of(&out)
        );
    }
}

#[test]
fn parse_quirks_agree_with_the_fortran_reader() {
    // A hand-crafted file exercising the legacy quirks that checked-in
    // inputs only show partially. Both parsers must agree on all of it.
    let dir = std::env::temp_dir().join(format!("snapwave_phase3_quirks_{}", std::process::id()));
    fs::create_dir_all(&dir).expect("create temp dir");
    let inp = dir.join("SnapWave.inp");
    fs::write(
        &inp,
        concat!(
            // First occurrence wins for duplicate keys.
            "hmin = 0.42\n",
            "hmin = 9.9\n",
            // Lines without '=' are ignored (used to "comment out" values).
            "map_file         ../output/ignored-because-there-is-no-equals-sign.nc\n",
            "gridfile  also_ignored.txt\n",
            // Keys are case-sensitive; leading blanks break matching.
            "TREF = 19990101 000000\n",
            "   tstart = 20240102 000000\n",
            // Unknown keywords are ignored.
            "wind = 0\n",
            "his_ee = 1\n",
            "some free-form prose without an equals sign\n",
            // The real grammar underneath.
            "tref = 20240417 000000\n",
            "tstart=20240418 010203\n",
            "tstop = 20240419 020304\n",
            "timestep = 61.5,\n",
            "niter = 33/44\n",
            "crit = .02\n",
            "encfile =   enclosure file with inner blanks.txt   \n",
            "u10 = 0.00\n", // not the literal '0.0': wind growth turns on
            "map_Cg = 1\n",
            "map_cg = 1\n", // unknown keyword (exact spelling matters)
        ),
    )
    .expect("write quirk input");

    let out = run_compare(&inp);
    assert_eq!(
        out.status.code(),
        Some(0),
        "--compare-input failed for the quirk input\n--- stdout ---\n{}\n--- stderr ---\n{}",
        stdout_of(&out),
        stderr_of(&out)
    );

    // And the non-positive interval is only an error with output enabled
    // (Fortran behaviour preserved by the Rust parser).
    let inp2 = dir.join("snapwave.inp");
    fs::write(&inp2, "map_interval = -1\n").expect("write interval input");
    let out = run_compare(&inp2);
    assert_eq!(
        out.status.code(),
        Some(0),
        "non-positive map_interval without map_file must be accepted by both parsers\n\
         --- stdout ---\n{}\n--- stderr ---\n{}",
        stdout_of(&out),
        stderr_of(&out)
    );

    let _ = fs::remove_dir_all(&dir);
}

#[test]
fn invalid_input_fails_as_a_wrapper_error_before_fortran_runs() {
    // Each case is a minimal SnapWave.inp that must fail in the *Rust*
    // parser (exit 2 with a wrapper error) before the Fortran core is
    // ever invoked — formerly Fortran `stop 1`/runtime-abort paths.
    let cases: &[(&str, &str, &str)] = &[
        ("unparseable real", "timestep = banana\n", "timestep"),
        ("decimal in integer", "niter = 3.5\n", "niter"),
        ("missing numeric value", "timestep =\n", "timestep"),
        ("integer overflow", "niter = 9999999999\n", "niter"),
        ("unparseable date", "tref = notadate\n", "tstart"),
        ("negative map interval with map output", "map_file = out.nc\nmap_interval = -1\n", "map_interval"),
        ("zero his interval with his output", "his_file = out.nc\nhis_interval = 0\n", "his_interval"),
    ];

    let dir = std::env::temp_dir().join(format!("snapwave_phase3_invalid_{}", std::process::id()));
    fs::create_dir_all(&dir).expect("create temp dir");

    for (label, content, needle) in cases {
        let inp = dir.join("SnapWave.inp");
        fs::write(&inp, content).expect("write invalid input");
        let out = wrapper().arg(&inp).output().expect("failed to run wrapper binary");

        let stderr = stderr_of(&out);
        assert_eq!(
            out.status.code(),
            Some(2),
            "{label}: expected wrapper exit code 2\n--- stdout ---\n{}\n--- stderr ---\n{stderr}",
            stdout_of(&out)
        );
        assert!(stderr.contains("error"), "{label}: stderr was:\n{stderr}");
        assert!(stderr.contains(needle), "{label}: error must name the offending keyword '{needle}':\n{stderr}");
        // The Fortran core never starts: its reader announces itself on
        // stdout before anything else.
        assert!(
            !stdout_of(&out).contains("Reading input file"),
            "{label}: the Fortran reader must not have been invoked, stdout was:\n{}",
            stdout_of(&out)
        );
    }

    let _ = fs::remove_dir_all(&dir);
}
