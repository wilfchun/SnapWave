//! Cargo-driven smoke test for the MWE (plan.md, Phase 5).
//!
//! Mirrors the reference flow used by the Nix flake `smoke-test` check:
//!   1. copy the coarse shoaling/refraction testcase to a temp directory;
//!   2. normalize Windows path separators in the copied `SnapWave.inp`;
//!   3. create the output directory the testcase writes to;
//!   4. run the Cargo-built Rust wrapper;
//!   5. verify the map/history NetCDF files exist;
//!   6. validate the NetCDF headers with `ncdump -h` (when available).
//!
//! Numeric comparisons against committed references are intentionally left
//! for a later phase; this test pins down structure and exit semantics.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

const TESTCASE_RUN_DIR: &str = "testcases/31_linear_shoaling_refraction/run";

#[test]
fn coarse_shoaling_refraction_runs_through_rust_wrapper() {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let src_run = manifest_dir.join(TESTCASE_RUN_DIR);
    assert!(src_run.is_dir(), "testcase run dir not found: {}", src_run.display());

    let work = std::env::temp_dir().join(format!("snapwave_mwe_{}", std::process::id()));
    if work.exists() {
        fs::remove_dir_all(&work).expect("clean stale work dir");
    }
    let dst_run = work.join("run");

    // 1. Copy the coarse run directory and the shared sibling input files
    //    that the (normalized) input references as ../<file>.
    copy_dir_recursive(&src_run.join("coarse"), &dst_run.join("coarse"));
    fs::create_dir_all(&dst_run).expect("create run dir");
    for entry in fs::read_dir(&src_run).expect("read run dir") {
        let entry = entry.expect("dir entry");
        if entry.file_type().expect("file type").is_file() {
            let to = dst_run.join(entry.file_name());
            fs::copy(entry.path(), &to).unwrap_or_else(|e| panic!("copy {}: {e}", entry.path().display()));
        }
    }

    // 2. Normalize Windows path separators (testcase is authored on Windows).
    let inp_path = dst_run.join("coarse").join("SnapWave.inp");
    let content = fs::read_to_string(&inp_path).expect("read SnapWave.inp");
    fs::write(&inp_path, content.replace('\\', "/")).expect("write normalized SnapWave.inp");

    // 3. The testcase writes to ../../output relative to the run directory.
    let output_dir = work.join("output");
    fs::create_dir_all(&output_dir).expect("create output dir");

    // 4. Run the Cargo-built wrapper binary.
    let bin = env!("CARGO_BIN_EXE_snapwave");
    let output = Command::new(bin)
        .arg(&inp_path)
        .output()
        .unwrap_or_else(|e| panic!("failed to run wrapper binary {bin}: {e}"));

    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        output.status.success(),
        "wrapper failed with status {:?}\n--- stdout ---\n{}\n--- stderr ---\n{}",
        output.status.code(),
        stdout,
        stderr
    );

    // 5. Verify generated files exist and are non-empty.
    let his = output_dir.join("shoalref_coarse_neu_his.nc");
    let map = output_dir.join("shoalref_coarse_neu_map.nc");
    for file in [&his, &map] {
        let len = fs::metadata(file)
            .unwrap_or_else(|e| panic!("expected output {} missing: {e}", file.display()))
            .len();
        assert!(len > 0, "output file {} is empty", file.display());
    }

    // 6. Schema checks against the captured baseline headers (Phase 0).
    if Command::new("ncdump").arg("-h").output().is_err() {
        eprintln!("note: ncdump not available, skipping NetCDF header validation");
    } else {
        check_ncdump_header(&his, &[
            "time = UNLIMITED",
            "stations = 201",
            "float time(time)",
            "float point_hm0(time, stations)",
            "float point_tp(time, stations)",
            "float point_wavdir(time, stations)",
        ]);
        check_ncdump_header(&map, &[
            "nmesh2d_node = ",
            "ntheta = ",
            "time = UNLIMITED",
            "float hm0(time, nmesh2d_node)",
            "float tp(time, nmesh2d_node)",
            "float ee(time, nmesh2d_node, ntheta)",
        ]);
    }

    // Keep artifacts when the assertions above fail; clean up on success.
    fs::remove_dir_all(&work).expect("clean up work dir");
}

#[test]
fn wrapper_rejects_missing_input_file() {
    let bin = env!("CARGO_BIN_EXE_snapwave");
    let output = Command::new(bin)
        .arg("/nonexistent/path/to/SnapWave.inp")
        .output()
        .expect("failed to run wrapper binary");
    assert!(
        !output.status.success(),
        "wrapper should fail on a missing input file, got status {:?}",
        output.status.code()
    );
    let stderr = String::from_utf8_lossy(&output.stderr);
    assert!(
        stderr.contains("SnapWave.inp"),
        "error message should mention the missing input path, got: {stderr}"
    );
}

fn check_ncdump_header(file: &Path, expected: &[&str]) {
    let output = Command::new("ncdump")
        .arg("-h")
        .arg(file)
        .output()
        .unwrap_or_else(|e| panic!("failed to run ncdump on {}: {e}", file.display()));
    assert!(
        output.status.success(),
        "ncdump -h failed on {}",
        file.display()
    );
    let header = String::from_utf8_lossy(&output.stdout);
    for needle in expected {
        assert!(
            header.contains(needle),
            "NetCDF header of {} is missing expected fragment {:?}\n--- header ---\n{}",
            file.display(),
            needle,
            header
        );
    }
}

fn copy_dir_recursive(src: &Path, dst: &Path) {
    fs::create_dir_all(dst).unwrap_or_else(|e| panic!("create {}: {e}", dst.display()));
    for entry in fs::read_dir(src).unwrap_or_else(|e| panic!("read {}: {e}", src.display())) {
        let entry = entry.expect("dir entry");
        let to = dst.join(entry.file_name());
        if entry.file_type().expect("file type").is_dir() {
            copy_dir_recursive(&entry.path(), &to);
        } else {
            fs::copy(entry.path(), &to)
                .unwrap_or_else(|e| panic!("copy {} -> {}: {e}", entry.path().display(), to.display()));
        }
    }
}
