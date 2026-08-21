//! Phase 8 route-parity tests (plan.md Phase 8) and Phase 10 time-loop
//! parity tests (plan.md Phase 10).
//!
//! The wrapper's default run path hands the Rust-owned domain state
//! (mesh, polylines, observation points, boundary series — see
//! `src_rust/state.rs`) to Fortran through `snapwave_run_capture_state_c`
//! instead of letting the Fortran readers re-read the files. The
//! `--legacy-mesh` flag forces the old Fortran-reading route, so running
//! one testcase both ways and comparing the outputs pins the handoff
//! against the unchanged readers *without* requiring the `make`-built
//! oracle (the regression suite adds that comparison when the oracle
//! exists).
//!
//! Phase 10 moves the time loop and output scheduling to Rust. The
//! `--fortran-time-loop` flag forces the old Fortran-owned time loop
//! (`snapwave_run_capture_c` / `snapwave_run_capture_state_c`), so
//! running one testcase both ways and comparing the outputs pins the
//! Rust-owned loop against the Fortran-owned loop.
//!
//! Unit tests for the layout/indexing conversions live with the code:
//! `src_rust/ffi_layout.rs` (one-based indices, column-major offsets)
//! and `src_rust/state.rs` (FFI widths, buffer layouts, deg→rad
//! recipes, ModelState scheduling logic) — plan.md Phase 8 acceptance:
//! "array shape and indexing conversions are covered by focused tests".

mod support;

use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use anyhow::{Context, Result};
use support::compare::compare_files;
use support::harness::{assert_success, CaseSpec, PreparedCase};
use support::ncdf::NcFile;

/// MWE regression target (same case as `tests/regression.rs`): NetCDF
/// mesh, boundary time series, enclosure + Neumann polylines,
/// observation points — everything the Phase 8 handoff covers.
const CASE_31_COARSE: CaseSpec = CaseSpec {
    name: "31_coarse",
    testcase_dir: "testcases/31_linear_shoaling_refraction",
    run_subdir: "run/coarse",
    map_file: "shoalref_coarse_neu_map.nc",
    his_file: "shoalref_coarse_neu_his.nc",
    expected_map_frames: Some(3),
    expected_his_frames: Some(3),
    use_map_baseline: false,
    use_his_baseline: false,
};

/// Single-point JONSWAP boundary variant (boundary_mode = 1 in the FFI
/// state) on a curvilinear mesh with per-iteration map output.
const CASE_32_CURVI: CaseSpec = CaseSpec {
    name: "32_curvi",
    testcase_dir: "testcases/32_curvi_island",
    run_subdir: "run",
    map_file: "snapwave_map.nc",
    his_file: "snapwave_his.nc",
    expected_map_frames: None,
    expected_his_frames: Some(1),
    use_map_baseline: false,
    use_his_baseline: false,
};

/// Run the Cargo-built wrapper deterministically (single OpenMP thread,
/// like the regression harness) with extra CLI flags before the input.
fn run_wrapper(wrapper: &Path, flags: &[&str], input: &Path, cwd: &Path) -> Result<Output> {
    let mut cmd = Command::new(wrapper);
    cmd.current_dir(cwd).env("OMP_NUM_THREADS", "1").args(flags).arg(input);
    cmd.output().with_context(|| format!("failed to spawn {}", wrapper.display()))
}

fn assert_output(prepared: &PreparedCase, file: &str) -> PathBuf {
    let path = prepared.output_dir.join(file);
    let meta = match std::fs::metadata(&path) {
        Ok(m) => m,
        Err(e) => panic!("expected output {} missing: {e} (run kept at {})", path.display(), prepared.work.display()),
    };
    assert!(meta.len() > 0, "output file {} is empty (run kept at {})", path.display(), prepared.work.display());
    path
}

/// Run `spec` through both input routes and require identical outputs
/// (strict comparison: same record counts, Phase-1 tolerances).
fn assert_route_parity(spec: &CaseSpec) {
    let wrapper = PathBuf::from(env!("CARGO_BIN_EXE_snapwave"));

    // --- default route: Rust-owned state handed to Fortran -------------
    let state = spec
        .prepare("state")
        .unwrap_or_else(|e| panic!("{}: preparing state-route copy: {e:#}", spec.name));
    let out = run_wrapper(&wrapper, &[], &state.inp, &state.run_dir)
        .unwrap_or_else(|e| panic!("{}: spawning wrapper (state route): {e:#}", spec.name));
    if let Err(e) = assert_success(&out, "wrapper (state route)") {
        panic!("{e}\n(state-route run kept at {})", state.work.display());
    }

    // --- legacy route: Fortran reads the mesh/text files itself --------
    let legacy = spec
        .prepare("legacymesh")
        .unwrap_or_else(|e| panic!("{}: preparing legacy-route copy: {e:#}", spec.name));
    let out = run_wrapper(&wrapper, &["--legacy-mesh"], &legacy.inp, &legacy.run_dir)
        .unwrap_or_else(|e| panic!("{}: spawning wrapper (legacy route): {e:#}", spec.name));
    if let Err(e) = assert_success(&out, "wrapper (legacy route)") {
        panic!("{e}\n(legacy-route run kept at {})", legacy.work.display());
    }

    for file in [spec.map_file, spec.his_file] {
        let state_file = assert_output(&state, file);
        let legacy_file = assert_output(&legacy, file);
        let state_nc = NcFile::open(&state_file)
            .with_context(|| format!("parsing state-route output {}", state_file.display()))
            .unwrap();
        let legacy_nc = NcFile::open(&legacy_file)
            .with_context(|| format!("parsing legacy-route output {}", legacy_file.display()))
            .unwrap();
        if let Err(report) = compare_files(
            &legacy_nc,
            &state_nc,
            "fortran-read route",
            "rust-state route",
            false,
        ) {
            panic!(
                "{}: {} output differs between the Fortran-read and Rust-state routes:\n{report}\n(runs kept at {} and {})",
                spec.name,
                file,
                state.work.display(),
                legacy.work.display()
            );
        }
    }

    state.cleanup();
    legacy.cleanup();
}

#[test]
fn state_route_matches_fortran_readers_31_coarse() {
    assert_route_parity(&CASE_31_COARSE);
}

#[test]
fn state_route_matches_fortran_readers_32_curvi_singlepoint() {
    assert_route_parity(&CASE_32_CURVI);
}

// ---------------------------------------------------------------------------
// Phase 10: Rust-owned time loop vs Fortran-owned time loop parity
// ---------------------------------------------------------------------------

/// Run `spec` through both time-loop routes and require identical outputs.
/// The default route uses the Phase 10 Rust-owned time loop; the
/// `--fortran-time-loop` route uses the old Fortran-owned time loop
/// (`snapwave_run_capture_state_c`). Both routes use the Rust-state
/// handoff for mesh/obs/boundary data (no `--legacy-mesh`).
fn assert_time_loop_parity(spec: &CaseSpec) {
    let wrapper = PathBuf::from(env!("CARGO_BIN_EXE_snapwave"));

    // --- Rust-owned time loop (default Phase 10 path) ------------------
    let rust_loop = spec
        .prepare("rustloop")
        .unwrap_or_else(|e| panic!("{}: preparing rust-loop copy: {e:#}", spec.name));
    let out = run_wrapper(&wrapper, &[], &rust_loop.inp, &rust_loop.run_dir)
        .unwrap_or_else(|e| panic!("{}: spawning wrapper (rust loop): {e:#}", spec.name));
    if let Err(e) = assert_success(&out, "wrapper (rust loop)") {
        panic!("{e}\n(rust-loop run kept at {})", rust_loop.work.display());
    }

    // --- Fortran-owned time loop (--fortran-time-loop) -----------------
    let f77_loop = spec
        .prepare("f77loop")
        .unwrap_or_else(|e| panic!("{}: preparing fortran-loop copy: {e:#}", spec.name));
    let out = run_wrapper(&wrapper, &["--fortran-time-loop"], &f77_loop.inp, &f77_loop.run_dir)
        .unwrap_or_else(|e| panic!("{}: spawning wrapper (fortran loop): {e:#}", spec.name));
    if let Err(e) = assert_success(&out, "wrapper (fortran loop)") {
        panic!("{e}\n(fortran-loop run kept at {})", f77_loop.work.display());
    }

    for file in [spec.map_file, spec.his_file] {
        let rust_file = assert_output(&rust_loop, file);
        let f77_file = assert_output(&f77_loop, file);
        let rust_nc = NcFile::open(&rust_file)
            .with_context(|| format!("parsing rust-loop output {}", rust_file.display()))
            .unwrap();
        let f77_nc = NcFile::open(&f77_file)
            .with_context(|| format!("parsing fortran-loop output {}", f77_file.display()))
            .unwrap();
        if let Err(report) = compare_files(
            &f77_nc,
            &rust_nc,
            "fortran time loop",
            "rust time loop",
            false,
        ) {
            panic!(
                "{}: {} output differs between the Fortran-owned and Rust-owned time loops:\n{report}\n(runs kept at {} and {})",
                spec.name,
                file,
                rust_loop.work.display(),
                f77_loop.work.display()
            );
        }
    }

    rust_loop.cleanup();
    f77_loop.cleanup();
}

#[test]
fn rust_time_loop_matches_fortran_time_loop_31_coarse() {
    assert_time_loop_parity(&CASE_31_COARSE);
}

#[test]
fn rust_time_loop_matches_fortran_time_loop_32_curvi_singlepoint() {
    assert_time_loop_parity(&CASE_32_CURVI);
}
