//! Phase-1 regression tests (plan.md, Phase 1: "Test Oracle And Baselines").
//!
//! For each case the harness:
//!   1. copies the testcase to a temp dir (checked-in inputs untouched,
//!      Windows separators normalized in the copy only);
//!   2. runs the Cargo-built Rust wrapper on the copy;
//!   3. compares the map/history NetCDF output against the **committed
//!      legacy baseline** under `testcases/<case>/output/` when one exists
//!      (lenient record counts, tolerance table in `support/compare.rs`);
//!   4. runs the **live Fortran oracle** (`make`-built binary, or
//!      `$SNAPWAVE_ORACLE`) on its own copy and compares wrapper vs oracle
//!      strictly — this covers families without committed baselines (e.g.
//!      the 31 map file) and double-checks everything else.
//!
//! Failing comparisons keep the temp directories and print their paths.
//! See `tests/README.md` for how to add cases.

mod support;

use std::path::{Path, PathBuf};

use support::compare::compare_files;
use support::harness::{assert_success, oracle_binary, run_command, CaseSpec, PreparedCase};
use support::ncdf::NcFile;

/// MWE regression target (AGENTS.md): linear shoaling/refraction, coarse
/// grid, 3 output frames (t = 0, 3600, 7200 s), 201 stations. A committed
/// history baseline exists; no committed map baseline does, so the map side
/// is pinned by the live oracle comparison.
const CASE_31_COARSE: CaseSpec = CaseSpec {
    name: "31_coarse",
    testcase_dir: "testcases/31_linear_shoaling_refraction",
    run_subdir: "run/coarse",
    map_file: "shoalref_coarse_neu_map.nc",
    his_file: "shoalref_coarse_neu_his.nc",
    expected_map_frames: Some(3),
    expected_his_frames: Some(3),
    use_map_baseline: true,
    use_his_baseline: true,
};

/// Broader validation case (plan.md Phase 1, step 4): curvilinear island
/// mesh, single-point JONSWAP boundary, enclosure file, per-iteration map
/// output (`ja_save_each_iter = 1`). The map frame count depends on solver
/// convergence, so it is not pinned; the history has a single frame
/// (tstart == tstop). The committed map baseline reflects current
/// behaviour; the committed his baseline does NOT — it was produced by an
/// older build that still wrote history output per solver iteration
/// (15 frames at t=1..15) — so it is disabled and the strict oracle
/// comparison pins the history family instead.
const CASE_32_CURVI: CaseSpec = CaseSpec {
    name: "32_curvi",
    testcase_dir: "testcases/32_curvi_island",
    run_subdir: "run",
    map_file: "snapwave_map.nc",
    his_file: "snapwave_his.nc",
    expected_map_frames: None,
    expected_his_frames: Some(1),
    use_map_baseline: true,
    use_his_baseline: false,
};

#[test]
fn mwe_31_linear_shoaling_refraction_coarse() {
    run_regression(&CASE_31_COARSE);
}

#[test]
fn broader_32_curvi_island() {
    run_regression(&CASE_32_CURVI);
}

struct Family<'a> {
    kind: &'static str,
    file: &'a str,
    expected_frames: Option<u64>,
}

fn families(spec: &CaseSpec) -> Vec<Family<'_>> {
    vec![
        Family { kind: "map", file: spec.map_file, expected_frames: spec.expected_map_frames },
        Family { kind: "his", file: spec.his_file, expected_frames: spec.expected_his_frames },
    ]
}

fn run_regression(spec: &CaseSpec) {
    let prepared = spec
        .prepare("rust")
        .unwrap_or_else(|e| panic!("preparing temp copy of {}: {e:#}", spec.name));
    let wrapper = PathBuf::from(env!("CARGO_BIN_EXE_snapwave"));

    let out = run_command(&wrapper, Some(prepared.inp.as_path()), &prepared.run_dir)
        .unwrap_or_else(|e| panic!("{}: spawning rust wrapper: {e:#}", spec.name));
    if let Err(e) = assert_success(&out, "rust wrapper") {
        panic!("{e}\n(output kept at {})", prepared.work.display());
    }

    // Wrapper output of every family, opened once and reused below.
    let mut wrapper_nc: Vec<(Family<'_>, NcFile)> = Vec::new();
    for fam in families(spec) {
        let produced = prepared.output_dir.join(fam.file);
        assert_nonempty(&produced, &prepared);
        let nc = NcFile::open(&produced)
            .unwrap_or_else(|e| panic!("{}: parsing wrapper output {}: {e:#}", spec.name, produced.display()));

        if let Some(expected) = fam.expected_frames {
            let t = nc
                .dim("time")
                .unwrap_or_else(|| panic!("{}: no time dimension in {}", spec.name, produced.display()));
            assert_eq!(
                t.len,
                expected,
                "{}: unexpected frame count in {} (run kept at {})",
                spec.name,
                produced.display(),
                prepared.work.display()
            );
        }
        wrapper_nc.push((fam, nc));
    }

    // --- committed baseline comparison ------------------------------------
    for (fam, act) in &wrapper_nc {
        let baseline = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join(spec.testcase_dir)
            .join("output")
            .join(fam.file);
        let use_baseline = match fam.kind {
            "map" => spec.use_map_baseline,
            _ => spec.use_his_baseline,
        };
        if baseline.is_file() && use_baseline {
            let base = NcFile::open(&baseline).unwrap_or_else(|e| {
                panic!("{}: parsing committed baseline {}: {e:#}", spec.name, baseline.display())
            });
            if let Err(report) =
                compare_files(&base, act, "committed baseline", "rust wrapper", true)
            {
                panic!(
                    "{} {} output does not match committed baseline:\n{report}\n(run kept at {})",
                    spec.name,
                    fam.kind,
                    prepared.work.display()
                );
            }
        } else if baseline.is_file() {
            // Reference exists but is known to predate current output
            // behaviour; comparing against it is meaningless. The strict
            // oracle comparison below pins this family instead.
            eprintln!(
                "note: {} {} baseline exists but is disabled (predates current output behaviour); oracle comparison pins this family",
                spec.name, fam.kind
            );
            assert!(
                act.var("time").is_some(),
                "{}: no time variable in {} output",
                spec.name,
                fam.kind
            );
        } else {
            // No committed reference: minimal structural pin only; the real
            // numeric pin for such families is the oracle comparison below.
            assert!(
                act.var("time").is_some(),
                "{}: no time variable in {} output",
                spec.name,
                fam.kind
            );
            eprintln!(
                "note: {} has no committed {} baseline; structural checks only (oracle comparison follows when available)",
                spec.name, fam.kind
            );
        }
    }

    // --- live Fortran oracle comparison ------------------------------------
    match oracle_binary() {
        Some(oracle) => {
            let oracle_prepared = spec
                .prepare("oracle")
                .unwrap_or_else(|e| panic!("{}: preparing oracle temp copy: {e:#}", spec.name));
            let out = run_command(&oracle, None, &oracle_prepared.run_dir)
                .unwrap_or_else(|e| panic!("{}: spawning fortran oracle {}: {e:#}", spec.name, oracle.display()));
            if let Err(e) = assert_success(&out, "fortran oracle") {
                panic!("{e}\n(oracle output kept at {})", oracle_prepared.work.display());
            }

            for (fam, act) in &wrapper_nc {
                let oracle_file = oracle_prepared.output_dir.join(fam.file);
                assert_nonempty(&oracle_file, &oracle_prepared);
                let onc = NcFile::open(&oracle_file).unwrap_or_else(|e| {
                    panic!("{}: parsing oracle output {}: {e:#}", spec.name, oracle_file.display())
                });
                if let Err(report) =
                    compare_files(&onc, act, "fortran oracle", "rust wrapper", false)
                {
                    panic!(
                        "{} {} output: rust wrapper diverges from fortran oracle:\n{report}\n(runs kept at {} and {})",
                        spec.name,
                        fam.kind,
                        prepared.work.display(),
                        oracle_prepared.work.display()
                    );
                }
            }
            oracle_prepared.cleanup();
        }
        None => eprintln!(
            "note: Fortran oracle not found for {} (run `make` or set SNAPWAVE_ORACLE); skipping wrapper-vs-oracle comparison",
            spec.name
        ),
    }

    prepared.cleanup();
}

fn assert_nonempty(file: &Path, prepared: &PreparedCase) {
    let meta = std::fs::metadata(file).unwrap_or_else(|e| {
        panic!("expected output {} missing: {e} (run kept at {})", file.display(), prepared.work.display())
    });
    assert!(
        meta.len() > 0,
        "output file {} is empty (run kept at {})",
        file.display(),
        prepared.work.display()
    );
}
