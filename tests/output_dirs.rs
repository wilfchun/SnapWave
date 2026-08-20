//! Phase 5 integration tests: filesystem and output-directory handling
//! (plan.md, Phase 5, step 5).
//!
//! All tests copy the coarse shoaling/refraction testcase (the MWE
//! regression target) to temp directories, normalize the Windows `\`
//! separators **only in the copies**, and never pre-create output
//! directories — the wrapper owns output-directory policy since this
//! phase (`src_rust/paths.rs`), so the tests pin:
//!
//! * outputs of a deeply nested run directory resolve relative to the
//!   input file's directory, and the missing `output/` parent is
//!   created by the wrapper;
//! * output directories configured through `map_file`/`his_file` are
//!   created when missing;
//! * an empty `map_file` disables map output (run succeeds, history
//!   file written, no map file);
//! * an empty `his_file` disables history output (run succeeds, map
//!   file written, no history file).

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

const TESTCASE_RUN_DIR: &str = "testcases/31_linear_shoaling_refraction/run";
const MAP_FILE: &str = "shoalref_coarse_neu_map.nc";
const HIS_FILE: &str = "shoalref_coarse_neu_his.nc";

/// Copy the coarse run directory and its shared sibling input files
/// into `run_root` (the `run/` directory of the copy), normalizing `\`
/// → `/` in the copied `SnapWave.inp` (testcases are authored on
/// Windows; AGENTS.md: normalize only temp copies). No output
/// directories are created — that is the wrapper's job now.
fn prepare_case(work: &Path) -> PathBuf {
    let src_run = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(TESTCASE_RUN_DIR);
    assert!(src_run.is_dir(), "testcase run dir not found: {}", src_run.display());

    let run_root = work.join("run");
    copy_dir_recursive(&src_run.join("coarse"), &run_root.join("coarse"));
    for entry in fs::read_dir(&src_run).expect("read run dir") {
        let entry = entry.expect("dir entry");
        if entry.file_type().expect("file type").is_file() {
            fs::copy(entry.path(), run_root.join(entry.file_name()))
                .unwrap_or_else(|e| panic!("copy {}: {e}", entry.path().display()));
        }
    }

    let inp = run_root.join("coarse").join("SnapWave.inp");
    let content = fs::read_to_string(&inp).expect("read SnapWave.inp");
    assert!(content.contains('\\'), "testcase must still use Windows separators");
    fs::write(&inp, content.replace('\\', "/")).expect("write normalized SnapWave.inp");
    inp
}

/// Replace (or empty) one `key = value` record of an input file,
/// preserving the first-occurrence-wins grammar (the replacement lands
/// on the original line position).
fn set_keyword(text: &str, key: &str, value: &str) -> String {
    let mut out = String::new();
    let mut replaced = false;
    for line in text.lines() {
        if !replaced
            && line.starts_with(key)
            && line[key.len()..].trim_start_matches(' ').starts_with('=')
        {
            out.push_str(&format!("{key} = {value}\n"));
            replaced = true;
        } else {
            out.push_str(line);
            out.push('\n');
        }
    }
    assert!(replaced, "keyword '{key}' not found in input");
    out
}

fn run_wrapper(args: &[&str], cwd: &Path) -> Output {
    Command::new(env!("CARGO_BIN_EXE_snapwave"))
        .args(args)
        .current_dir(cwd)
        .env("OMP_NUM_THREADS", "1")
        .output()
        .unwrap_or_else(|e| panic!("failed to run wrapper binary: {e}"))
}

fn assert_success(out: &Output, what: &str) {
    assert!(
        out.status.success(),
        "{what} failed with status {:?}\n--- stdout ---\n{}\n--- stderr ---\n{}",
        out.status.code(),
        String::from_utf8_lossy(&out.stdout),
        String::from_utf8_lossy(&out.stderr)
    );
}

fn fresh_work(name: &str) -> PathBuf {
    let work = std::env::temp_dir().join(format!("snapwave_phase5_{name}_{}", std::process::id()));
    if work.exists() {
        fs::remove_dir_all(&work).expect("clean stale work dir");
    }
    fs::create_dir_all(&work).expect("create work dir");
    work
}

/// plan.md Phase 5, step 5: a deeply nested run directory keeps output
/// resolution anchored at the input file's directory, and the wrapper
/// creates the missing `../../output` parent itself (the copy contains
/// no pre-created output directory).
#[test]
fn nested_run_directory_resolves_and_creates_output() {
    let work = fresh_work("nested");
    // Deep nesting below the process CWD; the wrapper is invoked with a
    // relative input path from `work`.
    let inp = prepare_case(&work.join("deep").join("ly").join("nested"));
    assert!(inp.is_file(), "input must exist at {}", inp.display());

    let rel = inp.strip_prefix(&work).unwrap().to_string_lossy().into_owned();
    assert!(Path::new(&rel).is_relative(), "input path must be relative: {rel}");
    let out = run_wrapper(&[&rel], &work);
    assert_success(&out, "nested run");

    // ../../output relative to run/coarse lands at deep/ly/nested/output
    // — and the wrapper must have created it (it did not exist).
    let output_dir = work.join("deep").join("ly").join("nested").join("output");
    for file in [MAP_FILE, HIS_FILE] {
        let path = output_dir.join(file);
        let len = fs::metadata(&path)
            .unwrap_or_else(|e| panic!("expected output {} missing: {e}", path.display()))
            .len();
        assert!(len > 0, "output file {} is empty", path.display());
    }

    fs::remove_dir_all(&work).expect("clean up work dir");
}

/// plan.md Phase 5, step 3/5: output directories named by `map_file` /
/// `his_file` that do not exist are created by the wrapper before the
/// Fortran core runs.
#[test]
fn missing_output_parents_are_created() {
    let work = fresh_work("mkdir");
    let inp = prepare_case(&work);

    let text = fs::read_to_string(&inp).expect("read copied SnapWave.inp");
    let text = set_keyword(&text, "map_file", "../made/up/m.nc");
    let text = set_keyword(&text, "his_file", "../made/up/h.nc");
    fs::write(&inp, text).expect("write modified SnapWave.inp");

    let made = work.join("run").join("made").join("up");
    assert!(!made.exists(), "precondition: output parent must not exist yet");

    let out = run_wrapper(&[inp.to_string_lossy().as_ref()], &work);
    assert_success(&out, "run with missing output parents");

    for file in ["m.nc", "h.nc"] {
        let path = made.join(file);
        let len = fs::metadata(&path)
            .unwrap_or_else(|e| panic!("expected output {} missing: {e}", path.display()))
            .len();
        assert!(len > 0, "output file {} is empty", path.display());
    }

    fs::remove_dir_all(&work).expect("clean up work dir");
}

/// plan.md Phase 5, step 5: an empty `map_file` disables map output
/// (legacy rule: only the empty string disables). The run must succeed
/// and write the history file only. Also pins the `--verbose` run
/// context, which reports the resolved (disabled) output paths.
#[test]
fn disabled_map_output_still_writes_history() {
    let work = fresh_work("nomap");
    let inp = prepare_case(&work);

    let text = fs::read_to_string(&inp).expect("read copied SnapWave.inp");
    let text = set_keyword(&text, "map_file", "");
    fs::write(&inp, text).expect("write modified SnapWave.inp");

    let out = run_wrapper(&["--verbose", inp.to_string_lossy().as_ref()], &work);
    assert_success(&out, "run with map output disabled");

    // Verbose run context names both output families (Phase 5).
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(stderr.contains("map output"), "--verbose stderr was:\n{stderr}");
    assert!(stderr.contains("(disabled)"), "--verbose stderr was:\n{stderr}");

    let output_dir = work.join("output");
    assert!(output_dir.join(HIS_FILE).is_file(), "history output must exist");
    assert!(!output_dir.join(MAP_FILE).exists(), "map output must not exist");
    // No stray map file next to the input either.
    assert!(!inp.parent().unwrap().join(MAP_FILE).exists());

    fs::remove_dir_all(&work).expect("clean up work dir");
}

/// plan.md Phase 5, step 5: an empty `his_file` disables history
/// output; the run succeeds and writes the map file only.
#[test]
fn disabled_history_output_still_writes_map() {
    let work = fresh_work("nohis");
    let inp = prepare_case(&work);

    let text = fs::read_to_string(&inp).expect("read copied SnapWave.inp");
    let text = set_keyword(&text, "his_file", "");
    fs::write(&inp, text).expect("write modified SnapWave.inp");

    let out = run_wrapper(&[inp.to_string_lossy().as_ref()], &work);
    assert_success(&out, "run with history output disabled");

    let output_dir = work.join("output");
    assert!(output_dir.join(MAP_FILE).is_file(), "map output must exist");
    assert!(!output_dir.join(HIS_FILE).exists(), "history output must not exist");

    fs::remove_dir_all(&work).expect("clean up work dir");
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
