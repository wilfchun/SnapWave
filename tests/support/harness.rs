//! Testcase preparation and process running for the Phase-1 harness
//! (plan.md Phase 1).
//!
//! Contract:
//! - Testcases are **copied** to a temp directory before anything runs, so
//!   checked-in inputs are never modified. Windows path separators are
//!   normalized (`\` → `/`) only in the copied `*.inp` files.
//! - The committed reference NetCDF files under `testcases/<case>/output/`
//!   are never copied into the run tree; they stay where they are and act as
//!   read-only baselines.
//! - The Rust wrapper is invoked with the absolute input path (it chdirs to
//!   the input directory itself, see `src_rust/main.rs`); the legacy Fortran
//!   oracle takes no arguments and is run with its CWD in the copied run
//!   directory, mirroring `runsnapwave.bat`.
//! - Both runs pin `OMP_NUM_THREADS=1`: OpenMP reduction order changes
//!   floating-point summation order, which would make numeric comparisons
//!   flaky for no scientific reason.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};

use anyhow::{bail, Context, Result};

/// Directories inside a testcase that never influence a model run.
/// `output` is excluded on purpose so committed baselines are neither copied
/// nor overwritten.
const SKIP_DIRS: [&str; 5] = [".svn", "scripts", "results", "analytical", "output"];

/// Input file names probed by `read_snapwave_input` (src/snapwave_input.f90).
const INPUT_NAMES: [&str; 6] =
    ["snapwave.inp", "SnapWave.inp", "SNAPWAVE.INP", "snapwave.INP", "Snapwave.INP", "SNAPWAVE.inp"];

pub struct CaseSpec {
    /// Short id used in temp dir names and log lines.
    pub name: &'static str,
    /// Testcase root relative to the crate manifest,
    /// e.g. `"testcases/31_linear_shoaling_refraction"`.
    pub testcase_dir: &'static str,
    /// Directory holding `SnapWave.inp`, relative to the testcase root; this
    /// becomes the process CWD (the Fortran reader resolves sibling paths
    /// relative to it).
    pub run_subdir: &'static str,
    /// Map output file name the testcase writes into `<root>/output/`.
    pub map_file: &'static str,
    /// History output file name the testcase writes into `<root>/output/`.
    pub his_file: &'static str,
    /// Pin the wrapper's map output frame count when known from the input
    /// times (tstop-tstart)/timestep; `None` leaves it unchecked (e.g. when
    /// the count depends on solver convergence, as with
    /// `ja_save_each_iter = 1`).
    pub expected_map_frames: Option<u64>,
    /// Pin the wrapper's history output frame count, same rules as above.
    pub expected_his_frames: Option<u64>,
    /// Compare map output against the committed baseline under
    /// `testcases/<case>/output/` when it exists. Disable only when the
    /// reference predates current output behaviour (documented per case).
    pub use_map_baseline: bool,
    /// Compare history output against the committed baseline; see
    /// `use_map_baseline`.
    pub use_his_baseline: bool,
}

pub struct PreparedCase {
    /// Temp work dir holding the copied testcase (removed by `cleanup`).
    pub work: PathBuf,
    /// Copied testcase root inside `work`.
    pub root: PathBuf,
    /// Copied run directory (process CWD for both binaries).
    pub run_dir: PathBuf,
    /// Output directory inside the copy (`<root>/output`).
    pub output_dir: PathBuf,
    /// Absolute path of the (normalized) `SnapWave.inp` in the copy.
    pub inp: PathBuf,
}

impl CaseSpec {
    /// Copy the testcase and prepare run/output directories under a fresh
    /// temp dir tagged `tag` (used to separate wrapper and oracle runs).
    pub fn prepare(&self, tag: &str) -> Result<PreparedCase> {
        let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        let src = manifest.join(self.testcase_dir);
        if !src.is_dir() {
            bail!("testcase directory not found: {}", src.display());
        }
        let case_root_name = match src.file_name() {
            Some(n) => n.to_os_string(),
            None => bail!("cannot derive testcase dir name from {}", src.display()),
        };
        let unique = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.as_micros())
            .unwrap_or(0);
        let work = std::env::temp_dir()
            .join(format!("snapwave_{}_{}_{tag}_{unique}", self.name, std::process::id()));
        if work.exists() {
            fs::remove_dir_all(&work)
                .with_context(|| format!("cleaning stale work dir {}", work.display()))?;
        }
        let root = work.join(&case_root_name);
        copy_case_inputs(&src, &root)?;
        // plan.md Phase 1, step 5: separator normalization happens on the
        // temp copies only.
        normalize_input_separators(&root)?;

        let run_dir = root.join(self.run_subdir);
        if !run_dir.is_dir() {
            bail!("run directory missing after copy: {}", run_dir.display());
        }
        let inp = find_input_file(&run_dir)?;
        let output_dir = root.join("output");
        fs::create_dir_all(&output_dir)
            .with_context(|| format!("creating output dir {}", output_dir.display()))?;

        Ok(PreparedCase { work, root, run_dir, output_dir, inp })
    }
}

impl PreparedCase {
    /// Remove the temp tree. Kept alive (with a notice) when
    /// `SNAPWAVE_TEST_KEEP=1` for debugging; failing tests intentionally do
    /// not call this.
    pub fn cleanup(&self) {
        if std::env::var("SNAPWAVE_TEST_KEEP").is_ok() {
            eprintln!("note: keeping temp run at {}", self.work.display());
            return;
        }
        let _ = fs::remove_dir_all(&self.work);
    }
}

/// Location of the legacy Fortran oracle: `$SNAPWAVE_ORACLE` if set, else the
/// default `make` target (`SnapWave/lnx64/bin/snapwave`). `None` disables
/// wrapper-vs-oracle comparison (committed baselines still apply).
pub fn oracle_binary() -> Option<PathBuf> {
    if let Ok(path) = std::env::var("SNAPWAVE_ORACLE") {
        if !path.is_empty() {
            return Some(PathBuf::from(path));
        }
    }
    let default =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("SnapWave").join("lnx64").join("bin").join("snapwave");
    if default.is_file() {
        Some(default)
    } else {
        None
    }
}

/// Run a SnapWave binary deterministically (single OpenMP thread) with the
/// given CWD; `arg` is the input path for the Rust wrapper, `None` for the
/// argument-less Fortran oracle.
pub fn run_command(program: &Path, arg: Option<&Path>, cwd: &Path) -> Result<Output> {
    let mut cmd = Command::new(program);
    cmd.current_dir(cwd).env("OMP_NUM_THREADS", "1");
    if let Some(a) = arg {
        cmd.arg(a);
    }
    cmd.output().with_context(|| format!("failed to spawn {}", program.display()))
}

pub fn assert_success(out: &Output, what: &str) -> Result<()> {
    if !out.status.success() {
        bail!(
            "{what} failed with status {:?}\n--- stdout ---\n{}\n--- stderr ---\n{}",
            out.status.code(),
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr)
        );
    }
    Ok(())
}

fn find_input_file(run_dir: &Path) -> Result<PathBuf> {
    for name in INPUT_NAMES {
        let p = run_dir.join(name);
        if p.is_file() {
            return Ok(p);
        }
    }
    bail!("no SnapWave.inp found in {}", run_dir.display());
}

fn copy_case_inputs(src: &Path, dst: &Path) -> Result<()> {
    fs::create_dir_all(dst).with_context(|| format!("creating {}", dst.display()))?;
    for entry in fs::read_dir(src).with_context(|| format!("reading {}", src.display()))? {
        let entry = entry.with_context(|| format!("reading dir entry in {}", src.display()))?;
        let file_type = entry.file_type()?;
        let name = entry.file_name();
        if file_type.is_dir() {
            if SKIP_DIRS.contains(&name.to_string_lossy().as_ref()) {
                continue;
            }
            copy_case_inputs(&entry.path(), &dst.join(&name))?;
        } else {
            fs::copy(entry.path(), dst.join(&name)).with_context(|| {
                format!("copying {} -> {}", entry.path().display(), dst.join(&name).display())
            })?;
        }
    }
    Ok(())
}

/// Replace `\` with `/` in every `*.inp` under the copied tree. Testcases are
/// authored on Windows; the Fortran readers on Linux only accept `/`.
fn normalize_input_separators(dir: &Path) -> Result<()> {
    for entry in fs::read_dir(dir).with_context(|| format!("reading {}", dir.display()))? {
        let entry = entry.with_context(|| format!("reading dir entry in {}", dir.display()))?;
        let path = entry.path();
        if entry.file_type()?.is_dir() {
            normalize_input_separators(&path)?;
        } else if entry.file_name().to_string_lossy().to_ascii_lowercase().ends_with(".inp") {
            let text = fs::read_to_string(&path)
                .with_context(|| format!("reading {}", path.display()))?;
            let normalized = text.replace('\\', "/");
            if normalized != text {
                fs::write(&path, normalized)
                    .with_context(|| format!("writing normalized {}", path.display()))?;
            }
        }
    }
    Ok(())
}
