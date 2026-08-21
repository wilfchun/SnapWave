//! Process-level run context owned by Rust (plan.md, Phase 2, step 2).
//!
//! `RunContext` captures everything the wrapper knows about *where* and
//! *how* a model runs before the model core takes over: the input paths,
//! the run directory, executable metadata and logging preferences. Path
//! resolution against the run directory and the output-directory policy
//! live in `crate::paths` since plan.md Phase 5; this struct remains the
//! process-level anchor (input path, run directory). Since Phase 12 the
//! whole model runs in Rust, so there is no `chdir` contract and no FFI
//! name conversion any more.

use std::ffi::OsString;
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};

use crate::paths::OutputPaths;

/// Executable metadata used for `--version` output and diagnostics.
pub struct ExeMeta {
    /// Program name as invoked (file name of `argv[0]`).
    pub name: String,
    /// Wrapper version (Cargo package version of this crate).
    pub version: &'static str,
}

impl ExeMeta {
    pub fn from_argv0(argv0: Option<&OsString>) -> Self {
        ExeMeta { name: crate::cli::program_name(argv0), version: env!("CARGO_PKG_VERSION") }
    }
}

/// Wrapper-side logging preferences.
pub struct LogPrefs {
    /// Print the run context before starting the model (`--verbose`).
    pub verbose: bool,
}

/// Validated, process-level description of one model run.
pub struct RunContext {
    /// Input path exactly as provided on the command line (used in
    /// user-facing messages).
    pub input_arg: PathBuf,
    /// Canonicalized absolute path of the validated input file.
    pub input_path: PathBuf,
    /// Directory of the input file; the model runs here. Every file
    /// reference of `SnapWave.inp` resolves against it (`crate::paths`).
    pub run_dir: PathBuf,
    /// Executable metadata (name, version).
    pub exe: ExeMeta,
    /// Logging preferences for the wrapper.
    pub log: LogPrefs,
}

impl RunContext {
    /// Validate the requested input path and resolve the run context.
    ///
    /// All wrapper-side user-facing input validation happens here, so the
    /// model core only ever sees an existing regular file (plan.md
    /// Phase 2 acceptance: Rust owns command-line validation and
    /// user-facing wrapper errors).
    pub fn new(input_arg: PathBuf, exe: ExeMeta, log: LogPrefs) -> Result<Self> {
        let input_path = input_arg
            .canonicalize()
            .with_context(|| format!("input file not found or not accessible: {}", input_arg.display()))?;
        if !input_path.is_file() {
            bail!("input path is not a file: {}", input_arg.display());
        }
        let run_dir = input_path
            .parent()
            .map(Path::to_path_buf)
            .with_context(|| format!("input path has no parent directory: {}", input_path.display()))?;

        Ok(RunContext { input_arg, input_path, run_dir, exe, log })
    }

    /// Human-readable summary printed for `--verbose`. Called after the
    /// input has been parsed so it can include the resolved output paths.
    pub fn describe(&self, outputs: &OutputPaths) -> String {
        let out = |p: &Option<PathBuf>| match p {
            Some(p) => p.display().to_string(),
            None => "(disabled)".to_string(),
        };
        [
            format!("snapwave {} run context", self.exe.version),
            format!("  input (as given) : {}", self.input_arg.display()),
            format!("  input (resolved) : {}", self.input_path.display()),
            format!("  run directory    : {}", self.run_dir.display()),
            format!("  map output       : {}", out(&outputs.map)),
            format!("  his output       : {}", out(&outputs.his)),
            "  model core       : Rust (plan.md Phase 12)".to_string(),
        ]
        .join("\n")
    }
}

