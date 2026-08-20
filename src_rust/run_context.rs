//! Process-level run context owned by Rust (plan.md, Phase 2, step 2).
//!
//! `RunContext` captures everything the wrapper knows about *where* and
//! *how* a model runs before the Fortran core takes over: the input paths,
//! the run directory, executable metadata and logging preferences. Path
//! resolution against the run directory and the output-directory policy
//! live in `crate::paths` since plan.md Phase 5; this struct remains the
//! process-level anchor (input path, run directory, legacy `chdir`
//! contract).

use std::ffi::{CString, OsStr, OsString};
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

/// Wrapper-side logging preferences. The Fortran core's own console output
/// is unchanged; these only control extra Rust-side diagnostics.
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
    /// Directory of the input file; the model runs here (legacy Fortran
    /// relative-path contract, see [`RunContext::enter_run_dir`]). Every
    /// file reference of `SnapWave.inp` resolves against it
    /// (`crate::paths`, plan.md Phase 5).
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
    /// Fortran facade only ever sees an existing regular file (plan.md
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
    /// input has been parsed so it can include the resolved output paths
    /// (plan.md Phase 5).
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
            "  model core       : Fortran via the C ABI facade".to_string(),
        ]
        .join("\n")
    }

    /// Legacy working-directory contract (plan.md, Phase 2, step 3).
    ///
    /// The Fortran readers resolve all sibling input and output paths
    /// relative to the current working directory, so until those readers
    /// and the NetCDF IO migrate to Rust (Phases 6 and 7) the wrapper
    /// must `chdir` into the input file's directory before calling the
    /// facade. (The input file itself has been Rust-selected and passed
    /// down explicitly since Phase 3, and the output *directories* are
    /// Rust-owned since Phase 5 — but the file *names* still cross to
    /// Fortran as CWD-relative text.) The contract is isolated in this
    /// single method so later phases can remove it cleanly.
    pub fn enter_run_dir(&self) -> Result<()> {
        std::env::set_current_dir(&self.run_dir)
            .with_context(|| format!("failed to change working directory to {}", self.run_dir.display()))
    }

    /// The input file name as a C string for the FFI facade
    /// (`snapwave_run_c`, see `src/snapwave_c_api.f90`; explicit length, no
    /// NUL termination).
    ///
    /// On Unix the raw OS bytes are used (file names need not be valid
    /// UTF-8); other platforms require valid Unicode for the conversion.
    /// Embedded NUL bytes cannot occur in real file names on the supported
    /// platforms, but the boundary conversion must still reject them
    /// instead of panicking (plan.md Phase 2, step 4).
    pub fn input_file_name_cstring(&self) -> Result<CString> {
        let name = self
            .input_path
            .file_name()
            .with_context(|| format!("cannot derive file name from {}", self.input_path.display()))?;
        os_str_to_cstring(name)
    }
}

/// Convert an arbitrary path to a C string for the FFI facade (same
/// platform rules as [`RunContext::input_file_name_cstring`]; used by the
/// Phase 3 comparison hook, which crosses the facade with full paths).
pub fn path_to_cstring(path: &Path) -> Result<CString> {
    os_str_to_cstring(path.as_os_str())
}

#[cfg(unix)]
fn os_str_to_cstring(name: &OsStr) -> Result<CString> {
    use std::os::unix::ffi::OsStrExt;
    CString::new(name.as_bytes()).with_context(|| {
        format!("input file name contains a NUL byte and cannot cross the FFI boundary: {}", name.to_string_lossy())
    })
}

#[cfg(not(unix))]
fn os_str_to_cstring(name: &OsStr) -> Result<CString> {
    let text = name.to_str().with_context(|| {
        format!("input file name is not valid Unicode and cannot cross the FFI boundary: {}", name.to_string_lossy())
    })?;
    CString::new(text).with_context(|| format!("input file name contains a NUL byte and cannot cross the FFI boundary: {text}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cstring_name_rejects_embedded_nul() {
        // A NUL byte is representable in an OsStr (just not in a real path
        // or C string); the boundary conversion must reject it cleanly.
        let err = os_str_to_cstring(OsStr::new("bad\0name")).expect_err("NUL must be rejected");
        assert!(format!("{err:#}").contains("NUL"), "error was: {err:#}");
    }

    #[test]
    #[cfg(unix)]
    fn cstring_name_accepts_non_utf8_bytes() {
        use std::os::unix::ffi::OsStrExt;
        let name = OsStr::from_bytes(b"SnapWave\xff.inp");
        let c = os_str_to_cstring(name).expect("non-UTF-8 names must cross FFI on unix");
        assert_eq!(c.as_bytes(), b"SnapWave\xff.inp");
    }
}
