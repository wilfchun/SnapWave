//! Process-level run context owned by Rust (plan.md, Phase 2, step 2).
//!
//! `RunContext` captures everything the wrapper knows about *where* and
//! *how* a model runs before the Fortran core takes over: the input paths,
//! the run directory, executable metadata and logging preferences. Later
//! phases extend or replace parts of it as input parsing (Phase 3) and
//! output-directory policy (Phase 5) migrate to Rust.

use std::ffi::{CString, OsStr, OsString};
use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};

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
    /// relative-path contract, see [`RunContext::enter_run_dir`]).
    pub run_dir: PathBuf,
    /// Expected output directory, if known. `None` for now: output file
    /// names and locations come from `SnapWave.inp` and are resolved by the
    /// Fortran core relative to the run directory. Rust takes ownership of
    /// output-directory policy in plan.md Phase 5.
    pub output_dir: Option<PathBuf>,
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

        Ok(RunContext { input_arg, input_path, run_dir, output_dir: None, exe, log })
    }

    /// Human-readable summary printed for `--verbose`.
    pub fn describe(&self) -> String {
        let output_dir = match &self.output_dir {
            Some(dir) => dir.display().to_string(),
            None => "defined by SnapWave.inp (Fortran-owned until plan.md Phase 5)".to_string(),
        };
        [
            format!("snapwave {} run context", self.exe.version),
            format!("  input (as given) : {}", self.input_arg.display()),
            format!("  input (resolved) : {}", self.input_path.display()),
            format!("  run directory    : {}", self.run_dir.display()),
            format!("  output directory : {output_dir}"),
            "  model core       : Fortran via the C ABI facade".to_string(),
        ]
        .join("\n")
    }

    /// Legacy working-directory contract (plan.md, Phase 2, step 3).
    ///
    /// The Fortran readers resolve all sibling input and output paths
    /// relative to the current working directory, so until those readers
    /// and the output handling migrate to Rust (Phases 5 and 6) the
    /// wrapper must `chdir` into the input file's directory before calling
    /// the facade. (The input file itself has been Rust-selected and
    /// passed down explicitly since Phase 3; everything else is still
    /// CWD-relative on the Fortran side.) The contract is isolated in this
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
