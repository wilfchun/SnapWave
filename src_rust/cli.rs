//! Command-line parsing for the SnapWave wrapper (plan.md, Phase 2, step 1).
//!
//! Hand-rolled on purpose: the surface is one positional input path plus a
//! handful of flags, and the repo style keeps runtime dependencies minimal
//! (`clap` is only warranted once argument parsing genuinely grows).
//!
//! Status-code contract (kept consistent with the pre-Phase-2 wrapper and
//! AGENTS.md):
//!   * `0` — success (`--help`, `--version`, or a model run reporting 0);
//!   * [`EXIT_USAGE`] (2) — wrapper-detected errors: bad arguments, an
//!     unusable input path, failure to enter the run directory;
//!   * any other non-zero code — a Fortran model status, passed through
//!     unchanged by `main`.
//!
//! Parsing works on OS strings, not `String`: arguments that are not valid
//! UTF-8 must produce a clean usage error, never a panic (plan.md Phase 2,
//! step 4). Embedded NUL bytes cannot reach this code through `execve` on
//! the supported platforms; the FFI-side name conversion still guards
//! against them (see `run_context.rs`).

use std::ffi::{OsStr, OsString};
use std::path::PathBuf;

/// Exit code for wrapper-detected (usage/validation) errors.
pub const EXIT_USAGE: i32 = 2;

/// The input file the Fortran core expects (used in help/error text).
const INPUT_FILE_HINT: &str = "SnapWave.inp";

/// What the user asked the wrapper to do.
#[derive(Debug)]
pub enum Invocation {
    /// Print help text and exit 0.
    Help(String),
    /// Print version information and exit 0.
    Version(String),
    /// Run a model.
    Run(RunCommand),
    /// Parse the input in Rust, compare against the legacy Fortran reader
    /// through the temporary Phase 3 hook, and exit without running the
    /// model (plan.md Phase 3, step 5; removed again in Phase 4+).
    CompareInput(RunCommand),
}

/// A requested model run, still unvalidated: input validation happens when
/// the [`crate::run_context::RunContext`] is built.
#[derive(Debug)]
pub struct RunCommand {
    /// Input path exactly as given on the command line.
    pub input: PathBuf,
    /// Logging preference: print the run context before starting the model.
    pub verbose: bool,
}

/// Parse `args` (including `argv[0]`) into an [`Invocation`].
///
/// Rules: `-h`/`--help` and `-V`/`--version` short-circuit wherever they
/// appear; `-v`/`--verbose` toggles verbose logging; exactly one positional
/// argument (the input path) is required; `--` ends option parsing so file
/// names starting with `-` stay usable.
pub fn parse(args: &[OsString]) -> Result<Invocation, String> {
    let prog = program_name(args.first());
    let mut input: Option<OsString> = None;
    let mut verbose = false;
    let mut compare_input = false;
    let mut only_positional = false;

    for arg in args.iter().skip(1) {
        if !only_positional && arg.as_os_str() == OsStr::new("--") {
            only_positional = true;
            continue;
        }
        // After `--` everything is positional, even `-`-prefixed names.
        let is_flag = !only_positional && is_flag_arg(arg.as_os_str());
        match (is_flag, arg.as_os_str().to_str()) {
            (true, Some("-h")) | (true, Some("--help")) => {
                return Ok(Invocation::Help(help_text(&prog)));
            }
            (true, Some("-V")) | (true, Some("--version")) => {
                return Ok(Invocation::Version(version_text()));
            }
            (true, Some("-v")) | (true, Some("--verbose")) => verbose = true,
            (true, Some("--compare-input")) => compare_input = true,
            (true, Some(other)) => return Err(format!("unknown option: {other}")),
            (true, None) => {
                return Err(format!("unknown option: {}", arg.to_string_lossy()));
            }
            (false, _) => {
                if input.replace(arg.clone()).is_some() {
                    return Err(format!("unexpected extra argument: {}", arg.to_string_lossy()));
                }
            }
        }
    }

    let command = |input: OsString| RunCommand { input: PathBuf::from(input), verbose };
    match input {
        Some(input) if compare_input => Ok(Invocation::CompareInput(command(input))),
        Some(input) => Ok(Invocation::Run(command(input))),
        None => Err(format!("missing input path (expected one {INPUT_FILE_HINT} file)")),
    }
}

/// Program name for usage/error lines: the file name of `argv[0]`, falling
/// back to the canonical binary name.
pub fn program_name(argv0: Option<&OsString>) -> String {
    argv0
        .and_then(|a| PathBuf::from(a).file_name().map(|f| f.to_string_lossy().into_owned()))
        .unwrap_or_else(|| "snapwave".to_string())
}

/// A lone `-` is not a flag (it is a legal, if odd, file name); anything
/// else starting with `-` is treated as an option.
fn is_flag_arg(arg: &OsStr) -> bool {
    let bytes = arg.as_encoded_bytes();
    bytes.len() > 1 && bytes[0] == b'-'
}

pub fn version_text() -> String {
    format!("snapwave {} (Rust wrapper around the Fortran core)\n", env!("CARGO_PKG_VERSION"))
}

pub fn help_text(prog: &str) -> String {
    format!(
        r#"SnapWave — fast, implicit, unstructured-grid short wave solver
(Rust wrapper around the Fortran core)

Usage: {prog} [OPTIONS] <INPUT>

Arguments:
  <INPUT>  Path to a {INPUT_FILE_HINT} file (absolute, or relative to the
           current working directory)

Options:
  -v, --verbose  Print the run context before starting the model
  --compare-input  Parse INPUT in Rust and compare the result against the
                   legacy Fortran reader AND the resolved-config handoff,
                   then exit without running the model (plan.md Phase 3/4)
  -h, --help     Print this help and exit
  -V, --version  Print version information and exit

The model runs from the input file's directory: the input and output files
named in SnapWave.inp are resolved by the Fortran core relative to it
(legacy contract, plan.md Phase 2).

Exit status:
  0  success
  2  invalid arguments or unusable input path (wrapper-detected)
  N  non-zero Fortran model status, passed through unchanged
"#
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `args` as seen by `parse`: argv[0] plus the given arguments.
    fn args(list: &[&str]) -> Vec<OsString> {
        std::iter::once(OsString::from("snapwave"))
            .chain(list.iter().map(OsString::from))
            .collect()
    }

    #[test]
    fn help_short_and_long() {
        for flag in ["-h", "--help"] {
            match parse(&args(&[flag])) {
                Ok(Invocation::Help(text)) => {
                    assert!(text.contains("Usage:"), "help must show usage for {flag}");
                    assert!(text.contains("--version"));
                    assert!(text.contains("--verbose"));
                }
                other => panic!("{flag} should request help, got {other:?}"),
            }
        }
    }

    #[test]
    fn version_short_and_long() {
        for flag in ["-V", "--version"] {
            match parse(&args(&[flag])) {
                Ok(Invocation::Version(text)) => {
                    assert!(text.contains(env!("CARGO_PKG_VERSION")), "version text: {text}");
                    assert!(text.contains("snapwave"));
                }
                other => panic!("{flag} should request version, got {other:?}"),
            }
        }
    }

    #[test]
    fn missing_input_is_an_error() {
        let err = parse(&args(&[])).expect_err("no input must be rejected");
        assert!(err.contains("missing input"), "error was: {err}");
    }

    #[test]
    fn unknown_flag_is_an_error() {
        let err = parse(&args(&["--frobnicate"])).expect_err("unknown flags must be rejected");
        assert!(err.contains("unknown option"), "error was: {err}");
    }

    #[test]
    fn extra_positional_argument_is_an_error() {
        let err = parse(&args(&["a.inp", "b.inp"])).expect_err("two inputs must be rejected");
        assert!(err.contains("unexpected extra argument"), "error was: {err}");
    }

    #[test]
    fn verbose_flag_is_parsed() {
        match parse(&args(&["--verbose", "SnapWave.inp"])) {
            Ok(Invocation::Run(cmd)) => {
                assert!(cmd.verbose);
                assert_eq!(cmd.input, PathBuf::from("SnapWave.inp"));
            }
            other => panic!("expected a run, got {other:?}"),
        }
        match parse(&args(&["-v", "SnapWave.inp"])) {
            Ok(Invocation::Run(cmd)) => assert!(cmd.verbose),
            other => panic!("expected a run, got {other:?}"),
        }
    }

    #[test]
    fn compare_input_flag_is_parsed() {
        match parse(&args(&["--compare-input", "SnapWave.inp"])) {
            Ok(Invocation::CompareInput(cmd)) => {
                assert_eq!(cmd.input, PathBuf::from("SnapWave.inp"));
                assert!(!cmd.verbose);
            }
            other => panic!("expected a compare-input invocation, got {other:?}"),
        }
        // Still needs the positional input path.
        let err = parse(&args(&["--compare-input"])).expect_err("input path still required");
        assert!(err.contains("missing input"), "error was: {err}");
        // Documented in the help text.
        assert!(help_text("snapwave").contains("--compare-input"));
    }

    #[test]
    fn double_dash_ends_option_parsing() {
        match parse(&args(&["--", "--help"])) {
            Ok(Invocation::Run(cmd)) => assert_eq!(cmd.input, PathBuf::from("--help")),
            other => panic!("`--` must make `--help` a file name, got {other:?}"),
        }
    }

    #[test]
    fn help_after_a_valid_input_still_wins() {
        assert!(matches!(parse(&args(&["SnapWave.inp", "--help"])), Ok(Invocation::Help(_))));
    }

    #[test]
    fn bare_dash_is_treated_as_input() {
        match parse(&args(&["-"])) {
            Ok(Invocation::Run(cmd)) => assert_eq!(cmd.input, PathBuf::from("-")),
            other => panic!("`-` must be a positional, got {other:?}"),
        }
    }

    #[test]
    #[cfg(unix)]
    fn non_utf8_flag_is_rejected_gracefully() {
        use std::os::unix::ffi::OsStringExt;
        let mut argv = args(&[]);
        // Starts with `-`, so it is parsed as an option, but is not valid
        // UTF-8 and therefore cannot be a known one.
        argv.push(OsString::from_vec(b"-frobnicate\xff".to_vec()));
        let err = parse(&argv).expect_err("non-UTF-8 flags cannot be known options");
        assert!(err.contains("unknown option"), "error was: {err}");
    }

    #[test]
    fn program_name_uses_argv0_file_name() {
        assert_eq!(program_name(Some(&OsString::from("/a/b/snapwave.exe"))), "snapwave.exe");
        assert_eq!(program_name(Some(&OsString::from("snapwave"))), "snapwave");
        assert_eq!(program_name(None), "snapwave");
    }
}
