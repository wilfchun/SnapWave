//! Thin Rust wrapper around the SnapWave Fortran core.
//!
//! Usage: snapwave [OPTIONS] <path/to/SnapWave.inp>
//!
//! Phase 2 (plan.md): Rust owns all process-level behaviour — argument
//! parsing (`cli`), input validation and the run context (`run_context`),
//! status-code semantics — while Fortran still performs all model work
//! through the coarse C ABI facade in `src/snapwave_c_api.f90`.
//!
//! Phase 3 (plan.md): the wrapper parses and validates `SnapWave.inp` in
//! Rust (`input`) before the Fortran core runs, so invalid input is a
//! wrapper error rather than a Fortran `stop`; the `--compare-input` mode
//! cross-checks the Rust parse against the legacy Fortran reader through
//! a temporary facade hook (`input_compare`).
//!
//! Phase 4 (plan.md): configuration defaults, validation and diagnostics
//! are fully Rust-owned. The wrapper resolves the entire configuration
//! (defaults, post-processing) and passes it to the Fortran facade as
//! canonical key=value text — Fortran no longer reads SnapWave.inp or
//! decides defaults on this route. The `--compare-input` mode now also
//! verifies the resolved-config handoff round-trips through Fortran.
//!
//! Phase 5 (plan.md): filesystem and output-directory policy is
//! Rust-owned. Every file reference of `SnapWave.inp` resolves against
//! the input file's directory in one place (`paths`), the output files
//! are `PathBuf`s, and the required output directories are created or
//! validated in Rust before the Fortran core runs. The legacy `chdir`
//! stays (isolated in `RunContext::enter_run_dir`) until the Phase 6-7
//! readers accept explicit paths.
//!
//! Status codes: 0 on success (including `--help`/`--version`), 2 on
//! wrapper-detected errors, and any non-zero Fortran status is passed
//! through unchanged.

mod cli;
mod diagnostics;
mod input;
mod input_compare;
mod paths;
mod run_context;

use anyhow::Context;
use std::ffi::{c_char, c_int, OsString};
use std::io::Write;

use anyhow::{bail, Result};

use cli::{Invocation, EXIT_USAGE};
use run_context::{path_to_cstring, ExeMeta, LogPrefs, RunContext};

// FFI boundary (AGENTS.md): signatures must match src/snapwave_c_api.f90
// exactly (`bind(C)`; explicit length, no reliance on NUL termination).
extern "C" {
    // plan.md Phase 4: the facade now receives the fully-resolved
    // configuration as canonical key=value text (Rust is the authority
    // for defaults, validation and post-processing).
    fn snapwave_run_c(config: *const c_char, config_len: c_int) -> c_int;

    // Phase 3 comparison hook: parse the input file with the legacy
    // Fortran reader and dump the resulting globals.
    fn snapwave_read_input_dump_c(
        input_path: *const c_char,
        input_path_len: c_int,
        dump_path: *const c_char,
        dump_path_len: c_int,
    ) -> c_int;

    // Phase 4 comparison hook: load the Rust-resolved configuration text
    // and dump the resulting globals, pinning the config handoff.
    fn snapwave_load_config_dump_c(
        config: *const c_char,
        config_len: c_int,
        dump_path: *const c_char,
        dump_path_len: c_int,
    ) -> c_int;
}

fn main() {
    // args_os (not args): invalid-UTF-8 arguments must produce a clean
    // usage error, not a panic (plan.md Phase 2, step 4).
    let argv: Vec<OsString> = std::env::args_os().collect();
    let exe = ExeMeta::from_argv0(argv.first());

    let code = match cli::parse(&argv) {
        Ok(Invocation::Help(text)) => {
            print!("{text}");
            let _ = std::io::stdout().flush();
            0
        }
        Ok(Invocation::Version(text)) => {
            print!("{text}");
            let _ = std::io::stdout().flush();
            0
        }
        Ok(Invocation::CompareInput(cmd)) => match compare_input_with_fortran(cmd, exe) {
            Ok(status) => status,
            Err(err) => {
                eprintln!("error: {err:#}");
                EXIT_USAGE
            }
        },
        Ok(Invocation::Run(cmd)) => match execute(cmd, exe) {
            Ok(status) => status,
            Err(err) => {
                eprintln!("error: {err:#}");
                EXIT_USAGE
            }
        },
        Err(message) => {
            eprintln!("error: {message}");
            eprintln!("usage: {} [options] <path/to/SnapWave.inp>", exe.name);
            eprintln!("try `{} --help` for more information", exe.name);
            EXIT_USAGE
        }
    };
    std::process::exit(code);
}

/// Run one model invocation; returns the process exit code (0, or a
/// Fortran status passed through unchanged).
fn execute(cmd: cli::RunCommand, exe: ExeMeta) -> Result<i32> {
    let ctx = RunContext::new(cmd.input, exe, LogPrefs { verbose: cmd.verbose })?;

    // plan.md Phase 4: parse, validate and resolve the entire configuration
    // in Rust (defaults, post-processing, diagnostics). The Fortran facade
    // receives the resolved config as canonical key=value text and no
    // longer reads SnapWave.inp or decides defaults on this route.
    let config = input::parse_file(&ctx.input_path)?;
    diagnostics::report_input_diagnostics(&config, ctx.log.verbose);

    // plan.md Phase 5: resolve every file reference against the input
    // file's directory and make the output-directory policy explicit:
    // missing output directories are created (and unusable ones
    // rejected) in Rust, before the Fortran core runs.
    let run_paths = paths::RunPaths::resolve(&ctx.run_dir, &config);

    // After parsing: the run context can now include the resolved
    // output paths (Phase 5).
    if ctx.log.verbose {
        eprintln!("{}", ctx.describe(&run_paths.outputs));
    }

    for dir in run_paths.outputs.prepare()? {
        // A filesystem action the legacy binary never took; report it
        // so users see where output directories come from.
        eprintln!("created output directory {}", dir.display());
    }

    let text = config.to_config_text();
    let c_text = std::ffi::CString::new(text)
        .with_context(|| "config text contains an embedded NUL byte")?;
    let config_len = c_text.as_bytes().len() as c_int;

    // Legacy chdir contract — isolated in RunContext until the plan.md
    // Phase 6-7 readers accept explicit paths.
    ctx.enter_run_dir()?;

    let status = unsafe { snapwave_run_c(c_text.as_ptr(), config_len) };

    // Preserve the Fortran exit status semantics: non-zero status fails
    // the process with the same code.
    Ok(status)
}

/// `--compare-input`: parse the input in Rust, run the legacy Fortran
/// reader through the temporary facade hook, and compare every resulting
/// global. Also verifies the Phase 4 resolved-config handoff (Rust ->
/// text -> Fortran globals) through a second hook. Both comparisons must
/// agree for the test to pass.
fn compare_input_with_fortran(cmd: cli::RunCommand, exe: ExeMeta) -> Result<i32> {
    let ctx = RunContext::new(cmd.input, exe, LogPrefs { verbose: cmd.verbose })?;
    let config = input::parse_file(&ctx.input_path)?;

    // Parse-only: the hooks open nothing but the input file / config text
    // itself, so neither the run-directory chdir nor a testcase copy is
    // needed here.

    // ---- (a) legacy reader comparison (Phase 3) --------------------------
    let c_path = path_to_cstring(&ctx.input_path)?;
    let dump_path =
        std::env::temp_dir().join(format!("snapwave-input-dump-{}.txt", std::process::id()));
    let c_dump = path_to_cstring(&dump_path)?;
    // The facade buffers are character(len=1024).
    for (what, len) in [("input", c_path.as_bytes().len()), ("dump", c_dump.as_bytes().len())] {
        if len > 1024 {
            bail!("{what} path is too long for the FFI facade (>1024 bytes)");
        }
    }

    let status = unsafe {
        snapwave_read_input_dump_c(
            c_path.as_ptr(),
            c_path.as_bytes().len() as c_int,
            c_dump.as_ptr(),
            c_dump.as_bytes().len() as c_int,
        )
    };
    if status != 0 {
        let _ = std::fs::remove_file(&dump_path);
        bail!("the Fortran input reader (read_snapwave_input) failed with status {status}");
    }

    let dump_text = std::fs::read_to_string(&dump_path)
        .with_context(|| format!("reading the Fortran input dump at {}", dump_path.display()))?;
    let _ = std::fs::remove_file(&dump_path);

    let count = input_compare::check(&config, &dump_text)
        .with_context(|| format!("comparing the Rust and Fortran parses of {}", ctx.input_path.display()))?;

    if ctx.log.verbose {
        eprintln!("legacy reader: input parse matches the Fortran globals ({count} values compared)");
    }

    // ---- (b) resolved-config handoff comparison (Phase 4) -----------------
    let text = config.to_config_text();
    let c_text = std::ffi::CString::new(text)
        .with_context(|| "config text contains an embedded NUL byte")?;
    let dump2_path =
        std::env::temp_dir().join(format!("snapwave-resolved-dump-{}.txt", std::process::id()));
    let c_dump2 = path_to_cstring(&dump2_path)?;
    // Only the dump *path* is limited to 1024 bytes (the facade uses
    // character(len=1024) for file paths); the config text itself is
    // dynamically allocated in Fortran and has no such limit.
    if c_dump2.as_bytes().len() > 1024 {
        bail!("dump path is too long for the FFI facade (>1024 bytes)");
    }

    let status = unsafe {
        snapwave_load_config_dump_c(
            c_text.as_ptr(),
            c_text.as_bytes().len() as c_int,
            c_dump2.as_ptr(),
            c_dump2.as_bytes().len() as c_int,
        )
    };
    if status != 0 {
        let _ = std::fs::remove_file(&dump2_path);
        bail!("the Fortran resolved-config reader (read_resolved_input) failed with status {status}");
    }

    let dump2_text = std::fs::read_to_string(&dump2_path)
        .with_context(|| format!("reading the Fortran resolved dump at {}", dump2_path.display()))?;
    let _ = std::fs::remove_file(&dump2_path);

    let count2 = input_compare::check(&config, &dump2_text).with_context(|| {
        format!("comparing the Rust config and the resolved Fortran globals of {}", ctx.input_path.display())
    })?;

    if ctx.log.verbose {
        eprintln!("resolved handoff: config round-trips through Fortran ({count2} values compared)");
    }
    Ok(0)
}
