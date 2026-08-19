//! Thin Rust wrapper around the SnapWave Fortran core.
//!
//! Usage: snapwave [OPTIONS] <path/to/SnapWave.inp>
//!
//! Phase 2 (plan.md): Rust owns all process-level behaviour — argument
//! parsing (`cli`), input validation and the run context (`run_context`),
//! status-code semantics — while Fortran still performs all model work
//! through the coarse C ABI facade in `src/snapwave_c_api.f90`.
//!
//! Phase 3 (plan.md): the wrapper additionally parses and validates
//! `SnapWave.inp` in Rust (`input`) before the Fortran core runs, so
//! invalid input is a wrapper error rather than a Fortran `stop`; the
//! `--compare-input` mode cross-checks the Rust parse against the legacy
//! Fortran reader through a temporary facade hook (`input_compare`).
//!
//! Status codes: 0 on success (including `--help`/`--version`), 2 on
//! wrapper-detected errors, and any non-zero Fortran status is passed
//! through unchanged.

mod cli;
mod input;
mod input_compare;
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
    fn snapwave_run_c(path: *const c_char, path_len: c_int) -> c_int;

    // Temporary Phase 3 comparison hook (plan.md Phase 3, step 5); removed
    // with the Fortran input reader once Phase 4+ retires it.
    fn snapwave_read_input_dump_c(
        input_path: *const c_char,
        input_path_len: c_int,
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

    if ctx.log.verbose {
        eprintln!("{}", ctx.describe());
    }

    // plan.md Phase 3: parse and validate SnapWave.inp in Rust before the
    // Fortran core runs, so invalid configuration is a wrapper error (and
    // never a Fortran `stop`/runtime abort). The Fortran reader remains
    // the numerical authority — it re-reads the same Rust-selected file —
    // until Phase 4 moves defaults and validation over.
    let _config = input::parse_file(&ctx.input_path)?;
    if ctx.log.verbose {
        eprintln!("input file parsed and validated in Rust (plan.md Phase 3)");
    }

    // Legacy chdir contract — isolated in RunContext until plan.md Phase 5
    // removes the need for a process-wide working directory.
    ctx.enter_run_dir()?;

    // Ownership: Fortran copies the characters immediately; the CString
    // simply outlives the call here.
    let c_path = ctx.input_file_name_cstring()?;
    let path_len = c_path.as_bytes().len() as c_int;

    let status = unsafe { snapwave_run_c(c_path.as_ptr(), path_len) };

    // Preserve the Fortran exit status semantics: non-zero status fails
    // the process with the same code (Fortran `stop 1` paths still
    // terminate the process directly inside the facade).
    Ok(status)
}

/// `--compare-input`: parse the input in Rust, run the legacy Fortran
/// reader through the temporary facade hook, and compare every resulting
/// global (plan.md Phase 3, step 5). Temporary scaffolding that goes away
/// again with the Phase 4+ input-reader migration.
fn compare_input_with_fortran(cmd: cli::RunCommand, exe: ExeMeta) -> Result<i32> {
    let ctx = RunContext::new(cmd.input, exe, LogPrefs { verbose: cmd.verbose })?;
    let config = input::parse_file(&ctx.input_path)?;

    // Parse-only: the hook opens nothing but the input file itself, so
    // neither the run-directory chdir nor a testcase copy is needed here.
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
        eprintln!("input parse matches the Fortran globals ({count} values compared)");
    }
    Ok(0)
}
