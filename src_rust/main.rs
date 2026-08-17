//! Thin Rust wrapper around the SnapWave Fortran core.
//!
//! Usage: snapwave <path/to/SnapWave.inp>
//!
//! The Fortran input reader (`read_snapwave_input`) probes the current
//! working directory for snapwave.inp / SnapWave.inp and resolves all
//! sibling input/output paths relative to it. The least invasive way to
//! keep the Fortran core unchanged is therefore to change the working
//! directory to the input file's parent and pass only the file name across
//! the FFI boundary (see plan.md, Phase 2).

use std::ffi::{c_char, c_int, CString};
use std::path::PathBuf;

use anyhow::{bail, Context, Result};

extern "C" {
    fn snapwave_run_c(path: *const c_char, path_len: c_int) -> c_int;
}

fn main() {
    if let Err(err) = run() {
        eprintln!("error: {err:#}");
        std::process::exit(2);
    }
}

fn run() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 2 {
        let prog = args.first().map(String::as_str).unwrap_or("snapwave");
        eprintln!("SnapWave (Rust wrapper around the Fortran core)");
        eprintln!("usage: {prog} <path/to/SnapWave.inp>");
        bail!("expected exactly one argument: path to a SnapWave.inp file");
    }

    // Resolve and validate the input file before touching Fortran.
    let input_arg = PathBuf::from(&args[1]);
    let input_path = input_arg
        .canonicalize()
        .with_context(|| format!("input file not found: {}", input_arg.display()))?;
    if !input_path.is_file() {
        bail!("input path is not a file: {}", input_path.display());
    }

    let file_name = input_path
        .file_name()
        .and_then(|n| n.to_str())
        .with_context(|| format!("cannot derive file name from {}", input_path.display()))?;
    let parent = input_path
        .parent()
        .with_context(|| format!("input path has no parent directory: {}", input_path.display()))?;

    // Legacy relative-path support: run from the input file's directory.
    std::env::set_current_dir(parent)
        .with_context(|| format!("failed to change working directory to {}", parent.display()))?;

    // Ownership: Fortran copies the characters immediately; the CString
    // simply outlives the call here.
    let c_path = CString::new(file_name)
        .with_context(|| format!("input file name is not a valid C string: {file_name}"))?;
    let path_len = c_path.as_bytes().len() as c_int;

    let status = unsafe { snapwave_run_c(c_path.as_ptr(), path_len) };

    // Preserve the Fortran exit status semantics: non-zero status fails
    // the process with the same code (Fortran `stop 1` paths still
    // terminate the process directly inside the facade).
    if status != 0 {
        std::process::exit(status);
    }

    Ok(())
}
