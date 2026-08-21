//! SnapWave — fast, implicit, unstructured-grid short wave solver, rewritten
//! in Rust (plan.md Phase 12: the Fortran core is retired from the Rust
//! build; the production executable is Rust-owned end to end).
//!
//! Usage: snapwave [OPTIONS] <path/to/SnapWave.inp>
//!
//! Phase 2 (plan.md): Rust owns all process-level behaviour — argument
//! parsing (`cli`), input validation and the run context (`run_context`),
//! status-code semantics.
//!
//! Phase 3/4 (plan.md): the wrapper parses, validates and resolves the
//! entire configuration in Rust (`input`) — defaults, diagnostics and the
//! non-positive-interval checks — before any model work starts.
//!
//! Phase 5 (plan.md): filesystem and output-directory policy is Rust-owned
//! (`paths`); file references resolve against the input file's directory and
//! the required output directories are created or validated in Rust.
//!
//! Phase 6 (plan.md): the auxiliary *text* inputs (observation points,
//! single-point JONSWAP and space/time-varying boundaries, wind, enclosure
//! and Neumann polylines) are parsed into Rust-owned structs (`text_input`).
//!
//! Phase 7 (plan.md): NetCDF input and output are Rust-owned — the mesh
//! reader (`mesh`) and the classic-format map/history writers (`netcdf`,
//! `output`).
//!
//! Phase 9 (plan.md): the derived geometry — surrounding points, upwind
//! neighbours, observation weights, boundary support-point mapping — is
//! ported to Rust (`geometry`, `interp`) and, since Phase 12, is wired into
//! the solver instead of being recomputed by Fortran.
//!
//! Phase 10 (plan.md): the time loop and output scheduling are Rust-owned
//! (`state::ModelState`).
//!
//! Phase 11 (plan.md): the numerical solver — celerities, dispersion,
//! breaking, wind input, vegetation, and the implicit 4-sweep
//! `solve_energy_balance2Dstat` — is ported to Rust (`solver`).
//!
//! Phase 12 (plan.md): the boundary-condition update, observation-point
//! update and the full model orchestration are ported to Rust (`model`),
//! and the Fortran/C/NetCDF sources are removed from the Cargo build. The
//! legacy `make` build remains available as the numerical oracle.
//!
//! Status codes: 0 on success (including `--help`/`--version`), 2 on
//! wrapper-detected errors.

mod capture;
mod cli;
mod date;
mod diagnostics;
mod geometry;
mod input;
mod interp;
mod mesh;
mod model;
mod netcdf;
mod output;
mod paths;
mod run_context;
mod solver;
mod state;
mod text_input;

use std::ffi::OsString;
use std::io::Write;

use anyhow::{bail, Context, Result};

use cli::{Invocation, EXIT_USAGE};
use run_context::{ExeMeta, LogPrefs, RunContext};

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

/// Run one model invocation; returns the process exit code (0 on success).
fn execute(cmd: cli::RunCommand, exe: ExeMeta) -> Result<i32> {
    let ctx = RunContext::new(cmd.input, exe, LogPrefs { verbose: cmd.verbose })?;

    // plan.md Phase 4: parse, validate and resolve the entire configuration
    // in Rust (defaults, post-processing, diagnostics).
    let config = input::parse_file(&ctx.input_path)?;
    diagnostics::report_input_diagnostics(&config, ctx.log.verbose);

    // plan.md Phase 5: resolve every file reference against the input
    // file's directory and make the output-directory policy explicit.
    let run_paths = paths::RunPaths::resolve(&ctx.run_dir, &config);

    // plan.md Phase 6: parse and validate every auxiliary text input in
    // Rust, so the data is Rust-owned before the timestep loop starts.
    let text_inputs = text_input::parse_all(&run_paths, &config)?;
    if ctx.log.verbose {
        diagnostics::report_text_input_diagnostics(&text_inputs);
    }

    if ctx.log.verbose {
        eprintln!("{}", ctx.describe(&run_paths.outputs));
    }

    for dir in run_paths.outputs.prepare()? {
        // A filesystem action the legacy binary never took; report it
        // so users see where output directories come from.
        eprintln!("created output directory {}", dir.display());
    }

    // plan.md Phase 12: read the mesh in Rust and run the whole model in
    // Rust. Only NetCDF meshes are supported (the structured/ASCII mesh
    // readers stay Fortran and are rejected with a clear error).
    if !state::rust_owns_gridfile(&config.grid.gridfile) {
        bail!(
            "gridfile '{}' is not a NetCDF mesh; the Rust build only supports NetCDF meshes",
            config.grid.gridfile
        );
    }
    let gridfile = run_paths
        .gridfile
        .as_ref()
        .context("the model requires a gridfile")?;
    let mesh = mesh::read_ugrid_netcdf(gridfile, config.grid.posdwn, config.grid.sferic)?;
    if ctx.log.verbose {
        diagnostics::report_model_diagnostics(mesh.no_nodes, mesh.no_faces, mesh.max_nodes, &text_inputs);
    }

    let mut model = model::Model::new(&config, &mesh, &text_inputs)?;
    model.run(run_paths.outputs.map.as_deref(), run_paths.outputs.his.as_deref())?;

    Ok(0)
}
