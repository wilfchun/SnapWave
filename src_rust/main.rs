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
//! Phase 6 (plan.md): the auxiliary *text* input readers are migrated to
//! Rust (`text_input`): observation points, single-point JONSWAP and
//! space/time-varying boundary files, wind, and the enclosure/neumann
//! polylines are parsed into Rust-owned structs and validated before the
//! model runs. The `--compare-text` mode pins the Rust parsers against the
//! unchanged Fortran readers through the temporary `snapwave_text_dump_c`
//! hook (`text_compare`).
//!
//! Phase 7 (plan.md): NetCDF input and output are Rust-owned. The mesh
//! reader (`mesh`, a port of `nc_read_net`) is pinned against the Fortran
//! oracle through `--compare-mesh`; the run path drives Fortran in capture
//! mode (`snapwave_run_capture_c` streams the output-time state to a temp
//! file) and the Rust writers (`output`, built on the hand-rolled classic
//! NetCDF writer in `netcdf`) write the map/history files.
//!
//! Phase 8 (plan.md): the non-solver domain data is explicit Rust state
//! (`state::DomainState`: config, mesh, boundary forcing, wind, obs
//! points, polylines, runtime scalars) and the run path hands it to the
//! coarse Fortran entry point `snapwave_run_capture_state_c` as one
//! `#[repr(C)]` struct of Rust-owned buffers — Fortran associates its
//! module globals with that memory instead of re-reading the files
//! (`ffi_layout` pins the one-based/column-major conversion facts).
//! Grid formats the Rust mesh reader does not cover (structured
//! index/mask, ASCII meshes) keep the Fortran readers, as does the
//! `--legacy-mesh` parity hook used by `tests/domain_state.rs`.
//!
//! Phase 9 (plan.md): the derived geometry — surrounding points, upwind
//! neighbours, observation interpolation weights, boundary support-point
//! mapping — is ported to Rust (`geometry`, `interp`, `date`) and pinned
//! against the unchanged Fortran routines by the `--compare-geometry` mode
//! (`geometry_compare`, via the `snapwave_geometry_dump_c` hook). Fortran
//! remains the runtime authority for the geometry; the ports exist to
//! prove parity.
//!
//! Phase 10 (plan.md): the time loop and output scheduling move to Rust.
//! Instead of Fortran's `run_time_loop`, Rust now drives the loop through
//! coarse entry points: `snapwave_init_capture_c` (or the legacy variant)
//! initialises the model and opens the capture stream, `snapwave_timestep_c`
//! computes one solver step, `snapwave_capture_{map,his}_c` capture output
//! when Rust's scheduling decides it is due, and
//! `snapwave_finalize_capture_c` closes the stream. The existing
//! `snapwave_run_capture_c` / `snapwave_run_capture_state_c` stay available
//! for the `--legacy-mesh` parity hook and the comparison hooks.
//!
//! Status codes: 0 on success (including `--help`/`--version`), 2 on
//! wrapper-detected errors, and any non-zero Fortran status is passed
//! through unchanged.

mod capture;
mod cli;
mod date;
mod diagnostics;
mod ffi_layout;
mod geometry;
mod geometry_compare;
mod input;
mod input_compare;
mod interp;
mod mesh;
mod netcdf;
mod output;
mod paths;
mod run_context;
mod state;
mod text_compare;
mod text_input;

use anyhow::Context;
use std::ffi::{c_char, c_int, OsString};
use std::io::Write;

use anyhow::{bail, Result};

use cli::{Invocation, EXIT_USAGE};
use run_context::{path_to_cstring, ExeMeta, LogPrefs, RunContext};

// FFI boundary (AGENTS.md): signatures must match src/snapwave_c_api.f90
// exactly (`bind(C)`; explicit length, no reliance on NUL termination).
extern "C" {
    // plan.md Phase 7: run the model in capture mode — Fortran streams the
    // output-time state to `capture_path` instead of writing NetCDF, and the
    // Rust writer (src_rust/output.rs) replays it into the real files.
    fn snapwave_run_capture_c(
        config: *const c_char,
        config_len: c_int,
        capture_path: *const c_char,
        capture_path_len: c_int,
    ) -> c_int;

    // plan.md Phase 8, step 5: the coarse Fortran entry point that
    // consumes Rust-prepared state. `state` points at a
    // `state::SnapWaveStateC` (the #[repr(C)] mirror of Fortran's
    // snapwave_state_t) whose buffers Rust owns and keeps alive for the
    // duration of the call.
    fn snapwave_run_capture_state_c(
        config: *const c_char,
        config_len: c_int,
        capture_path: *const c_char,
        capture_path_len: c_int,
        state: *const state::SnapWaveStateC,
    ) -> c_int;

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

    // Phase 6 comparison hook: read the auxiliary text inputs with the
    // unchanged Fortran readers and dump the resulting globals, pinning the
    // Rust parsers in `text_input` against the numerical oracle.
    fn snapwave_text_dump_c(
        config: *const c_char,
        config_len: c_int,
        dump_path: *const c_char,
        dump_path_len: c_int,
    ) -> c_int;

    // Phase 7 step 2 comparison hook: read the mesh NetCDF with the
    // unchanged nc_read_net reader and dump the resulting globals, pinning
    // the Rust mesh reader in `mesh` against the numerical oracle.
    fn snapwave_mesh_dump_c(
        config: *const c_char,
        config_len: c_int,
        dump_path: *const c_char,
        dump_path_len: c_int,
    ) -> c_int;

    // Phase 9 comparison hook: compute the derived geometry (surrounding
    // points, upwind neighbours, observation weights, boundary mapping) with
    // the unchanged Fortran routines and dump the resulting globals, pinning
    // the Rust ports in `geometry`/`interp` against the numerical oracle.
    fn snapwave_geometry_dump_c(
        config: *const c_char,
        config_len: c_int,
        dump_path: *const c_char,
        dump_path_len: c_int,
    ) -> c_int;

    // ---- plan.md Phase 10: Rust-owned time loop entry points ----------

    // Initialise the model with Rust-owned state, open the capture stream
    // and write static output data, but do NOT run the time loop.
    fn snapwave_init_capture_c(
        config: *const c_char,
        config_len: c_int,
        capture_path: *const c_char,
        capture_path_len: c_int,
        state: *const state::SnapWaveStateC,
    ) -> c_int;

    // Initialise the model by reading files (legacy route), open the
    // capture stream and write static output data, but do NOT run the
    // time loop.
    fn snapwave_init_legacy_capture_c(
        config: *const c_char,
        config_len: c_int,
        capture_path: *const c_char,
        capture_path_len: c_int,
    ) -> c_int;

    // One solver timestep: update_boundary_conditions(t) + compute_wave_field(t).
    fn snapwave_timestep_c(
        t: f64,
        it: c_int,
    ) -> c_int;

    // Capture map output at the current time (writes to the capture stream).
    fn snapwave_capture_map_c(
        t: f64,
        ntmapout: c_int,
    ) -> c_int;

    // Capture history output at the current time: update_obs_points() +
    // ncoutput_update_his (writes to the capture stream).
    fn snapwave_capture_his_c(
        t: f64,
        nthisout: c_int,
    ) -> c_int;

    // Finalize the capture run: ncoutput_finalize() + ncoutput_capture_end().
    fn snapwave_finalize_capture_c() -> c_int;
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
        Ok(Invocation::CompareText(cmd)) => match compare_text_with_fortran(cmd, exe) {
            Ok(status) => status,
            Err(err) => {
                eprintln!("error: {err:#}");
                EXIT_USAGE
            }
        },
        Ok(Invocation::CompareMesh(cmd)) => match compare_mesh_with_fortran(cmd, exe) {
            Ok(status) => status,
            Err(err) => {
                eprintln!("error: {err:#}");
                EXIT_USAGE
            }
        },
        Ok(Invocation::CompareGeometry(cmd)) => match compare_geometry_with_fortran(cmd, exe) {
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

    // plan.md Phase 6: parse and validate every auxiliary text input in
    // Rust, so the data is Rust-owned before the timestep loop starts.
    // Since Phase 8 the obs points, polylines and boundary series also
    // cross to Fortran as Rust state (below); the wind input stays
    // Fortran-read until the Phase 9 value-or-file interpolation moves,
    // and the Rust parse is validated against the Fortran readers by
    // `--compare-text`.
    let text_inputs = text_input::parse_all(&run_paths, &config)?;
    if ctx.log.verbose {
        diagnostics::report_text_input_diagnostics(&text_inputs);
    }

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

    // plan.md Phase 7: run the model in capture mode. Fortran streams the
    // output-time state to a temp file (instead of writing NetCDF) and the
    // Rust writers below replay it into the real map/history files.
    let capture_path =
        std::env::temp_dir().join(format!("snapwave-capture-{}.bin", std::process::id()));
    let c_capture = path_to_cstring(&capture_path)?;
    // The facade buffers file paths as character(len=1024).
    if c_capture.as_bytes().len() > 1024 {
        bail!("capture path is too long for the FFI facade (>1024 bytes)");
    }

    // plan.md Phase 8: when the Rust mesh reader owns this gridfile,
    // hand the whole non-solver domain state (mesh, polylines, obs
    // points, boundary series) to Fortran as Rust-owned buffers instead
    // of letting the Fortran readers re-read the files. Other grid
    // formats (structured index/mask, ASCII meshes) keep the Fortran
    // readers, as does the `--legacy-mesh` parity hook. The dispatch
    // mirrors the `ext == 'nc'` check of initialize_snapwave_domain
    // verbatim (see state::rust_owns_gridfile).
    let use_rust_state =
        !cmd.legacy_mesh && state::rust_owns_gridfile(&config.grid.gridfile);

    let ffi_state = if use_rust_state {
        let gridfile = run_paths
            .gridfile
            .as_ref()
            .context("the NetCDF mesh route requires a gridfile")?;
        let mesh =
            mesh::read_ugrid_netcdf(gridfile, config.grid.posdwn, config.grid.sferic)?;
        let domain = state::DomainState::new(&config, &mesh, &text_inputs);
        if ctx.log.verbose {
            diagnostics::report_domain_state_diagnostics(&domain);
        }
        Some(state::FfiState::build(&domain)?)
    } else {
        None
    };

    // Legacy chdir contract — isolated in RunContext until the plan.md
    // Phase 9 readers (wind/fw value-or-file inputs) accept explicit
    // paths. The Phase 8 handoff removed the mesh/text-file reads on
    // the Rust-state route, but the remaining Fortran-side file opens
    // still resolve CWD-relative.
    ctx.enter_run_dir()?;

    // plan.md Phase 10: Rust-owned time loop and output scheduling.
    // When `--fortran-time-loop` is set (a hidden test hook), use the
    // old Fortran-owned time loop path for parity comparison.
    if cmd.fortran_time_loop {
        let status = unsafe {
            match &ffi_state {
                Some(ffi) => {
                    let c_state = ffi.c_state();
                    snapwave_run_capture_state_c(
                        c_text.as_ptr(),
                        c_text.as_bytes().len() as c_int,
                        c_capture.as_ptr(),
                        c_capture.as_bytes().len() as c_int,
                        &c_state as *const state::SnapWaveStateC,
                    )
                }
                None => snapwave_run_capture_c(
                    c_text.as_ptr(),
                    c_text.as_bytes().len() as c_int,
                    c_capture.as_ptr(),
                    c_capture.as_bytes().len() as c_int,
                ),
            }
        };

        if status != 0 {
            let _ = std::fs::remove_file(&capture_path);
            return Ok(status);
        }
    } else {
        // ---- Phase 10: Rust-owned time loop -------------------------------
        let init_status = unsafe {
            match &ffi_state {
                Some(ffi) => {
                    let c_state = ffi.c_state();
                    snapwave_init_capture_c(
                        c_text.as_ptr(),
                        c_text.as_bytes().len() as c_int,
                        c_capture.as_ptr(),
                        c_capture.as_bytes().len() as c_int,
                        &c_state as *const state::SnapWaveStateC,
                    )
                }
                None => snapwave_init_legacy_capture_c(
                    c_text.as_ptr(),
                    c_text.as_bytes().len() as c_int,
                    c_capture.as_ptr(),
                    c_capture.as_bytes().len() as c_int,
                ),
            }
        };

        if init_status != 0 {
            let _ = std::fs::remove_file(&capture_path);
            return Ok(init_status);
        }

        let mut model = state::ModelState::new(&config);
        let nobs = text_inputs.obs.as_ref().map_or(0, |o| o.points.len() as i32);
        let map_file_nonempty = !config.output.map_file.is_empty();
        let his_file_nonempty = !config.output.his_file.is_empty();

        while model.is_running() {
            model.advance_iteration();

            let step_status = unsafe { snapwave_timestep_c(model.t, model.it) };
            if step_status != 0 {
                unsafe { snapwave_finalize_capture_c(); }
                let _ = std::fs::remove_file(&capture_path);
                return Ok(step_status);
            }

            // History output check (mirrors the Fortran order: his before map).
            if model.should_output_his(his_file_nonempty, nobs) {
                model.record_his_output(config.output.his_interval);
                let his_status =
                    unsafe { snapwave_capture_his_c(model.t, model.his_output_count) };
                if his_status != 0 {
                    unsafe { snapwave_finalize_capture_c(); }
                    let _ = std::fs::remove_file(&capture_path);
                    return Ok(his_status);
                }
            }

            // Map output check.
            if model.should_output_map(config.output.ja_save_each_iter, map_file_nonempty) {
                model.record_map_output(config.output.map_interval);
                let map_status =
                    unsafe { snapwave_capture_map_c(model.t, model.map_output_count) };
                if map_status != 0 {
                    unsafe { snapwave_finalize_capture_c(); }
                    let _ = std::fs::remove_file(&capture_path);
                    return Ok(map_status);
                }
            }

            model.advance_time();
        }

        let final_status = unsafe { snapwave_finalize_capture_c() };
        if final_status != 0 {
            let _ = std::fs::remove_file(&capture_path);
            return Ok(final_status);
        }
    }

    let capture = capture::read_capture(&capture_path, &config)?;
    let _ = std::fs::remove_file(&capture_path);

    if let Some(sm) = &capture.static_map {
        if let Some(map_path) = &run_paths.outputs.map {
            output::write_map(map_path, &config, sm, &capture.map_records)?;
        }
    }
    if let Some(sh) = &capture.static_his {
        if let Some(his_path) = &run_paths.outputs.his {
            output::write_his(his_path, &config, sh, &capture.his_records)?;
        }
    }

    Ok(0)
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

/// `--compare-text`: parse the auxiliary text inputs in Rust, run the
/// unchanged Fortran readers through the temporary Phase 6 hook, and
/// compare the resulting globals. Exits 0 on agreement, without running
/// the model (plan.md Phase 6, step 2).
fn compare_text_with_fortran(cmd: cli::RunCommand, exe: ExeMeta) -> Result<i32> {
    let ctx = RunContext::new(cmd.input, exe, LogPrefs { verbose: cmd.verbose })?;
    let config = input::parse_file(&ctx.input_path)?;
    let run_paths = paths::RunPaths::resolve(&ctx.run_dir, &config);

    // Rust-parse the auxiliary text inputs (family 1-5 of plan.md Phase 6).
    let rust = text_input::parse_all(&run_paths, &config)?;
    if ctx.log.verbose {
        diagnostics::report_text_input_diagnostics(&rust);
    }

    // The Fortran hook reads the mesh and the sibling text files relative
    // to the run directory (same chdir contract as a real run).
    ctx.enter_run_dir()?;

    let text = config.to_config_text();
    let c_text = std::ffi::CString::new(text)
        .with_context(|| "config text contains an embedded NUL byte")?;

    let dump_path =
        std::env::temp_dir().join(format!("snapwave-text-dump-{}.txt", std::process::id()));
    let c_dump = path_to_cstring(&dump_path)?;
    // Only the dump *path* is limited to 1024 bytes (the facade uses
    // character(len=1024) for file paths); the config text is dynamically
    // allocated in Fortran.
    if c_dump.as_bytes().len() > 1024 {
        bail!("dump path is too long for the FFI facade (>1024 bytes)");
    }

    let status = unsafe {
        snapwave_text_dump_c(
            c_text.as_ptr(),
            c_text.as_bytes().len() as c_int,
            c_dump.as_ptr(),
            c_dump.as_bytes().len() as c_int,
        )
    };
    if status != 0 {
        let _ = std::fs::remove_file(&dump_path);
        bail!("the Fortran text readers failed with status {status}");
    }

    let dump_text = std::fs::read_to_string(&dump_path)
        .with_context(|| format!("reading the Fortran text dump at {}", dump_path.display()))?;
    let _ = std::fs::remove_file(&dump_path);

    let count = text_compare::check(&rust, &dump_text).with_context(|| {
        format!("comparing the Rust and Fortran text-input parses of {}", ctx.input_path.display())
    })?;

    if ctx.log.verbose {
        eprintln!("text inputs: Rust parse matches the Fortran globals ({count} values compared)");
    }
    Ok(0)
}

/// `--compare-mesh`: read the mesh NetCDF in Rust, run the unchanged
/// Fortran `nc_read_net` reader through the temporary Phase 7 hook, and
/// compare the resulting globals. Exits 0 on agreement, without running the
/// model (plan.md Phase 7, step 2).
fn compare_mesh_with_fortran(cmd: cli::RunCommand, exe: ExeMeta) -> Result<i32> {
    let ctx = RunContext::new(cmd.input, exe, LogPrefs { verbose: cmd.verbose })?;
    let config = input::parse_file(&ctx.input_path)?;
    let run_paths = paths::RunPaths::resolve(&ctx.run_dir, &config);

    let gridfile = run_paths
        .gridfile
        .as_ref()
        .context("mesh comparison requires a gridfile")?;
    let rust_mesh = mesh::read_ugrid_netcdf(gridfile, config.grid.posdwn, config.grid.sferic)?;
    if ctx.log.verbose {
        eprintln!(
            "mesh: read {} nodes, {} faces (max {} nodes/face) in Rust",
            rust_mesh.no_nodes, rust_mesh.no_faces, rust_mesh.max_nodes
        );
    }

    // The Fortran hook reads the gridfile (and any sibling inputs) relative
    // to the run directory (same chdir contract as a real run).
    ctx.enter_run_dir()?;

    let text = config.to_config_text();
    let c_text = std::ffi::CString::new(text)
        .with_context(|| "config text contains an embedded NUL byte")?;

    let dump_path =
        std::env::temp_dir().join(format!("snapwave-mesh-dump-{}.txt", std::process::id()));
    let c_dump = path_to_cstring(&dump_path)?;
    // The facade buffers file paths as character(len=1024).
    if c_dump.as_bytes().len() > 1024 {
        bail!("dump path is too long for the FFI facade (>1024 bytes)");
    }

    let status = unsafe {
        snapwave_mesh_dump_c(
            c_text.as_ptr(),
            c_text.as_bytes().len() as c_int,
            c_dump.as_ptr(),
            c_dump.as_bytes().len() as c_int,
        )
    };
    if status != 0 {
        let _ = std::fs::remove_file(&dump_path);
        bail!("the Fortran mesh reader (nc_read_net) failed with status {status}");
    }

    let dump_text = std::fs::read_to_string(&dump_path)
        .with_context(|| format!("reading the Fortran mesh dump at {}", dump_path.display()))?;
    let _ = std::fs::remove_file(&dump_path);

    let count = mesh::check(&rust_mesh, &dump_text).with_context(|| {
        format!("comparing the Rust and Fortran mesh reads of {}", gridfile.display())
    })?;

    if ctx.log.verbose {
        eprintln!("mesh: Rust read matches the Fortran nc_read_net ({count} values compared)");
    }
    Ok(0)
}

/// `--compare-geometry`: compute the derived geometry (surrounding points,
/// upwind neighbours, observation interpolation weights, boundary
/// support-point mapping) in Rust, run the unchanged Fortran routines
/// through the temporary Phase 9 hook, and compare the resulting globals.
/// Exits 0 on agreement, without running the model (plan.md Phase 9).
fn compare_geometry_with_fortran(cmd: cli::RunCommand, exe: ExeMeta) -> Result<i32> {
    let ctx = RunContext::new(cmd.input, exe, LogPrefs { verbose: cmd.verbose })?;
    let config = input::parse_file(&ctx.input_path)?;
    let run_paths = paths::RunPaths::resolve(&ctx.run_dir, &config);

    let gridfile = run_paths
        .gridfile
        .as_ref()
        .context("geometry comparison requires a gridfile")?;
    if !state::rust_owns_gridfile(&config.grid.gridfile) {
        bail!(
            "geometry comparison only supports NetCDF meshes (gridfile '{}' is not one); \
             the ASCII/structured mesh geometry is not yet Rust-owned",
            config.grid.gridfile
        );
    }

    // Rust-side geometry from the Rust-owned mesh and text inputs.
    let mesh = mesh::read_ugrid_netcdf(gridfile, config.grid.posdwn, config.grid.sferic)?;
    let text = text_input::parse_all(&run_paths, &config)?;
    let geometry = geometry_compare::compute_geometry(&mesh, &config, &text);
    if ctx.log.verbose {
        eprintln!(
            "geometry: Rust computed {} surrounding-point lists, {} upwind directions, {} boundary nodes",
            geometry.domain.no_nodes, geometry.domain.ntheta360, geometry.domain.nb
        );
    }

    // The Fortran hook reads the mesh and sibling boundary files relative to
    // the run directory (same chdir contract as a real run).
    ctx.enter_run_dir()?;

    let text_config = config.to_config_text();
    let c_text = std::ffi::CString::new(text_config)
        .with_context(|| "config text contains an embedded NUL byte")?;

    let dump_path =
        std::env::temp_dir().join(format!("snapwave-geometry-dump-{}.txt", std::process::id()));
    let c_dump = path_to_cstring(&dump_path)?;
    // The facade buffers file paths as character(len=1024).
    if c_dump.as_bytes().len() > 1024 {
        bail!("dump path is too long for the FFI facade (>1024 bytes)");
    }

    let status = unsafe {
        snapwave_geometry_dump_c(
            c_text.as_ptr(),
            c_text.as_bytes().len() as c_int,
            c_dump.as_ptr(),
            c_dump.as_bytes().len() as c_int,
        )
    };
    if status != 0 {
        let _ = std::fs::remove_file(&dump_path);
        bail!("the Fortran geometry routines (initialize_snapwave_domain/read_obs_points/read_boundary_data) failed with status {status}");
    }

    let dump_text = std::fs::read_to_string(&dump_path)
        .with_context(|| format!("reading the Fortran geometry dump at {}", dump_path.display()))?;
    let _ = std::fs::remove_file(&dump_path);

    let count = geometry_compare::check(&geometry, &dump_text).with_context(|| {
        format!("comparing the Rust and Fortran geometry of {}", gridfile.display())
    })?;

    if ctx.log.verbose {
        eprintln!("geometry: Rust computation matches the Fortran routines ({count} values compared)");
    }
    Ok(0)
}
