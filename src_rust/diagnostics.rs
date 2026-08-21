//! Structured wrapper diagnostics (plan.md Phase 4, step 3).
//!
//! The legacy Fortran reader (`read_snapwave_input`) printed informational
//! messages to stdout ("Reading input file ...", "Wind growth turned on.",
//! etc.). On the Rust path those messages are now owned by the wrapper and
//! emitted through this module, so the Fortran facade no longer produces
//! configuration-related console chatter.

use crate::input::SnapWaveInput;
use crate::state::DomainState;
use crate::text_input::{BoundaryInput, ParsedTextInputs, WindInput};

/// Emit the informational messages that the legacy Fortran reader used to
/// print, now structured and Rust-owned. Called after parsing and before
/// the model run.
pub fn report_input_diagnostics(cfg: &SnapWaveInput, verbose: bool) {
    if verbose {
        eprintln!("input: parsed and validated in Rust (plan.md Phase 3/4)");
    }
    // The wind-switch message is the one semantically-meaningful
    // informational line from the legacy reader; preserve it.
    if cfg.wind.enabled {
        eprintln!("   Wind growth turned on.");
    } else {
        eprintln!("   Uniform wave period in entire domain.");
    }
}

/// One-line verbose summary of the Rust-parsed auxiliary text inputs
/// (plan.md Phase 6). Emitted only in verbose mode.
pub fn report_text_input_diagnostics(text: &ParsedTextInputs) {
    let obs = match &text.obs {
        Some(o) => format!("{} points", o.len()),
        None => "none".to_string(),
    };
    let boundary = match &text.boundary {
        BoundaryInput::None => "none".to_string(),
        BoundaryInput::Single(j) => format!("single-point JONSWAP ({} records)", j.len()),
        BoundaryInput::Timeseries(s) => {
            format!("{} support points x {} time steps", s.nwbnd, s.ntwbnd)
        }
    };
    let wind = match &text.wind {
        WindInput::Uniform(u) => match (u.u10, u.u10dir_deg) {
            (Some(mag), Some(dir)) => format!("uniform ({mag} m/s @ {dir} deg)"),
            _ => "uniform (file-backed)".to_string(),
        },
        WindInput::List(l) => format!("list ({} records)", l.len()),
    };
    let enc = match &text.enclosure {
        Some(p) => format!("{} points", p.len()),
        None => "none".to_string(),
    };
    let neu = match &text.neumann {
        Some(p) => format!("{} points", p.len()),
        None => "none".to_string(),
    };
    eprintln!("text inputs: obs {obs}; boundary {boundary}; wind {wind}; enclosure {enc}; neumann {neu}");
}

/// Verbose summary of the Phase 8 state handoff: what the Rust-owned
/// domain state contains before the buffers cross to Fortran. Emitted
/// only in verbose mode, just before the facade call.
pub fn report_domain_state_diagnostics(domain: &DomainState) {
    let mode = match &domain.text.boundary {
        BoundaryInput::None => "none".to_string(),
        BoundaryInput::Single(j) => format!("single-point, {} records", j.len()),
        BoundaryInput::Timeseries(s) => {
            format!("timeseries, {} points x {} times", s.nwbnd, s.ntwbnd)
        }
    };
    let obs = domain.text.obs.as_ref().map_or(0, |o| o.len());
    eprintln!(
        "domain state: Rust-owned mesh ({} nodes, {} faces, max {} nodes/face) + boundary ({mode}) + {obs} obs points handed to the Fortran core (plan.md Phase 8)",
        domain.mesh.no_nodes, domain.mesh.no_faces, domain.mesh.max_nodes
    );
    eprintln!(
        "runtime state: tstart {} s, tstop {} s, timestep {} s, map interval {} s, his interval {} s",
        domain.runtime.tstart,
        domain.runtime.tstop,
        domain.runtime.timestep,
        domain.runtime.map_interval,
        domain.runtime.his_interval
    );
}