//! Structured wrapper diagnostics (plan.md Phase 4, step 3).
//!
//! The legacy Fortran reader (`read_snapwave_input`) printed informational
//! messages to stdout ("Reading input file ...", "Wind growth turned on.",
//! etc.). On the Rust path those messages are now owned by the wrapper and
//! emitted through this module, so the Fortran facade no longer produces
//! configuration-related console chatter.

use crate::input::SnapWaveInput;

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