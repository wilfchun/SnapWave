//! Retired (plan.md Phase 12).
//!
//! These tests pinned the Phase 8 state handoff (`--legacy-mesh`) and the
//! Phase 10 Rust-owned time loop (`--fortran-time-loop`) against the legacy
//! Fortran routes. The Fortran core is gone from the Rust build, so those
//! parity routes and their tests are retired. The pure-Rust model run
//! (`src_rust/model.rs`) is pinned end-to-end by `tests/regression.rs`
//! against the committed baselines and the live Fortran oracle.
