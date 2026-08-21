//! Retired (plan.md Phase 12).
//!
//! These tests pinned the Rust solver ports against the legacy Fortran solver
//! through the temporary `--compare-solver` hook. The Fortran core is gone
//! from the Rust build, so the hook and these tests are retired. The solver
//! ports are driven by the model run (`src_rust/model.rs`) and unit-tested in
//! `src_rust/solver.rs`; their end-to-end correctness is pinned by
//! `tests/regression.rs`.
