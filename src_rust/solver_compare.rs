//! Retired (plan.md Phase 12).
//!
//! This module drove the temporary `--compare-solver` hook, which pinned
//! the Rust solver ports against the legacy Fortran solver. The Fortran core
//! is gone from the Rust build, so the hook and its comparison module are no
//! longer part of the crate. The solver ports are now driven directly by
//! `crate::model` with real geometry, and unit-tested in `crate::solver`.
