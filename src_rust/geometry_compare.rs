//! Retired (plan.md Phase 12).
//!
//! This module drove the temporary `--compare-geometry` hook, which pinned
//! the Rust geometry ports against the legacy Fortran routines. The Fortran
//! core is gone from the Rust build, so the hook and its comparison module
//! are no longer part of the crate. The geometry ports are now wired into
//! the model run (`crate::model`) and unit-tested in `crate::geometry` and
//! `crate::interp`.
