//! Retired (plan.md Phase 12).
//!
//! These tests pinned the Rust geometry ports against the legacy Fortran
//! routines through the temporary `--compare-geometry` hook. The Fortran core
//! is gone from the Rust build, so the hook and these tests are retired. The
//! geometry ports are wired into the model run (`src_rust/model.rs`) and
//! unit-tested in `src_rust/geometry.rs` and `src_rust/interp.rs`; their
//! end-to-end correctness is pinned by `tests/regression.rs`.
