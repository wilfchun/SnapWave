//! Retired (plan.md Phase 12).
//!
//! These tests pinned the Rust mesh reader against the legacy Fortran
//! `nc_read_net` reader through the temporary `--compare-mesh` hook. The
//! Fortran core is gone from the Rust build, so the hook and these tests are
//! retired. The Rust mesh reader is the sole authority and is unit-tested in
//! `src_rust/mesh.rs`; the map/history writers are pinned end-to-end by
//! `tests/regression.rs`.
