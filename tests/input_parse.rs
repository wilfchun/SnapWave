//! Retired (plan.md Phase 12).
//!
//! These tests pinned the Rust input parser against the legacy Fortran
//! reader through the temporary `--compare-input` hook. The Fortran core is
//! gone from the Rust build, so the hook and these tests are retired. The
//! Rust parser is the sole authority and is unit-tested in
//! `src_rust/input.rs`; invalid-input wrapper behaviour is still covered by
//! `tests/cli.rs` and `tests/mwe.rs`.
