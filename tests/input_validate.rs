//! Retired (plan.md Phase 12).
//!
//! These tests pinned the Phase 4 validation (bad output intervals, the
//! resolved-config handoff) against the Fortran reader via the temporary
//! `--compare-input` hook. The Fortran core is gone from the Rust build, so
//! the hook-dependent tests are retired. Validation lives in
//! `src_rust/input.rs` (unit-tested) and `src_rust/paths.rs`; missing-input
//! behaviour is covered by `tests/cli.rs` and `tests/mwe.rs`.
