//! Retired (plan.md Phase 12).
//!
//! These tests pinned the Rust auxiliary-text parsers against the legacy
//! Fortran readers through the temporary `--compare-text` hook. The Fortran
//! core is gone from the Rust build, so the hook and these tests are retired.
//! The Rust parsers are the sole authority and are unit-tested in
//! `src_rust/text_input.rs`.
