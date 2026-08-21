//! Retired (plan.md Phase 12).
//!
//! This module drove the temporary `--compare-input` hook, which pinned the
//! Rust input parser against the legacy Fortran reader. The Fortran core is
//! gone from the Rust build, so the hook and its comparison module are no
//! longer part of the crate. The Rust parser is now the sole authority and
//! is unit-tested in `crate::input`.
