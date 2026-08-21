//! Retired (plan.md Phase 12).
//!
//! This module drove the temporary `--compare-text` hook, which pinned the
//! Rust auxiliary-text parsers against the legacy Fortran readers. The
//! Fortran core is gone from the Rust build, so the hook and its comparison
//! module are no longer part of the crate. The Rust parsers are now the sole
//! authority and are unit-tested in `crate::text_input`.
