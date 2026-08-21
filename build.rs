//! Retired (plan.md Phase 12).
//!
//! This file used to be the Cargo build orchestrator: it compiled the
//! Fortran core, the bundled Triangle C code, and emitted the NetCDF/OpenMP/
//! Fortran-runtime link directives. Since Phase 12 the model is 100% Rust and
//! `Cargo.toml` no longer sets `build = "build.rs"`, so nothing compiles or
//! links any Fortran/C/NetCDF source. Kept as a stub for history; the legacy
//! Fortran oracle is built by the `Makefile` instead.
