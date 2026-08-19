//! Shared support code for the Phase-1 regression harness (plan.md Phase 1).
//!
//! Split into three modules:
//! - [`ncdf`]: dependency-free reader for NetCDF classic (CDF-1/CDF-2) files,
//!   so the tests do not depend on `ncdump` being installed.
//! - [`harness`]: testcase copying, Windows separator normalization, and
//!   running the Rust wrapper / legacy Fortran oracle on temp copies.
//! - [`compare`]: schema pinning plus numeric comparison with per-variable
//!   tolerances.
//!
//! See `tests/README.md` for how to add new cases.
//!
//! `allow(dead_code)`: this module is compiled into several test crates
//! (regression.rs, ncdf_parser.rs, cli.rs) that each use only part of the
//! API.

#![allow(dead_code)]

pub mod compare;
pub mod harness;
pub mod ncdf;
