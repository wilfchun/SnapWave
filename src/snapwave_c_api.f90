! Retired (plan.md Phase 12).
!
! This module was the coarse C ABI facade between the Rust wrapper and the
! Fortran core. Since Phase 12 the model runs entirely in Rust, so no Rust
! runtime path calls Fortran any more and the facade is removed from the
! Cargo build. The unchanged Fortran solver lives on in the legacy `make`
! build (`src/snapwave.f90`), which remains the numerical oracle.
