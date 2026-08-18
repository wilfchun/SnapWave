# SnapWave Rust Migration Plan

This plan starts from the current state: a Cargo-built Rust binary owns
`main`, links the existing Fortran/C implementation, and calls the model through
the coarse facade in `src/snapwave_c_api.f90`. The completed wrapper/bootstrap
work is not repeated here. The remaining goal is to move SnapWave to 100% Rust
while keeping the Fortran implementation available as the numerical oracle
until each migrated path is proven equivalent.

## Migration Rules

- Preserve scientific behaviour. Each migration either matches the Fortran
  oracle within documented tolerances or is explicitly tracked as a physics
  change outside this rewrite plan.
- Migrate outside-in. Move command-line handling, configuration, filesystem
  behaviour, text readers, and NetCDF IO before data structures and solver
  numerics.
- Keep the FFI boundary coarse until a phase deliberately replaces a subsystem.
  Do not bind individual numerical routines as a shortcut.
- Keep `cargo build`, `cargo test`, and the legacy `make` build green while
  Fortran remains in the tree.
- Add regression coverage before replacing behaviour. Structural checks are
  enough for early IO work; solver-facing work needs numeric comparisons with
  tolerances.
- Treat Fortran module globals as shared legacy state. Prefer Rust-owned
  structs for new code and add narrow transfer points only where needed.

## Current Implementation Review

The current wrapper is a good strangler entry point: Rust validates the input
path, changes to the input directory to preserve legacy relative path semantics,
and calls one C ABI facade. The facade mirrors `src/snapwave.f90` closely, which
keeps the initial behavioural surface small.

Important follow-up items:

- `src/snapwave_c_api.f90` accepts an explicit input path but still calls
  `read_snapwave_input()`, which probes hard-coded file names. That is fine for
  the bootstrap but should be removed when input parsing moves to Rust.
- `read_snapwave_input()` still contains `stop 1` paths. Any error path made
  reachable from Rust-owned orchestration should become a status/error return.
- `flake.nix` currently packages the legacy Makefile executable, while Cargo
  builds the Rust wrapper. Align Nix checks/package outputs with the Cargo
  wrapper before relying on Nix as migration CI.
- The smoke test validates execution and NetCDF shape, not numeric equivalence.
  Numeric baselines are required before replacing domain, boundary, or solver
  behaviour.

## Phase 1: Test Oracle And Baselines

Goal: make the Fortran behaviour cheap to compare against before migrating more
code.

Steps:

1. Add a test harness that can run the Fortran oracle and the Rust-wrapped path
   on copied testcases without modifying checked-in inputs.
2. Store or generate baseline summaries for the MWE map and history NetCDF
   outputs: dimensions, variables, attributes, time coordinates, and selected
   numeric arrays.
3. Define numeric tolerances per output family. Start with pragmatic absolute
   and relative tolerances for `hm0`, `tp`, `wavdir`, `ee`, `depth`, and station
   outputs.
4. Add at least one broader validation case beyond the current coarse
   shoaling/refraction MWE before migrating any solver-adjacent code.
5. Normalize Windows path separators only in temporary testcase copies.

Acceptance:

- `cargo test` compares wrapper output to the Fortran oracle or committed
  baseline data.
- Test output identifies which variable, index/time, and tolerance failed.
- The test harness is documented enough to add new cases without reading the
  Fortran build internals.

Status: implemented (2026-08-18). The harness lives in `tests/support/` with
entry points in `tests/regression.rs` and is documented in `tests/README.md`.
Committed legacy outputs under `testcases/<case>/output/` serve as the
numeric baselines (case 31: history only — no committed map file exists;
case 32 curvi island — the broader validation case with JONSWAP boundary and
curvilinear mesh —: map only, because its committed history file predates
the removal of per-iteration history output). Wrapper-versus-oracle
comparison activates automatically when the legacy `make` binary exists or
`SNAPWAVE_ORACLE` points to one; it is strict, while committed-baseline
comparison allows record-count drift across compilers plus a documented
allowlist of schema additions the legacy files predate (`point_dirspr`).
Tolerances are defined per variable in `tests/support/compare.rs`; a
dependency-free classic-NetCDF reader (`tests/support/ncdf.rs`) provides
schema and numeric access so tests do not depend on `ncdump` availability.
Acceptance verified: `cargo test` passes end to end with the `make`-built
oracle present (wrapper output matches the oracle exactly and the committed
baselines within tolerance).

## Phase 2: Rust CLI And Run Context

Goal: move all process-level behaviour into Rust while Fortran still performs
model work.

Steps:

1. Replace manual argument handling with a small Rust CLI module. Add `--help`,
   `--version`, explicit input path handling, and consistent status code
   semantics.
2. Introduce a Rust `RunContext` containing the original input path, run
   directory, output directory expectations, executable metadata, and logging
   preferences.
3. Keep the legacy `chdir` contract until input parsing is migrated, but isolate
   it so later phases can remove it cleanly.
4. Add tests for missing input, directory input, invalid UTF-8 or embedded NUL
   file names where supported, and relative path handling.
5. Align `flake.nix` package and smoke test with the Cargo-built wrapper.

Acceptance:

- Rust owns all command-line validation and user-facing wrapper errors.
- The same testcase runs through `cargo run`, `cargo test`, and Nix checks.
- Fortran still receives only the minimum information needed to run.

## Phase 3: SnapWave.inp Parsing In Rust

Goal: parse and validate the main configuration file in Rust, while still
feeding equivalent values to the Fortran model.

Steps:

1. Document the exact current `SnapWave.inp` grammar: keyword matching,
   comments, whitespace, duplicate keys, defaults, case sensitivity, and value
   parsing quirks.
2. Add a Rust parser that preserves those semantics first. Do not clean up or
   modernize the format during the initial migration.
3. Represent configuration as typed Rust structs grouped by concern:
   time/control, grid/domain, boundary forcing, wind, output, vegetation, and
   diagnostics.
4. Add parser unit tests covering defaults and each known keyword in
   `src/snapwave_input.f90`.
5. Add a temporary Fortran comparison hook that reads the legacy globals after
   `read_snapwave_input()` and compares them to the Rust parse result.
6. Replace hard-coded Fortran filename probing with Rust-selected input paths.

Acceptance:

- Rust can parse every checked-in `SnapWave.inp`.
- Rust parse results match Fortran globals for representative cases.
- Wrapper failures for invalid input are Rust errors, not Fortran `stop`.

## Phase 4: Config Defaults, Validation, And Diagnostics

Goal: move configuration defaults and non-numerical diagnostics out of Fortran.

Steps:

1. Move all default values from `read_snapwave_input()` into Rust constants or
   typed defaults.
2. Implement validation for non-physical or unsupported settings already
   rejected by Fortran, such as non-positive output intervals.
3. Preserve legacy warning and informational behaviour where tests depend on it,
   but use structured Rust diagnostics internally.
4. Add validation tests for bad intervals, missing required files, unsupported
   combinations, and optional output settings.
5. Convert remaining input-reader `stop` paths used by the facade into status
   returns or remove those paths from the facade route.

Acceptance:

- Fortran no longer decides main config defaults for the Rust path.
- Invalid configuration fails before domain initialization.
- Legacy and Rust paths agree on accepted testcase configuration.

## Phase 5: Filesystem And Output Directory Handling

Goal: make all path resolution and output directory behaviour explicit in Rust.

Steps:

1. Centralize path resolution relative to the input file directory.
2. Model output paths as Rust `PathBuf`s rather than fixed-width Fortran
   strings on the Rust-owned path.
3. Create or validate required output directories in Rust.
4. Preserve legacy handling of `none`, empty strings, and relative paths until a
   format change is intentionally introduced.
5. Add tests for nested run directories, missing output parents, disabled map
   output, disabled history output, and Windows-style separators in temp copies.

Acceptance:

- Rust owns run directory and output directory policy.
- The wrapper no longer needs global `chdir` after downstream readers accept
  explicit paths.
- Existing testcases continue to run unchanged from the user's perspective.

## Phase 6: Text Input Readers

Goal: migrate auxiliary text readers before NetCDF and numerical domain logic.

Order:

1. Observation points from `src/snapwave_obspoints.f90`.
2. Single-point JONSWAP boundary files.
3. Boundary location and time-series files: `bndfile`, `bhsfile`, `btpfile`,
   `bwdfile`, `bdsfile`, `bzsfile`.
4. Wind list files and uniform/file-backed wind inputs.
5. Boundary enclosure and Neumann boundary files.
6. Plain text mesh/sample readers currently embedded in `snapwave_domain.f90`.

Steps:

1. For each reader, write Rust parsing tests against checked-in testcase files.
2. Compare Rust parsed data to Fortran globals through a temporary oracle
   facade.
3. Introduce Rust-owned data structs and pass only required data to Fortran at
   a coarse boundary.
4. Remove or bypass the matching Fortran reader only after comparison tests
   pass.

Acceptance:

- Auxiliary text input data is Rust-owned before the timestep loop starts.
- Fortran reader code remains available until each reader is replaced and
  covered.
- No parser migration changes numerical output beyond defined tolerances.

## Phase 7: NetCDF Input And Output

Goal: move NetCDF schema handling and file IO to Rust while leaving solver state
updates in Fortran until later phases.

Steps:

1. Choose a Rust NetCDF strategy that works under the existing Nix and system
   toolchains. Prefer bindings to the installed NetCDF library over vendoring.
2. Port `nc_read_net()` behaviour first for mesh NetCDF input and compare node,
   face, bathymetry, mask, and metadata arrays with the Fortran path.
3. Build Rust map/history NetCDF writers that reproduce dimensions, variable
   names, attribute strings, fill values, ordering, and time indexing.
4. Add schema regression tests using `ncdump -h` plus numeric reads for selected
   arrays.
5. Switch output writing to Rust using snapshots of Fortran state at output
   times.
6. Remove Fortran NetCDF output from the Rust path once map/history parity is
   proven.

Acceptance:

- Rust writes map and history files accepted by existing downstream tooling.
- Headers and selected numeric variables match Fortran outputs within
  tolerance.
- NetCDF errors are surfaced as Rust errors with filename and operation context.

## Phase 8: Rust Domain And Mesh Data Structures

Goal: replace the global Fortran data module with explicit Rust state for
non-solver model data.

Steps:

1. Design Rust structs for config, mesh, boundary forcing, wind forcing,
   observation points, output selection, and runtime state.
2. Preserve Fortran-compatible scalar widths and indexing conventions where
   data still crosses FFI.
3. Introduce conversion helpers for Fortran one-based indices and column-major
   arrays. Keep these conversions localized and heavily tested.
4. Move allocation ownership for non-solver arrays to Rust first.
5. Keep a coarse Fortran entry point that consumes Rust-prepared state and runs
   the old numerical loop.

Acceptance:

- New Rust code no longer depends on `snapwave_data` globals for migrated
  subsystems.
- The Fortran path can still run as oracle.
- Array shape and indexing conversions are covered by focused tests.

## Phase 9: Geometry, Interpolation, And Lookup Utilities

Goal: migrate supporting numerical utilities before the wave solver itself.

Scope:

- Date/time conversion utilities.
- Generic interpolation helpers from `interp.F90`.
- Boundary and observation point interpolation.
- Surrounding-point and upwind-neighbour preprocessing.
- K-d tree usage and Triangle integration decisions.

Steps:

1. Split pure utility routines from file-reading and global-state routines.
2. Port small deterministic routines first with unit tests using hand-checked
   fixtures.
3. For larger geometry routines, compare Rust outputs to Fortran for real
   testcase meshes.
4. Decide whether to keep C Triangle and a Rust k-d tree crate temporarily, or
   replace both with Rust-native libraries after parity tests exist.
5. Only after parity is proven, remove corresponding Fortran utility modules
   from the Cargo build list.

Acceptance:

- Mesh preprocessing, interpolation weights, and boundary/observation mappings
  match the Fortran oracle on representative meshes.
- Any third-party replacement is justified by tests and licensing review.

## Phase 10: Solver State Boundary

Goal: prepare for solver migration without rewriting solver physics yet.

Steps:

1. Define a Rust `ModelState` that contains all arrays and scalars needed by
   one timestep.
2. Add a coarse Fortran function that computes one timestep or one full run from
   explicit state rather than implicit module globals, if this can be done
   without destabilizing the Fortran oracle.
3. Add snapshot tests around pre-step and post-step state for small cases.
4. Make output scheduling Rust-owned so solver code only updates model state.
5. Document all numerical invariants discovered during state extraction.

Acceptance:

- Rust owns orchestration of the time loop and output scheduling.
- Fortran solver can still be called as one coarse numerical kernel.
- State snapshots make solver rewrites reviewable in small pieces.

## Phase 11: Solver Internals In Rust

Goal: migrate the numerical solver last, preserving validated behaviour.

Suggested order:

1. Small pure routines: tridiagonal solve, sorting helpers, date-independent
   math helpers.
2. Wave celerity, group velocity, and dispersion-related calculations.
3. Dissipation/source terms: breaking, bottom friction, vegetation, wind input.
4. Directional spreading and boundary spectrum construction.
5. The implicit sweep solver and convergence loop.
6. Infragravity-specific paths.
7. OpenMP parallel sections mapped to Rust parallelism only after scalar parity
   is achieved.

Steps:

1. Port one routine at a time with direct Fortran oracle tests.
2. Use deterministic fixtures before full testcase comparisons.
3. Keep floating-point operation order as close as practical until parity is
   established.
4. Define tolerances per routine and per full-output variable.
5. Benchmark only after correctness is pinned.

Acceptance:

- Full testcase outputs match the Fortran oracle within documented tolerances.
- Rust solver performance is measured against the Fortran baseline.
- The Rust path no longer links migrated Fortran solver modules.

## Phase 12: Retire Fortran From The Rust Build

Goal: complete the strangler rewrite and remove Fortran as a runtime
dependency.

Steps:

1. Remove migrated Fortran sources from `build.rs` in dependency order.
2. Retain the legacy Makefile/oracle path until final sign-off, or archive it
   behind a clearly named compatibility target.
3. Remove `src/snapwave_c_api.f90` once no Rust runtime path calls Fortran.
4. Simplify Cargo and Nix builds to Rust plus any remaining C/third-party
   dependencies.
5. Update README, AGENTS.md, license notes, and developer setup instructions.
6. Run all regression cases and produce a final parity report.

Acceptance:

- `cargo build` and `cargo test` do not compile or link SnapWave Fortran
  sources.
- The production executable is Rust-owned end to end.
- The final parity report records tolerances, tested cases, and any intentional
  deviations.

## Ongoing Workstreams

Regression coverage:

- Expand from the MWE to cases covering single-point boundaries,
  space/time-varying boundaries, wind, vegetation, Neumann boundaries,
  observation points, map-only output, history-only output, and NetCDF mesh
  variants.
- Keep generated NetCDF files out of version control unless deliberately stored
  as compact fixtures or baseline extracts.

Build and packaging:

- Keep `Makefile` and `build.rs` source ordering synchronized while Fortran is
  compiled by both systems.
- Make Nix build the same executable path used by Cargo migration tests.
- Avoid changing bundled third-party code except for build integration.

Documentation:

- Update this plan when a phase reaches acceptance or when a constraint changes.
- Record migration-specific tolerances and known behavioural quirks near the
  tests that enforce them.
