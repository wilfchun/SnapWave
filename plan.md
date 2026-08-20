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
- Resolved in Phase 2: the Nix default package and the smoke-test check now
  build and run the Cargo wrapper; the legacy Makefile executable stays
  available as `packages.snapwave-legacy` (the Fortran oracle).
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

Status: implemented (2026-08-18). The CLI lives in `src_rust/cli.rs`
(hand-rolled, no new dependencies): `--help`/`-h`, `--version`/`-V`,
`--verbose`, `--` separator and exactly one positional input path, parsed
from OS strings so invalid UTF-8 is a clean error rather than a panic.
Status codes: 0 for success (including help/version), 2 for
wrapper-detected errors, Fortran statuses passed through unchanged. The
process-level run state lives in `src_rust/run_context.rs` as `RunContext`
(input path as given plus resolved, run directory, output-directory
expectation — `None` until Phase 5 owns output policy —, executable
metadata, logging preferences). The legacy `chdir` contract is isolated in
`RunContext::enter_run_dir()` and the FFI file-name conversion in
`RunContext::input_file_name_cstring()` (raw bytes on Unix, NUL-rejecting)
so Phase 5 can remove them cleanly. Integration tests in `tests/cli.rs`
cover the flags, usage errors, missing/directory/invalid-UTF-8 input and a
relative-path run of the coarse testcase (embedded NUL cannot cross
`execve`, so it is unit-tested at the conversion boundary).
`flake.nix` now builds and smoke-tests the Cargo wrapper as the default
package (`packages.snapwave`, via `rustPlatform.buildRustPackage` with the
committed `Cargo.lock`); the legacy `make` build remains available as
`packages.snapwave-legacy` for use as the `SNAPWAVE_ORACLE` oracle.

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

Status: implemented (2026-08-19). The parser lives in `src_rust/input.rs`
with the grammar documented in its module docs (step 1) and configuration
as typed structs grouped by concern (step 3). Legacy semantics are
preserved verbatim, including the quirks: per-key first-match-wins,
case-sensitive exact key matching (leading blanks break matching), records
without `=` and unknown keywords silently ignored, Fortran
list-directed value semantics (including `D` exponents and first-token
extraction), the 256-character line buffer, per-field character
truncation widths (15/232/256), `map_interval`/`his_interval` defaulting
to the parsed `timestep` (with the non-positive-interval `stop 1` checks
converted to Rust errors), the string-equality `u10 == '0.0'` wind
switch, `mmax`/`nmax` plus two dummy rows, and `tstart`/`tstop` computed
with the same Fliegel & Van Flandern Julian-day arithmetic. The
temporary comparison hook (step 5) is `snapwave_read_input_dump_c` in
`src/snapwave_c_api.f90`, exposed through the wrapper's `--compare-input`
mode (`src_rust/input_compare.rs`): it runs the legacy Fortran reader
with the Rust-selected file name and dumps every global it sets
(reals as exact IEEE bit patterns); `tests/input_parse.rs` asserts
agreement for every checked-in testcase input plus quirk and invalid-input
cases. Step 6 is done via the optional `input_file` argument of
`read_snapwave_input()` (the stand-alone `make` build keeps the
name-probing behaviour), and the wrapper now parses and validates the
input in Rust before calling the facade, so invalid input is a wrapper
error (exit 2) with the Fortran core never invoked. Defaults, broader
validation and the remaining input-reader `stop` paths moved to Rust in
Phase 4.

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

Status: implemented (2026-08-19). Defaults live in `src_rust/input.rs::defaults`
as named constants (step 1). The Rust parser resolves the entire configuration
(defaults, validation, post-processing) and passes it to the Fortran facade as
canonical `key=value` text through `snapwave_run_c`; Fortran stores it via the
new `read_resolved_input` subroutine in `src/snapwave_input.f90` and no longer
reads `SnapWave.inp` or decides defaults on the Rust route (acceptance:
"Fortran no longer decides main config defaults"). Validation (step 2) covers
the non-positive output-interval checks (already in Rust since Phase 3) plus
the input-file existence check in `RunContext`; all invalid-configuration
errors are wrapper exit 2 before the Fortran core runs (acceptance: "Invalid
configuration fails before domain initialization"). Structured diagnostics
live in `src_rust/diagnostics.rs` (step 3): the wind-switch informational
message is preserved, and verbose mode reports the parse/validation status.
Validation tests in `tests/input_validate.rs` cover bad intervals, optional
output settings, missing input file, and the resolved-config handoff (step 4).
The five `stop 1` paths in `read_snapwave_input` are converted to status
returns via an optional `status` argument (step 5); the legacy stand-alone
program omits the argument and keeps the original `stop` behaviour. The
`--compare-input` mode now verifies BOTH the legacy-reader equivalence
(Phase 3) and the resolved-config handoff (Phase 4) through the new
`snapwave_load_config_dump_c` hook, so "Legacy and Rust paths agree on
accepted testcase configuration" is pinned by every `cargo test` run.

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

Status: implemented (2026-08-20). Path resolution is centralized in
`src_rust/paths.rs` (step 1): every file reference of `SnapWave.inp`
resolves against the input file's directory through one module
(`RunPaths::resolve`), ready for the Phase 6 readers to consume as
Rust-owned paths. Output paths are `PathBuf`s (step 2) and the required
output directories are created — or, when the path already is a
directory or the parent cannot be created, rejected with a clean
wrapper error — in Rust before the Fortran core runs (step 3). This is
a deliberate wrapper improvement over the legacy binary, which fails
inside the NetCDF writer when the directory is missing; it changes no
scientific behaviour. Legacy semantics are preserved verbatim (step 4):
empty strings mean "not configured", `bndfile`/`encfile`/`neumannfile`/
`obsfile` disable on any value whose first four characters are `none`
(mirroring the `name(1:4) /= 'none'` guards), map/history output
disables only on the empty string (`none` would be a literal file
name), relative paths (including `..`) join verbatim, and Windows `\`
separators are deliberately not normalized (that stays a temp-copy
concern of the test harness). Tests (step 5) live in
`tests/output_dirs.rs` (nested run directories, missing output parents
both created and uncreatable, disabled map output, disabled history
output, Windows-style separators in temp copies) plus unit tests in
`src_rust/paths.rs`; `tests/mwe.rs` and the flake smoke test no longer
pre-create `output/`, pinning wrapper ownership (the Phase-1 harness
still pre-creates it for the legacy oracle run, which cannot create
directories itself). Acceptance: Rust owns run/output directory policy;
the `chdir` removal is deliberately deferred — the downstream readers
accept explicit paths only from Phase 6 on, so the contract stays
isolated in `RunContext::enter_run_dir()`; all existing testcases run
unchanged (regression suite green, including the live-oracle
comparison).

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

Status: implemented (2026-08-20). The parsers live in `src_rust/text_input.rs`
with the Fortran list-directed semantics preserved verbatim (blank-line
skipping, blank/comma/slash separators, quoted literals skipped by numeric
reads — so a dangling `'` after the last number, as in one checked-in obs
file, still parses — `D` exponents, `character*32` name truncation, the
`station_%04d` default name, the per-file `t_bwv` overwrite of the boundary
time-series reader, and the `wd`/`ds` degree→radian conversions using the
same `deg2rad = 4*atan(1)/180d0` recipe as `snapwave_data.f90`). Rust-owned
structs (step 3) group the data by reader: `ObsPoints`, `JonswapSeries`,
`BoundarySeries`, `WindInput` (`Uniform`/`List`), `Polyline` (shared by the
enclosure and Neumann readers), plus `AsciiMesh` and `SamplePoints` for the
plain-text mesh/sample readers. The comparison hook (step 2) is
`snapwave_text_dump_c` in `src/snapwave_c_api.f90`: it loads the resolved
config, runs `initialize_snapwave_domain`, `read_obs_points`,
`read_boundary_data` and `read_wind_data`, then dumps the resulting globals
(reals as IEEE-754 bit patterns) for comparison by
`src_rust/text_compare.rs` (reals compared within the 1e-6/1e-9 relative
tolerances already used by the Phase 3/4 comparison; integers/names exact).
The wrapper's `--compare-text` mode (`src_rust/main.rs`) drives that
comparison; `tests/text_input.rs` runs it on the checked-in cases (31 coarse:
obs+boundary time series+enclosure+Neumann; 32/33: single-point JONSWAP +
enclosure; 45 haringvliet: quoted obs names) and asserts agreement. The run
path (`execute` in `src_rust/main.rs`) now parses and validates every
auxiliary text input in Rust before the model runs (acceptance: data is
Rust-owned before the timestep loop), reporting a verbose summary through
`diagnostics::report_text_input_diagnostics`. The Fortran readers remain the
runtime authority (acceptance: reader code remains available) — the
coarse-boundary handoff of the parsed data and the bypass of the Fortran
readers are deferred to Phase 8 (data structures) / Phase 9 (interpolation),
because the readers also compute `make_map_fm`/`find_boundary_indices`
interpolation weights that are out of scope for text parsing. Family 6 has no
checked-in testcase (every mesh is NetCDF; `fw60.xyz` sample interpolation is
Phase 9), so `AsciiMesh`/`SamplePoints` are unit-tested only and not part of
the oracle comparison yet.

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

Status: implemented (2026-08-20). The NetCDF strategy (step 1) is a
hand-rolled classic-format (CDF-1) writer + reader in `src_rust/netcdf.rs`,
chosen over bindings to the installed library to keep the build
dependency-free under Nix and identical on every toolchain — the same
reason the read side was hand-rolled in Phase 1 (`tests/support/ncdf.rs`);
the rationale is documented in the module docs. Step 2 ports `nc_read_net`
into `src_rust/mesh.rs` (`read_ugrid_netcdf`: old/new dimension-name
detection, coordinate widening to real*8, `zb = -posdwn*zb`, the
`-1 -> 0 -> -999` fourth-node chain, and the `|y(1)|>90` sferic fix) and
pins it against the unchanged Fortran reader through the temporary
`snapwave_mesh_dump_c` hook, driven by the wrapper's `--compare-mesh` mode
and `tests/netcdf_io.rs`. Steps 3–6 switch output writing to Rust:
`snapwave_ncoutput` gains a capture mode (off by default, so the legacy
`make` oracle and the retained `snapwave_run_c` facade path are unaffected)
that streams the exact buffers `nf90_put_var` would receive into a
little-endian stream file; the wrapper runs the model through the new
`snapwave_run_capture_c` facade, reads the stream (`src_rust/capture.rs`)
and writes map/history with `src_rust/output.rs`, reproducing the Fortran
schema (dimension/variable names, attribute strings, fill values, ordering,
time indexing) verbatim — including the per-solver-iteration map output of
`ja_save_each_iter`, because the capture intercepts `ncoutput_update_map`
itself. Values Fortran never writes (`mesh2d`, `crs`, `station_id`,
`point_zb`, `total_runtime`, `average_dt`) are emitted with the NetCDF
default fill values. Acceptance is pinned by the existing regression suite:
the wrapper's map/history files are now Rust-written and
`tests/regression.rs` compares them against the committed baselines and the
live Fortran oracle, while `tests/mwe.rs` and the flake smoke test keep the
structural (`ncdump -h`) checks. NetCDF errors surface as Rust errors with
filename and operation context (`anyhow` context in `netcdf.rs`/`output.rs`).
The Fortran mesh readers remain the runtime authority inside a model run —
handing Rust-owned mesh data to Fortran is the Phase 8 data-structure
handoff.

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
