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

Status: implemented (2026-08-20). Step 1: `DomainState` in
`src_rust/state.rs` composes the already-Rust-owned parses — config
(`input`), mesh (`mesh`), boundary forcing, wind forcing, observation
points and polylines (`text_input`) — plus a `RuntimeState` holding the
scheduling scalars `run_time_loop` initialises (tstart/tstop/timestep,
intervals, next-output times, counters, output tolerance; the documented
Phase 10 seam where output scheduling becomes Rust-owned). Steps 2-3:
`src_rust/ffi_layout.rs` is the localized conversion layer — one-based
index helpers and a `ColMajor` layout type whose tests pin the Fortran
memory formula `(i1-1) + d1*(i2-1) + ...` and the two facts that make
the SnapWave handoffs copy-free: the node-major `face_nodes` buffer *is*
the column-major flattening of `face_nodes(4, no_faces)`, and the
time-major series layout *is* the column-major flattening of
`hs_bwv(nwbnd, ntwbnd)`. Steps 4-5: the `snapwave_data` globals for the
mesh, the enclosure/Neumann polylines, the observation-point coordinates
and the boundary series switched from `allocatable` to `pointer`
(`nameobs` stays allocatable — character arrays do not associate
portably and are copied from 32-byte blank-padded records); the new
coarse entry point `snapwave_run_capture_state_c` receives one
`#[repr(C)]` struct (`SnapWaveStateC` ↔ the `bind(C)` `snapwave_state_t`
mirror; Fortran-compatible widths: `real*8`/`real*4`/`integer*1`/
`integer`) and associates the globals with the Rust-owned buffers via
`c_f_pointer` — allocation ownership of those non-solver arrays is
Rust's, and Fortran no longer reads those files on this route. Reading
was split from derived computation so nothing else moved:
`initialize_snapwave_domain` gained an optional `mesh_from_rust` flag
(skipping the readers and the post-processing the Rust reader already
applies — the `zb = -posdwn*zb` flip must not run twice), and
`init_obs_points_from_state` / `init_boundary_from_state` perform the
derived tails of `read_obs_points` / `read_boundary_data` (weights via
`make_map_fm` / `find_boundary_indices` stay Fortran — Phase 9).
Wind and the `fw`/`fwig` value-or-file inputs remain Fortran-read: their
file-backed branch needs `triintfast` mesh interpolation (Phase 9); the
uniform wind data already crosses as resolved config text. The wrapper
dispatches exactly like `initialize_snapwave_domain` did
(`rust_owns_gridfile` mirrors the two-characters-after-the-first-dot
extension check, quirks included): NetCDF grids take the state route,
other formats keep `snapwave_run_capture_c`, as does the hidden
`--legacy-mesh` parity hook. Acceptance: the run path hands migrated
subsystem data from Rust state (no Fortran re-read of mesh, polylines,
obs or boundary files); the oracle is intact — `make` and the unchanged
`snapwave_run_c`/reader route still work, all `--compare-*` hooks are
untouched, and the regression suite passes with the live-oracle
comparison active (outputs identical); conversions are covered by
focused unit tests in `ffi_layout.rs`/`state.rs` plus the route-parity
tests `tests/domain_state.rs` (case 31 timeseries mode, case 32
single-point JONSWAP mode: state route vs `--legacy-mesh`, strict
comparison).

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

Status: implemented (2026-08-20). Step 1: the date/time utilities are split
into `src_rust/date.rs` (`julian_date`, `parse_date15`, `seconds_between` =
`time_difference`, `date_to_iso8601`), and `src_rust/input.rs` now delegates
to them instead of carrying private copies. `convert_fewsdate` is deliberately
not ported: it has no callers and reads an uninitialised local `trefstr`, so
there is no oracle behaviour to preserve (documented in `date.rs`). Step 2:
`src_rust/interp.rs` ports the small deterministic helpers with hand-checked
fixture tests — `binary_search`, `linear_interp`, `linear_interp_2d`, `hunt`,
`indexx`/`sort`, `ipon`, `bilin5`, `triangle_intp`, the trapezoidal family
(`trapezoidal`, `trapezoidal_cyclic`, `interp_using_trapez_rule`,
`interp_in_cyclic_function`) and the curvilinear-grid mapping routines
(`make_map`, `mkmap_step`, `grmap`, `grmap2`, `grmap_sg`) — plus
`make_map_fm`. Every `real*4`/`real*8` width and single-precision literal
(the Fortran default real kind) is preserved so the ports match bit-for-bit.
Step 3: `src_rust/geometry.rs` ports the runtime geometry — `fm_surrounding_points`
(surrounding-point ring + `real*4` plane fit via `plane_fit`/`solve_linear_system`),
`find_upwind_neighbours`/`intersect_angle`, `neuboundaries` and
`find_boundary_indices` — and the new facade hook `snapwave_geometry_dump_c`
(`src/snapwave_c_api.f90`) runs the unchanged Fortran routines
(`initialize_snapwave_domain` + `read_obs_points` + `read_boundary_data`) and
dumps the resulting globals (`kp`, `dhdx`/`dhdy`, `w360`/`prev360`/`ds360`,
`msk`, `neumannconnected`, `nmindbnd`/`neubnd`, `wobs`/`irefobs`/`nrefobs`,
`ind1/ind2/fac_bwv_cst`) in the canonical sectioned format. The wrapper's
`--compare-geometry` mode (`src_rust/geometry_compare.rs`) computes the same
geometry in Rust and compares (integers exact; reals bit-exact first, then the
1e-6/1e-9 relative tolerances already used by the Phase 3/4/6 comparisons,
because the geometry involves `hypot`/`tan`/`atan2` libm results that may
drift one ulp between runtimes). `tests/geometry.rs` runs it on the checked-in
cases (31 coarse triangles, 32 curvilinear quads, 33 circle reef, 45
haringvliet), pinning acceptance "mesh preprocessing, interpolation weights
and boundary/observation mappings match the Fortran oracle on representative
meshes". Step 4 (decision): the sample-point path `triintfast`/`findtri_kdtree`/
`dlaun` — backed by the bundled C Triangle and the Fortran `kdtree2` wrapper —
stays Fortran, because it is only reachable through the value-or-file
`fw`/`fwig`/`u10`/`u10dir` inputs when those name a file, and no checked-in
testcase exercises that branch; replacing read-only third-party code with
Rust-native Delaunay/k-d-tree crates would add runtime dependencies with no
oracle testcase to justify the licensing/parity review (recorded in
`src_rust/geometry.rs`). Step 5 is deliberately not taken: parity is proven
through the hook, but the Fortran routines remain the runtime authority and
stay in the build until a later phase wires the Rust results into the state
handoff.

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

Status: implemented (2026-08-21). Step 1: `ModelState` in `src_rust/state.rs`
owns the time-loop scheduling state — current time `t`, iteration counter
`it`, next-output times, output counters, and the `output_tol` tolerance —
with scheduling predicates (`should_output_map`, `should_output_his`) and
advancement methods (`record_map_output`, `record_his_output`,
`advance_time`) that mirror the Fortran `run_time_loop` logic exactly,
including the `do while` advancement that handles timesteps larger than the
output interval. Unit tests in `state.rs` pin the initialisation, the
`output_tol` floor, the predicates (including the `ja_save_each_iter`,
empty-filename and `nobs==0` suppression cases), the interval-advancement
`do while` loop, and the `t <= tstop` termination condition.

Step 2: six new coarse Fortran entry points in `src/snapwave_c_api.f90`:
`snapwave_init_capture_c` (Rust-state route: config load, state association,
domain/boundary/wind init, capture stream open, static output write — stops
before the time loop), `snapwave_init_legacy_capture_c` (legacy route: same
but reads files instead of consuming Rust state), `snapwave_timestep_c(t, it)`
(one solver step: `update_boundary_conditions(t)` + `compute_wave_field(t)`),
`snapwave_capture_map_c(t, ntmapout)` (capture map output via
`ncoutput_update_map`), `snapwave_capture_his_c(t, nthisout)` (capture
history output via `update_obs_points` + `ncoutput_update_his`), and
`snapwave_finalize_capture_c` (close capture stream + reset capture mode).
The existing `snapwave_run_capture_c` / `snapwave_run_capture_state_c` stay
available for the `--fortran-time-loop` parity hook and the comparison hooks;
they are unchanged. The shared init tail (`ncoutput_capture_begin` +
`ncoutput_init`) is factored into `run_init_tail`.

Step 3: route-parity tests in `tests/domain_state.rs`
(`rust_time_loop_matches_fortran_time_loop_31_coarse`,
`rust_time_loop_matches_fortran_time_loop_32_curvi_singlepoint`) run the
MWE (boundary timeseries, enclosure+Neumann, obs points) and the curvilinear
single-point-JONSWAP case through both the Rust-owned time loop (default
Phase 10 path) and the Fortran-owned time loop (`--fortran-time-loop` hidden
flag, which forces the old `snapwave_run_capture_state_c` path), then compare
the map and history NetCDF outputs with strict Phase-1 tolerances. Both
comparisons are bit-identical — the Rust scheduling logic reproduces the
Fortran loop's output schedule exactly.

Step 4: the `execute()` function in `src_rust/main.rs` now drives the time
loop from Rust. After initialisation (Phase 10 init entry points), a
`while model.is_running()` loop calls `snapwave_timestep_c` for each
iteration, checks the Rust-owned scheduling predicates, and calls
`snapwave_capture_{map,his}_c` when output is due. On any Fortran error
the capture stream is finalised (closed) before bailing out. The
`--fortran-time-loop` flag provides an escape hatch that uses the old
Fortran-owned loop for parity comparison. Output scheduling is fully
Rust-owned: Fortran no longer decides when to write output on the default
route.

Step 5: numerical invariants documented in `ModelState`'s module docs and
in the Fortran entry-point comments:
- Both output schedules start at `tstart` (first iteration outputs at t=0).
- Output fires when `t >= next_*_output - output_tol` where
  `output_tol = max(1d-6, |timestep|*1d-6)`.
- After output, `next_*_output` advances by the interval in a `do while`
  loop (handles timestep > interval).
- History output additionally requires `nobs > 0` and a non-empty
  `his_filename`.
- Map output is suppressed when `ja_save_each_iter != 0`.
- Time advances by `timestep` after output checks (not before).
- The loop runs while `t <= tstop` (the last iteration is at the largest
  `t` not exceeding `tstop`).
- The iteration counter `it` is 1-based and incremented at the top of each
  iteration, before the solver step.
- The output counters (`map_output_count`, `his_output_count`) are 1-based
  and passed to `ncoutput_update_{map,his}` as the NetCDF record index.

Acceptance verified: `cargo test` passes all 148 tests (including the 2 new
Phase 10 time-loop parity tests and the 2 existing Phase 8 route-parity
tests); `cargo build` succeeds; the legacy `make` build is unaffected (the
new entry points are additive and the existing facades are unchanged). The
regression suite confirms the Phase 10 output matches both the committed
baselines and the live Fortran oracle.

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
