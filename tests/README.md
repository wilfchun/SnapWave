# SnapWave test harness

Cargo integration tests live in this directory. Everything needed for the
Phase-1 regression machinery (plan.md, Phase 1: "Test Oracle And Baselines")
is described here; adding a new case does not require reading the Fortran
build internals.

## Files

| File | Purpose |
| --- | --- |
| `mwe.rs` | Original MWE smoke test (structure only, `ncdump`-optional) |
| `regression.rs` | Phase-1 regression: runs cases and compares NetCDF output against committed baselines and the live Fortran oracle |
| `ncdf_parser.rs` | Self-tests of the NetCDF reader against committed fixtures (no model run) |
| `support/ncdf.rs` | Dependency-free reader for NetCDF classic (CDF-1/CDF-2) files |
| `support/harness.rs` | Testcase copying, separator normalization, wrapper/oracle execution |
| `support/compare.rs` | Schema pinning + numeric comparison, tolerance table |

The reader is hand-rolled on purpose: SnapWave writes classic-format files,
and numeric comparison must work anywhere `cargo test` runs (there is no
guarantee `ncdump` is on PATH) while reporting variable, index and tolerance
for every failure.

## Running

```sh
nix develop        # toolchain: gfortran, netcdf, cargo
cargo test
```

Each regression case:

1. copies the testcase to a temp dir (checked-in inputs are never modified;
   Windows `\` separators are normalized **only in the copy**);
2. runs the Cargo-built wrapper (`CARGO_BIN_EXE_snapwave`) there;
3. compares map/history output to the committed baseline under
   `testcases/<case>/output/` when that file exists;
4. runs the Fortran oracle on its own copy and compares wrapper vs oracle.

Both binaries run with `OMP_NUM_THREADS=1`: OpenMP reduction order changes
floating-point summation order, which would make numeric comparisons flaky
without any scientific meaning.

## The Fortran oracle

The oracle is the legacy `make` build, `SnapWave/lnx64/bin/snapwave`
(argument-less, reads `SnapWave.inp` from the CWD). Comparison against it
turns on automatically when that file exists; `SNAPWAVE_ORACLE=/path/to/binary`
overrides the location (e.g. pointing at a Nix-built one). Without an oracle,
committed-baseline comparisons still run; families without a committed
baseline fall back to structural checks and print a notice.

```sh
make               # enable the oracle comparison
cargo test
```

Wrapper-vs-oracle comparison is strict (identical code, compiler and
threads, so record counts must match). Committed baselines were produced by
a different platform/compiler, so differing record counts (solver iteration
counts can drift) only warn there and the common prefix is compared, and a
small documented allowlist of schema additions
(`LEGACY_BASELINE_SCHEMA_ADDITIONS` in `support/compare.rs`, currently
`point_dirspr`) tolerates variables the legacy reference files predate —
those are pinned strictly by the oracle comparison instead.

## Baselines and tolerances

Baselines are the committed legacy outputs under `testcases/<case>/output/`
(case 31: history only — the committed map file does not exist; case 32: map
only — the committed history file predates the removal of per-iteration
history output and is disabled via `use_his_baseline: false` in
`regression.rs`). They are read in place and never copied or overwritten.
Families without an active committed baseline (absent or disabled) are pinned
by the strict live-oracle comparison instead.

Numeric tolerances are defined per variable in
`support/compare.rs::tolerance_for` — pragmatic absolute + relative
tolerances for `hm0`, `tp`, directions (circular, guarded on `Hm0 > 0.01 m`),
`ee`, `depth`, dissipation terms, mesh geometry, and time coordinates.
Additional rules:

- `_FillValue` (`-999999`) positions must match exactly; values are compared
  only where both sides are non-fill.
- Directional quantities are compared on the circle (0/360 wrap) and only
  where wave height is non-negligible in both files (the guard rule above).
- `total_runtime`, `average_dt` and `station_id` are excluded from numeric
  comparison: SnapWave never writes them (see `snapwave_ncoutput.F90`).
- The `Build-Revision-Date-Netcdf-library` attribute is excluded from schema
  comparison (it embeds the netcdf library version).

Adjust tolerances only with a stated reason; they are the documented
Phase-1 tolerance definition.

## Adding a case

1. Add a `CaseSpec` constant in `regression.rs` — the field documentation on
   the struct explains each entry (testcase dir, run subdir, output file
   names, optional pinned frame counts).
2. Add a one-line `#[test]` function calling `run_regression(&YOUR_CASE)`.
3. If reference output exists, commit it under
   `testcases/<case>/output/<map|his file name>`; it is picked up by file
   name automatically.

## Debugging failures

Failing tests keep their temp directories and print the paths; rerun the
binary inside them to reproduce. Set `SNAPWAVE_TEST_KEEP=1` to always keep
temp runs. Failure reports name the variable, the decoded index (e.g.
`[time=1, stations=142]`), both values, the difference and the tolerance
that was exceeded.

## Environment variables

| Variable | Effect |
| --- | --- |
| `SNAPWAVE_ORACLE` | Path to the Fortran oracle binary (default `SnapWave/lnx64/bin/snapwave`) |
| `SNAPWAVE_TEST_KEEP` | Keep temp run directories even on success |
