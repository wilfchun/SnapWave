# AGENTS.md

Working conventions for AI coding agents (and humans) in this repository.
Read this before making changes. `plan.md` is the authoritative migration
plan; read it before any structural change.

## What this project is

SnapWave is a fast, implicit, unstructured-grid short wave solver, originally
Fortran (with bundled C: Triangle; Fortran: kdtree2). It has been rewritten
in Rust using a **strangler fig** approach:

- subsystems moved to Rust **outside-in**, one at a time, only after
  regression tests pinned their behaviour;
- the Fortran implementation stayed available as the oracle until each
  migrated Rust path passed regression.

Current state: the migration is complete. Phases 1–11 moved the CLI, run
context, input parsing, config defaults/validation, filesystem handling,
text readers, NetCDF IO, data structures, geometry/interpolation, time
loop and solver internals to Rust. **Phase 12 is done**: the Cargo build
no longer compiles or links any Fortran/C/NetCDF source — the production
executable is 100% Rust (`src_rust/`, orchestrated by `src_rust/model.rs`).
The legacy `make` build (`src/`, `utils_lgpl/`, `third_party_open/`) is
retained unchanged as the numerical oracle (`flake.nix` `snapwave-legacy`,
or `SNAPWAVE_ORACLE`), used only by the regression suite.

## Non-negotiable rules

1. **The committed regression baselines are the numerical authority.** Do
   not "improve" scientific or numerical behaviour while changing the Rust
   solver. A change either preserves outputs (within the tolerances in
   `tests/support/compare.rs`) or intentionally alters physics — never both
   in one step.
2. **The Fortran oracle is the reference.** Any change to `src_rust/model.rs`
   or `src_rust/solver.rs` must be validated against the `make`-built oracle
   through `cargo test` before it is trusted. (The FFI facade is gone; the
   oracle runs as a separate binary, not inside the Rust process.)
3. **Keep both builds green.** The legacy `make` build (Fortran oracle) and
   `cargo build` (pure Rust) must both keep working. The Rust build needs no
   Fortran/C/NetCDF toolchain.
4. **`src/snapwave.f90` is deliberately not compiled into the Cargo binary**
   (Rust provides `main`). It is retained as the oracle's main program; keep
   it working with the `make` build.
5. **Follow the plan.md phase order** when picking what to move next
   (currently: Phases 1–10 are done; Phase 11 — solver internals — is
   next). Do not jump ahead to the solver.
6. **Third-party code is read-only in practice.** Do not modify
   `third_party_open/` or `utils_lgpl/` except for build integration fixes.
   Keep licensing (LGPL, see `LICENSE`) intact.

## Repository layout

| Path | Contents |
| --- | --- |
| `src/` | Fortran solver sources (compile order matters) — retained for the legacy `make` oracle only |
| `src_rust/` | 100%-Rust production binary sources (`[[bin]] path` in `Cargo.toml`) |
| `Makefile` | Legacy stand-alone Fortran build → `SnapWave/lnx64/bin/snapwave` (the numerical oracle) |
| `Cargo.toml` | Crate manifest; pure-Rust binary also named `snapwave` (no build script) |
| `flake.nix` | Nix dev shell, pure-Rust package, legacy-oracle package and smoke-test check |
| `third_party_open/` | Bundled Triangle (C) and kdtree2 (Fortran) — used only by the legacy `make` build |
| `utils_lgpl/` | Deltares common utils + kdtree wrapper — used only by the legacy `make` build |
| `tests/` | Cargo integration tests: `mwe.rs` smoke test; Phase-1 harness (`regression.rs`, `support/`, see `tests/README.md`); Phase-2 CLI tests (`cli.rs`) |
| `testcases/` | Validation cases; `31_linear_shoaling_refraction/run/coarse` is the MWE regression target |
| `plan.md` | Migration plan (phases, constraints, risks) |
| `doc/` | Reference manuals (physics; consult before touching numerics) |
| `SnapWave.sln`, `scripts/postbuild.ps1`, `SnapWave.zip` | Legacy Windows/VS artifacts; leave alone |

## Build, run, test

Toolchain prerequisites for the production build: a Rust toolchain
(`rustc`, `cargo`) only — **no** Fortran/C/NetCDF compiler or library is
needed. Building the legacy oracle additionally needs `gfortran`, `cc`,
`nf-config` (netCDF-Fortran); `ncdump` (netCDF tools) is optional but used
by schema checks. A Nix dev shell provides everything: `nix develop` (or
direnv).

```sh
cargo build                        # pure Rust; no Fortran/C/NetCDF
cargo test                         # smoke, regression and CLI tests
cargo run -- path/to/SnapWave.inp  # run the model
cargo run -- --help                # CLI usage (Phase 2)

make                # legacy stand-alone Fortran build (the numerical oracle)
make DEBUG=1
make clean

nix build                    # default package: pure-Rust binary
nix build .#snapwave-legacy  # legacy Makefile build (the Fortran oracle)
nix flake check              # Nix build + smoke test of the Rust binary
```

## Regression testing rules

- The MWE regression target is
  `testcases/31_linear_shoaling_refraction/run/coarse/SnapWave.inp`.
- Testcases are authored on Windows: **normalize `\` → `/` only on temp
  copies** (`cargo test` and the flake check already do this). Never commit
  modified testcase inputs.
- Since Phase 12 all file references resolve through `src_rust/paths.rs`
  against the input file's directory; there is no `chdir` and no FFI name
  conversion any more.
- Existing checks are structural: exit status, output file presence, NetCDF
  headers via `ncdump -h`. Exact floating-point reproducibility is not
  expected — tolerances live in `tests/support/compare.rs`.
- Phase-1 numeric regression lives in `tests/regression.rs`: the pure-Rust
  binary's output is compared against committed testcase reference NetCDF
  files and, when the legacy `make` build exists (or `SNAPWAVE_ORACLE` is
  set), against the live Fortran oracle. Per-variable tolerances live in
  `tests/support/compare.rs`; see `tests/README.md` before adding cases.

## Known pitfalls

- **The Fortran oracle is a separate process.** It runs the legacy `make`
  binary with its CWD in the run directory and reads `SnapWave.inp` itself.
  The Rust binary and the oracle never share a process.
- **Array ordering**: Fortran is column-major, Rust is row-major. The Rust
  solver keeps column-major flat buffers (`ee`, `w`, `ds`, `prev`, `ctheta`)
  in exactly the Fortran memory order; see the module docs in `src_rust/`.
- Build artefacts (`SnapWave/lnx64/`, `target/`) are gitignored — never
  commit them or generated NetCDF outputs.

## Code style

### Rust (`src_rust/`, `tests/`)

- Edition 2021; keep dependencies minimal (`anyhow`; add `clap` only when
  argument parsing genuinely grows).
- Document the **why**, not the what; reference the relevant plan.md phase
  in comments where the design is non-obvious.
- Use `anyhow` with `.with_context(...)` for user-facing failures; the
  wrapper exits with status 2 on its own errors.

### Fortran (`src/`)

- Retained as the numerical oracle; keep it unchanged except for genuine
  build-integration fixes, and keep `make` green.

### Definition of done for any change

1. `cargo build` and `cargo test` pass.
2. If Fortran sources changed: `make` also builds.
3. No unintended scientific behaviour change; regression coverage extended
   if outputs or a migrated subsystem were touched.
4. plan.md updated when a phase milestone is reached or constraints change.
