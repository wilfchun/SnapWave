# AGENTS.md

Working conventions for AI coding agents (and humans) in this repository.
Read this before making changes. `plan.md` is the authoritative migration
plan; read it before any structural change.

## What this project is

SnapWave is a fast, implicit, unstructured-grid short wave solver, originally
Fortran (with bundled C: Triangle; Fortran: kdtree2). It is being rewritten in
Rust using a **strangler fig** approach:

- a Rust binary (Cargo-orchestrated) provides `main` and calls the unchanged
  Fortran solver through one coarse C ABI facade (`src/snapwave_c_api.f90`);
- subsystems move to Rust **outside-in**, one at a time, only after regression
  tests pin their current behaviour;
- the Fortran implementation stays available as the oracle until each migrated
  Rust path passes regression.

Current state: the Cargo-orchestrated wrapper/bootstrap is complete, and
plan.md Phases 1 (test oracle and baselines), 2 (Rust CLI and run
context), 3 (`SnapWave.inp` parsing in Rust, validated against the
Fortran reader through the temporary `--compare-input` hook), 4 (config
defaults, validation, and diagnostics) and 5 (filesystem and output
directory handling in Rust) are done. Phase 6 (text input readers) is
the active frontier.

## Non-negotiable rules

1. **Fortran is the numerical authority.** Do not "improve" scientific or
   numerical behaviour while migrating it. A change either preserves outputs
   (within defined tolerances) or intentionally alters physics — never both
   in one step.
2. **Keep the FFI boundary coarse.** All crossing goes through the facade in
   `src/snapwave_c_api.f90`. Do not bind individual solver, interpolation or
   NetCDF routines until the corresponding plan.md phase calls for it.
3. **Keep both builds green.** The legacy `make` build and `cargo build` must
   both keep working unless a phase explicitly retires one. When adding or
   removing Fortran sources, update the compile-order lists in **both** the
   `Makefile` and `build.rs` (module order matters).
4. **`src/snapwave.f90` is deliberately not compiled into the Cargo binary**
   (Rust provides `main`). When changing the model lifecycle, keep
   `src/snapwave.f90` and `src/snapwave_c_api.f90` in sync.
5. **Follow the plan.md phase order** when picking what to move next
   (currently: CLI args → `SnapWave.inp` parsing → config
   defaults/diagnostics → output directory handling → text input readers →
   NetCDF output → mesh/domain data structures → solver internals). Do not
   jump ahead to the solver.
6. **Third-party code is read-only in practice.** Do not modify
   `third_party_open/` or `utils_lgpl/` except for build integration fixes.
   Keep licensing (LGPL, see `LICENSE`) intact.

## Repository layout

| Path | Contents |
| --- | --- |
| `src/` | Fortran solver sources (compile order matters) |
| `src_rust/` | Rust wrapper sources (`[[bin]] path` in `Cargo.toml`) |
| `build.rs` | Cargo orchestrator: Fortran + Triangle compilation, NetCDF/OMP/runtime linking |
| `Makefile` | Legacy stand-alone Fortran build → `SnapWave/lnx64/bin/snapwave` |
| `Cargo.toml` | Crate manifest; Rust binary is also named `snapwave` |
| `flake.nix` | Nix dev shell, package and smoke-test check |
| `third_party_open/` | Bundled Triangle (C) and kdtree2 (Fortran) |
| `utils_lgpl/` | Deltares common utils + kdtree wrapper |
| `tests/` | Cargo integration tests: `mwe.rs` smoke test; Phase-1 harness (`regression.rs`, `support/`, see `tests/README.md`); Phase-2 CLI tests (`cli.rs`) |
| `testcases/` | Validation cases; `31_linear_shoaling_refraction/run/coarse` is the MWE regression target |
| `plan.md` | Migration plan (phases, constraints, risks) |
| `doc/` | Reference manuals (physics; consult before touching numerics) |
| `SnapWave.sln`, `scripts/postbuild.ps1`, `SnapWave.zip` | Legacy Windows/VS artifacts; leave alone |

## Build, run, test

Toolchain prerequisites: `gfortran`, `cc`, `nf-config` (netCDF-Fortran);
`ncdump` (netCDF tools) is optional but used by schema checks. A Nix dev
shell provides everything: `nix develop` (or direnv).

```sh
cargo build                        # Rust + all Fortran/C objects + link
DEBUG=1 cargo build                # mirrors `make DEBUG=1` (g, O0, checks)
cargo test                         # smoke, regression and CLI tests
cargo run -- path/to/SnapWave.inp  # run the model through the wrapper
cargo run -- --help                # wrapper CLI usage (Phase 2)

make                # legacy stand-alone Fortran build (must keep working)
make DEBUG=1
make clean

nix build                    # default package: Cargo-built Rust wrapper
nix build .#snapwave-legacy  # legacy Makefile build (the Fortran oracle)
nix flake check              # Nix build + smoke test of the Rust wrapper
```

Environment variables respected by `build.rs`: `FC`, `CC`, `NF_CONFIG`,
`DEBUG`. Defaults mirror the `Makefile` (`gfortran`, `cc`, `nf-config`).

## Regression testing rules

- The MWE regression target is
  `testcases/31_linear_shoaling_refraction/run/coarse/SnapWave.inp`.
- Testcases are authored on Windows: **normalize `\` → `/` only on temp
  copies** (`cargo test` and the flake check already do this). Never commit
  modified testcase inputs.
- The wrapper `chdir`s to the input file's parent before calling the
  facade, because the Fortran readers resolve sibling input and output
  file names relative to the CWD (the configuration itself has crossed
  as resolved text since Phase 4). Output *directories* are Rust-owned
  since Phase 5 (`src_rust/paths.rs` creates/validates them before the
  core runs), but the file *names* still resolve CWD-relative in
  Fortran. Preserve the chdir contract until the remaining readers and
  the NetCDF IO move to Rust (Phases 6-7). The chdir is isolated in
  `RunContext::enter_run_dir()` (`src_rust/run_context.rs`) so those
  phases can remove it cleanly.
- Existing checks are structural: exit status, output file presence, NetCDF
  headers via `ncdump -h`. **Any change touching output or migrated
  subsystems must extend `tests/mwe.rs`** (or add tests) to cover it.
  Exact floating-point reproducibility is not expected — define tolerances
  when adding numeric checks.
- Phase-1 numeric regression lives in `tests/regression.rs`: wrapper output
  is compared against committed testcase reference NetCDF files and, when the
  legacy `make` build exists (or `SNAPWAVE_ORACLE` is set), against the live
  Fortran oracle. Per-variable tolerances live in
  `tests/support/compare.rs`; see `tests/README.md` before adding cases.
- Each newly migrated subsystem keeps the Fortran path callable as the
  oracle and is validated against it before the Fortran side is retired.

## Known pitfalls

- **Fortran module global state**: one model run per process. Repeated calls
  into the facade from the same process are unsafe until reset/finalize
  logic exists. Don't "fix" this casually.
- **Fortran `stop` kills the whole Rust process.** When a `stop` becomes
  reachable through the facade, prefer converting it to a status return —
  but only on paths the facade actually hits.
- **Array ordering**: Fortran is column-major, Rust is row-major. Avoid
  passing large multidimensional arrays across FFI (a plan.md constraint).
- **Static link order matters**: the Fortran archive must precede the
  Triangle archive (`build.rs` handles this; don't reorder).
- **`nf-config --flibs` tokens** need careful translation in `build.rs`
  (`-L`/`-l`/absolute paths/`-Wl,`).
- **OpenMP + libgfortran linking** through rustc is toolchain-specific;
  `build.rs` queries the compiler for paths (works under Nix too).
- Build artefacts (`SnapWave/lnx64/`, `target/`) are gitignored — never
  commit them or generated NetCDF outputs.

## Code style

### Rust (`src_rust/`, `tests/`, `build.rs`)

- Edition 2021; keep dependencies minimal (`anyhow`; add `clap` only when
  argument parsing genuinely grows).
- Document the **why**, not the what; reference the relevant plan.md phase
  in comments where the design is non-obvious.
- Use `anyhow` with `.with_context(...)` for user-facing failures; the
  wrapper exits with status 2 on its own errors and passes through non-zero
  Fortran status codes unchanged.
- FFI signatures live in `src_rust/main.rs` and must match
  `src/snapwave_c_api.f90` exactly (`bind(C)`; explicit length, no reliance
  on NUL termination).

### Fortran (`src/`)

- Match the existing style of the file you are editing.
- New boundary/facade code goes in `src/snapwave_c_api.f90`; keep solver
  internals untouched except where a plan.md phase explicitly opens them.
- Prefer status returns over `stop`/`abort` on newly exposed paths.
- Any new Fortran source must be added to the `Makefile` and `build.rs`
  lists in correct dependency order, and must not break `make`.

### Definition of done for any change

1. `cargo build` and `cargo test` pass.
2. If Fortran sources changed: `make` also builds.
3. No unintended scientific behaviour change; regression coverage extended
   if outputs or a migrated subsystem were touched.
4. plan.md updated when a phase milestone is reached or constraints change.
