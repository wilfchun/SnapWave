# SnapWave
 Fast, implicit, unstructured-grid short wave solver

## Compiling and Installation

Executables are installed under `./snapwave/lnx64/bin` or `.\snapwave\x64\<configuration>\bin`.

### Windows
Compiling is supported under VS with the Intel ifx compiler. Set the ONEAPI_ROOT environment variable before building `snapwave.sln`. The necessary dll's are collected as a post-build step based on the value of this variable.

### Linux — pure-Rust build (default)

The production binary is 100% Rust (`src_rust/`) and needs only a Rust
toolchain — no Fortran, C or NetCDF toolchain:

``` bash
cargo build
cargo test
cargo run -- path/to/SnapWave.inp
```

### Linux — legacy Fortran oracle

The original Fortran implementation is retained unchanged as the numerical
oracle, built with the `makefile` (tested with `ifort`, `ifx` and `gfortran`).
Make sure your system netcdf libraries are compatible with the compiler used.

``` bash
sudo add-apt-repository universe
sudo apt update
sudo apt install libnetcdff-dev gfortran gcc
make clean
make
```

To automatically remove module and object files, run with
``` bash
make clean
make STRIP_BUILD=1
```

## Notes

- The checked-in `testcases/31_linear_shoaling_refraction/run/coarse/SnapWave.inp`
  uses Windows path separators (`..\`); on Linux, copy the testcase and
  normalize them first (`sed -i 's|\\|/|g' SnapWave.inp`). `cargo test` does
  this automatically on a temp copy.
- All input and output file references in `SnapWave.inp` resolve against the
  input file's directory (the Rust model runs there; there is no `chdir`
  and no FFI).
- The `make`-based stand-alone Fortran build remains available as the
  numerical oracle (`SNAPWAVE_ORACLE`); the Rust regression suite compares
  against it when it exists.
