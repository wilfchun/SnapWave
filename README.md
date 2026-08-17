# SnapWave
 Fast, implicit, unstructured-grid short wave solver
 
## Compiling and Installation

Executables are installed under `./snapwave/lnx64/bin` or `.\snapwave\x64\<configuration>\bin`.  

### Windows
Compiling is supported under VS with the Intel ifx compiler. Set the ONEAPI_ROOT environment variable before building `snapwave.sln`. The necessary dll's are collected as a post-build step based on the value of this variable.

### Linux
A `makefile` is provided with the repository. Compiling was tested with `ifort`, `ifx` and `gfortran`. Make sure your system netcdf libraries are compatible with the compiler used.
 
For example, to build under Ubuntu:
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

## Rust wrapper (Cargo)

A minimal Rust binary wraps the unchanged Fortran solver through a coarse C ABI
(`src/snapwave_c_api.f90`) and Cargo orchestrates the full build: Rust, Fortran
objects, the bundled Triangle C code, NetCDF flags and final linking.

``` bash
# build wrapper + all Fortran/C objects (requires gfortran, cc, nf-config)
cargo build

# smoke test: runs the coarse linear shoaling/refraction testcase and
# verifies the generated NetCDF map/history outputs (uses ncdump if present)
cargo test

# run the model through the Rust wrapper
cargo run -- path/to/SnapWave.inp
```

Notes:
- The Fortran input reader resolves all paths relative to the working
  directory, so the wrapper changes to the input file's directory and passes
  only the file name across the FFI boundary.
- The checked-in `testcases/31_linear_shoaling_refraction/run/coarse/SnapWave.inp`
  uses Windows path separators (`..\`); on Linux, copy the testcase and
  normalize them first (`sed -i 's|\\|/|g' SnapWave.inp`). `cargo test` does
  this automatically on a temp copy.
- `DEBUG=1 cargo build` mirrors `make DEBUG=1` (g, O0, checks, backtrace).
- The `make`-based stand-alone Fortran build remains available unchanged.
