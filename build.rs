//! Cargo build script for the SnapWave Rust wrapper.
//!
//! Cargo is the single orchestrator:
//!   * compiles the bundled Triangle C sources (via the `cc` crate);
//!   * compiles all Fortran sources in the dependency order used by the
//!     Makefile (plus the new src/snapwave_c_api.f90 facade);
//!   * archives the Fortran objects and emits the link directives needed
//!     for NetCDF (via `nf-config`), OpenMP and the Fortran runtime.
//!
//! The stand-alone `src/snapwave.f90` program is deliberately NOT compiled:
//! the Rust binary provides `main`.

use std::env;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

fn run_cmd(cmd: &mut Command) -> String {
    let output = cmd
        .output()
        .unwrap_or_else(|e| panic!("failed to run {:?}: {}", cmd, e));
    let stdout = String::from_utf8_lossy(&output.stdout).trim().to_string();
    let stderr = String::from_utf8_lossy(&output.stderr).trim().to_string();
    if !output.status.success() {
        panic!(
            "command failed ({:?}):\n--- stdout ---\n{}\n--- stderr ---\n{}",
            cmd, stdout, stderr
        );
    }
    stdout
}

fn tool_file(tool: &str, file: &str) -> Option<PathBuf> {
    // NB: the option and its value must be joined; gfortran rejects the
    // separated form ("-print-file-name" "<file>").
    let out = run_cmd(Command::new(tool).arg(format!("-print-file-name={file}")));
    let path = PathBuf::from(&out);
    // The compiler drivers echo the bare name when the file cannot be found.
    if path.is_absolute() && path.is_file() {
        Some(path)
    } else {
        None
    }
}

/// Compile a single Fortran source, mirroring the flags of the Makefile:
/// OMP + module dir + nf-config includes + preprocessor + free line length.
fn compile_fortran(fc: &str, base_flags: &[String], src: &Path, obj: &Path) {
    let mut cmd = Command::new(fc);
    cmd.arg("-c")
        .args(base_flags)
        .arg(src)
        .arg("-o")
        .arg(obj);
    println!("compiling fortran: {}", src.display());
    run_cmd(&mut cmd);
}

fn main() {
    let manifest_dir = PathBuf::from(env::var("CARGO_MANIFEST_DIR").expect("CARGO_MANIFEST_DIR"));
    let out_dir = PathBuf::from(env::var("OUT_DIR").expect("OUT_DIR"));
    let obj_dir = out_dir.join("obj");
    let mod_dir = out_dir.join("mod");
    fs::create_dir_all(&obj_dir).unwrap();
    fs::create_dir_all(&mod_dir).unwrap();

    // ------------------------------------------------------------------
    // Toolchain selection (mirrors the Makefile defaults)
    // ------------------------------------------------------------------
    let fc = env::var("FC").unwrap_or_else(|_| "gfortran".to_string());
    let fc_base = Path::new(&fc)
        .file_name()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_default();
    let is_intel = fc_base == "ifx" || fc_base == "ifort";
    let debug = env::var("DEBUG").map(|v| v == "1").unwrap_or(false);

    let (omp_flag, mod_flag, pp_flag, line_flag) = if is_intel {
        ("-qopenmp".to_string(), format!("-module{}", mod_dir.display()), "-fpp".to_string(), String::new())
    } else {
        (
            "-fopenmp".to_string(),
            format!("-J{}", mod_dir.display()),
            "-cpp".to_string(),
            "-ffree-line-length-none".to_string(),
        )
    };

    // NetCDF Fortran flags
    let nf_config = env::var("NF_CONFIG").unwrap_or_else(|_| "nf-config".to_string());
    let nc_fflags = run_cmd(Command::new(&nf_config).arg("--fflags"));
    let nc_flibs = run_cmd(Command::new(&nf_config).arg("--flibs"));

    let mut base_flags: Vec<String> = Vec::new();
    if debug {
        base_flags.extend(["-g".to_string(), "-O0".to_string()]);
        if !is_intel {
            base_flags.extend(["-fcheck=all".to_string(), "-fbacktrace".to_string()]);
        }
    } else {
        base_flags.push("-O2".to_string());
    }
    base_flags.push(omp_flag.clone());
    base_flags.push(mod_flag);
    base_flags.push(format!("-I{}", mod_dir.display()));
    base_flags.extend(nc_fflags.split_whitespace().map(String::from));
    base_flags.push(pp_flag);
    if !line_flag.is_empty() {
        base_flags.push(line_flag);
    }

    // ------------------------------------------------------------------
    // Fortran sources in the dependency order of the Makefile
    // ------------------------------------------------------------------
    let fortran_sources: Vec<PathBuf> = [
        "third_party_open/kdtree2/src-f90/kdtree2.f90",
        "utils_lgpl/deltares_common/src/deltares_common_modules.f90",
        "utils_lgpl/deltares_common/src/malloc.f90",
        "utils_lgpl/deltares_common/src/m_ec_triangle.f90",
        "utils_lgpl/kdtree_wrapper/src/kdtreeWrapper.f90",
        "src/snapwave_data.f90",
        "src/snapwave_date.f90",
        "src/snapwave_results.f90",
        "src/interp.F90",
        "src/snapwave_input.f90",
        "src/snapwave_windsource.f90",
        "src/snapwave_ncoutput.F90",
        "src/snapwave_domain.f90",
        "src/snapwave_boundaries.f90",
        "src/snapwave_obspoints.f90",
        "src/snapwave_solver.f90",
        // New C ABI facade; replaces src/snapwave.f90 (Rust provides main).
        "src/snapwave_c_api.f90",
    ]
    .iter()
    .map(|rel| manifest_dir.join(rel))
    .collect();

    // Triangle C sources
    let triangle_sources: Vec<PathBuf> = [
        "third_party_open/triangle/triangle.c",
        "third_party_open/triangle/tricall2.c",
    ]
    .iter()
    .map(|rel| manifest_dir.join(rel))
    .collect();

    // Watch all inputs so Cargo re-runs this script only when needed.
    println!("cargo:rerun-if-changed=build.rs");
    for src in fortran_sources.iter().chain(triangle_sources.iter()) {
        println!("cargo:rerun-if-changed={}", src.display());
    }
    for key in ["FC", "CC", "NF_CONFIG", "DEBUG"] {
        println!("cargo:rerun-if-env-changed={}", key);
    }

    // ------------------------------------------------------------------
    // Compile Fortran. Module (.mod) files make fine-grained caching
    // error-prone, so if any source changed we rebuild all of them
    // (a full build is fast and always dependency-ordered).
    // ------------------------------------------------------------------
    let objects: Vec<PathBuf> = fortran_sources
        .iter()
        .map(|src| obj_dir.join(format!("{}.o", file_stem(src))))
        .collect();

    let needs_rebuild = objects.iter().zip(&fortran_sources).any(|(obj, src)| {
        !obj.is_file()
            || fs::metadata(src)
                .and_then(|m| m.modified())
                .ok()
                .and_then(|src_t| fs::metadata(obj).and_then(|m| m.modified()).ok().map(|obj_t| src_t > obj_t))
                .unwrap_or(true)
    });

    if needs_rebuild {
        for (src, obj) in fortran_sources.iter().zip(&objects) {
            compile_fortran(&fc, &base_flags, src, obj);
        }
    } else {
        println!("cargo:warning=fortran objects up to date, skipping recompile");
    }

    // Archive the Fortran objects into one static library (recreate it so
    // members from an older source list can never linger).
    let fortran_lib = out_dir.join("libsnapwave_fortran.a");
    let _ = fs::remove_file(&fortran_lib);
    let mut ar = Command::new(env::var("AR").unwrap_or_else(|_| "ar".to_string()));
    ar.arg("crs").arg(&fortran_lib).args(&objects);
    run_cmd(&mut ar);
    println!("cargo:rustc-link-search=native={}", out_dir.display());
    // Note: emitted before the Triangle archive so that the linker resolves
    // Fortran -> triangle symbol references correctly (static lib ordering).
    println!("cargo:rustc-link-lib=static=snapwave_fortran");

    // ------------------------------------------------------------------
    // NetCDF libraries from `nf-config --flibs`
    // ------------------------------------------------------------------
    for tok in nc_flibs.split_whitespace() {
        if let Some(dir) = tok.strip_prefix("-L") {
            println!("cargo:rustc-link-search=native={}", dir);
        } else if let Some(lib) = tok.strip_prefix("-l") {
            println!("cargo:rustc-link-lib=dylib={}", lib);
        } else if (tok.ends_with(".a") || tok.ends_with(".so")) && tok.starts_with('/') {
            // Absolute path to a library file.
            let path = PathBuf::from(tok);
            if let Some(dir) = path.parent() {
                println!("cargo:rustc-link-search=native={}", dir.display());
            }
            if let Some(stem) = path.file_name().and_then(|s| s.to_str()) {
                let stem = stem
                    .trim_end_matches(".a")
                    .trim_end_matches(".so")
                    .trim_start_matches("lib");
                let kind = if tok.ends_with(".a") { "static" } else { "dylib" };
                println!("cargo:rustc-link-lib={}={}", kind, stem);
            }
        } else if tok.starts_with("-Wl,") {
            println!("cargo:rustc-link-arg={}", tok);
        } else {
            println!("cargo:warning=ignoring unrecognized nf-config flib token: {}", tok);
        }
    }

    // ------------------------------------------------------------------
    // Fortran runtime + OpenMP: rustc links with the C driver, which does
    // not automatically pull in libgfortran/libgomp, so add them explicitly
    // (with the compiler-reported search paths so this also works in Nix).
    // ------------------------------------------------------------------
    if is_intel {
        println!("cargo:rustc-link-lib=dylib=ifcore");
        println!("cargo:rustc-link-lib=dylib=iomp5");
    } else {
        if let Some(path) = tool_file(&fc, "libgfortran.so") {
            if let Some(dir) = path.parent() {
                println!("cargo:rustc-link-search=native={}", dir.display());
            }
        }
        println!("cargo:rustc-link-lib=dylib=gfortran");

        let cc = env::var("CC").unwrap_or_else(|_| "cc".to_string());
        if let Some(path) = tool_file(&cc, "libgomp.so") {
            if let Some(dir) = path.parent() {
                println!("cargo:rustc-link-search=native={}", dir.display());
            }
        }
        println!("cargo:rustc-link-lib=dylib=gomp");

        // Only needed by some gfortran configurations; link it when present.
        if let Some(path) = tool_file(&fc, "libquadmath.so") {
            if let Some(dir) = path.parent() {
                println!("cargo:rustc-link-search=native={}", dir.display());
            }
            println!("cargo:rustc-link-lib=dylib=quadmath");
        }
    }

    // ------------------------------------------------------------------
    // Triangle (C), compiled last so its archive follows the Fortran
    // archive in link order (Fortran code calls triangle routines).
    // ------------------------------------------------------------------
    let mut triangle = cc::Build::new();
    triangle
        .define("ANSI_DECLARATORS", None)
        .opt_level(2)
        .files(triangle_sources.iter().map(|p| p.to_string_lossy().to_string()));
    if debug {
        triangle.debug(true).opt_level(0);
    }
    triangle.compile("snapwave_triangle");
}

fn file_stem(path: &Path) -> String {
    path.file_stem()
        .map(|s| s.to_string_lossy().to_string())
        .expect("source file without stem")
}
