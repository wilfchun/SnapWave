//! Run-relative path resolution and output-directory policy
//! (plan.md, Phase 5).
//!
//! # What moves to Rust in this phase
//!
//! The Fortran core resolves every file name from `SnapWave.inp`
//! relative to the process working directory; the wrapper `chdir`s to
//! the input file's directory first (`RunContext::enter_run_dir`), so
//! effectively everything resolves relative to that directory. This
//! module makes the same resolution explicit on the Rust side:
//!
//! * every file reference of the configuration resolves against the
//!   input file's directory through this one module (Phase 5, step 1),
//!   so the Phase 6+ readers can consume Rust-owned paths instead of
//!   re-deriving the legacy rules;
//! * the map/history output files are `PathBuf`s, not fixed-width
//!   Fortran strings (step 2), and their directories are created or
//!   validated in Rust before the Fortran core runs (step 3).
//!
//! Creating missing output directories is a deliberate wrapper
//! behaviour improvement over the legacy binary, which fails inside
//! the NetCDF writer when the directory is missing (test harnesses
//! used to pre-create `output/` for exactly that reason). It changes
//! no scientific behaviour.
//!
//! The downstream readers still open the *raw* file names relative to
//! the CWD, so until they migrate (Phases 6-7) the resolved paths are
//! authoritative Rust-side knowledge: they drive the output-directory
//! policy and `--verbose` diagnostics now, and will be handed to the
//! Rust readers later. The `chdir` contract itself therefore stays
//! (isolated in `RunContext::enter_run_dir`).
//!
//! # Legacy semantics preserved verbatim (step 4)
//!
//! * **Empty string** — a file reference of `''` is "not configured".
//! * **`none`** — `bndfile`, `encfile`, `neumannfile` and `obsfile`
//!   disable on *any* value whose first four characters are `none`,
//!   mirroring the `if (<name>(1:4) /= 'none')` guards in
//!   `src/snapwave_boundaries.f90`, `src/snapwave_domain.f90` and
//!   `src/snapwave_obspoints.f90` (so a value like `nonevil.txt` is
//!   also "disabled" to the legacy readers, and stays so here).
//! * **Output files** — only the *empty string* disables map/history
//!   output (`if (map_filename /= '')` in `src/snapwave_c_api.f90`);
//!   a value of `none` would be a literal output file name.
//! * **Relative paths** — including `..` segments — join verbatim
//!   against the input directory. Windows `\` separators are NOT
//!   normalized here: on Linux a `\` is an ordinary file-name
//!   character, so an un-normalized testcase resolves to a differently
//!   named file exactly like the Fortran readers do (normalization is
//!   a temp-copy concern of the test harness, see AGENTS.md).
//! * **Value-or-file strings** — `fw`, `fwig`, `u10`, `u10dir` hold
//!   either a uniform value or a file name; the Fortran readers decide
//!   with `inquire(file=..., exist=...)` relative to the CWD (e.g.
//!   `src/snapwave_domain.f90` for `fw`, `src/snapwave_boundaries.f90`
//!   for `u10`). The resolved candidates below only record where such
//!   a file *would* live; the disambiguation stays Fortran-side until
//!   the readers migrate (plan.md Phase 6).

use std::path::{Path, PathBuf};

use anyhow::{bail, Context, Result};

use crate::input::SnapWaveInput;

/// Every file reference of `SnapWave.inp`, resolved against the input
/// file's directory. `None` = not configured (per-field legacy
/// empty/`none` semantics, see the module docs).
#[derive(Debug)]
// The input-side fields have no consumer yet: they exist so the
// Phase 6+ readers take Rust-owned paths instead of re-deriving the
// legacy resolution rules (plan.md Phase 5, step 1). Only `outputs`
// drives wrapper behaviour today.
#[allow(dead_code)]
pub struct RunPaths {
    // ---- grid / domain
    /// Mesh/grid file (`.grd`, UGRID NetCDF, ...); read unconditionally
    /// by `initialize_snapwave_domain` (default name `.txt`).
    pub gridfile: Option<PathBuf>,
    pub depfile: Option<PathBuf>,
    pub mskfile: Option<PathBuf>,
    pub indfile: Option<PathBuf>,
    pub upwfile: Option<PathBuf>,
    // ---- boundary forcing
    pub jonswapfile: Option<PathBuf>,
    pub bndfile: Option<PathBuf>,
    pub encfile: Option<PathBuf>,
    pub neumannfile: Option<PathBuf>,
    pub bhsfile: Option<PathBuf>,
    pub btpfile: Option<PathBuf>,
    pub bwdfile: Option<PathBuf>,
    pub bdsfile: Option<PathBuf>,
    pub bzsfile: Option<PathBuf>,
    pub obsfile: Option<PathBuf>,
    // ---- wind
    pub windlistfile: Option<PathBuf>,
    /// Where a `u10` *file* would live (value-or-file, see module docs).
    pub u10_candidate: Option<PathBuf>,
    /// Where a `u10dir` file would live (value-or-file).
    pub u10dir_candidate: Option<PathBuf>,
    // ---- vegetation
    pub vegmapfile: Option<PathBuf>,
    // ---- physics
    /// Where a `fw` file would live (value-or-file).
    pub fw_candidate: Option<PathBuf>,
    /// Where a `fwig` file would live (value-or-file).
    pub fw_ig_candidate: Option<PathBuf>,
    // ---- output (Phase 5, step 2: PathBufs, not Fortran strings)
    pub outputs: OutputPaths,
}

/// Resolved map/history output paths. `None` = that output family is
/// disabled (the legacy empty-string rule; `none` is *not* special
/// here and would be a literal file name).
#[derive(Debug)]
pub struct OutputPaths {
    pub map: Option<PathBuf>,
    pub his: Option<PathBuf>,
}

impl OutputPaths {
    /// Create or validate the required output directories
    /// (plan.md Phase 5, step 3). For every enabled output family:
    ///
    /// * reject an output path that already exists as a directory
    ///   (the NetCDF writer would fail on it) with a clean wrapper
    ///   error;
    /// * create missing parent directories (including `..`-relative
    ///   ones, joined verbatim);
    /// * reject a parent that exists but is not a directory.
    ///
    /// Returns the directories that were created (for reporting); an
    /// already-present directory is not created (or reported) again,
    /// so a shared map/his directory is handled once.
    pub fn prepare(&self) -> Result<Vec<PathBuf>> {
        let mut created = Vec::new();
        for (kind, path) in [("map", &self.map), ("his", &self.his)] {
            let Some(path) = path else { continue };
            if path.is_dir() {
                bail!("{kind} output path is an existing directory: {}", path.display());
            }
            // Resolved paths always hang off the (absolute) run
            // directory, so a parent exists; guard anyway rather than
            // unwrap on a future call site with a bare relative path.
            let Some(parent) = path.parent() else {
                bail!("{kind} output path has no parent directory: {}", path.display());
            };
            if parent.is_dir() {
                continue;
            }
            if parent.exists() {
                bail!("{kind} output parent exists but is not a directory: {}", parent.display());
            }
            std::fs::create_dir_all(parent)
                .with_context(|| format!("cannot create {kind} output directory {}", parent.display()))?;
            created.push(parent.to_path_buf());
        }
        Ok(created)
    }
}

impl RunPaths {
    /// Resolve every file reference of `cfg` against `run_dir` (the
    /// input file's directory — the same anchor the Fortran readers
    /// effectively use, see the module docs).
    pub fn resolve(run_dir: &Path, cfg: &SnapWaveInput) -> Self {
        RunPaths {
            gridfile: file_ref(run_dir, &cfg.grid.gridfile),
            depfile: file_ref(run_dir, &cfg.grid.depfile),
            mskfile: file_ref(run_dir, &cfg.grid.mskfile),
            indfile: file_ref(run_dir, &cfg.grid.indfile),
            upwfile: file_ref(run_dir, &cfg.grid.upwfile),
            jonswapfile: file_ref(run_dir, &cfg.boundary.jonswapfile),
            bndfile: none4_ref(run_dir, &cfg.boundary.bndfile),
            encfile: none4_ref(run_dir, &cfg.boundary.encfile),
            neumannfile: none4_ref(run_dir, &cfg.boundary.neumannfile),
            bhsfile: file_ref(run_dir, &cfg.boundary.bhsfile),
            btpfile: file_ref(run_dir, &cfg.boundary.btpfile),
            bwdfile: file_ref(run_dir, &cfg.boundary.bwdfile),
            bdsfile: file_ref(run_dir, &cfg.boundary.bdsfile),
            bzsfile: file_ref(run_dir, &cfg.boundary.bzsfile),
            obsfile: none4_ref(run_dir, &cfg.boundary.obsfile),
            windlistfile: file_ref(run_dir, &cfg.wind.windlistfile),
            u10_candidate: file_ref(run_dir, &cfg.wind.u10),
            u10dir_candidate: file_ref(run_dir, &cfg.wind.u10dir),
            vegmapfile: file_ref(run_dir, &cfg.vegetation.vegmapfile),
            fw_candidate: file_ref(run_dir, &cfg.physics.fw),
            fw_ig_candidate: file_ref(run_dir, &cfg.physics.fw_ig),
            outputs: OutputPaths {
                map: file_ref(run_dir, &cfg.output.map_file),
                his: file_ref(run_dir, &cfg.output.his_file),
            },
        }
    }
}

/// Plain file reference: the empty string means "not configured" (all
/// readers guard with `/= ''`).
fn file_ref(run_dir: &Path, name: &str) -> Option<PathBuf> {
    (!name.is_empty()).then(|| run_dir.join(name))
}

/// File reference whose reader disables on a leading `none`: the
/// Fortran guards compare `name(1:4) == 'none'`, so any value starting
/// with those four characters is disabled (values shorter than or
/// different in the first four characters stay configured).
fn none4_ref(run_dir: &Path, name: &str) -> Option<PathBuf> {
    (!name.is_empty() && !name.starts_with("none")).then(|| run_dir.join(name))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::input::parse_str;

    // A scratch dir per test (create_dir_all/is_dir probes need a real
    // filesystem; paths never leave the temp dir).
    fn scratch(name: &str) -> PathBuf {
        let dir = std::env::temp_dir()
            .join(format!("snapwave_paths_{}_{}", name, std::process::id()))
            .join("rundir");
        let _ = std::fs::remove_dir_all(dir.parent().unwrap());
        std::fs::create_dir_all(&dir).expect("create scratch run dir");
        dir
    }

    #[test]
    fn empty_and_plain_references() {
        let run_dir = Path::new("/run");
        assert_eq!(file_ref(run_dir, ""), None);
        assert_eq!(file_ref(run_dir, "mesh.grd"), Some(PathBuf::from("/run/mesh.grd")));
        // Relative paths with .. join verbatim, no normalization.
        assert_eq!(file_ref(run_dir, "../../output/m.nc"), Some(PathBuf::from("/run/../../output/m.nc")));
        // Windows separators are ordinary characters on this side of
        // the FFI (normalization is a test-copy concern).
        assert_eq!(file_ref(run_dir, r"..\m.nc"), Some(PathBuf::from(r"/run/..\m.nc")));
    }

    #[test]
    fn none4_references_mirror_the_fortran_guards() {
        let run_dir = Path::new("/run");
        assert_eq!(none4_ref(run_dir, "none"), None, "'none' disables");
        assert_eq!(none4_ref(run_dir, "nonevil.txt"), None, "first four chars 'none' disable (legacy quirk)");
        assert_eq!(none4_ref(run_dir, ""), None, "empty is not configured");
        assert_eq!(none4_ref(run_dir, "no"), Some(PathBuf::from("/run/no")), "shorter values stay configured");
        assert_eq!(none4_ref(run_dir, "non"), Some(PathBuf::from("/run/non")));
        assert_eq!(none4_ref(run_dir, "bnd.txt"), Some(PathBuf::from("/run/bnd.txt")));
    }

    #[test]
    fn output_references_disable_only_on_the_empty_string() {
        let run_dir = Path::new("/run");
        let cfg = parse_str("map_file = none\n").unwrap();
        let paths = RunPaths::resolve(run_dir, &cfg);
        assert_eq!(paths.outputs.map, Some(PathBuf::from("/run/none")), "'none' is a literal output file name");
        assert_eq!(paths.outputs.his, None, "empty his_file disables history output");
    }

    #[test]
    fn resolve_maps_every_file_reference_family() {
        let cfg = parse_str(
            "gridfile = g.grd\ndepfile = d.dep\nbndfile = b.txt\nencfile = none\nobsfile = none\n\
             neumannfile = n.txt\njonswapfile = j.txt\nwindlistfile = w.txt\nu10 = 10.5\nfw = fw.xyz\n\
             vegmapfile = v.nc\nmap_file = ..\\..\\output\\m.nc\n",
        )
        .unwrap();
        let paths = RunPaths::resolve(Path::new("/run"), &cfg);

        assert_eq!(paths.gridfile, Some(PathBuf::from("/run/g.grd")));
        assert_eq!(paths.depfile, Some(PathBuf::from("/run/d.dep")));
        assert_eq!(paths.bndfile, Some(PathBuf::from("/run/b.txt")));
        assert_eq!(paths.encfile, None);
        assert_eq!(paths.obsfile, None);
        assert_eq!(paths.neumannfile, Some(PathBuf::from("/run/n.txt")));
        assert_eq!(paths.jonswapfile, Some(PathBuf::from("/run/j.txt")));
        assert_eq!(paths.windlistfile, Some(PathBuf::from("/run/w.txt")));
        // Value-or-file candidates resolve even when they hold values;
        // the disambiguation stays Fortran-side (see module docs).
        assert_eq!(paths.u10_candidate, Some(PathBuf::from("/run/10.5")));
        assert_eq!(paths.fw_candidate, Some(PathBuf::from("/run/fw.xyz")));
        assert_eq!(paths.vegmapfile, Some(PathBuf::from("/run/v.nc")));
        // The Windows-authored map path joins verbatim (not normalized).
        assert_eq!(paths.outputs.map, Some(PathBuf::from(r"/run/..\..\output\m.nc")));
    }

    #[test]
    fn prepare_creates_missing_parents_once_and_reports_them() {
        let run_dir = scratch("prepare_creates");
        let out = OutputPaths {
            map: Some(run_dir.join("../made/deeper/m.nc")),
            his: Some(run_dir.join("../made/deeper/h.nc")),
        };
        let created = out.prepare().expect("missing parents must be created");
        // One shared parent, reported once, actually on disk.
        let expected = run_dir.join("../made/deeper");
        assert!(expected.is_dir(), "parent must exist after prepare");
        assert_eq!(created, vec![expected]);

        // Second call is a no-op (existing directories are not re-created).
        assert!(out.prepare().expect("prepare on existing dirs must succeed").is_empty());
        let _ = std::fs::remove_dir_all(run_dir.parent().unwrap());
    }

    #[test]
    fn prepare_accepts_existing_parents() {
        let run_dir = scratch("prepare_existing");
        let out = OutputPaths { map: Some(run_dir.join("m.nc")), his: None };
        assert!(out.prepare().expect("run dir itself must suffice").is_empty());
        let _ = std::fs::remove_dir_all(run_dir.parent().unwrap());
    }

    #[test]
    fn prepare_rejects_an_output_path_that_is_a_directory() {
        let run_dir = scratch("prepare_dirout");
        let blocker = run_dir.join("iamdir");
        std::fs::create_dir(&blocker).expect("create blocker dir");
        let out = OutputPaths { map: Some(blocker), his: None };
        let err = out.prepare().expect_err("directory as output file must fail");
        let msg = format!("{err:#}");
        assert!(msg.contains("map") && msg.contains("directory"), "error was: {msg}");
        let _ = std::fs::remove_dir_all(run_dir.parent().unwrap());
    }

    #[test]
    fn prepare_rejects_a_parent_blocked_by_a_file() {
        let run_dir = scratch("prepare_blocked");
        let blocker = run_dir.join("blocker");
        std::fs::write(&blocker, b"").expect("create blocker file");
        let out = OutputPaths { map: Some(run_dir.join("blocker/sub/m.nc")), his: None };
        let err = out.prepare().expect_err("parent under a file must fail");
        assert!(format!("{err:#}").contains("map"), "error must name the output family");
        let _ = std::fs::remove_dir_all(run_dir.parent().unwrap());
    }

    #[test]
    fn prepare_skips_disabled_outputs() {
        let out = OutputPaths { map: None, his: None };
        assert!(out.prepare().expect("nothing to do").is_empty());
    }
}
