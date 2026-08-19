//! Rust parser for the SnapWave model input file (`SnapWave.inp`),
//! plan.md Phase 3, steps 1-3.
//!
//! # Grammar of `SnapWave.inp` (documenting `read_snapwave_input` in
//! `src/snapwave_input.f90`, the authority this parser mirrors)
//!
//! The file is a flat list of `key = value` records. There is no comment
//! syntax, no section structure and no quoting. `read_snapwave_input`
//! performs a *separate* rewind-and-scan for every keyword, which is
//! observationally equivalent to these rules:
//!
//! * **Records**: lines are read with `read '(a)'` into a
//!   `character(len=256)` buffer; longer records are truncated at 256
//!   characters. Records are split on `\n` only (gfortran keeps a `\r` as
//!   an ordinary character; here `\r` is additionally accepted as a
//!   numeric value separator, which gfortran list-directed reads also
//!   tolerate in practice).
//! * **Key matching**: the key is everything before the *first* `=` of the
//!   record with trailing blanks trimmed (`trim`). Matching is exact and
//!   **case-sensitive** (`c_dispT` and `map_Cg` only match in that exact
//!   spelling); leading blanks are *not* trimmed, so an indented key never
//!   matches. A record without `=` — or with `=` as its first character —
//!   is silently ignored; several checked-in inputs rely on this to carry
//!   commented-out lines such as `gridfile  some_old_mesh.nc`.
//! * **Duplicate keys**: for every keyword the *first* matching record in
//!   file order wins (the Fortran scan `exit`s on the first hit).
//! * **Unknown keywords** are silently ignored (`wind`, `zs0`, `fw0`,
//!   `nHrel`, `his_ee`, `uuupwindfile`, ... appear in checked-in inputs).
//! * **Character values** (`read_char_input`): everything after the first
//!   `=`, right-trimmed of blanks, then left-adjusted (`adjustl` strips
//!   leading blanks, *not* tabs). Inner blanks are preserved
//!   (`tref = 20240417 000000` keeps the space). The value is then
//!   truncated to the width of the receiving global: `character*15` for
//!   the three date strings, `character*232` for file names,
//!   `character*256` for `fw`/`fwig`/`u10`/`u10dir`.
//! * **Numeric values** (`read_int_input`/`read_real_input`): Fortran
//!   list-directed read of the first item of the record remainder: the
//!   first token delimited by blank/tab/comma/slash. Integers are 32-bit;
//!   reals are `real*4` and accept the Fortran forms seen in the wild,
//!   including `D` exponents (`1d-5`), leading dots (`.02`) and trailing
//!   dots (`3600.`). A missing or malformed value is a hard Fortran
//!   runtime error; here it is a clean Rust error (Phase 3 acceptance:
//!   wrapper failures for invalid input are Rust errors, not Fortran
//!   `stop`/runtime aborts).
//! * **Defaults and post-processing** mirror `read_snapwave_input`
//!   exactly, including the quirks:
//!   - `map_interval`/`his_interval` default to the *parsed* `timestep`
//!     value; a non-positive interval is an error only when the
//!     corresponding output file name is non-empty (a Fortran `stop 1`
//!     path, made a Rust error here).
//!   - the wind switch is off only when the *parsed* `u10` value equals
//!     the literal string `0.0` **and** `windlistfile` is empty — so
//!     `u10 = 0.00` turns wind growth on.
//!   - `mmax`/`nmax` are the input values **plus two** dummy rows/columns
//!     added by `read_snapwave_input`; the struct stores the model-facing
//!     values.
//!   - `tstart`/`tstop` are the seconds between `tref` and the respective
//!     date strings, computed with the same Fliegel & Van Flandern Julian
//!     day formula (truncating integer division) as `snapwave_date.f90`.
//!   - `sigmin`/`sigmax` default to `8*atan(1)/25` resp. `8*atan(1)/1`
//!     evaluated in single-precision arithmetic.
//!   - `restart` and `writetestfiles` are the corresponding integer
//!     keywords interpreted as booleans (non-zero → true).
//!
//! Known deviations (documented; none reachable for checked-in inputs):
//! non-UTF-8 file content is a Rust error where Fortran reads raw bytes;
//! list-directed repeat counts (`3*5.0`) and radix constants are rejected.
//!
//! Field names deliberately stay close to the Fortran globals so the
//! comparison against `snapwave_read_input_dump_c` (see
//! `crate::input_compare`) and later phases remain greppable against
//! `src/snapwave_data.f90`.

use std::path::Path;

use anyhow::{bail, Context, Result};

/// Line buffer width of the Fortran reader (`character(len=256)`).
const FORTRAN_LINE_LEN: usize = 256;

/// Widths of the receiving Fortran globals (character kinds in
/// `src/snapwave_data.f90`); longer values are truncated on assignment.
const WIDTH_DATE: usize = 15;
const WIDTH_FILENAME: usize = 232;
const WIDTH_VALUE_STR: usize = 256;

/// Characters that terminate the first list-directed value of a record
/// (blank, tab, comma, slash — plus carriage return, see module docs).
const LIST_SEPARATORS: [char; 5] = [' ', '\t', ',', '/', '\r'];

/// Time and run-control settings (globals `trefstr`…`tstop`, `timestep`,
/// `dt`, `niter`, `crit`, `restart`).
#[derive(Debug)]
pub struct TimeControl {
    /// Reference date string `yyyymmdd hhmmss` (global `trefstr`).
    pub tref: String,
    /// Start date string (global `tstartstr`).
    pub tstart_str: String,
    /// Stop date string (global `tstopstr`).
    pub tstop_str: String,
    /// Seconds between `tref` and `tstart` (real*8 global `tstart`).
    pub tstart: f64,
    /// Seconds between `tref` and `tstop` (real*8 global `tstop`).
    pub tstop: f64,
    /// Time step [s] (global `timestep`).
    pub timestep: f32,
    /// Time step used by the solver internals [s] (global `dt`).
    pub dt: f32,
    /// Maximum number of iterations (global `niter`).
    pub niter: i32,
    /// Relative accuracy stopping criterion (global `crit`).
    pub crit: f32,
    /// `restart` keyword, non-zero → true (logical global `restart`).
    pub restart: bool,
}

/// Grid / domain description (globals `nmax`…`sector`, `gridfile`,
/// `depfile`, `mskfile`, `indfile`, `upwfile`).
#[derive(Debug)]
pub struct GridConfig {
    /// Number of cells in the first grid direction **plus two dummy rows**
    /// (as left behind in the global by `read_snapwave_input`).
    pub mmax: i32,
    /// Number of cells in the second grid direction plus two dummy rows.
    pub nmax: i32,
    pub dx: f32,
    pub dy: f32,
    pub x0: f32,
    pub y0: f32,
    pub rotation: f32,
    /// Multiplier for bed level values (global `posdwn`).
    pub posdwn: f32,
    /// Spherical (1) or cartesian (0) coordinates (global `sferic`).
    pub sferic: i32,
    /// Directional grid resolution [deg] (global `dtheta`).
    pub dtheta: f32,
    /// Directional sector width [deg] (global `sector`).
    pub sector: f32,
    /// Mesh/grid file name (Delft3D `.grd`, UGRID NetCDF, …).
    pub gridfile: String,
    /// Bathymetry file name (optional).
    pub depfile: String,
    /// Mask file name (optional).
    pub mskfile: String,
    /// Index file name (optional).
    pub indfile: String,
    /// Pre-computed upwind-neighbour file name (optional).
    pub upwfile: String,
}

/// Boundary forcing file names (globals `jonswapfile`…`bzsfile`, `obsfile`,
/// `tol`).
#[derive(Debug)]
pub struct BoundaryConfig {
    /// Single-point JONSWAP spectrum file (optional).
    pub jonswapfile: String,
    /// Space/time-varying boundary support points file (or `none`).
    pub bndfile: String,
    /// Boundary enclosure polyline file (or `none`).
    pub encfile: String,
    /// Neumann boundary polyline file (or `none`).
    pub neumannfile: String,
    pub bhsfile: String,
    pub btpfile: String,
    pub bwdfile: String,
    pub bdsfile: String,
    pub bzsfile: String,
    /// Observation points file (or `none`).
    pub obsfile: String,
    /// Tolerance [m] for boundary points (global `tol`).
    pub tol: f32,
}

/// Wind input (globals `u10str`, `u10dirstr`, `windlistfile`, `mwind`,
/// `wind`). `u10`/`u10dir` stay strings: they hold either a uniform value
/// or a file name, disambiguated later by the Fortran readers.
#[derive(Debug)]
pub struct WindConfig {
    /// Wind speed: uniform value or file name (global `u10str`).
    pub u10: String,
    /// Wind direction [deg nautical]: uniform value or file name.
    pub u10dir: String,
    /// Wind list file name (optional).
    pub windlistfile: String,
    /// Wind input formulation switch (global `mwind`).
    pub mwind: i32,
    /// Wind growth enabled (global `wind`): false only when the parsed
    /// `u10` equals `"0.0"` and `windlistfile` is empty.
    pub enabled: bool,
}

/// Output selection (globals `map_filename`…`ja_save_each_iter`).
#[derive(Debug)]
pub struct OutputConfig {
    pub map_file: String,
    pub his_file: String,
    /// Map output interval [s]; defaults to the parsed `timestep`.
    pub map_interval: f32,
    /// History output interval [s]; defaults to the parsed `timestep`.
    pub his_interval: f32,
    pub map_depth: i32,
    pub map_Hm0: i32,
    pub map_Hig: i32,
    pub map_Tp: i32,
    pub map_dir: i32,
    pub map_dirspr: i32,
    pub map_cg: i32,
    pub map_Dw: i32,
    pub map_Df: i32,
    pub map_SwE: i32,
    pub map_SwA: i32,
    pub map_sig: i32,
    pub map_u10: i32,
    pub map_Dveg: i32,
    pub map_ee: i32,
    pub map_ctheta: i32,
    /// Save map output after each solver iteration (0 = final only).
    pub ja_save_each_iter: i32,
}

/// Diagnostics switches (global `writetestfiles`).
#[derive(Debug)]
pub struct Diagnostics {
    /// `writetestfiles` keyword, non-zero → true.
    pub writetestfiles: bool,
}

/// Vegetation settings (globals `ja_vegetation`, `vegmapfile`).
#[derive(Debug)]
pub struct VegetationConfig {
    pub ja_vegetation: i32,
    pub vegmapfile: String,
}

/// Solver physics knobs (globals `gamma`…`fw_igstr`, `Tpini`, `zsini`,
/// `sigmin`, `sigmax`, `jadcgdx`, `c_dispT`, `ig`, `upwindref`).
#[derive(Debug)]
pub struct PhysicsConfig {
    pub gamma: f32,
    pub alpha: f32,
    pub gammax: f32,
    pub hmin: f32,
    pub fwcutoff: f32,
    /// Friction coefficient: uniform value or file name (global `fwstr`).
    pub fw: String,
    /// IG-wave friction: uniform value or file name (global `fw_igstr`).
    pub fw_ig: String,
    /// Initial wave period [s] (global `Tpini`).
    pub Tpini: f32,
    /// Initial water level [m] (global `zsini`).
    pub zsini: f32,
    /// Minimum frequency [rad/s] (global `sigmin`).
    pub sigmin: f32,
    /// Maximum frequency [rad/s] (global `sigmax`).
    pub sigmax: f32,
    pub jadcgdx: i32,
    pub c_dispT: f32,
    /// Infragravity waves switch (global `ig`).
    pub ig: i32,
    pub upwindref: i32,
}

/// Parsed and validated `SnapWave.inp` configuration, grouped by concern
/// (plan.md Phase 3, step 3). Values are the *model-facing* globals after
/// `read_snapwave_input` post-processing (see the quirks in the module
/// docs).
#[derive(Debug)]
pub struct SnapWaveInput {
    pub time: TimeControl,
    pub grid: GridConfig,
    pub boundary: BoundaryConfig,
    pub wind: WindConfig,
    pub output: OutputConfig,
    pub diagnostics: Diagnostics,
    pub vegetation: VegetationConfig,
    pub physics: PhysicsConfig,
}

/// One `key = value` record after Fortran's line handling.
struct Entry {
    line_no: usize,
    key: String,
    /// Raw record text after the first `=` (blanks not yet trimmed).
    value: String,
}

/// Keyword table with first-occurrence-wins lookup, equivalent to the
/// per-keyword rewind-and-scan of `read_snapwave_input`.
struct Keywords {
    entries: Vec<Entry>,
}

impl Keywords {
    fn new(text: &str) -> Self {
        let mut entries = Vec::new();
        for (idx, record) in text.split('\n').enumerate() {
            let line_no = idx + 1;
            // Truncate to the 256-character line buffer; cut on a char
            // boundary (Fortran truncates bytes, a multi-byte character
            // straddling the limit is beyond pathological for this file
            // format).
            let mut end = FORTRAN_LINE_LEN.min(record.len());
            while end > 0 && !record.is_char_boundary(end) {
                end -= 1;
            }
            let line = &record[..end];
            let Some(eq) = line.find('=') else {
                continue; // no '=': record silently ignored
            };
            // Trailing blanks trimmed, leading blanks kept (Fortran
            // compares trim(line(1:j-1)) == keyword).
            let key = line[..eq].trim_end_matches(' ');
            if key.is_empty() {
                continue; // '=' at position 1: empty key never matches
            }
            entries.push(Entry { line_no, key: key.to_string(), value: line[eq + 1..].to_string() });
        }
        Keywords { entries }
    }

    /// First matching record for `key`, or `None` when absent.
    fn first(&self, key: &str) -> Option<&Entry> {
        self.entries.iter().find(|e| e.key == key)
    }
}

/// Parse `text` (the full content of a `SnapWave.inp`) into a
/// [`SnapWaveInput`]. Errors carry the offending line and keyword.
pub fn parse_str(text: &str) -> Result<SnapWaveInput> {
    let kw = Keywords::new(text);

    // Reads follow the exact call order of read_snapwave_input(); the
    // order matters only where a default references an earlier keyword
    // (map_interval/his_interval default to the parsed timestep).
    let nmax = int_of(&kw, "nmax", 0)?;
    let mmax = int_of(&kw, "mmax", 0)?;
    let dx = real_of(&kw, "dx", 0.0)?;
    let dy = real_of(&kw, "dy", 0.0)?;
    let x0 = real_of(&kw, "x0", 0.0)?;
    let y0 = real_of(&kw, "y0", 0.0)?;
    let rotation = real_of(&kw, "rotation", 0.0)?;
    let posdwn = real_of(&kw, "posdwn", -1.0)?;
    let trefstr = char_of(&kw, "tref", "20000101 000000", WIDTH_DATE)?;
    let tstartstr = char_of(&kw, "tstart", "20000101 000000", WIDTH_DATE)?;
    let tstopstr = char_of(&kw, "tstop", "20000101 000000", WIDTH_DATE)?;
    let timestep = real_of(&kw, "timestep", 3600.0)?;
    let niter = int_of(&kw, "niter", 10)?;
    let crit = real_of(&kw, "crit", 0.00001)?;
    let dt = real_of(&kw, "dt", 36000.0)?;
    let gamma = real_of(&kw, "gamma", 0.7)?;
    let alpha = real_of(&kw, "alpha", 1.0)?;
    let hmin = real_of(&kw, "hmin", 0.1)?;
    let gammax = real_of(&kw, "gammax", 0.6)?;
    let gridfile = char_of(&kw, "gridfile", ".txt", WIDTH_FILENAME)?;
    let sferic = int_of(&kw, "sferic", 0)?;
    let fwstr = char_of(&kw, "fw", "0.01", WIDTH_VALUE_STR)?;
    let fw_igstr = char_of(&kw, "fwig", "0.015", WIDTH_VALUE_STR)?;
    let fwcutoff = real_of(&kw, "fwcutoff", 200.0)?;
    let tol = real_of(&kw, "tol", 10.0)?;
    let dtheta = real_of(&kw, "dtheta", 10.0)?;
    let sector = real_of(&kw, "sector", 180.0)?;
    let jonswapfile = char_of(&kw, "jonswapfile", "", WIDTH_FILENAME)?;
    let windlistfile = char_of(&kw, "windlistfile", "", WIDTH_FILENAME)?;
    let bndfile = char_of(&kw, "bndfile", "none", WIDTH_FILENAME)?;
    let encfile = char_of(&kw, "encfile", "none", WIDTH_FILENAME)?;
    let neumannfile = char_of(&kw, "neumannfile", "none", WIDTH_FILENAME)?;
    let bhsfile = char_of(&kw, "bhsfile", "", WIDTH_FILENAME)?;
    let btpfile = char_of(&kw, "btpfile", "", WIDTH_FILENAME)?;
    let bwdfile = char_of(&kw, "bwdfile", "", WIDTH_FILENAME)?;
    let bdsfile = char_of(&kw, "bdsfile", "", WIDTH_FILENAME)?;
    let bzsfile = char_of(&kw, "bzsfile", "", WIDTH_FILENAME)?;
    let upwfile = char_of(&kw, "upwfile", "", WIDTH_FILENAME)?;
    let mskfile = char_of(&kw, "mskfile", "", WIDTH_FILENAME)?;
    let indfile = char_of(&kw, "indfile", "", WIDTH_FILENAME)?;
    let depfile = char_of(&kw, "depfile", "", WIDTH_FILENAME)?;
    let obsfile = char_of(&kw, "obsfile", "none", WIDTH_FILENAME)?;
    let map_filename = char_of(&kw, "map_file", "", WIDTH_FILENAME)?;
    let his_filename = char_of(&kw, "his_file", "", WIDTH_FILENAME)?;
    let map_interval = real_of(&kw, "map_interval", timestep)?;
    let his_interval = real_of(&kw, "his_interval", timestep)?;
    // The two `stop 1` checks of read_snapwave_input, made Rust errors so
    // the wrapper fails before the Fortran core is ever invoked.
    if !map_filename.is_empty() && map_interval <= 0.0 {
        bail!("map_interval must be positive. (keyword 'map_interval' with map_file = '{map_filename}')");
    }
    if !his_filename.is_empty() && his_interval <= 0.0 {
        bail!("his_interval must be positive. (keyword 'his_interval' with his_file = '{his_filename}')");
    }
    let map_dep = int_of(&kw, "map_depth", 1)?;
    let map_Hm0 = int_of(&kw, "map_Hm0", 1)?;
    let map_Hig = int_of(&kw, "map_Hig", 0)?;
    let map_Tp = int_of(&kw, "map_Tp", 1)?;
    let map_dir = int_of(&kw, "map_dir", 1)?;
    let map_dirspr = int_of(&kw, "map_dirspr", 0)?;
    let map_cg = int_of(&kw, "map_Cg", 0)?;
    let map_Dw = int_of(&kw, "map_Dw", 0)?;
    let map_Df = int_of(&kw, "map_Df", 0)?;
    let map_SwE = int_of(&kw, "map_SwE", 0)?;
    let map_SwA = int_of(&kw, "map_SwA", 0)?;
    let map_sig = int_of(&kw, "map_sig", 0)?;
    let map_u10 = int_of(&kw, "map_u10", 0)?;
    let map_Dveg = int_of(&kw, "map_Dveg", 0)?;
    let writetestfiles_kw = int_of(&kw, "writetestfiles", 0)?;
    let ja_save_each_iter = int_of(&kw, "ja_save_each_iter", 0)?;
    let map_ee = int_of(&kw, "map_ee", 0)?;
    let map_ctheta = int_of(&kw, "map_ctheta", 0)?;
    let irestart = int_of(&kw, "restart", 0)?;
    let u10str = char_of(&kw, "u10", "0.0", WIDTH_VALUE_STR)?;
    let u10dirstr = char_of(&kw, "u10dir", "270.0", WIDTH_VALUE_STR)?;
    let Tpini = real_of(&kw, "Tpini", 1.0)?;
    let mwind = int_of(&kw, "mwind", 2)?;
    // Fortran default expressions: 8.0*atan(1.0)/25.0 and 8.0*atan(1.0)/1.0,
    // evaluated in single precision (real*4) arithmetic.
    let two_pi_f32 = 8.0f32 * 1.0f32.atan();
    let sigmin = real_of(&kw, "sigmin", two_pi_f32 / 25.0)?;
    let sigmax = real_of(&kw, "sigmax", two_pi_f32 / 1.0)?;
    let jadcgdx = int_of(&kw, "jadcgdx", 1)?;
    let c_dispT = real_of(&kw, "c_dispT", 1.0)?;
    let zsini = real_of(&kw, "zsini", 0.0)?;
    let ig = int_of(&kw, "ig", 0)?;
    let upwindref = int_of(&kw, "upwindref", 0)?;
    let ja_vegetation = int_of(&kw, "ja_vegetation", 0)?;
    let vegmapfile = char_of(&kw, "vegmapfile", ".txt", WIDTH_FILENAME)?;

    // Wind switch: string comparison against the *parsed* u10 value, so
    // '0.00' (or a file name) turns wind growth on.
    let wind_enabled = !(u10str == "0.0" && windlistfile.is_empty());

    // Reference-relative times via the snapwave_date conversion.
    let tstart_sec = seconds_between(&trefstr, &tstartstr)
        .with_context(|| format!("keyword 'tstart': date '{}'", tstartstr))?;
    let tstop_sec =
        seconds_between(&trefstr, &tstopstr).with_context(|| format!("keyword 'tstop': date '{}'", tstopstr))?;

    Ok(SnapWaveInput {
        time: TimeControl {
            tref: trefstr,
            tstart_str: tstartstr,
            tstop_str: tstopstr,
            tstart: tstart_sec,
            tstop: tstop_sec,
            timestep,
            dt,
            niter,
            crit,
            restart: irestart != 0,
        },
        grid: GridConfig {
            // read_snapwave_input adds two dummy rows/columns to the input
            // values before storing the globals.
            mmax: mmax.wrapping_add(2),
            nmax: nmax.wrapping_add(2),
            dx,
            dy,
            x0,
            y0,
            rotation,
            posdwn,
            sferic,
            dtheta,
            sector,
            gridfile,
            depfile,
            mskfile,
            indfile,
            upwfile,
        },
        boundary: BoundaryConfig {
            jonswapfile,
            bndfile,
            encfile,
            neumannfile,
            bhsfile,
            btpfile,
            bwdfile,
            bdsfile,
            bzsfile,
            obsfile,
            tol,
        },
        wind: WindConfig { u10: u10str, u10dir: u10dirstr, windlistfile, mwind, enabled: wind_enabled },
        output: OutputConfig {
            map_file: map_filename,
            his_file: his_filename,
            map_interval,
            his_interval,
            map_depth: map_dep,
            map_Hm0,
            map_Hig,
            map_Tp,
            map_dir,
            map_dirspr,
            map_cg,
            map_Dw,
            map_Df,
            map_SwE,
            map_SwA,
            map_sig,
            map_u10,
            map_Dveg,
            map_ee,
            map_ctheta,
            ja_save_each_iter,
        },
        diagnostics: Diagnostics { writetestfiles: writetestfiles_kw != 0 },
        vegetation: VegetationConfig { ja_vegetation, vegmapfile },
        physics: PhysicsConfig {
            gamma,
            alpha,
            gammax,
            hmin,
            fwcutoff,
            fw: fwstr,
            fw_ig: fw_igstr,
            Tpini,
            zsini,
            sigmin,
            sigmax,
            jadcgdx,
            c_dispT,
            ig,
            upwindref,
        },
    })
}

/// Parse the input file at `path` (adds file context to every error).
pub fn parse_file(path: &Path) -> Result<SnapWaveInput> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("cannot read input file {}", path.display()))?;
    parse_str(&text).with_context(|| format!("invalid input file {}", path.display()))
}

// ----------------------------------------------------------------------
// Value readers mirroring read_int_input / read_real_input /
// read_char_input
// ----------------------------------------------------------------------

fn int_of(kw: &Keywords, key: &str, default: i32) -> Result<i32> {
    match kw.first(key) {
        Some(e) => {
            let tok = first_token(&e.value);
            if tok.is_empty() {
                bail!("keyword '{key}' (line {}): missing integer value", e.line_no);
            }
            tok.parse::<i32>()
                .with_context(|| format!("keyword '{key}' (line {}): invalid integer value '{tok}'", e.line_no))
        }
        None => Ok(default),
    }
}

fn real_of(kw: &Keywords, key: &str, default: f32) -> Result<f32> {
    match kw.first(key) {
        Some(e) => {
            let tok = first_token(&e.value);
            if tok.is_empty() {
                bail!("keyword '{key}' (line {}): missing real value", e.line_no);
            }
            // Fortran list-directed reals accept D exponents ('1d-5');
            // Rust does not, so normalize before parsing.
            let normalized: String =
                tok.chars().map(|c| if c == 'd' || c == 'D' { 'e' } else { c }).collect();
            normalized.parse::<f32>().with_context(|| {
                format!("keyword '{key}' (line {}): invalid real value '{tok}'", e.line_no)
            })
        }
        None => Ok(default),
    }
}

fn char_of(kw: &Keywords, key: &str, default: &str, width: usize) -> Result<String> {
    Ok(match kw.first(key) {
        // adjustl(trim(...)): trailing then leading blanks removed; then
        // assignment truncates to the width of the receiving global.
        Some(e) => truncate_chars(e.value.trim_end_matches(' ').trim_start_matches(' '), width),
        None => default.to_string(),
    })
}

/// First list-directed value of a record remainder: leading separators are
/// skipped, the value ends at the first separator.
fn first_token(raw: &str) -> &str {
    let is_sep = |c: char| LIST_SEPARATORS.contains(&c);
    let start = raw.find(|c| !is_sep(c)).unwrap_or(raw.len());
    let rest = &raw[start..];
    let end = rest.find(is_sep).unwrap_or(rest.len());
    &rest[..end]
}

/// Fortran character assignment truncation (to `width` bytes).
fn truncate_chars(s: &str, width: usize) -> String {
    if s.len() <= width {
        return s.to_string();
    }
    let mut end = width;
    while end > 0 && !s.is_char_boundary(end) {
        end -= 1;
    }
    s[..end].to_string()
}

// ----------------------------------------------------------------------
// Date handling mirroring snapwave_date.f90
// ----------------------------------------------------------------------

/// Calendar date and time fields (yyyy, mm, dd, hh, mn, ss).
type DateFields = (i32, i32, i32, i32, i32, i32);

/// Parse the fixed-position date read by `snapwave_date` with the format
/// `'(I4,2I2,1X,3I2)'` from a `character*15` value: `yyyymmdd hhmmss`.
/// The character at position 9 is skipped, blanks inside the numeric
/// fields are ignored (Fortran BLANK='NULL' default), blank fields read as
/// zero, and any other non-digit character is an error (a Fortran
/// formatted-read runtime error).
fn parse_date15(s: &str) -> Result<DateFields> {
    let padded = format!("{s:<width$}", width = WIDTH_DATE);
    let b = padded.as_bytes();
    let field = |range: std::ops::Range<usize>| -> Result<i32> {
        let mut digits = String::new();
        for &c in &b[range] {
            match c {
                b' ' => continue,
                b'0'..=b'9' => digits.push(c as char),
                _ => bail!("invalid character '{}' in date '{s}' (expected yyyymmdd hhmmss)", c as char),
            }
        }
        if digits.is_empty() {
            Ok(0)
        } else {
            Ok(digits.parse::<i32>()?)
        }
    };
    Ok((
        field(0..4)?,   // yyyy
        field(4..6)?,   // mm
        field(6..8)?,   // dd
        field(9..11)?,  // hh (position 8 is the skipped separator)
        field(11..13)?, // mn
        field(13..15)?, // ss
    ))
}

/// Fliegel & Van Flandern Julian day number, identical to `julian_date`
/// in `src/snapwave_date.f90`. Both Fortran and Rust integer division
/// truncate toward zero, which this formula relies on for months < 3.
/// (Computed in i64 so pathological date ranges cannot overflow; Fortran
/// 32-bit integers would wrap only beyond ~68-year spans.)
fn julian_date(yyyy: i32, mm: i32, dd: i32) -> i64 {
    let (yyyy, mm, dd) = (yyyy as i64, mm as i64, dd as i64);
    dd - 32075 + 1461 * (yyyy + 4800 + (mm - 14) / 12) / 4
        + 367 * (mm - 2 - ((mm - 14) / 12) * 12) / 12
        - 3 * ((yyyy + 4900 + (mm - 14) / 12) / 100) / 4
}

/// Seconds between two `yyyymmdd hhmmss` strings (date2 - date1), as
/// `time_difference` in `src/snapwave_date.f90` computes for the globals
/// `tstart`/`tstop` (seconds relative to `tref`).
fn seconds_between(date1: &str, date2: &str) -> Result<f64> {
    let (y1, m1, d1, h1, n1, s1) = parse_date15(date1)?;
    let (y2, m2, d2, h2, n2, s2) = parse_date15(date2)?;
    let jul1 = julian_date(y1, m1, d1);
    let jul2 = julian_date(y2, m2, d2);
    let sec1 = (h1 as i64) * 3600 + (n1 as i64) * 60 + s1 as i64;
    let sec2 = (h2 as i64) * 3600 + (n2 as i64) * 60 + s2 as i64;
    Ok(((jul2 - jul1) * 86400 + sec2 - sec1) as f64)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `8.0 * atan(1.0)` in single-precision arithmetic, so the test
    /// assertions compute exactly the same value as the parser's runtime
    /// default expression (not a pre-computed constant, which may differ
    /// by 1 ulp across libm implementations).
    fn two_pi_f32() -> f32 {
        8.0f32 * 1.0f32.atan()
    }

    #[test]
    fn empty_input_yields_the_documented_defaults() {
        let cfg = parse_str("").expect("empty input must parse to defaults");

        assert_eq!(cfg.time.tref, "20000101 000000");
        assert_eq!(cfg.time.tstart_str, "20000101 000000");
        assert_eq!(cfg.time.tstop_str, "20000101 000000");
        assert_eq!(cfg.time.tstart, 0.0);
        assert_eq!(cfg.time.tstop, 0.0);
        assert_eq!(cfg.time.timestep, 3600.0);
        assert_eq!(cfg.time.dt, 36000.0);
        assert_eq!(cfg.time.niter, 10);
        assert_eq!(cfg.time.crit, 0.00001f32);
        assert!(!cfg.time.restart);

        assert_eq!(cfg.grid.mmax, 2); // 0 input + 2 dummy rows
        assert_eq!(cfg.grid.nmax, 2);
        assert_eq!(cfg.grid.dx, 0.0);
        assert_eq!(cfg.grid.dy, 0.0);
        assert_eq!(cfg.grid.x0, 0.0);
        assert_eq!(cfg.grid.y0, 0.0);
        assert_eq!(cfg.grid.rotation, 0.0);
        assert_eq!(cfg.grid.posdwn, -1.0);
        assert_eq!(cfg.grid.sferic, 0);
        assert_eq!(cfg.grid.dtheta, 10.0);
        assert_eq!(cfg.grid.sector, 180.0);
        assert_eq!(cfg.grid.gridfile, ".txt");
        assert_eq!(cfg.grid.depfile, "");
        assert_eq!(cfg.grid.mskfile, "");
        assert_eq!(cfg.grid.indfile, "");
        assert_eq!(cfg.grid.upwfile, "");

        assert_eq!(cfg.boundary.jonswapfile, "");
        assert_eq!(cfg.boundary.bndfile, "none");
        assert_eq!(cfg.boundary.encfile, "none");
        assert_eq!(cfg.boundary.neumannfile, "none");
        assert_eq!(cfg.boundary.bhsfile, "");
        assert_eq!(cfg.boundary.btpfile, "");
        assert_eq!(cfg.boundary.bwdfile, "");
        assert_eq!(cfg.boundary.bdsfile, "");
        assert_eq!(cfg.boundary.bzsfile, "");
        assert_eq!(cfg.boundary.obsfile, "none");
        assert_eq!(cfg.boundary.tol, 10.0);

        assert_eq!(cfg.wind.u10, "0.0");
        assert_eq!(cfg.wind.u10dir, "270.0");
        assert_eq!(cfg.wind.windlistfile, "");
        assert_eq!(cfg.wind.mwind, 2);
        assert!(!cfg.wind.enabled, "default u10 '0.0' disables wind growth");

        assert_eq!(cfg.output.map_file, "");
        assert_eq!(cfg.output.his_file, "");
        assert_eq!(cfg.output.map_interval, 3600.0); // defaults to timestep
        assert_eq!(cfg.output.his_interval, 3600.0);
        assert_eq!(cfg.output.map_depth, 1);
        assert_eq!(cfg.output.map_Hm0, 1);
        assert_eq!(cfg.output.map_Hig, 0);
        assert_eq!(cfg.output.map_Tp, 1);
        assert_eq!(cfg.output.map_dir, 1);
        assert_eq!(cfg.output.map_dirspr, 0);
        assert_eq!(cfg.output.map_cg, 0);
        assert_eq!(cfg.output.map_Dw, 0);
        assert_eq!(cfg.output.map_Df, 0);
        assert_eq!(cfg.output.map_SwE, 0);
        assert_eq!(cfg.output.map_SwA, 0);
        assert_eq!(cfg.output.map_sig, 0);
        assert_eq!(cfg.output.map_u10, 0);
        assert_eq!(cfg.output.map_Dveg, 0);
        assert_eq!(cfg.output.map_ee, 0);
        assert_eq!(cfg.output.map_ctheta, 0);
        assert_eq!(cfg.output.ja_save_each_iter, 0);

        assert!(!cfg.diagnostics.writetestfiles);

        assert_eq!(cfg.vegetation.ja_vegetation, 0);
        assert_eq!(cfg.vegetation.vegmapfile, ".txt");

        assert_eq!(cfg.physics.gamma, 0.7);
        assert_eq!(cfg.physics.alpha, 1.0);
        assert_eq!(cfg.physics.gammax, 0.6);
        assert_eq!(cfg.physics.hmin, 0.1);
        assert_eq!(cfg.physics.fwcutoff, 200.0);
        assert_eq!(cfg.physics.fw, "0.01");
        assert_eq!(cfg.physics.fw_ig, "0.015");
        assert_eq!(cfg.physics.Tpini, 1.0);
        assert_eq!(cfg.physics.zsini, 0.0);
        assert_eq!(cfg.physics.sigmin, two_pi_f32() / 25.0);
        assert_eq!(cfg.physics.sigmax, two_pi_f32() / 1.0);
        assert_eq!(cfg.physics.jadcgdx, 1);
        assert_eq!(cfg.physics.c_dispT, 1.0);
        assert_eq!(cfg.physics.ig, 0);
        assert_eq!(cfg.physics.upwindref, 0);
    }

    /// One record per keyword read by read_snapwave_input, each with a
    /// distinctive value; pins that every known keyword is recognized.
    #[test]
    fn every_keyword_is_recognized() {
        let text = "\
nmax = 101
mmax = 102
dx = 1.5
dy = 2.5
x0 = 100.25
y0 = 200.5
rotation = 12.5
posdwn = 1.0
tref = 20210101 010203
tstart = 20210203 040506
tstop = 20220304 070809
timestep = 123.5
niter = 77
crit = 0.001
dt = 999.0
gamma = 0.55
alpha = 0.65
hmin = 0.42
gammax = 0.62
gridfile = grid.grd
sferic = 1
fw = fwfile.xyz
fwig = fwigfile.xyz
fwcutoff = 150.0
tol = 7.5
dtheta = 2.5
sector = 90.0
jonswapfile = jon.txt
windlistfile = wind.txt
bndfile = bnd.txt
encfile = enc.txt
neumannfile = neu.txt
bhsfile = hs.txt
btpfile = tp.txt
bwdfile = wd.txt
bdsfile = ds.txt
bzsfile = zs.txt
upwfile = upw.upw
mskfile = msk.msk
indfile = ind.ind
depfile = dep.dep
obsfile = obs.obs
map_file = mapout.nc
his_file = hisout.nc
map_interval = 11.5
his_interval = 12.5
map_depth = 0
map_Hm0 = 0
map_Hig = 1
map_Tp = 0
map_dir = 0
map_dirspr = 1
map_Cg = 1
map_Dw = 1
map_Df = 1
map_SwE = 1
map_SwA = 1
map_sig = 1
map_u10 = 1
map_Dveg = 1
writetestfiles = 1
ja_save_each_iter = 1
map_ee = 1
map_ctheta = 1
restart = 1
u10 = 10.5
u10dir = 250.0
Tpini = 12.25
mwind = 3
sigmin = 0.1
sigmax = 3.3
jadcgdx = 0
c_dispT = 0.5
zsini = 0.75
ig = 1
upwindref = 1
ja_vegetation = 1
vegmapfile = veg.nc
";
        let cfg = parse_str(text).expect("all-keyword input must parse");

        assert_eq!(cfg.grid.nmax, 103);
        assert_eq!(cfg.grid.mmax, 104);
        assert_eq!(cfg.grid.dx, 1.5);
        assert_eq!(cfg.grid.dy, 2.5);
        assert_eq!(cfg.grid.x0, 100.25);
        assert_eq!(cfg.grid.y0, 200.5);
        assert_eq!(cfg.grid.rotation, 12.5);
        assert_eq!(cfg.grid.posdwn, 1.0);
        assert_eq!(cfg.grid.gridfile, "grid.grd");
        assert_eq!(cfg.grid.sferic, 1);
        assert_eq!(cfg.grid.dtheta, 2.5);
        assert_eq!(cfg.grid.sector, 90.0);
        assert_eq!(cfg.grid.depfile, "dep.dep");
        assert_eq!(cfg.grid.mskfile, "msk.msk");
        assert_eq!(cfg.grid.indfile, "ind.ind");
        assert_eq!(cfg.grid.upwfile, "upw.upw");

        assert_eq!(cfg.time.tref, "20210101 010203");
        assert_eq!(cfg.time.tstart_str, "20210203 040506");
        assert_eq!(cfg.time.tstop_str, "20220304 070809");
        // Hand-computed with the snapwave_date formula: 2021-01-01 to
        // 2021-02-03 is 33 days; time-of-day difference 3h3m3s.
        assert_eq!(cfg.time.tstart, 2_862_183.0);
        // 2021-01-01 to 2022-03-04 is 427 days (2021 not a leap year);
        // time-of-day difference 6h6m6s.
        assert_eq!(cfg.time.tstop, 36_914_766.0);
        assert_eq!(cfg.time.timestep, 123.5);
        assert_eq!(cfg.time.dt, 999.0);
        assert_eq!(cfg.time.niter, 77);
        assert_eq!(cfg.time.crit, 0.001);
        assert!(cfg.time.restart);

        assert_eq!(cfg.boundary.jonswapfile, "jon.txt");
        assert_eq!(cfg.boundary.bndfile, "bnd.txt");
        assert_eq!(cfg.boundary.encfile, "enc.txt");
        assert_eq!(cfg.boundary.neumannfile, "neu.txt");
        assert_eq!(cfg.boundary.bhsfile, "hs.txt");
        assert_eq!(cfg.boundary.btpfile, "tp.txt");
        assert_eq!(cfg.boundary.bwdfile, "wd.txt");
        assert_eq!(cfg.boundary.bdsfile, "ds.txt");
        assert_eq!(cfg.boundary.bzsfile, "zs.txt");
        assert_eq!(cfg.boundary.obsfile, "obs.obs");
        assert_eq!(cfg.boundary.tol, 7.5);

        assert_eq!(cfg.wind.u10, "10.5");
        assert_eq!(cfg.wind.u10dir, "250.0");
        assert_eq!(cfg.wind.windlistfile, "wind.txt");
        assert_eq!(cfg.wind.mwind, 3);
        assert!(cfg.wind.enabled, "non-'0.0' u10 enables wind growth");

        assert_eq!(cfg.output.map_file, "mapout.nc");
        assert_eq!(cfg.output.his_file, "hisout.nc");
        assert_eq!(cfg.output.map_interval, 11.5);
        assert_eq!(cfg.output.his_interval, 12.5);
        assert_eq!(cfg.output.map_depth, 0);
        assert_eq!(cfg.output.map_Hm0, 0);
        assert_eq!(cfg.output.map_Hig, 1);
        assert_eq!(cfg.output.map_Tp, 0);
        assert_eq!(cfg.output.map_dir, 0);
        assert_eq!(cfg.output.map_dirspr, 1);
        assert_eq!(cfg.output.map_cg, 1);
        assert_eq!(cfg.output.map_Dw, 1);
        assert_eq!(cfg.output.map_Df, 1);
        assert_eq!(cfg.output.map_SwE, 1);
        assert_eq!(cfg.output.map_SwA, 1);
        assert_eq!(cfg.output.map_sig, 1);
        assert_eq!(cfg.output.map_u10, 1);
        assert_eq!(cfg.output.map_Dveg, 1);
        assert_eq!(cfg.output.map_ee, 1);
        assert_eq!(cfg.output.map_ctheta, 1);
        assert_eq!(cfg.output.ja_save_each_iter, 1);

        assert!(cfg.diagnostics.writetestfiles);

        assert_eq!(cfg.vegetation.ja_vegetation, 1);
        assert_eq!(cfg.vegetation.vegmapfile, "veg.nc");

        assert_eq!(cfg.physics.gamma, 0.55);
        assert_eq!(cfg.physics.alpha, 0.65);
        assert_eq!(cfg.physics.gammax, 0.62);
        assert_eq!(cfg.physics.hmin, 0.42);
        assert_eq!(cfg.physics.fwcutoff, 150.0);
        assert_eq!(cfg.physics.fw, "fwfile.xyz");
        assert_eq!(cfg.physics.fw_ig, "fwigfile.xyz");
        assert_eq!(cfg.physics.Tpini, 12.25);
        assert_eq!(cfg.physics.zsini, 0.75);
        assert_eq!(cfg.physics.sigmin, 0.1);
        assert_eq!(cfg.physics.sigmax, 3.3);
        assert_eq!(cfg.physics.jadcgdx, 0);
        assert_eq!(cfg.physics.c_dispT, 0.5);
        assert_eq!(cfg.physics.ig, 1);
        assert_eq!(cfg.physics.upwindref, 1);
    }

    #[test]
    fn first_occurrence_of_duplicate_keys_wins() {
        let cfg = parse_str("hmin = 0.42\nhmin = 9.9\n").unwrap();
        assert_eq!(cfg.physics.hmin, 0.42);
    }

    #[test]
    fn unknown_keywords_and_records_without_equals_are_ignored() {
        // `his_ee`, `wind`, `uuupwindfile` appear in checked-in inputs;
        // `map_file` without '=' is how testcases "comment out" a line.
        let cfg = parse_str("\
his_ee = 1
wind = 0
uuupwindfile =
zs0=0.
map_file         ../output/should_be_ignored.nc
just some prose line
")
        .unwrap();
        assert_eq!(cfg.output.map_file, "");
        assert_eq!(cfg.output.map_ee, 0);
        assert!(!cfg.wind.enabled);
    }

    #[test]
    fn keys_are_case_sensitive() {
        // Exact-spelling keywords only: TREF/c_dispt/map_cg do not match.
        let cfg = parse_str("\
TREF = 20240417 000000
tref = 20200101 000000
c_dispt = 9.0
c_dispT = 0.5
map_cg = 1
map_Cg = 1
")
        .unwrap();
        assert_eq!(cfg.time.tref, "20200101 000000");
        assert_eq!(cfg.physics.c_dispT, 0.5);
        assert_eq!(cfg.output.map_cg, 1);
    }

    #[test]
    fn leading_blanks_break_key_matching() {
        // Fortran compares trim(line(1:j-1)) == keyword without stripping
        // leading blanks, so indented keys never match.
        let cfg = parse_str("   tref = 20240417 000000\ntstop =20200101 000000\n").unwrap();
        assert_eq!(cfg.time.tref, "20000101 000000"); // default kept
        assert_eq!(cfg.time.tstop_str, "20200101 000000"); // no space needed around '='
    }

    #[test]
    fn char_values_strip_outer_blanks_but_keep_inner_blanks() {
        let cfg = parse_str("tref =    20240417 000000   \nencfile =\tkeeps-leading-tab\n").unwrap();
        assert_eq!(cfg.time.tref, "20240417 000000");
        // adjustl strips blanks only, not tabs.
        assert_eq!(cfg.boundary.encfile, "\tkeeps-leading-tab");
    }

    #[test]
    fn char_values_truncate_to_fortran_field_widths() {
        let long_date = "20240417 000000999"; // 19 chars -> character*15
        let long_file = "f".repeat(300); // -> character*232
        let long_val = "v".repeat(300); // field is character*256 but line buf is 256 chars
        let cfg = parse_str(&format!(
            "tref = {long_date}\nencfile = {long_file}\nfw = {long_val}\n"
        ))
        .unwrap();
        assert_eq!(cfg.time.tref, "20240417 000000");
        assert_eq!(cfg.boundary.encfile, "f".repeat(232));
        assert_eq!(cfg.physics.fw, "v".repeat(251));
    }

    #[test]
    fn real_values_accept_fortran_list_directed_forms() {
        let cfg = parse_str("\
crit = .02
timestep = 3600.
hmin = 1d-2
gamma = 2.5D-1
dt = 1e3
niter = 10
")
        .unwrap();
        assert_eq!(cfg.time.crit, 0.02);
        assert_eq!(cfg.time.timestep, 3600.0);
        assert_eq!(cfg.physics.hmin, 0.01);
        assert_eq!(cfg.physics.gamma, 0.25);
        assert_eq!(cfg.time.dt, 1000.0);
    }

    #[test]
    fn numeric_values_stop_at_the_first_separator() {
        // List-directed reads take the first item; the rest is ignored.
        let cfg = parse_str("timestep = 60, 120 240\nniter = 200/300\n").unwrap();
        assert_eq!(cfg.time.timestep, 60.0);
        assert_eq!(cfg.time.niter, 200);
    }

    #[test]
    fn line_buffer_truncates_at_256_characters() {
        // A '=' beyond position 256 is invisible to the Fortran reader.
        let key = "k".repeat(250);
        let line = format!("{key} {}", "x".repeat(10)); // '=' would be at 251
        assert!(!line.contains('='));
        let line_with_eq = format!("{}={}", "k".repeat(255), "v");
        let cfg = parse_str(&format!("hmin = 0.42\n{line_with_eq}\n")).unwrap();
        assert_eq!(cfg.physics.hmin, 0.42); // the overlong line is ignored

        // 256-character truncation also applies to values.
        let long = "a".repeat(260);
        let cfg = parse_str(&format!("obsfile = {long}\n")).unwrap();
        assert_eq!(cfg.boundary.obsfile, "a".repeat(232)); // field width is tighter
    }

    #[test]
    fn map_interval_defaults_to_the_parsed_timestep() {
        let cfg = parse_str("timestep = 300\nmap_file = m.nc\nhis_file = h.nc\n").unwrap();
        assert_eq!(cfg.output.map_interval, 300.0);
        assert_eq!(cfg.output.his_interval, 300.0);
    }

    #[test]
    fn wind_switch_quirk_depends_on_the_exact_u10_string() {
        assert!(!parse_str("u10 = 0.0\n").unwrap().wind.enabled);
        assert!(parse_str("u10 = 0.00\n").unwrap().wind.enabled, "'0.00' != '0.0'");
        assert!(!parse_str("").unwrap().wind.enabled);
        assert!(parse_str("windlistfile = w.txt\n").unwrap().wind.enabled);
        assert!(parse_str("u10 = 12.3\n").unwrap().wind.enabled);
    }

    #[test]
    fn mmax_and_nmax_include_two_dummy_rows() {
        let cfg = parse_str("mmax = 10\nnmax = 20\n").unwrap();
        assert_eq!(cfg.grid.mmax, 12);
        assert_eq!(cfg.grid.nmax, 22);
    }

    #[test]
    fn julian_day_matches_reference_values() {
        // Reference values from the snapwave_date.f90 header and the
        // Fliegel & Van Flandern paper: 1970-01-01 -> 2440588,
        // 2000-01-01 -> 2451545.
        assert_eq!(julian_date(1970, 1, 1), 2440588);
        assert_eq!(julian_date(2000, 1, 1), 2451545);
        // Consistency across the March boundary where the formula's
        // truncating division matters (Jan/Feb belong to the previous
        // "Roman" year).
        assert_eq!(julian_date(2000, 2, 28) - julian_date(2000, 1, 31), 28);
        assert_eq!(julian_date(2000, 3, 1) - julian_date(2000, 2, 28), 2); // leap day
        assert_eq!(julian_date(2001, 3, 1) - julian_date(2001, 2, 28), 1);
    }

    #[test]
    fn date_parsing_mirrors_the_fortran_format() {
        // Position 9 is skipped entirely (the '1X' in '(I4,2I2,1X,3I2)').
        assert_eq!(parse_date15("20240417 000000").unwrap(), (2024, 4, 17, 0, 0, 0));
        assert_eq!(parse_date15("20240417T010203").unwrap(), (2024, 4, 17, 1, 2, 3));
        // Blanks are ignored inside numeric fields; blank fields read 0,
        // so a date-only value means midnight.
        assert_eq!(parse_date15("20240417").unwrap(), (2024, 4, 17, 0, 0, 0));
        assert_eq!(parse_date15("  240417 000000").unwrap(), (24, 4, 17, 0, 0, 0));
        // Non-digits are an error (a Fortran formatted-read abort).
        assert!(parse_date15("notadate1234567").is_err());
        assert!(parse_date15("2024041X 000000").is_err());
    }

    #[test]
    fn seconds_between_computes_signed_differences() {
        assert_eq!(seconds_between("20240417 000000", "20240417 000000").unwrap(), 0.0);
        assert_eq!(seconds_between("20240417 000000", "20240418 000000").unwrap(), 86400.0);
        assert_eq!(seconds_between("20240417 000000", "20240416 235959").unwrap(), -1.0);
        // Cross month and year boundaries.
        assert_eq!(seconds_between("20240131 235959", "20240201 000000").unwrap(), 1.0);
        assert_eq!(seconds_between("20231231 235959", "20240101 000000").unwrap(), 1.0);
    }

    #[test]
    fn missing_numeric_value_is_an_error() {
        let err = parse_str("timestep =\n").expect_err("empty numeric value must fail");
        assert!(format!("{err:#}").contains("timestep"), "error was: {err:#}");
        let err = parse_str("timestep =   \n").expect_err("blank numeric value must fail");
        assert!(format!("{err:#}").contains("missing"), "error was: {err:#}");
    }

    #[test]
    fn invalid_numeric_values_are_errors() {
        assert!(parse_str("timestep = banana\n").is_err());
        assert!(parse_str("niter = 3.5\n").is_err(), "decimal in integer field is rejected");
        assert!(parse_str("gamma = 0.7x\n").is_err(), "trailing garbage is rejected");
        assert!(parse_str("niter = 9999999999\n").is_err(), "32-bit overflow is rejected");
    }

    #[test]
    fn nonpositive_interval_with_output_enabled_is_an_error() {
        let err = parse_str("map_file = out.nc\nmap_interval = -1\n").expect_err("must fail");
        assert!(format!("{err:#}").contains("map_interval"), "error was: {err:#}");
        let err = parse_str("his_file = out.nc\nhis_interval = 0\n").expect_err("must fail");
        assert!(format!("{err:#}").contains("his_interval"), "error was: {err:#}");
        // Without the corresponding output file, Fortran accepts it.
        assert!(parse_str("map_interval = -1\n").is_ok());
        assert!(parse_str("his_interval = 0\n").is_ok());
    }

    #[test]
    fn invalid_dates_are_errors() {
        let err = parse_str("tref = notadate\n").expect_err("bad date must fail");
        assert!(format!("{err:#}").contains("notadate"), "error was: {err:#}");
    }
}
