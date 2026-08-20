//! Rust parsers for the auxiliary *text* input files of SnapWave
//! (plan.md, Phase 6).
//!
//! # What this module covers
//!
//! The Fortran core reads several plain-text files in addition to
//! `SnapWave.inp`. This module moves their *parsing* into Rust-owned data
//! structs so that the auxiliary input data is Rust-owned before the
//! timestep loop starts (Phase 6 acceptance). The readers mirrored here:
//!
//! | # | File(s)                              | Fortran reader                                       |
//! |---|--------------------------------------|------------------------------------------------------|
//! | 1 | observation points (`obsfile`)       | `read_obs_points` in `src/snapwave_obspoints.f90`    |
//! | 2 | single-point JONSWAP (`jonswapfile`) | `read_boundary_data_singlepoint`                     |
//! | 3 | boundary time series (`bndfile`, `bhsfile`, `btpfile`, `bwdfile`, `bdsfile`, `bzsfile`) | `read_boundary_data_timeseries` |
//! | 4 | wind list + uniform wind (`windlistfile`, `u10`, `u10dir`) | `read_wind_data` / `read_wind_data_from_list` |
//! | 5 | boundary enclosure + Neumann polylines (`encfile`, `neumannfile`) | `read_boundary_enclosure` / `read_neumann_boundary` |
//! | 6 | plain-text mesh / sample files        | `initialize_snapwave_domain` (ASCII branch) + `read_interpolate_map_input` |
//!
//! Families 1–5 are wired into the run path ([`parse_all`]) and pinned
//! against the Fortran oracle through the temporary `snapwave_text_dump_c`
//! hook (see `crate::text_compare`). Family 6 has no checked-in testcase
//! (every checked-in mesh is NetCDF; the sample-file interpolation
//! `triintfast` belongs to Phase 9), so its parsers are provided here with
//! unit tests but are not part of the oracle comparison yet.
//!
//! # Fortran list-directed semantics preserved
//!
//! * Values are blank/comma separated; `/` terminates a record; a quoted
//!   `'…'` or `"…"` literal is skipped by a numeric read (so a trailing
//!   `'` after the last number — present in one checked-in obs file — does
//!   not break the numeric read).
//! * Blank lines are skipped by the record-counting loops (`read (…, *,
//!   iostat=…) dummy`), so the parsed length is the number of non-blank
//!   lines.
//! * Reals accept Fortran `D`/`d` exponents; integer-ish tokens read into
//!   real fields (`-999`, `270.`) are accepted.
//! * The boundary time-series reader overwrites `t_bwv` once per file, so
//!   the *final* time column is the one from `bzsfile` ("Times in btp and
//!   bwd files must be the same as in bhs file!"). This module reproduces
//!   that overwrite.
//! * Wave direction / spreading are stored as read (degrees) and converted
//!   exactly the way `snapwave_data.f90` does (`wd = (270 - wd) * deg2rad`,
//!   `ds = ds * deg2rad`), available via `wd_rad()` / `ds_rad()`.
//!
//! # Deferred to later phases
//!
//! * Feeding this data back to Fortran and bypassing the Fortran readers is
//!   the Phase 8 (data structures) / Phase 9 (interpolation) handoff: the
//!   readers also compute interpolation weights (`make_map_fm`,
//!   `find_boundary_indices`) or sample interpolation (`triintfast`) that
//!   are not part of text parsing. The Fortran readers therefore stay the
//!   runtime authority for now and are pinned by the oracle comparison.

use std::path::Path;

use anyhow::{bail, Context, Result};

use crate::input::SnapWaveInput;
use crate::paths::RunPaths;

/// Width of the `nameobs` global (`character*32` in `snapwave_data.f90`);
/// longer names are truncated on assignment.
const WIDTH_NAME: usize = 32;

// ----------------------------------------------------------------------
// Conversions shared with the Fortran data module
// ----------------------------------------------------------------------

/// `pi = 4.*atan(1.)` evaluated in single precision (the `pi` parameter of
/// `snapwave_data.f90`).
pub fn pi_f32() -> f32 {
    4.0f32 * 1.0f32.atan()
}

/// `deg2rad = pi / 180d0` evaluated the way `snapwave_data.f90` does it:
/// single-precision `pi` divided by a double-precision `180d0`, the
/// mixed-precision result truncated back to `real*4` on assignment.
pub fn deg2rad_f32() -> f32 {
    (pi_f32() as f64 / 180.0f64) as f32
}

// ----------------------------------------------------------------------
// Data structs (Rust-owned; plan.md Phase 6, step 3)
// ----------------------------------------------------------------------

/// One observation point as read from the observation file.
#[derive(Debug, Clone, PartialEq)]
pub struct ObsPoint {
    pub x: f64,
    pub y: f64,
    /// `nameobs`: the quoted name, or the generated `station_%04d` default
    /// (truncated to `character*32`, matching the Fortran global).
    pub name: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ObsPoints {
    pub points: Vec<ObsPoint>,
}

impl ObsPoints {
    pub fn len(&self) -> usize {
        self.points.len()
    }
}

/// Single-point JONSWAP boundary time series (`jonswapfile`), columns
/// `t Hs Tp dir ds zs`. Directions/spreading are kept as read (degrees);
/// see [`JonswapSeries::wd_rad`] / [`JonswapSeries::ds_rad`] for the
/// Fortran conversion.
#[derive(Debug, Clone, PartialEq)]
pub struct JonswapSeries {
    pub t: Vec<f32>,
    pub hs: Vec<f32>,
    pub tp: Vec<f32>,
    pub wd_deg: Vec<f32>,
    pub ds_deg: Vec<f32>,
    pub zs: Vec<f32>,
}

impl JonswapSeries {
    pub fn len(&self) -> usize {
        self.t.len()
    }

    /// `wd_bwv = (270.0 - wd) * deg2rad` (going-to cartesian radians).
    pub fn wd_rad(&self) -> Vec<f32> {
        self.wd_deg.iter().map(|&d| (270.0f32 - d) * deg2rad_f32()).collect()
    }

    /// `ds_bwv = ds * deg2rad`.
    pub fn ds_rad(&self) -> Vec<f32> {
        self.ds_deg.iter().map(|&d| d * deg2rad_f32()).collect()
    }
}

/// Space- and time-varying boundary data (`bndfile` + `bhsfile`/`btpfile`/
/// `bwdfile`/`bdsfile`/`bzsfile`).
///
/// The 2-D fields are flattened in *time-major* order (the file's row
/// order): index `itb * nwbnd + ib` for time step `itb`, boundary point
/// `ib` — matching how the five files are read row-by-row.
#[derive(Debug, Clone, PartialEq)]
pub struct BoundarySeries {
    pub nwbnd: usize,
    pub ntwbnd: usize,
    /// Boundary support point coordinates (`real*8`, from `bndfile`).
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    /// Final time column (from `bzsfile`, after the per-file overwrite).
    pub t: Vec<f32>,
    pub hs: Vec<f32>,
    pub tp: Vec<f32>,
    pub wd_deg: Vec<f32>,
    pub ds_deg: Vec<f32>,
    pub zs: Vec<f32>,
}

impl BoundarySeries {
    /// `wd_bwv = (270.0 - wd) * deg2rad`.
    pub fn wd_rad(&self) -> Vec<f32> {
        self.wd_deg.iter().map(|&d| (270.0f32 - d) * deg2rad_f32()).collect()
    }

    /// `ds_bwv = ds * deg2rad`.
    pub fn ds_rad(&self) -> Vec<f32> {
        self.ds_deg.iter().map(|&d| d * deg2rad_f32()).collect()
    }
}

/// A wind-list record (`read(11,*) t_u10_bwv, u10str, u10dirstr`): the two
/// trailing columns are either uniform values or file names, decided later
/// by `inquire` in `read_wind_data_from_list`.
#[derive(Debug, Clone, PartialEq)]
pub struct WindListRecord {
    pub t: f32,
    pub u10: String,
    pub u10dir: String,
}

#[derive(Debug, Clone, PartialEq)]
pub struct WindList {
    pub records: Vec<WindListRecord>,
}

impl WindList {
    pub fn len(&self) -> usize {
        self.records.len()
    }

    pub fn t(&self) -> Vec<f32> {
        self.records.iter().map(|r| r.t).collect()
    }
}

/// Uniform wind field (`windlistfile` empty). `u10`/`u10dir` are `None`
/// when the configured string is a *file name* rather than a number — the
/// file-backed interpolation is Phase 9 and not part of text parsing.
#[derive(Debug, Clone, PartialEq)]
pub struct UniformWind {
    pub u10: Option<f32>,
    pub u10dir_deg: Option<f32>,
}

/// Boundary input, mirroring the `jonswapfile` vs `bndfile` dispatch of
/// `read_boundary_data`.
#[derive(Debug, Clone, PartialEq)]
pub enum BoundaryInput {
    /// Neither `jonswapfile` nor `bndfile` configured.
    None,
    /// Single-point JONSWAP file.
    Single(JonswapSeries),
    /// Space/time-varying boundary files.
    Timeseries(BoundarySeries),
}

/// Wind input, mirroring the `windlistfile` dispatch of `read_wind_data`.
#[derive(Debug, Clone, PartialEq)]
pub enum WindInput {
    /// `windlistfile` empty: one uniform (or file-backed) wind field.
    Uniform(UniformWind),
    /// `windlistfile` present: a time series of uniform/file-backed winds.
    List(WindList),
}

/// A polyline of `x,y` points (`real*8`), shared by the boundary-enclosure
/// and Neumann-boundary readers. `-999` separator points are preserved
/// verbatim.
#[derive(Debug, Clone, PartialEq)]
pub struct Polyline {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
}

impl Polyline {
    pub fn len(&self) -> usize {
        self.x.len()
    }
}

/// Everything `parse_all` reads, grouped by reader (plan.md Phase 6).
#[derive(Debug, Clone, PartialEq)]
pub struct ParsedTextInputs {
    pub obs: Option<ObsPoints>,
    pub boundary: BoundaryInput,
    pub wind: WindInput,
    pub enclosure: Option<Polyline>,
    pub neumann: Option<Polyline>,
}

/// ASCII unstructured mesh, from the `else` (non-NetCDF, non-index) branch
/// of `initialize_snapwave_domain` (plan.md Phase 6, family 6).
// `allow(dead_code)`: family 6 has no checked-in testcase (every mesh is
// NetCDF), so the parser is only exercised by unit tests until such a case
// exists or the mesh reader migrates (plan.md Phase 8).
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq)]
pub struct AsciiMesh {
    pub no_nodes: usize,
    pub no_faces: usize,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub zb: Vec<f32>,
    pub msk: Vec<i32>,
    /// `face_nodes(4, no_faces)`, flattened `[face][node]` (node-major).
    pub face_nodes: Vec<i32>,
}

/// Sample points (x, y, z triples), the input of `read_interpolate_map_input`
/// (plan.md Phase 6, family 6).
// `allow(dead_code)`: no checked-in case uses the plain-text sample reader
// (sample interpolation is Phase 9); exercised by unit tests only.
#[allow(dead_code)]
#[derive(Debug, Clone, PartialEq)]
pub struct SamplePoints {
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub z: Vec<f64>,
}

// ----------------------------------------------------------------------
// Low-level list-directed helpers
// ----------------------------------------------------------------------

/// Non-blank lines, mirroring the record-counting loops that read a single
/// list-directed value per record (`read(500, *, iostat=stat) dummy` skips
/// blank records).
fn non_blank_lines(text: &str) -> Vec<&str> {
    text.lines().map(|l| l.trim_end_matches('\r')).filter(|l| !l.trim().is_empty()).collect()
}

/// List-directed tokens of one line: blank/comma separated; `/` ends the
/// record; `'…'`/`"…"` literals are skipped (a numeric read never consumes
/// them). A trailing quote is therefore not part of a number token.
fn list_tokens(line: &str) -> Vec<&str> {
    let b = line.as_bytes();
    let mut out = Vec::new();
    let mut i = 0;
    while i < b.len() {
        let c = b[i] as char;
        if c.is_ascii_whitespace() || c == ',' {
            i += 1;
        } else if c == '/' {
            break;
        } else if c == '\'' || c == '"' {
            i += 1;
            while i < b.len() && (b[i] as char) != c {
                i += 1;
            }
            if i < b.len() {
                i += 1; // closing quote
            }
        } else {
            let start = i;
            while i < b.len() {
                let c = b[i] as char;
                if c.is_ascii_whitespace() || c == ',' || c == '/' || c == '\'' || c == '"' {
                    break;
                }
                i += 1;
            }
            out.push(&line[start..i]);
        }
    }
    out
}

/// Parse a Fortran list-directed real (accepts `D`/`d` exponents).
fn parse_f64(tok: &str) -> Result<f64> {
    let norm: String = tok.chars().map(|c| if c == 'd' || c == 'D' { 'e' } else { c }).collect();
    norm.parse::<f64>().with_context(|| format!("invalid real value '{tok}'"))
}

fn parse_f32(tok: &str) -> Result<f32> {
    let norm: String = tok.chars().map(|c| if c == 'd' || c == 'D' { 'e' } else { c }).collect();
    norm.parse::<f32>().with_context(|| format!("invalid real value '{tok}'"))
}

// `parse_i32` is only reached from `parse_ascii_mesh` (family 6), which has
// no checked-in testcase yet; see the `#[allow(dead_code)]` notes below.
#[allow(dead_code)]
fn parse_i32(tok: &str) -> Result<i32> {
    tok.parse::<i32>().with_context(|| format!("invalid integer value '{tok}'"))
}

/// First `n` list-directed real values of a line (real*8 reads).
fn first_n_reals_f64(line: &str, n: usize) -> Result<Vec<f64>> {
    list_tokens(line).into_iter().take(n).map(parse_f64).collect()
}

/// First `n` list-directed real values of a line (real*4 reads).
fn first_n_reals_f32(line: &str, n: usize) -> Result<Vec<f32>> {
    list_tokens(line).into_iter().take(n).map(parse_f32).collect()
}

/// Fortran character-assignment truncation (to `width` bytes).
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

fn read_file(path: &Path) -> Result<String> {
    std::fs::read_to_string(path).with_context(|| format!("cannot read text input {}", path.display()))
}

// ----------------------------------------------------------------------
// Family 1: observation points
// ----------------------------------------------------------------------

/// Parse an observation-points file (plan.md Phase 6, family 1).
///
/// Mirrors `read_obs_points`: each non-blank line is one point; the first
/// two list-directed values are `x`,`y` (`real*8`); an optional quoted
/// name follows. Without a quote the name defaults to `station_%04d`.
pub fn parse_obs_points(text: &str) -> Result<ObsPoints> {
    let lines = non_blank_lines(text);
    let mut points = Vec::with_capacity(lines.len());
    for (idx, line) in lines.iter().enumerate() {
        let n = idx + 1;
        let vals = first_n_reals_f64(line, 2)?;
        if vals.len() < 2 {
            bail!("observation point {}: expected 'x y' coordinates in line '{}'", n, line.trim());
        }
        points.push(ObsPoint { x: vals[0], y: vals[1], name: obs_name(line, n) });
    }
    Ok(ObsPoints { points })
}

/// Name extraction mirroring the quote logic of `read_obs_points`.
fn obs_name(line: &str, n: usize) -> String {
    let sq = line.find('\'');
    let dq = line.find('"');
    let name = match (sq, dq) {
        (None, None) => format!("station_{n:04}"),
        (Some(j1), _) => quoted_name(&line[j1 + 1..], '\''),
        (None, Some(j1)) => quoted_name(&line[j1 + 1..], '"'),
    };
    truncate_chars(&name, WIDTH_NAME)
}

/// Text between the opening quote (already consumed) and the next matching
/// quote, trimmed of blanks (`adjustl(trim(...))`); an unmatched closing
/// quote yields an empty name (the Fortran substring `line2(1:-1)`).
fn quoted_name(after: &str, q: char) -> String {
    let after = after.trim_matches(' ');
    match after.find(q) {
        Some(j2) => after[..j2].trim_matches(' ').to_string(),
        None => String::new(),
    }
}

// ----------------------------------------------------------------------
// Family 2: single-point JONSWAP boundary files
// ----------------------------------------------------------------------

/// Parse a single-point JONSWAP boundary file (plan.md Phase 6, family 2).
pub fn parse_jonswap(text: &str) -> Result<JonswapSeries> {
    let lines = non_blank_lines(text);
    let mut t = Vec::with_capacity(lines.len());
    let mut hs = Vec::with_capacity(lines.len());
    let mut tp = Vec::with_capacity(lines.len());
    let mut wd = Vec::with_capacity(lines.len());
    let mut ds = Vec::with_capacity(lines.len());
    let mut zs = Vec::with_capacity(lines.len());
    for (idx, line) in lines.iter().enumerate() {
        let v = first_n_reals_f32(line, 6)?;
        if v.len() < 6 {
            bail!(
                "jonswap record {}: expected 't Hm0 Tp dir ds zs' in line '{}'",
                idx + 1,
                line.trim()
            );
        }
        t.push(v[0]);
        hs.push(v[1]);
        tp.push(v[2]);
        wd.push(v[3]);
        ds.push(v[4]);
        zs.push(v[5]);
    }
    Ok(JonswapSeries { t, hs, tp, wd_deg: wd, ds_deg: ds, zs })
}

// ----------------------------------------------------------------------
// Family 3: boundary location and time-series files
// ----------------------------------------------------------------------

/// Parse the space/time-varying boundary files (plan.md Phase 6, family 3).
///
/// `nwbnd` comes from the number of points in `bnd`; `ntwbnd` from the
/// number of time rows in `bhs`. The time column is read (and overwritten)
/// by each of the five series files, so the final `t` is the `bzs` column.
pub fn parse_boundary_timeseries(
    bnd: &str,
    bhs: &str,
    btp: &str,
    bwd: &str,
    bds: &str,
    bzs: &str,
) -> Result<BoundarySeries> {
    let bnd_lines = non_blank_lines(bnd);
    let nwbnd = bnd_lines.len();
    if nwbnd == 0 {
        bail!("boundary points file has no non-blank lines");
    }
    let mut x = Vec::with_capacity(nwbnd);
    let mut y = Vec::with_capacity(nwbnd);
    for (idx, line) in bnd_lines.iter().enumerate() {
        let v = first_n_reals_f64(line, 2)?;
        if v.len() < 2 {
            bail!("boundary point {}: expected 'x y' in line '{}'", idx + 1, line.trim());
        }
        x.push(v[0]);
        y.push(v[1]);
    }

    let ntwbnd = non_blank_lines(bhs).len();
    if ntwbnd == 0 {
        bail!("boundary Hs file has no non-blank lines");
    }

    // t is overwritten by each subsequent read; final value is the bzs column.
    let mut t = vec![0.0f32; ntwbnd];
    let hs = read_series(bhs, ntwbnd, nwbnd, &mut t, "Hs")?;
    let tp = read_series(btp, ntwbnd, nwbnd, &mut t, "Tp")?;
    let wd = read_series(bwd, ntwbnd, nwbnd, &mut t, "dir")?;
    let ds = read_series(bds, ntwbnd, nwbnd, &mut t, "dspr")?;
    let zs = read_series(bzs, ntwbnd, nwbnd, &mut t, "water level")?;

    Ok(BoundarySeries { nwbnd, ntwbnd, x, y, t, hs, tp, wd_deg: wd, ds_deg: ds, zs })
}

/// Read one series file: `ntwbnd` rows of `time + nwbnd` values, storing the
/// time column into `t` (overwriting) and returning the flattened values in
/// time-major order (`itb * nwbnd + ib`).
fn read_series(text: &str, ntwbnd: usize, nwbnd: usize, t: &mut [f32], what: &str) -> Result<Vec<f32>> {
    let lines = non_blank_lines(text);
    if lines.len() < ntwbnd {
        bail!(
            "{what} file has {} non-blank lines but the Hs file defines {} time steps",
            lines.len(),
            ntwbnd
        );
    }
    let mut out = vec![0.0f32; ntwbnd * nwbnd];
    for itb in 0..ntwbnd {
        let v = first_n_reals_f32(lines[itb], 1 + nwbnd)?;
        if v.len() < 1 + nwbnd {
            bail!(
                "{what} file line {}: expected time + {} values in '{}'",
                itb + 1,
                nwbnd,
                lines[itb].trim()
            );
        }
        t[itb] = v[0];
        for ib in 0..nwbnd {
            out[itb * nwbnd + ib] = v[1 + ib];
        }
    }
    Ok(out)
}

// ----------------------------------------------------------------------
// Family 4: wind list files and uniform wind
// ----------------------------------------------------------------------

/// Parse a wind-list file (plan.md Phase 6, family 4): one `time u10 u10dir`
/// record per line, the trailing columns being either uniform values or
/// file names.
pub fn parse_wind_list(text: &str) -> Result<WindList> {
    let lines = non_blank_lines(text);
    let mut records = Vec::with_capacity(lines.len());
    for (idx, line) in lines.iter().enumerate() {
        let toks = list_tokens(line);
        if toks.len() < 3 {
            bail!("wind list line {}: expected 'time u10 u10dir' in '{}'", idx + 1, line.trim());
        }
        records.push(WindListRecord {
            t: parse_f32(toks[0])?,
            u10: toks[1].to_string(),
            u10dir: toks[2].to_string(),
        });
    }
    Ok(WindList { records })
}

/// Uniform wind from the `u10`/`u10dir` configuration strings. `None` means
/// the string is a file name (interpolation is Phase 9).
pub fn parse_uniform_wind(u10str: &str, u10dirstr: &str) -> UniformWind {
    UniformWind { u10: parse_f32(u10str).ok(), u10dir_deg: parse_f32(u10dirstr).ok() }
}

// ----------------------------------------------------------------------
// Family 5: boundary enclosure and Neumann polylines
// ----------------------------------------------------------------------

/// Parse a polyline file (`encfile` or `neumannfile`): `x y` per non-blank
/// line, read as `real*8`. `-999` separators are preserved.
pub fn parse_polyline(text: &str) -> Result<Polyline> {
    let lines = non_blank_lines(text);
    let mut x = Vec::with_capacity(lines.len());
    let mut y = Vec::with_capacity(lines.len());
    for (idx, line) in lines.iter().enumerate() {
        let v = first_n_reals_f64(line, 2)?;
        if v.len() < 2 {
            bail!("polyline point {}: expected 'x y' in line '{}'", idx + 1, line.trim());
        }
        x.push(v[0]);
        y.push(v[1]);
    }
    Ok(Polyline { x, y })
}

// ----------------------------------------------------------------------
// Family 6: plain-text mesh / sample readers
// ----------------------------------------------------------------------

/// Parse a plain-text unstructured mesh (plan.md Phase 6, family 6): the
/// `else` branch of `initialize_snapwave_domain` — `no_nodes no_faces`,
/// then `x y zb msk` per node, then `face_nodes(1..4)` per face.
///
/// Note: this reader has no checked-in testcase (all meshes are NetCDF);
/// post-processing (`zb = -posdwn * max(zb, -200)`, `face_nodes(4) == 0 ->
/// -999`) happens outside the read and is not reproduced here.
// `allow(dead_code)`: see the note on [`AsciiMesh`].
#[allow(dead_code)]
pub fn parse_ascii_mesh(text: &str) -> Result<AsciiMesh> {
    let lines = non_blank_lines(text);
    if lines.len() < 1 {
        bail!("empty mesh file");
    }
    let header = list_tokens(lines[0]);
    if header.len() < 2 {
        bail!("mesh header: expected 'no_nodes no_faces' in '{}'", lines[0].trim());
    }
    let no_nodes = parse_i32(header[0])? as usize;
    let no_faces = parse_i32(header[1])? as usize;
    if lines.len() < 1 + no_nodes + no_faces {
        bail!(
            "mesh file has {} lines but {} nodes + {} faces + header are required",
            lines.len(),
            no_nodes,
            no_faces
        );
    }

    let mut x = Vec::with_capacity(no_nodes);
    let mut y = Vec::with_capacity(no_nodes);
    let mut zb = Vec::with_capacity(no_nodes);
    let mut msk = Vec::with_capacity(no_nodes);
    for k in 0..no_nodes {
        let v = lines[1 + k];
        let toks = list_tokens(v);
        if toks.len() < 4 {
            bail!("mesh node {}: expected 'x y zb msk' in '{}'", k + 1, v.trim());
        }
        x.push(parse_f64(toks[0])?);
        y.push(parse_f64(toks[1])?);
        zb.push(parse_f32(toks[2])?);
        msk.push(parse_i32(toks[3])?);
    }

    let mut face_nodes = vec![0i32; no_faces * 4];
    for f in 0..no_faces {
        let v = lines[1 + no_nodes + f];
        let toks = list_tokens(v);
        if toks.len() < 4 {
            bail!("mesh face {}: expected 4 node indices in '{}'", f + 1, v.trim());
        }
        for j in 0..4 {
            face_nodes[f * 4 + j] = parse_i32(toks[j])?;
        }
    }

    Ok(AsciiMesh { no_nodes, no_faces, x, y, zb, msk, face_nodes })
}

/// Parse a sample-points file (plan.md Phase 6, family 6): `x y z` triples
/// (`real*8`), the input of `read_interpolate_map_input` (the interpolation
/// itself is Phase 9).
// `allow(dead_code)`: see the note on [`SamplePoints`].
#[allow(dead_code)]
pub fn parse_samples(text: &str) -> Result<SamplePoints> {
    let lines = non_blank_lines(text);
    let mut x = Vec::with_capacity(lines.len());
    let mut y = Vec::with_capacity(lines.len());
    let mut z = Vec::with_capacity(lines.len());
    for (idx, line) in lines.iter().enumerate() {
        let v = first_n_reals_f64(line, 3)?;
        if v.len() < 3 {
            bail!("sample {}: expected 'x y z' in '{}'", idx + 1, line.trim());
        }
        x.push(v[0]);
        y.push(v[1]);
        z.push(v[2]);
    }
    Ok(SamplePoints { x, y, z })
}

// ----------------------------------------------------------------------
// Orchestration (plan.md Phase 6, step 3)
// ----------------------------------------------------------------------

/// Parse every auxiliary text input referenced by `cfg`, resolving paths
/// through `paths` (plan.md Phase 5). Families 1–5; family 6 (mesh/samples)
/// is not part of this run-path validation because every checked-in mesh is
/// NetCDF and sample interpolation is Phase 9.
pub fn parse_all(paths: &RunPaths, cfg: &SnapWaveInput) -> Result<ParsedTextInputs> {
    let obs = match &paths.obsfile {
        Some(p) => Some(parse_obs_points(&read_file(p)?)?),
        None => None,
    };

    let boundary = if let Some(p) = &paths.jonswapfile {
        BoundaryInput::Single(parse_jonswap(&read_file(p)?)?)
    } else if let Some(p) = &paths.bndfile {
        let bhs = paths
            .bhsfile
            .as_ref()
            .context("bndfile configured without bhsfile")?;
        let btp = paths.btpfile.as_ref().context("bndfile configured without btpfile")?;
        let bwd = paths.bwdfile.as_ref().context("bndfile configured without bwdfile")?;
        let bds = paths.bdsfile.as_ref().context("bndfile configured without bdsfile")?;
        let bzs = paths.bzsfile.as_ref().context("bndfile configured without bzsfile")?;
        BoundaryInput::Timeseries(parse_boundary_timeseries(
            &read_file(p)?,
            &read_file(bhs)?,
            &read_file(btp)?,
            &read_file(bwd)?,
            &read_file(bds)?,
            &read_file(bzs)?,
        )?)
    } else {
        BoundaryInput::None
    };

    let wind = if let Some(p) = &paths.windlistfile {
        WindInput::List(parse_wind_list(&read_file(p)?)?)
    } else {
        WindInput::Uniform(parse_uniform_wind(&cfg.wind.u10, &cfg.wind.u10dir))
    };

    let enclosure = match &paths.encfile {
        Some(p) => Some(parse_polyline(&read_file(p)?)?),
        None => None,
    };
    let neumann = match &paths.neumannfile {
        Some(p) => Some(parse_polyline(&read_file(p)?)?),
        None => None,
    };

    Ok(ParsedTextInputs { obs, boundary, wind, enclosure, neumann })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deg2rad_matches_the_fortran_parameter_recipe() {
        // pi = 4*atan(1) (f32), deg2rad = pi/180d0 truncated to f32.
        let deg2rad = deg2rad_f32();
        let expected = (pi_f32() as f64 / 180.0f64) as f32;
        assert_eq!(deg2rad, expected);
        // Roughly pi/180.
        assert!((deg2rad - 0.017_453_292).abs() < 1e-6);
    }

    #[test]
    fn list_tokens_handles_separators_quotes_and_slash() {
        assert_eq!(list_tokens(" 1.0  2.0\t3.0 "), vec!["1.0", "2.0", "3.0"]);
        assert_eq!(list_tokens("1.0,2.0"), vec!["1.0", "2.0"]);
        assert_eq!(list_tokens("1.0 / 2.0"), vec!["1.0"]); // slash ends the record
        assert_eq!(list_tokens("-64.629997  17.684000'"), vec!["-64.629997", "17.684000"]);
        assert_eq!(list_tokens("97861 514017 '8'"), vec!["97861", "514017"]);
    }

    #[test]
    fn obs_points_without_names_get_station_defaults() {
        let text = "2.0e+03 5.0e+03\n1.9e+03 5.0e+03\n";
        let obs = parse_obs_points(text).unwrap();
        assert_eq!(obs.len(), 2);
        assert_eq!(obs.points[0].x, 2000.0);
        assert_eq!(obs.points[0].y, 5000.0);
        assert_eq!(obs.points[0].name, "station_0001");
        assert_eq!(obs.points[1].name, "station_0002");
    }

    #[test]
    fn obs_points_with_single_quoted_names() {
        let text = "6960.2  11691.1    'obs1'\n9920.0  9964.8     'obs2'\n";
        let obs = parse_obs_points(text).unwrap();
        assert_eq!(obs.points[0].name, "obs1");
        assert_eq!(obs.points[1].name, "obs2");
    }

    #[test]
    fn obs_point_with_trailing_quote_reads_empty_name() {
        // The checked-in `obspoints_stcroix.txt` second line has a dangling
        // `'`: Fortran reads the two reals, then extracts an empty name.
        let text = "-64.719345  17.767447\n-64.629997  17.684000'\n";
        let obs = parse_obs_points(text).unwrap();
        assert_eq!(obs.len(), 2);
        assert_eq!(obs.points[0].name, "station_0001");
        assert_eq!(obs.points[1].x, -64.629997);
        assert_eq!(obs.points[1].y, 17.684000);
        assert_eq!(obs.points[1].name, "");
    }

    #[test]
    fn obs_names_truncate_to_character32() {
        let text = "1.0 2.0 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJK'\n";
        let obs = parse_obs_points(text).unwrap();
        assert_eq!(obs.points[0].name.len(), 32);
        assert_eq!(obs.points[0].name, "abcdefghijklmnopqrstuvwxyzABCDEF");
    }

    #[test]
    fn blank_lines_are_skipped() {
        // reef_enc.txt leads with a blank line: 5 points, not 6.
        let text = "\n-460 -460\n5 -460\n5 460\n-460 460\n-460 -460\n";
        let poly = parse_polyline(text).unwrap();
        assert_eq!(poly.len(), 5);
        assert_eq!(poly.x[0], -460.0);
        assert_eq!(poly.y[4], -460.0);
    }

    #[test]
    fn jonswap_parses_columns_and_converts_directions() {
        let text = "0 .1000 12.0 270. 5.0 1.0\n3600 .1000 12.0 90. 5.0 1.0\n";
        let j = parse_jonswap(text).unwrap();
        assert_eq!(j.len(), 2);
        assert_eq!(j.t, vec![0.0, 3600.0]);
        assert_eq!(j.hs, vec![0.1, 0.1]);
        assert_eq!(j.tp, vec![12.0, 12.0]);
        // (270 - 270) * deg2rad = 0 ; (270 - 90) * deg2rad = 180*deg2rad = pi
        let wd = j.wd_rad();
        assert_eq!(wd[0], 0.0);
        assert!((wd[1] - pi_f32()).abs() < 1e-4);
        let ds = j.ds_rad();
        assert!((ds[0] - 5.0 * deg2rad_f32()).abs() < 1e-6);
    }

    #[test]
    fn boundary_timeseries_single_point() {
        // Case 31: one boundary point, three time steps.
        let bnd = "2000 5000\n";
        let bhs = "0.0 1.00\n3600.0 1.00\n7200.0 1.00\n";
        let btp = "0.0 5.00\n3600.0 5.00\n7200.0 5.00\n";
        let bwd = "0.0 90.0\n3600.0 45.0\n7200.0 30.0\n";
        let bds = "0.0 5\n3600.0 5\n7200.0 5\n";
        let bzs = "0.0 0.00\n3600.0 0.00\n7200.0 0.00\n";
        let s = parse_boundary_timeseries(bnd, bhs, btp, bwd, bds, bzs).unwrap();
        assert_eq!(s.nwbnd, 1);
        assert_eq!(s.ntwbnd, 3);
        assert_eq!(s.x, vec![2000.0]);
        assert_eq!(s.y, vec![5000.0]);
        assert_eq!(s.t, vec![0.0, 3600.0, 7200.0]);
        assert_eq!(s.hs, vec![1.0, 1.0, 1.0]);
        assert_eq!(s.tp, vec![5.0, 5.0, 5.0]);
        // wd: (270-90), (270-45), (270-30) degrees -> radians
        let wd = s.wd_rad();
        assert!((wd[0] - 180.0 * deg2rad_f32()).abs() < 1e-5);
        assert!((wd[1] - 225.0 * deg2rad_f32()).abs() < 1e-5);
        assert!((wd[2] - 240.0 * deg2rad_f32()).abs() < 1e-5);
    }

    #[test]
    fn boundary_timeseries_multi_point_row_major() {
        let bnd = "10 20\n30 40\n";
        let bhs = "0.0 1.0 2.0\n100.0 3.0 4.0\n";
        let btp = "0.0 5.0 6.0\n100.0 7.0 8.0\n";
        let bwd = "0.0 270.0 270.0\n100.0 270.0 270.0\n";
        let bds = "0.0 5.0 5.0\n100.0 5.0 5.0\n";
        let bzs = "0.0 0.0 0.0\n100.0 0.0 0.0\n";
        let s = parse_boundary_timeseries(bnd, bhs, btp, bwd, bds, bzs).unwrap();
        assert_eq!(s.nwbnd, 2);
        assert_eq!(s.ntwbnd, 2);
        assert_eq!(s.hs, vec![1.0, 2.0, 3.0, 4.0]); // itb=0: 1,2 ; itb=1: 3,4
        assert_eq!(s.tp, vec![5.0, 6.0, 7.0, 8.0]);
    }

    #[test]
    fn wind_list_parses_time_and_spec_columns() {
        let text = "0 10.5 windfile_1.txt\n3600 11.0 windfile_2.txt\n";
        let w = parse_wind_list(text).unwrap();
        assert_eq!(w.len(), 2);
        assert_eq!(w.t(), vec![0.0, 3600.0]);
        assert_eq!(w.records[0].u10, "10.5");
        assert_eq!(w.records[0].u10dir, "windfile_1.txt");
    }

    #[test]
    fn uniform_wind_parses_values_and_files() {
        let u = parse_uniform_wind("0.0", "270.0");
        assert_eq!(u.u10, Some(0.0));
        assert_eq!(u.u10dir_deg, Some(270.0));
        let f = parse_uniform_wind("windfile.txt", "270.0");
        assert_eq!(f.u10, None);
        assert_eq!(f.u10dir_deg, Some(270.0));
    }

    #[test]
    fn neumann_polyline_keeps_separator_points() {
        let text = "0 0\n1990 0\n-999 -999\n0 10000\n1990 10000\n";
        let poly = parse_polyline(text).unwrap();
        assert_eq!(poly.len(), 5);
        assert_eq!(poly.x[2], -999.0);
        assert_eq!(poly.y[2], -999.0);
    }

    #[test]
    fn ascii_mesh_parses_header_nodes_and_faces() {
        let text = "3 2\n\
                    0.0 0.0 -5.0 2\n\
                    10.0 0.0 -4.0 1\n\
                    0.0 10.0 -3.0 1\n\
                    1 2 3 0\n\
                    1 3 3 2\n";
        let m = parse_ascii_mesh(text).unwrap();
        assert_eq!(m.no_nodes, 3);
        assert_eq!(m.no_faces, 2);
        assert_eq!(m.x, vec![0.0, 10.0, 0.0]);
        assert_eq!(m.zb, vec![-5.0, -4.0, -3.0]);
        assert_eq!(m.msk, vec![2, 1, 1]);
        assert_eq!(m.face_nodes, vec![1, 2, 3, 0, 1, 3, 3, 2]);
    }

    #[test]
    fn samples_parse_xyz_triples() {
        let text = "113.5 -22.5 0.05\n113.5 -22.495 0.05\n";
        let s = parse_samples(text).unwrap();
        assert_eq!(s.x.len(), 2);
        assert_eq!(s.x[0], 113.5);
        assert_eq!(s.z[1], 0.05);
    }
}
