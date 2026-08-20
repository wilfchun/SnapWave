//! Reader for the Fortran output-state capture stream (plan.md, Phase 7,
//! step 5: "switch output writing to Rust using snapshots of Fortran state
//! at output times").
//!
//! # Why a capture stream
//!
//! The Fortran solver remains the numerical authority through Phase 11, so
//! the Rust NetCDF writers need the state the Fortran code would have
//! written. Instead of copying arrays across the FFI boundary (which
//! AGENTS.md discourages), the facade runs the model in *capture mode*:
//! `snapwave_ncoutput` writes the exact buffers it would have handed to
//! `nf90_put_var` into a stream file, and the Rust writer replays them into
//! real NetCDF files. This keeps the FFI boundary coarse (one path in, one
//! file out) and preserves the Fortran computation bit-for-bit — the buffers
//! are the same `real*4` values, including the `where (depth < hmin)` fill
//! masking and the `modulo(270 - …*rad2deg, 360.)` direction wrapping.
//!
//! # Format (little-endian; the supported platforms are x86_64)
//!
//! A header, then tagged blocks until EOF. Every multi-byte value is
//! written by gfortran `access='stream'` (native byte order, which on the
//! supported platforms is little-endian, the same convention as the existing
//! `snapwave.upw` file).
//!
//! ```text
//! header:  magic "SWCA" (4 bytes)  +  version u32 = 1
//! block:   tag u32  +  payload
//! tag 1 = STATIC_MAP, 2 = STATIC_HIS, 3 = MAP_RECORD, 4 = HIS_RECORD
//! string: len u32  +  len bytes
//! ```
//!
//! The per-record field presence mirrors the Fortran `ncoutput_update_map` /
//! `ncoutput_update_his` conditionals and is therefore derived from the
//! resolved configuration, which Rust already owns.

use std::path::Path;

use anyhow::{bail, Context, Result};

use crate::input::SnapWaveInput;

const TAG_STATIC_MAP: u32 = 1;
const TAG_STATIC_HIS: u32 = 2;
const TAG_MAP_RECORD: u32 = 3;
const TAG_HIS_RECORD: u32 = 4;

/// Static map-side state (what `ncoutput_map_init` would have written).
#[derive(Debug, Clone)]
pub struct StaticMap {
    pub no_nodes: usize,
    pub no_faces: usize,
    pub max_nodes: usize,
    pub ntheta: usize,
    pub sferic: i32,
    pub tref_iso8601: String,
    pub libvers: String,
    pub x: Vec<f64>,
    pub y: Vec<f64>,
    pub zb: Vec<f32>,
    /// `face_nodes(1:max_nodes,:)` section, node-major within a face
    /// (`[face*max_nodes + node]`).
    pub face_nodes: Vec<i32>,
    pub fw: Vec<f32>,
    pub fw_ig: Vec<f32>,
    pub veg: Option<VegStatic>,
}

#[derive(Debug, Clone)]
pub struct VegStatic {
    pub ah: Vec<f32>,
    pub bstems: Vec<f32>,
    pub nstems: Vec<f32>,
}

/// Static history-side state (observation points). `tref_iso8601` and
/// `libvers` are duplicated here so a history-only run (no map file) still
/// has them for the `time` variable and the library-version attribute.
#[derive(Debug, Clone)]
pub struct StaticHis {
    pub tref_iso8601: String,
    pub libvers: String,
    pub nobs: usize,
    pub xobs: Vec<f64>,
    pub yobs: Vec<f64>,
    /// `nameobs`, trimmed of trailing blanks.
    pub names: Vec<String>,
}

#[derive(Debug, Clone, Default)]
pub struct MapRecord {
    pub t: f64,
    pub ntmapout: i32,
    pub depth: Option<Vec<f32>>,
    pub hm0: Option<Vec<f32>>,
    pub hm0_ig: Option<Vec<f32>>,
    pub tp: Option<Vec<f32>>,
    pub wd: Option<Vec<f32>>,
    pub wdspr: Option<Vec<f32>>,
    pub cg: Option<Vec<f32>>,
    pub dw: Option<Vec<f32>>,
    pub df: Option<Vec<f32>>,
    pub sw: Option<Vec<f32>>,
    pub st: Option<Vec<f32>>,
    pub sig: Option<Vec<f32>>,
    pub u10: Option<Vec<f32>>,
    pub u10dir: Option<Vec<f32>>,
    pub dveg: Option<Vec<f32>>,
    pub ee: Option<Vec<f32>>,
    pub ctheta: Option<Vec<f32>>,
    pub theta_deg: Option<Vec<f32>>,
}

#[derive(Debug, Clone)]
pub struct HisRecord {
    pub t: f64,
    pub nthisout: i32,
    pub zs: Vec<f32>,
    pub hm0: Vec<f32>,
    pub tp: Vec<f32>,
    pub wavdir: Vec<f32>,
    pub dirspr: Vec<f32>,
    pub hm0ig: Option<Vec<f32>>,
    pub dw: Vec<f32>,
    pub df: Vec<f32>,
    pub sw: Option<Vec<f32>>,
    pub st: Option<Vec<f32>>,
}

/// Everything captured from one model run.
#[derive(Debug, Clone, Default)]
pub struct Capture {
    pub static_map: Option<StaticMap>,
    pub static_his: Option<StaticHis>,
    pub map_records: Vec<MapRecord>,
    pub his_records: Vec<HisRecord>,
}

/// Read and parse a capture stream. Field presence is derived from `cfg`.
pub fn read_capture(path: &Path, cfg: &SnapWaveInput) -> Result<Capture> {
    let data = std::fs::read(path).with_context(|| format!("reading capture stream {}", path.display()))?;
    let mut r = Cursor { data: &data, pos: 0 };

    if r.take(4)? != b"SWCA" {
        bail!("capture stream {}: bad magic (not a SnapWave capture)", path.display());
    }
    if r.u32()? != 1 {
        bail!("capture stream {}: unsupported version", path.display());
    }

    let wind = cfg.wind.enabled;
    let ig = cfg.physics.ig == 1;
    let veg = cfg.vegetation.ja_vegetation == 1;

    let mut cap = Capture::default();
    while r.pos < data.len() {
        let tag = r.u32()?;
        match tag {
            TAG_STATIC_MAP => {
                cap.static_map = Some(read_static_map(&mut r, veg)?);
            }
            TAG_STATIC_HIS => {
                cap.static_his = Some(read_static_his(&mut r)?);
            }
            TAG_MAP_RECORD => {
                let ntheta = cap.static_map.as_ref().map(|m| m.ntheta).unwrap_or(0);
                let no_nodes = cap.static_map.as_ref().map(|m| m.no_nodes).unwrap_or(0);
                cap.map_records.push(read_map_record(&mut r, cfg, ntheta, no_nodes)?);
            }
            TAG_HIS_RECORD => {
                let nobs = cap.static_his.as_ref().map(|h| h.nobs).unwrap_or(0);
                cap.his_records.push(read_his_record(&mut r, ig, wind, nobs)?);
            }
            other => bail!("capture stream {}: unknown block tag {other}", path.display()),
        }
    }
    Ok(cap)
}

fn read_static_map(r: &mut Cursor, veg: bool) -> Result<StaticMap> {
    let no_nodes = r.u32()? as usize;
    let no_faces = r.u32()? as usize;
    let max_nodes = r.u32()? as usize;
    let ntheta = r.u32()? as usize;
    let sferic = r.i32()?;
    let tref_iso8601 = r.string()?;
    let libvers = r.string()?;
    let x = r.f64_array(no_nodes)?;
    let y = r.f64_array(no_nodes)?;
    let zb = r.f32_array(no_nodes)?;
    let face_nodes = r.i32_array(max_nodes * no_faces)?;
    let fw = r.f32_array(no_nodes)?;
    let fw_ig = r.f32_array(no_nodes)?;
    let ja_vegetation = r.i32()?;
    let veg = if ja_vegetation == 1 && veg {
        Some(VegStatic {
            ah: r.f32_array(no_nodes)?,
            bstems: r.f32_array(no_nodes)?,
            nstems: r.f32_array(no_nodes)?,
        })
    } else {
        None
    };
    Ok(StaticMap {
        no_nodes,
        no_faces,
        max_nodes,
        ntheta,
        sferic,
        tref_iso8601,
        libvers,
        x,
        y,
        zb,
        face_nodes,
        fw,
        fw_ig,
        veg,
    })
}

fn read_static_his(r: &mut Cursor) -> Result<StaticHis> {
    let tref_iso8601 = r.string()?;
    let libvers = r.string()?;
    let nobs = r.u32()? as usize;
    let xobs = r.f64_array(nobs)?;
    let yobs = r.f64_array(nobs)?;
    let mut names = Vec::with_capacity(nobs);
    for _ in 0..nobs {
        let raw = r.take(32)?;
        let s = String::from_utf8_lossy(raw).trim_end_matches(' ').to_string();
        names.push(s);
    }
    Ok(StaticHis { tref_iso8601, libvers, nobs, xobs, yobs, names })
}

#[allow(clippy::too_many_arguments)]
fn read_map_record(
    r: &mut Cursor,
    cfg: &SnapWaveInput,
    ntheta: usize,
    no_nodes: usize,
) -> Result<MapRecord> {
    let o = &cfg.output;
    let wind = cfg.wind.enabled;
    let ig = cfg.physics.ig == 1;
    let veg = cfg.vegetation.ja_vegetation == 1;

    let t = r.f64()?;
    let ntmapout = r.i32()?;

    let opt = |cond: bool, r: &mut Cursor, n: usize| -> Result<Option<Vec<f32>>> {
        if cond { Ok(Some(r.f32_array(n)?)) } else { Ok(None) }
    };

    let depth = opt(o.map_depth == 1, r, no_nodes)?;
    let hm0 = opt(o.map_Hm0 == 1, r, no_nodes)?;
    let hm0_ig = opt(ig && o.map_Hig == 1, r, no_nodes)?;
    let tp = opt(o.map_Tp == 1, r, no_nodes)?;
    let wd = opt(o.map_dir == 1, r, no_nodes)?;
    let wdspr = opt(o.map_dirspr == 1, r, no_nodes)?;
    let cg = opt(o.map_cg == 1, r, no_nodes)?;
    let dw = opt(o.map_Dw == 1, r, no_nodes)?;
    let df = opt(o.map_Df == 1, r, no_nodes)?;
    let sw = opt(wind && o.map_SwE == 1, r, no_nodes)?;
    let st = opt(wind && o.map_SwA == 1, r, no_nodes)?;
    let sig = opt(wind && o.map_sig == 1, r, no_nodes)?;
    let u10 = opt(wind && o.map_u10 == 1, r, no_nodes)?;
    let u10dir = opt(wind && o.map_u10 == 1, r, no_nodes)?;
    let dveg = opt(veg && o.map_Dveg == 1, r, no_nodes)?;
    let ee = opt(o.map_ee == 1, r, ntheta * no_nodes)?;
    let ctheta = opt(o.map_ctheta == 1, r, ntheta * no_nodes)?;
    let theta_deg = opt(o.map_ee == 1 || o.map_ctheta == 1, r, ntheta)?;

    Ok(MapRecord {
        t,
        ntmapout,
        depth,
        hm0,
        hm0_ig,
        tp,
        wd,
        wdspr,
        cg,
        dw,
        df,
        sw,
        st,
        sig,
        u10,
        u10dir,
        dveg,
        ee,
        ctheta,
        theta_deg,
    })
}

fn read_his_record(r: &mut Cursor, ig: bool, wind: bool, nobs: usize) -> Result<HisRecord> {
    let t = r.f64()?;
    let nthisout = r.i32()?;
    let zs = r.f32_array(nobs)?;
    let hm0 = r.f32_array(nobs)?;
    let tp = r.f32_array(nobs)?;
    let wavdir = r.f32_array(nobs)?;
    let dirspr = r.f32_array(nobs)?;
    let hm0ig = if ig { Some(r.f32_array(nobs)?) } else { None };
    let dw = r.f32_array(nobs)?;
    let df = r.f32_array(nobs)?;
    let sw = if wind { Some(r.f32_array(nobs)?) } else { None };
    let st = if wind { Some(r.f32_array(nobs)?) } else { None };
    Ok(HisRecord { t, nthisout, zs, hm0, tp, wavdir, dirspr, hm0ig, dw, df, sw, st })
}

// ----------------------------------------------------------------------
// Little-endian byte cursor
// ----------------------------------------------------------------------

struct Cursor<'a> {
    data: &'a [u8],
    pos: usize,
}

impl<'a> Cursor<'a> {
    fn take(&mut self, n: usize) -> Result<&'a [u8]> {
        if self.pos.saturating_add(n) > self.data.len() {
            bail!("capture stream truncated at offset {} (need {n} bytes)", self.pos);
        }
        let s = &self.data[self.pos..self.pos + n];
        self.pos += n;
        Ok(s)
    }

    fn u32(&mut self) -> Result<u32> {
        Ok(u32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }

    fn i32(&mut self) -> Result<i32> {
        Ok(i32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }

    fn f32(&mut self) -> Result<f32> {
        Ok(f32::from_le_bytes(self.take(4)?.try_into().unwrap()))
    }

    fn f64(&mut self) -> Result<f64> {
        Ok(f64::from_le_bytes(self.take(8)?.try_into().unwrap()))
    }

    fn string(&mut self) -> Result<String> {
        let n = self.u32()? as usize;
        let bytes = self.take(n)?;
        Ok(String::from_utf8_lossy(bytes).into_owned())
    }

    fn f32_array(&mut self, n: usize) -> Result<Vec<f32>> {
        let mut out = Vec::with_capacity(n);
        for _ in 0..n {
            out.push(self.f32()?);
        }
        Ok(out)
    }

    fn f64_array(&mut self, n: usize) -> Result<Vec<f64>> {
        let mut out = Vec::with_capacity(n);
        for _ in 0..n {
            out.push(self.f64()?);
        }
        Ok(out)
    }

    fn i32_array(&mut self, n: usize) -> Result<Vec<i32>> {
        let mut out = Vec::with_capacity(n);
        for _ in 0..n {
            out.push(self.i32()?);
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // The capture reader is exercised end-to-end by tests/netcdf_io.rs
    // (through a real model run); here we only pin the stream-format helpers
    // that don't require a Fortran binary.

    #[test]
    fn cursor_reads_little_endian_values() {
        let mut data = Vec::new();
        data.extend_from_slice(b"SWCA");
        data.extend_from_slice(&1u32.to_le_bytes());
        let mut r = Cursor { data: &data, pos: 0 };
        assert_eq!(r.take(4).unwrap(), b"SWCA");
        assert_eq!(r.u32().unwrap(), 1);
        assert_eq!(r.pos, data.len());
    }
}
