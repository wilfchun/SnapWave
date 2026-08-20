//! Classic-format (CDF-1) NetCDF reader and writer (plan.md, Phase 7,
//! step 1: "choose a Rust NetCDF strategy").
//!
//! # Why hand-rolled
//!
//! The plan prefers bindings to the installed NetCDF library, but SnapWave
//! only ever reads/writes the *classic* format (CDF-1/2, no NETCDF4 flag),
//! and the project already hand-rolls its read side for exactly this reason
//! (`tests/support/ncdf.rs`, plan.md Phase 1: numeric comparison must work
//! anywhere `cargo test` runs and must report variable/index/tolerance for
//! every failure). Phase 7 needs the *write* side too, and it must produce
//! files the existing dependency-free reader — and `ncdump` — accept, with
//! no new crate dependency (which would also complicate the Nix build's
//! `cargoLock` vendoring). A ~300-line classic-format writer keeps the whole
//! NetCDF surface dependency-free and identical on every toolchain.
//!
//! # Format notes
//!
//! * All integers in the file are **big-endian** (the classic format is
//!   network byte order); floats likewise.
//! * Variables store data in C order (row-major). The Fortran writer
//!   declares dimensions Fortran-style (fastest first) and netCDF-Fortran
//!   reverses them, so a Fortran variable `(/nmesh2d_node, time/)` becomes
//!   `(time, nmesh2d_node)` on disk — the caller is responsible for handing
//!   this writer C-ordered dimensions and C-ordered data.
//! * Record variables (those whose leading dimension is the unlimited one)
//!   are interleaved per record: record 0 of every record variable, then
//!   record 1, … The stride between a variable's slabs is the sum of all
//!   record-variable slab sizes (padded to a 4-byte boundary).
//! * Variable data begins on 4-byte boundaries; `vsize`/`begin` in the
//!   header reflect that padding.

use std::path::Path;

use anyhow::{anyhow, bail, Context, Result};

// NetCDF classic format constants (from the NetCDF Classic Format Spec).
const NC_DIMENSION: u32 = 0x0A;
const NC_VARIABLE: u32 = 0x0B;
const NC_ATTRIBUTE: u32 = 0x0C;

const NC_BYTE: u32 = 1;
const NC_CHAR: u32 = 2;
const NC_SHORT: u32 = 3;
const NC_INT: u32 = 4;
const NC_FLOAT: u32 = 5;
const NC_DOUBLE: u32 = 6;

/// NetCDF default fill values (used to match what the C library writes for
/// variables that are defined but never `put`).
pub const NC_FILL_INT: i32 = -2_147_483_647;
pub const NC_FILL_FLOAT: f32 = 9.969_209_968_386_869e36;

// ----------------------------------------------------------------------
// Shared value types
// ----------------------------------------------------------------------

/// Variable/attribute types of the classic format. SnapWave only *writes*
/// char/int/float, but the mesh reader must tolerate files that also carry
/// byte/short/double variables (it simply never reads them).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum NcType {
    Byte,
    Char,
    Short,
    Int,
    Float,
    Double,
}

impl NcType {
    fn code(self) -> u32 {
        match self {
            NcType::Byte => NC_BYTE,
            NcType::Char => NC_CHAR,
            NcType::Short => NC_SHORT,
            NcType::Int => NC_INT,
            NcType::Float => NC_FLOAT,
            NcType::Double => NC_DOUBLE,
        }
    }

    fn from_code(code: u32) -> Result<Self> {
        match code {
            NC_BYTE => Ok(NcType::Byte),
            NC_CHAR => Ok(NcType::Char),
            NC_SHORT => Ok(NcType::Short),
            NC_INT => Ok(NcType::Int),
            NC_FLOAT => Ok(NcType::Float),
            NC_DOUBLE => Ok(NcType::Double),
            other => bail!("unknown netcdf type code {other}"),
        }
    }

    fn size(self) -> usize {
        match self {
            NcType::Byte | NcType::Char => 1,
            NcType::Short => 2,
            NcType::Int | NcType::Float => 4,
            NcType::Double => 8,
        }
    }

    /// Human-readable type name for error messages.
    fn name(self) -> &'static str {
        match self {
            NcType::Byte => "byte",
            NcType::Char => "char",
            NcType::Short => "short",
            NcType::Int => "int",
            NcType::Float => "float",
            NcType::Double => "double",
        }
    }
}

/// Attribute value (only the kinds SnapWave writes).
#[derive(Clone, Debug)]
pub enum AttrValue {
    Text(String),
    Int(i32),
    Float(f32),
}

#[derive(Clone, Debug)]
pub struct Attr {
    pub name: String,
    pub value: AttrValue,
}

impl Attr {
    pub fn text(name: &str, value: &str) -> Self {
        Attr { name: name.to_string(), value: AttrValue::Text(value.to_string()) }
    }

    pub fn int(name: &str, value: i32) -> Self {
        Attr { name: name.to_string(), value: AttrValue::Int(value) }
    }

    pub fn float(name: &str, value: f32) -> Self {
        Attr { name: name.to_string(), value: AttrValue::Float(value) }
    }

    fn typecode(&self) -> u32 {
        match self.value {
            AttrValue::Text(_) => NC_CHAR,
            AttrValue::Int(_) => NC_INT,
            AttrValue::Float(_) => NC_FLOAT,
        }
    }

    fn elem_size(&self) -> usize {
        match self.value {
            // One nc_char per byte of the string.
            AttrValue::Text(_) => 1,
            AttrValue::Int(_) | AttrValue::Float(_) => 4,
        }
    }

    fn elems(&self) -> usize {
        match self.value {
            // The classic format stores a character attribute with
            // `nelems` = the string length (not 1, unlike numeric scalars).
            AttrValue::Text(ref s) => s.len(),
            AttrValue::Int(_) | AttrValue::Float(_) => 1,
        }
    }
}

// ----------------------------------------------------------------------
// Writer
// ----------------------------------------------------------------------

/// Data of one variable. Fixed variables hold their complete (C-ordered)
/// byte payload; record variables hold one byte payload per record.
#[derive(Clone, Debug)]
pub enum VarData {
    Fixed(Vec<u8>),
    Record(Vec<Vec<u8>>),
}

#[derive(Clone, Debug)]
pub struct Var {
    pub name: String,
    /// Indices into [`Writer::dims`], **C order** (the unlimited dimension,
    /// if any, must be first).
    pub dims: Vec<usize>,
    pub typ: NcType,
    pub attrs: Vec<Attr>,
    pub data: VarData,
}

/// In-memory classic-format NetCDF file builder. Everything is assembled in
/// memory and written out by [`Writer::build`]; SnapWave's map/history files
/// are small enough that streaming is not worth the extra complexity (a
/// future optimization, see the module docs).
#[derive(Clone, Debug, Default)]
pub struct Writer {
    /// Dimension name and length; `None` marks the (single) unlimited
    /// dimension.
    pub dims: Vec<(String, Option<u64>)>,
    pub global_attrs: Vec<Attr>,
    pub vars: Vec<Var>,
}

impl Writer {
    pub fn new() -> Self {
        Writer::default()
    }

    /// Add a dimension, returning its index. `len == None` marks the
    /// unlimited (record) dimension.
    pub fn dim(&mut self, name: &str, len: Option<u64>) -> usize {
        self.dims.push((name.to_string(), len));
        self.dims.len() - 1
    }

    /// Add a global attribute.
    pub fn global_attr(&mut self, attr: Attr) {
        self.global_attrs.push(attr);
    }

    /// Add a variable (with its complete data).
    pub fn var(&mut self, var: Var) {
        self.vars.push(var);
    }

    /// Serialize the file. Returns the complete classic-format byte image.
    pub fn build(&self) -> Result<Vec<u8>> {
        // Exactly one unlimited dimension, and it must be referenced only as
        // the leading dimension of record variables.
        let unlimited: Vec<usize> =
            self.dims.iter().enumerate().filter(|(_, d)| d.1.is_none()).map(|(i, _)| i).collect();
        if unlimited.len() > 1 {
            bail!("more than one unlimited dimension");
        }

        // Determine record counts up front (all record variables must agree),
        // and validate that each variable's data size matches its dimensions.
        let mut numrecs: Option<u64> = None;
        for v in &self.vars {
            let is_record = v.is_record(self);
            match (&v.data, is_record) {
                (VarData::Record(slabs), true) => {
                    let n = slabs.len() as u64;
                    match numrecs {
                        Some(m) if m != n => bail!(
                            "record variable '{}' has {} records, expected {m}",
                            v.name,
                            n
                        ),
                        _ => numrecs = Some(n),
                    }
                    let expected = v.elems_per_slice(self)? * v.typ.size() as u64;
                    for slab in slabs {
                        if slab.len() as u64 != expected {
                            bail!(
                                "record variable '{}': slab has {} bytes, expected {expected}",
                                v.name,
                                slab.len()
                            );
                        }
                    }
                }
                (VarData::Fixed(data), false) => {
                    let expected = v.data_size(self)? as u64;
                    if data.len() as u64 != expected {
                        bail!(
                            "fixed variable '{}': data has {} bytes, expected {expected}",
                            v.name,
                            data.len()
                        );
                    }
                }
                (VarData::Fixed(_), true) => bail!("record variable '{}' has no records", v.name),
                (VarData::Record(_), false) => {
                    bail!("variable '{}' has records but no unlimited dimension", v.name)
                }
            }
        }
        let numrecs = numrecs.unwrap_or(0);

        // Sizes of the header pieces (independent of the begin values, which
        // occupy a fixed 4 bytes each).
        let dim_list_size = list_size(self.dims.len(), |i| {
            string_size(&self.dims[i].0) + 4
        });
        let gattr_size = attr_list_size(&self.global_attrs);
        let var_list_size = list_size(self.vars.len(), |i| {
            string_size(&self.vars[i].name) + 4 + 4 * self.vars[i].dims.len()
                + attr_list_size(&self.vars[i].attrs)
                + 4
                + 4
                + 4 // type, vsize, begin
        });
        let header_size = 4 + 4 + dim_list_size + gattr_size + var_list_size;

        // Fixed-variable layout, then record-variable layout.
        let mut begins = vec![0u32; self.vars.len()];
        let mut vsizes = vec![0u32; self.vars.len()];
        let mut fixed_total: u64 = 0;
        for (i, v) in self.vars.iter().enumerate() {
            if v.is_record(self) {
                continue;
            }
            begins[i] = header_size as u32 + fixed_total as u32;
            let size = v.data_size(self)? as u64;
            vsizes[i] = pad4(size) as u32;
            fixed_total += pad4(size);
        }
        let mut record_offset: u64 = 0;
        for (i, v) in self.vars.iter().enumerate() {
            if !v.is_record(self) {
                continue;
            }
            begins[i] = header_size as u32 + fixed_total as u32 + record_offset as u32;
            let size = v.elems_per_slice(self)? * v.typ.size() as u64;
            vsizes[i] = pad4(size) as u32;
            record_offset += pad4(size);
        }

        let mut out = Vec::with_capacity(header_size as usize + fixed_total as usize + 64);
        // Header
        out.extend_from_slice(b"CDF\x01");
        put_u32(&mut out, numrecs as u32);
        put_list(&mut out, NC_DIMENSION, self.dims.len(), |buf, i| {
            let (name, len) = &self.dims[i];
            put_string(buf, name);
            put_u32(buf, len.unwrap_or(0) as u32);
        });
        put_attr_list(&mut out, &self.global_attrs);
        put_list(&mut out, NC_VARIABLE, self.vars.len(), |buf, i| {
            let v = &self.vars[i];
            put_string(buf, &v.name);
            put_u32(buf, v.dims.len() as u32);
            for &d in &v.dims {
                put_u32(buf, d as u32);
            }
            put_attr_list(buf, &v.attrs);
            put_u32(buf, v.typ.code());
            put_u32(buf, vsizes[i]);
            put_u32(buf, begins[i]);
        });

        // Fixed-variable data.
        for (i, v) in self.vars.iter().enumerate() {
            if v.is_record(self) {
                continue;
            }
            let VarData::Fixed(data) = &v.data else { unreachable!() };
            out.extend_from_slice(data);
            let pad = pad4(data.len() as u64) - data.len() as u64;
            out.extend(std::iter::repeat(0u8).take(pad as usize));
        }

        // Record-variable data, interleaved per record.
        for rec in 0..numrecs {
            for (i, v) in self.vars.iter().enumerate() {
                if !v.is_record(self) {
                    continue;
                }
                let VarData::Record(slabs) = &v.data else { unreachable!() };
                let slab = &slabs[rec as usize];
                out.extend_from_slice(slab);
                let pad = pad4(slab.len() as u64) - slab.len() as u64;
                out.extend(std::iter::repeat(0u8).take(pad as usize));
            }
        }

        Ok(out)
    }

    /// Write the built file to disk (parent directory must already exist;
    /// directory creation is `crate::paths`'s job, plan.md Phase 5).
    pub fn write_to(&self, path: &Path) -> Result<()> {
        let bytes = self.build().with_context(|| format!("building NetCDF file {}", path.display()))?;
        std::fs::write(path, bytes)
            .with_context(|| format!("writing NetCDF output {}", path.display()))
    }
}

impl Var {
    fn is_record(&self, w: &Writer) -> bool {
        match self.dims.first() {
            Some(&id) => w.dims[id].1.is_none(),
            None => false,
        }
    }

    /// Bytes of a fixed variable's complete data.
    fn data_size(&self, w: &Writer) -> Result<usize> {
        let elems: u64 = self
            .dims
            .iter()
            .map(|&id| w.dims[id].1.ok_or_else(|| anyhow!("fixed variable '{}' has an unlimited dimension", self.name)))
            .collect::<Result<Vec<_>>>()?
            .iter()
            .product();
        Ok(elems as usize * self.typ.size())
    }

    /// Elements of one record slab (all non-leading dimensions).
    fn elems_per_slice(&self, w: &Writer) -> Result<u64> {
        let ids: &[usize] = if self.is_record(w) { &self.dims[1..] } else { &self.dims[..] };
        let mut product: u64 = 1;
        for &id in ids {
            let len = w
                .dims
                .get(id)
                .map(|d| d.1.unwrap_or(0))
                .ok_or_else(|| anyhow!("variable '{}' references unknown dimension {id}", self.name))?;
            product *= len;
        }
        Ok(product)
    }
}

// ----------------------------------------------------------------------
// Big-endian serialization helpers
// ----------------------------------------------------------------------

fn put_u32(buf: &mut Vec<u8>, v: u32) {
    buf.extend_from_slice(&v.to_be_bytes());
}

/// Zero-padded-to-4-byte-boundary string, length prefixed.
fn put_string(buf: &mut Vec<u8>, s: &str) {
    put_u32(buf, s.len() as u32);
    buf.extend_from_slice(s.as_bytes());
    pad_to(buf, s.len());
}

/// Pad `buf` so its total length is a multiple of 4 (classic-format
/// string/attribute padding).
fn pad_to(buf: &mut Vec<u8>, elem_len: usize) {
    let pad = (4 - elem_len % 4) % 4;
    buf.extend(std::iter::repeat(0u8).take(pad));
}

fn pad4(n: u64) -> u64 {
    (n + 3) & !3
}

/// Size of one length-prefixed, padded string.
fn string_size(s: &str) -> usize {
    4 + s.len() + (4 - s.len() % 4) % 4
}

/// Size of a classic-format list header (`tag` + `count` + element headers).
fn list_size(count: usize, elem: impl Fn(usize) -> usize) -> usize {
    if count == 0 {
        8 // absent list: tag + count
    } else {
        8 + (0..count).map(elem).sum::<usize>()
    }
}

fn attr_size(a: &Attr) -> usize {
    let nbytes = a.elems() * a.elem_size();
    string_size(&a.name) + 4 + 4 + nbytes + (4 - nbytes % 4) % 4
}

fn attr_list_size(attrs: &[Attr]) -> usize {
    list_size(attrs.len(), |i| attr_size(&attrs[i]))
}

/// Write a classic-format list (`tag`, `count`, then one callback per element).
fn put_list(buf: &mut Vec<u8>, tag: u32, count: usize, elem: impl Fn(&mut Vec<u8>, usize)) {
    if count == 0 {
        put_u32(buf, 0);
        put_u32(buf, 0);
        return;
    }
    put_u32(buf, tag);
    put_u32(buf, count as u32);
    for i in 0..count {
        elem(buf, i);
    }
}

fn put_attr_list(buf: &mut Vec<u8>, attrs: &[Attr]) {
    put_list(buf, NC_ATTRIBUTE, attrs.len(), |b, i| {
        let a = &attrs[i];
        put_string(b, &a.name);
        put_u32(b, a.typecode());
        put_u32(b, a.elems() as u32);
        match &a.value {
            AttrValue::Text(s) => {
                b.extend_from_slice(s.as_bytes());
                pad_to(b, s.len());
            }
            AttrValue::Int(v) => b.extend_from_slice(&v.to_be_bytes()),
            AttrValue::Float(v) => b.extend_from_slice(&v.to_be_bytes()),
        }
    });
}

// ----------------------------------------------------------------------
// Reader (focused: enough for mesh input and round-trip unit tests)
// ----------------------------------------------------------------------

struct Rdr<'a> {
    b: &'a [u8],
    pos: usize,
}

impl<'a> Rdr<'a> {
    fn need(&self, n: usize) -> Result<()> {
        if self.pos.saturating_add(n) > self.b.len() {
            bail!("truncated netcdf file: need {n} bytes at offset {}", self.pos);
        }
        Ok(())
    }

    fn u32(&mut self) -> Result<u32> {
        self.need(4)?;
        let v = u32::from_be_bytes(self.b[self.pos..self.pos + 4].try_into().unwrap());
        self.pos += 4;
        Ok(v)
    }

    fn u64(&mut self) -> Result<u64> {
        self.need(8)?;
        let v = u64::from_be_bytes(self.b[self.pos..self.pos + 8].try_into().unwrap());
        self.pos += 8;
        Ok(v)
    }

    fn take(&mut self, n: usize) -> Result<&'a [u8]> {
        self.need(n)?;
        let s = &self.b[self.pos..self.pos + n];
        self.pos += n;
        Ok(s)
    }

    fn pad4(&mut self, n: usize) {
        self.pos += (4 - n % 4) % 4;
    }

    fn string(&mut self) -> Result<String> {
        let n = self.u32()? as usize;
        let bytes = self.take(n)?;
        let s = String::from_utf8_lossy(bytes).into_owned();
        self.pad4(n);
        Ok(s)
    }
}

/// A parsed classic-format file, trimmed to what the mesh reader needs.
pub struct NetcdfFile {
    data: Vec<u8>,
    numrecs: u64,
    dims: Vec<(String, u64, bool)>, // name, len, unlimited
    vars: Vec<RVar>,
    recsize: u64,
}

pub struct RVar {
    pub name: String,
    dim_ids: Vec<usize>,
    typ: NcType,
    vsize: u64,
    begin: u64,
    is_record: bool,
}

impl NetcdfFile {
    pub fn open(path: &Path) -> Result<Self> {
        let data = std::fs::read(path).with_context(|| format!("reading netcdf file {}", path.display()))?;
        NetcdfFile::from_bytes(data).with_context(|| format!("parsing {}", path.display()))
    }

    pub fn from_bytes(data: Vec<u8>) -> Result<Self> {
        if data.len() < 8 {
            bail!("file too short for a netcdf classic header");
        }
        if &data[0..3] != b"CDF" {
            bail!("not a netcdf classic file (bad magic)");
        }
        let version = data[3];
        if version != 1 && version != 2 {
            bail!("unsupported netcdf classic version {version} (expected CDF-1/2)");
        }
        let offset_size = if version == 1 { 4 } else { 8 };

        let mut r = Rdr { b: &data, pos: 4 };
        let numrecs_raw = r.u32()?;
        if numrecs_raw == u32::MAX {
            bail!("streaming numrecs is not supported");
        }
        let numrecs = numrecs_raw as u64;

        // dim_list
        let mut dims: Vec<(String, u64, bool)> = Vec::new();
        let tag = r.u32()?;
        if tag == 0 {
            let nelems = r.u32()?;
            if nelems != 0 {
                bail!("malformed absent dimension list (count {nelems})");
            }
        } else if tag == NC_DIMENSION {
            let count = r.u32()? as usize;
            for _ in 0..count {
                let name = r.string()?;
                let len = r.u32()? as u64;
                dims.push((name, len, len == 0));
            }
        } else {
            bail!("expected NC_DIMENSION tag, got {tag}");
        }
        if dims.iter().filter(|d| d.2).count() > 1 {
            bail!("more than one unlimited dimension");
        }

        // global attrs (parsed and discarded: the mesh reader does not need them)
        let _gattrs = read_attrs(&mut r)?;

        // var_list
        let mut vars: Vec<RVar> = Vec::new();
        let tag = r.u32()?;
        if tag == 0 {
            let nelems = r.u32()?;
            if nelems != 0 {
                bail!("malformed absent variable list (count {nelems})");
            }
        } else if tag == NC_VARIABLE {
            let count = r.u32()? as usize;
            for _ in 0..count {
                let name = r.string()?;
                let ndims = r.u32()? as usize;
                let mut dim_ids = Vec::with_capacity(ndims);
                for _ in 0..ndims {
                    let id = r.u32()? as usize;
                    if id >= dims.len() {
                        bail!("variable {name}: dimension id {id} out of range");
                    }
                    dim_ids.push(id);
                }
                let _vattrs = read_attrs(&mut r)?;
                let typ = NcType::from_code(r.u32()?)?;
                let vsize = r.u32()? as u64;
                let begin = if offset_size == 4 { r.u32()? as u64 } else { r.u64()? };
                let is_record = !dim_ids.is_empty() && dims[dim_ids[0]].2;
                if dim_ids.iter().skip(1).any(|&id| dims[id].2) {
                    bail!("variable {name}: unlimited dimension not in leading position");
                }
                vars.push(RVar { name, dim_ids, typ, vsize, begin, is_record });
            }
        } else {
            bail!("expected NC_VARIABLE tag, got {tag}");
        }

        let recsize: u64 = vars.iter().filter(|v| v.is_record).map(|v| v.vsize).sum();

        Ok(NetcdfFile { data, numrecs, dims, vars, recsize })
    }

    /// Length of a dimension by name (`None` if absent).
    pub fn dim(&self, name: &str) -> Option<u64> {
        self.dims.iter().find(|d| d.0 == name).map(|d| if d.2 { self.numrecs } else { d.1 })
    }

    fn var(&self, name: &str) -> Option<&RVar> {
        self.vars.iter().find(|v| v.name == name)
    }

    fn elems_per_slice(&self, var: &RVar) -> u64 {
        let ids = if var.is_record { &var.dim_ids[1..] } else { &var.dim_ids[..] };
        ids.iter().map(|&id| self.dims[id].1.max(1)).product()
    }

    /// Base byte offset of one record slab (`s`), or the fixed offset.
    fn base(&self, var: &RVar, s: u64) -> u64 {
        if var.is_record { var.begin + s * self.recsize } else { var.begin }
    }

    /// Read one numeric variable as `f64`, accepting either a `float` or a
    /// `double` source (netCDF-Fortran reads mesh coordinates into `real*8`
    /// regardless of the on-disk type, so both are valid mesh files).
    pub fn read_f64(&self, name: &str) -> Result<Vec<f64>> {
        let var = self.var(name).ok_or_else(|| anyhow!("variable {name:?} not found"))?;
        let per_slice = self.elems_per_slice(var) as usize;
        let slices = if var.is_record { self.numrecs } else { 1 };
        let width = var.typ.size();
        let mut out = Vec::with_capacity(per_slice * slices as usize);
        for s in 0..slices {
            let base = self.base(var, s);
            for i in 0..per_slice {
                let pos = base as usize + i * width;
                let bytes = self
                    .data
                    .get(pos..pos + width)
                    .ok_or_else(|| anyhow!("variable {name:?}: data beyond end of file"))?;
                let v = match var.typ {
                    NcType::Float => f32::from_be_bytes(bytes.try_into().unwrap()) as f64,
                    NcType::Double => f64::from_be_bytes(bytes.try_into().unwrap()),
                    other => bail!("variable {name:?} is {} but a float/double was requested", other.name()),
                };
                out.push(v);
            }
        }
        Ok(out)
    }

    /// Read a numeric variable as `f32` (C order, record-major), accepting a
    /// `float` source (and a `double` source, converted the same way
    /// netCDF-Fortran converts into a `real*4` target).
    pub fn read_f32(&self, name: &str) -> Result<Vec<f32>> {
        let var = self.var(name).ok_or_else(|| anyhow!("variable {name:?} not found"))?;
        let per_slice = self.elems_per_slice(var) as usize;
        let slices = if var.is_record { self.numrecs } else { 1 };
        let width = var.typ.size();
        let mut out = Vec::with_capacity(per_slice * slices as usize);
        for s in 0..slices {
            let base = self.base(var, s);
            for i in 0..per_slice {
                let pos = base as usize + i * width;
                let bytes = self
                    .data
                    .get(pos..pos + width)
                    .ok_or_else(|| anyhow!("variable {name:?}: data beyond end of file"))?;
                let v = match var.typ {
                    NcType::Float => f32::from_be_bytes(bytes.try_into().unwrap()),
                    NcType::Double => f64::from_be_bytes(bytes.try_into().unwrap()) as f32,
                    other => bail!("variable {name:?} is {} but a float was requested", other.name()),
                };
                out.push(v);
            }
        }
        Ok(out)
    }

    /// Read an int variable (same layout rules as [`NetcdfFile::read_f32`]),
    /// accepting an `int` or `short` source (widened like netCDF-Fortran).
    pub fn read_i32(&self, name: &str) -> Result<Vec<i32>> {
        let var = self.var(name).ok_or_else(|| anyhow!("variable {name:?} not found"))?;
        let per_slice = self.elems_per_slice(var) as usize;
        let slices = if var.is_record { self.numrecs } else { 1 };
        let width = var.typ.size();
        let mut out = Vec::with_capacity(per_slice * slices as usize);
        for s in 0..slices {
            let base = self.base(var, s);
            for i in 0..per_slice {
                let pos = base as usize + i * width;
                let bytes = self
                    .data
                    .get(pos..pos + width)
                    .ok_or_else(|| anyhow!("variable {name:?}: data beyond end of file"))?;
                let v = match var.typ {
                    NcType::Int => i32::from_be_bytes(bytes.try_into().unwrap()),
                    NcType::Short => {
                        i16::from_be_bytes(bytes.try_into().unwrap()) as i32
                    }
                    other => bail!("variable {name:?} is {} but an int was requested", other.name()),
                };
                out.push(v);
            }
        }
        Ok(out)
    }
}

/// Read and discard an attribute list (the mesh reader needs none of them).
fn read_attrs(r: &mut Rdr) -> Result<()> {
    let tag = r.u32()?;
    if tag == 0 {
        let nelems = r.u32()?;
        if nelems != 0 {
            bail!("malformed absent attribute list (count {nelems})");
        }
        return Ok(());
    }
    if tag != NC_ATTRIBUTE {
        bail!("expected NC_ATTRIBUTE tag, got {tag}");
    }
    let count = r.u32()? as usize;
    for _ in 0..count {
        let _name = r.string()?;
        let typ = NcType::from_code(r.u32()?)?;
        let nelems = r.u32()? as usize;
        let nbytes = nelems * typ.size();
        r.take(nbytes)?;
        r.pad4(nbytes);
    }
    Ok(())
}

// ----------------------------------------------------------------------
// Helpers for building typed payloads
// ----------------------------------------------------------------------

/// Big-endian byte image of a `&[f32]` (C order).
pub fn f32_payload(values: &[f32]) -> Vec<u8> {
    let mut out = Vec::with_capacity(values.len() * 4);
    for v in values {
        out.extend_from_slice(&v.to_be_bytes());
    }
    out
}

/// Big-endian byte image of a `&[i32]` (C order).
pub fn i32_payload(values: &[i32]) -> Vec<u8> {
    let mut out = Vec::with_capacity(values.len() * 4);
    for v in values {
        out.extend_from_slice(&v.to_be_bytes());
    }
    out
}

/// Big-endian byte image of a single `i32`.
pub fn i32_single(value: i32) -> Vec<u8> {
    value.to_be_bytes().to_vec()
}

/// Big-endian byte image of a single `f32`.
pub fn f32_single(value: f32) -> Vec<u8> {
    value.to_be_bytes().to_vec()
}

/// A `&[f64]` converted element-wise to `f32` (round-to-nearest, matching
/// netCDF-Fortran's real*8 -> real*4 conversion), then to a big-endian
/// payload. Used for `mesh2d_node_x/y`, `station_x/y` and `time`, which the
/// Fortran writer holds as real*8 but writes to NF90_FLOAT variables.
pub fn f64_as_f32_payload(values: &[f64]) -> Vec<u8> {
    let mut out = Vec::with_capacity(values.len() * 4);
    for v in values {
        out.extend_from_slice(&(*v as f32).to_be_bytes());
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a tiny file and read it back through [`NetcdfFile`], checking
    /// the round-trip of a fixed float var, a record float var and a fixed
    /// int var (including the interleaved record layout).
    #[test]
    fn round_trip_fixed_and_record_variables() {
        let mut w = Writer::new();
        let node = w.dim("node", Some(3));
        let time = w.dim("time", None);
        w.var(Var {
            name: "x".into(),
            dims: vec![node],
            typ: NcType::Float,
            attrs: vec![Attr::text("units", "m")],
            data: VarData::Fixed(f32_payload(&[1.0, 2.5, -3.25])),
        });
        w.var(Var {
            name: "face_nodes".into(),
            dims: vec![node],
            typ: NcType::Int,
            attrs: vec![Attr::int("_FillValue", -999)],
            data: VarData::Fixed(i32_payload(&[1, 2, 3])),
        });
        w.var(Var {
            name: "hm0".into(),
            dims: vec![time, node],
            typ: NcType::Float,
            attrs: vec![],
            data: VarData::Record(vec![
                f32_payload(&[0.0, 1.0, 2.0]),
                f32_payload(&[3.0, 4.0, 5.0]),
            ]),
        });

        let bytes = w.build().expect("build");
        let f = NetcdfFile::from_bytes(bytes).expect("parse");

        assert_eq!(f.dim("node"), Some(3));
        assert_eq!(f.dim("time"), Some(2));
        assert_eq!(f.read_f32("x").expect("x"), vec![1.0, 2.5, -3.25]);
        assert_eq!(f.read_i32("face_nodes").expect("face_nodes"), vec![1, 2, 3]);
        // Record variable is record-major: record 0 then record 1.
        assert_eq!(f.read_f32("hm0").expect("hm0"), vec![0.0, 1.0, 2.0, 3.0, 4.0, 5.0]);
    }

    #[test]
    fn record_count_mismatch_is_rejected() {
        let mut w = Writer::new();
        let node = w.dim("node", Some(2));
        let time = w.dim("time", None);
        w.var(Var {
            name: "a".into(),
            dims: vec![time, node],
            typ: NcType::Float,
            attrs: vec![],
            data: VarData::Record(vec![f32_payload(&[1.0, 2.0])]),
        });
        w.var(Var {
            name: "b".into(),
            dims: vec![time, node],
            typ: NcType::Float,
            attrs: vec![],
            data: VarData::Record(vec![f32_payload(&[1.0, 2.0]), f32_payload(&[3.0, 4.0])]),
        });
        assert!(w.build().is_err(), "mismatched record counts must fail");
    }

    #[test]
    fn fill_constants_match_the_netcdf_definitions() {
        assert_eq!(NC_FILL_INT, i32::MIN + 1);
        assert_eq!(NC_FILL_FLOAT.to_bits(), 0x7CF0_0000);
    }
}
