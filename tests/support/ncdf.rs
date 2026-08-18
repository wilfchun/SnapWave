//! Minimal read-only parser for NetCDF classic files (CDF-1 / CDF-2).
//!
//! Why hand-rolled (plan.md Phase 1): SnapWave writes classic-format files
//! (`NF90_CLOBBER`, no NETCDF4 flag), and the regression tests must read
//! schema *and* numbers wherever `cargo test` runs — including environments
//! without `ncdump` on PATH. Depending on an external tool for numeric
//! comparison would also make failure reporting coarser than the Phase-1
//! acceptance criteria require (variable, index, tolerance).
//!
//! Layout notes that matter to callers:
//! - All integers in the file are big-endian.
//! - Variables are stored in C order. The Fortran writer declares dimensions
//!   Fortran-style (fastest first), and the netCDF-Fortran layer reverses
//!   them, so a variable declared `(/nmesh2d_node, time/)` in
//!   `snapwave_ncoutput.F90` appears as `(time, nmesh2d_node)` here.
//! - For record variables (those using the unlimited dimension, which is
//!   always the leading dimension in C order) one record slab per variable is
//!   interleaved in the data section; the stride between a variable's slabs
//!   is the sum of all record-variable slab sizes.
//! - [`NcFile::read_f32`] returns data record-major and row-major within each
//!   slab, e.g. `hm0[t * nmesh2d_node + k]`.

use std::fs;
use std::path::Path;

use anyhow::{anyhow, bail, Context, Result};

// Format tags and type codes from the NetCDF Classic Format Specification.
const NC_DIMENSION: u32 = 0x0A;
const NC_VARIABLE: u32 = 0x0B;
const NC_ATTRIBUTE: u32 = 0x0C;

const NC_BYTE: u32 = 1;
const NC_CHAR: u32 = 2;
const NC_SHORT: u32 = 3;
const NC_INT: u32 = 4;
const NC_FLOAT: u32 = 5;
const NC_DOUBLE: u32 = 6;

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

    pub fn name(self) -> &'static str {
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

#[derive(Clone, Debug)]
pub enum NcAttrValue {
    Text(String),
    Bytes(Vec<i8>),
    Shorts(Vec<i16>),
    Ints(Vec<i32>),
    Floats(Vec<f32>),
    Doubles(Vec<f64>),
}

#[derive(Clone, Debug)]
pub struct NcDim {
    pub name: String,
    /// Declared size; 0 marks the unlimited (record) dimension.
    pub declared_len: u64,
    /// Effective size; the record count for the unlimited dimension.
    pub len: u64,
    pub unlimited: bool,
}

#[derive(Clone, Debug)]
pub struct NcAttr {
    pub name: String,
    pub value: NcAttrValue,
}

#[derive(Clone, Debug)]
pub struct NcVar {
    pub name: String,
    pub dim_ids: Vec<usize>,
    pub typ: NcType,
    pub attrs: Vec<NcAttr>,
    /// Padded on-disk size in bytes of one record slab (record variables only).
    vsize: u64,
    /// Absolute byte offset of the first slab.
    begin: u64,
    pub is_record: bool,
}

pub struct NcFile {
    data: Vec<u8>,
    pub numrecs: u64,
    pub dims: Vec<NcDim>,
    pub global_attrs: Vec<NcAttr>,
    pub vars: Vec<NcVar>,
    /// Bytes occupied by one interleaved record (sum of record slab sizes).
    recsize: u64,
}

struct Reader<'a> {
    b: &'a [u8],
    pos: usize,
}

impl<'a> Reader<'a> {
    fn need(&self, n: usize) -> Result<()> {
        if self.pos.saturating_add(n) > self.b.len() {
            bail!("truncated file: need {n} bytes at offset {}", self.pos);
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

    /// Classic-format strings and attribute values are zero-padded to a
    /// 4-byte boundary; skip the padding bytes.
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

fn attr_list(r: &mut Reader) -> Result<Vec<NcAttr>> {
    let tag = r.u32()?;
    if tag == 0 {
        // ABSENT list: two zero words, no elements.
        let nelems = r.u32()?;
        if nelems != 0 {
            bail!("malformed absent attribute list (count {nelems})");
        }
        return Ok(Vec::new());
    }
    if tag != NC_ATTRIBUTE {
        bail!("expected NC_ATTRIBUTE tag, got {tag}");
    }
    let count = r.u32()? as usize;
    let mut out = Vec::with_capacity(count);
    for _ in 0..count {
        let name = r.string()?;
        let typ = NcType::from_code(r.u32()?)?;
        let nelems = r.u32()? as usize;
        let nbytes = nelems * typ.size();
        let bytes = r.take(nbytes)?;
        let value = match typ {
            NcType::Char => NcAttrValue::Text(String::from_utf8_lossy(bytes).into_owned()),
            NcType::Byte => NcAttrValue::Bytes(bytes.iter().map(|&b| b as i8).collect()),
            NcType::Short => NcAttrValue::Shorts(
                bytes
                    .chunks_exact(2)
                    .map(|c| i16::from_be_bytes([c[0], c[1]]))
                    .collect(),
            ),
            NcType::Int => NcAttrValue::Ints(
                bytes
                    .chunks_exact(4)
                    .map(|c| i32::from_be_bytes(c.try_into().unwrap()))
                    .collect(),
            ),
            NcType::Float => NcAttrValue::Floats(
                bytes
                    .chunks_exact(4)
                    .map(|c| f32::from_be_bytes(c.try_into().unwrap()))
                    .collect(),
            ),
            NcType::Double => NcAttrValue::Doubles(
                bytes
                    .chunks_exact(8)
                    .map(|c| f64::from_be_bytes(c.try_into().unwrap()))
                    .collect(),
            ),
        };
        r.pad4(nbytes);
        out.push(NcAttr { name, value });
    }
    Ok(out)
}

impl NcFile {
    pub fn open(path: &Path) -> Result<NcFile> {
        let data =
            fs::read(path).with_context(|| format!("reading netcdf file {}", path.display()))?;
        NcFile::from_bytes(data).with_context(|| format!("parsing {}", path.display()))
    }

    pub fn from_bytes(data: Vec<u8>) -> Result<NcFile> {
        if data.len() < 8 {
            bail!("file too short for a netcdf classic header");
        }
        if &data[0..3] != b"CDF" {
            bail!("not a netcdf classic file (bad magic)");
        }
        let version = data[3];
        if version == 5 {
            bail!("CDF-5 is not supported (SnapWave writes classic CDF-1/2)");
        }
        if version != 1 && version != 2 {
            bail!("unknown netcdf classic version {version}");
        }
        let offset_size = if version == 1 { 4 } else { 8 };

        let mut r = Reader { b: &data, pos: 4 };
        let numrecs_raw = r.u32()?;
        if numrecs_raw == u32::MAX {
            bail!("streaming numrecs is not supported");
        }
        let numrecs = numrecs_raw as u64;

        // dim_list
        let mut dims: Vec<NcDim> = Vec::new();
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
                dims.push(NcDim { name, declared_len: len, len, unlimited: len == 0 });
            }
        } else {
            bail!("expected NC_DIMENSION tag, got {tag}");
        }
        if dims.iter().filter(|d| d.unlimited).count() > 1 {
            bail!("more than one unlimited dimension");
        }
        for d in dims.iter_mut() {
            if d.unlimited {
                d.len = numrecs;
            }
        }

        let global_attrs = attr_list(&mut r)?;

        // var_list: name, ndims, dimids, vatt_list, type, vsize, begin
        let mut vars: Vec<NcVar> = Vec::new();
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
                let attrs = attr_list(&mut r)?;
                let typ = NcType::from_code(r.u32()?)?;
                let vsize = r.u32()? as u64;
                let begin = if offset_size == 4 { r.u32()? as u64 } else { r.u64()? };
                let is_record = !dim_ids.is_empty() && dims[dim_ids[0]].unlimited;
                if dim_ids.iter().skip(1).any(|&id| dims[id].unlimited) {
                    bail!("variable {name}: unlimited dimension not in leading (record) position");
                }
                vars.push(NcVar { name, dim_ids, typ, attrs, vsize, begin, is_record });
            }
        } else {
            bail!("expected NC_VARIABLE tag, got {tag}");
        }

        let recsize: u64 = vars.iter().filter(|v| v.is_record).map(|v| v.vsize).sum();

        Ok(NcFile { data, numrecs, dims, global_attrs, vars, recsize })
    }

    pub fn dim(&self, name: &str) -> Option<&NcDim> {
        self.dims.iter().find(|d| d.name == name)
    }

    pub fn var(&self, name: &str) -> Option<&NcVar> {
        self.vars.iter().find(|v| v.name == name)
    }

    pub fn var_dim_names(&self, name: &str) -> Option<Vec<&str>> {
        self.var(name).map(|v| v.dim_ids.iter().map(|&id| self.dims[id].name.as_str()).collect())
    }

    /// Number of records (frames) in the file.
    pub fn record_count(&self) -> u64 {
        self.numrecs
    }

    /// Elements of one record slab (record vars) or of the whole variable.
    fn elems_per_slice(&self, var: &NcVar) -> u64 {
        let ids = if var.is_record { &var.dim_ids[1..] } else { &var.dim_ids[..] };
        ids.iter().map(|&id| self.dims[id].len).product()
    }

    /// Total element count (all records included).
    pub fn total_elems(&self, var: &NcVar) -> u64 {
        if var.is_record {
            self.numrecs * self.elems_per_slice(var)
        } else {
            self.elems_per_slice(var)
        }
    }

    fn read_be_at(&self, name: &str, offset: u64, i: usize, width: usize) -> Result<&[u8]> {
        let pos = offset as usize + i * width;
        if pos.saturating_add(width) > self.data.len() {
            bail!(
                "variable {name:?}: data at offset {pos} (+{width}) lies beyond end of file ({} bytes)",
                self.data.len()
            );
        }
        Ok(&self.data[pos..pos + width])
    }

    /// Read a float variable, record-major and row-major (C order) within each
    /// slab. See the module docs for the Fortran/C dimension-ordering note.
    pub fn read_f32(&self, name: &str) -> Result<Vec<f32>> {
        let var = self.var(name).ok_or_else(|| anyhow!("variable {name:?} not found"))?;
        if var.typ != NcType::Float {
            bail!("variable {name:?} is {} but float was requested", var.typ.name());
        }
        let per_slice = self.elems_per_slice(var) as usize;
        let slices = if var.is_record { self.numrecs } else { 1 };
        let mut out = Vec::with_capacity(self.total_elems(var) as usize);
        for s in 0..slices {
            let base = if var.is_record { var.begin + s * self.recsize } else { var.begin };
            for i in 0..per_slice {
                let bytes = self.read_be_at(name, base, i, 4)?;
                out.push(f32::from_be_bytes(bytes.try_into().unwrap()));
            }
        }
        Ok(out)
    }

    /// Read an int variable (e.g. `mesh2d_face_nodes`), same layout rules as
    /// [`NcFile::read_f32`]. Compared exactly in the regression tests.
    pub fn read_i32(&self, name: &str) -> Result<Vec<i32>> {
        let var = self.var(name).ok_or_else(|| anyhow!("variable {name:?} not found"))?;
        if var.typ != NcType::Int {
            bail!("variable {name:?} is {} but int was requested", var.typ.name());
        }
        let per_slice = self.elems_per_slice(var) as usize;
        let slices = if var.is_record { self.numrecs } else { 1 };
        let mut out = Vec::with_capacity(self.total_elems(var) as usize);
        for s in 0..slices {
            let base = if var.is_record { var.begin + s * self.recsize } else { var.begin };
            for i in 0..per_slice {
                let bytes = self.read_be_at(name, base, i, 4)?;
                out.push(i32::from_be_bytes(bytes.try_into().unwrap()));
            }
        }
        Ok(out)
    }
}
