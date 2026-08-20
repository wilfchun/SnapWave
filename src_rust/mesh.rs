//! Rust port of the UGRID mesh NetCDF reader `nc_read_net`
//! (plan.md, Phase 7, step 2).
//!
//! `nc_read_net` (in `src/snapwave_ncoutput.F90`) opens the `gridfile`,
//! detects the old (`mesh2d_nNodes`/`mesh2d_nFaces`/`mesh2d_nMax_face_nodes`)
//! versus new (`nmesh2d_node`/`nmesh2d_face`/`max_nmesh2d_face_nodes`)
//! dimension names, and reads node coordinates (`mesh2d_node_x/y/z`) and
//! face connectivity (`mesh2d_face_nodes`). `initialize_snapwave_domain`
//! then post-processes: `zb = -posdwn * zb`, and a missing fourth node of a
//! face (`face_nodes(4,:) == 0`) becomes the `-999` sentinel.
//!
//! [`read_ugrid_netcdf`] reproduces exactly that behaviour (minus the parts
//! that never affect the output: the commented-out `standard_name` check and
//! the `writetestfiles` ASCII dump). It is validated against the unchanged
//! Fortran reader through the temporary `snapwave_mesh_dump_c` hook (see
//! [`check`]), driven by the wrapper's `--compare-mesh` mode. The mesh data
//! itself is not yet handed to Fortran — that is the Phase 8 data-structure
//! handoff; for now the run path still lets Fortran read the mesh, and the
//! Rust reader exists to pin the port.

use std::path::Path;

use anyhow::{bail, Context, Result};

use crate::netcdf::NetcdfFile;

/// A UGRID mesh as read by `nc_read_net` + the domain post-processing.
#[derive(Debug, Clone, PartialEq)]
pub struct Mesh {
    pub no_nodes: usize,
    pub no_faces: usize,
    /// `max_nmesh2d_face_nodes` (3 for pure triangles, 4 for quads/mixed).
    pub max_nodes: usize,
    /// `sferic` after the `abs(y(1)) > 90` fix.
    pub sferic: i32,
    /// Node x coordinates (real*8; widened from the file's float values).
    pub x: Vec<f64>,
    /// Node y coordinates.
    pub y: Vec<f64>,
    /// Bed level after `zb = -posdwn * zb` (real*4).
    pub zb: Vec<f32>,
    /// Mask (always 1 after `nc_read_net`).
    pub msk: Vec<i32>,
    /// `face_nodes(4, no_faces)`, flattened node-major
    /// (`face_nodes[f*4 + j]` is node `j+1` of face `f+1`). A missing fourth
    /// node is `-999`.
    pub face_nodes: Vec<i32>,
}

/// Read a UGRID mesh NetCDF file, mirroring `nc_read_net` and the
/// `zb = -posdwn*zb` / fourth-node post-processing of
/// `initialize_snapwave_domain`.
pub fn read_ugrid_netcdf(path: &Path, posdwn: f32, sferic: i32) -> Result<Mesh> {
    let f = NetcdfFile::open(path)
        .with_context(|| format!("reading mesh NetCDF {}", path.display()))?;

    // Dimension-name detection: `mesh2d_nNodes` marks the old naming scheme.
    let (no_nodes, max_nodes, no_faces) = match f.dim("mesh2d_nNodes") {
        Some(n) => {
            let max = f.dim("mesh2d_nMax_face_nodes").ok_or_else(|| {
                anyhow_missing_dim(path, "mesh2d_nMax_face_nodes")
            })?;
            let nf = f.dim("mesh2d_nFaces").ok_or_else(|| anyhow_missing_dim(path, "mesh2d_nFaces"))?;
            (n, max, nf)
        }
        None => {
            let n = f.dim("nmesh2d_node").ok_or_else(|| anyhow_missing_dim(path, "nmesh2d_node"))?;
            let max = f
                .dim("max_nmesh2d_face_nodes")
                .ok_or_else(|| anyhow_missing_dim(path, "max_nmesh2d_face_nodes"))?;
            let nf = f.dim("nmesh2d_face").ok_or_else(|| anyhow_missing_dim(path, "nmesh2d_face"))?;
            (n, max, nf)
        }
    };
    let (no_nodes, max_nodes, no_faces) = (no_nodes as usize, max_nodes as usize, no_faces as usize);

    let x = f.read_f64("mesh2d_node_x").with_context(|| format!("{}: mesh2d_node_x", path.display()))?;
    let y = f.read_f64("mesh2d_node_y").with_context(|| format!("{}: mesh2d_node_y", path.display()))?;
    let zb_raw = f.read_f32("mesh2d_node_z").with_context(|| format!("{}: mesh2d_node_z", path.display()))?;
    let face_nodes_temp =
        f.read_i32("mesh2d_face_nodes").with_context(|| format!("{}: mesh2d_face_nodes", path.display()))?;

    if x.len() != no_nodes || y.len() != no_nodes || zb_raw.len() != no_nodes {
        bail!(
            "mesh {}: node variable lengths (x={}, y={}, z={}) do not match nmesh2d_node={}",
            path.display(),
            x.len(),
            y.len(),
            zb_raw.len(),
            no_nodes
        );
    }
    if face_nodes_temp.len() != max_nodes * no_faces {
        bail!(
            "mesh {}: mesh2d_face_nodes has {} elements, expected max_nmesh2d_face_nodes({}) * nmesh2d_face({}) = {}",
            path.display(),
            face_nodes_temp.len(),
            max_nodes,
            no_faces,
            max_nodes * no_faces
        );
    }

    // `x`/`y` are real*8 in snapwave_data; read_f64 already returns f64
    // (widening a float source, which is lossless).

    // `zb = -posdwn * zb` (both real*4).
    let zb: Vec<f32> = zb_raw.into_iter().map(|z| -posdwn * z).collect();

    // face_nodes(1:max_nodes,:) = face_nodes_temp; then -1 -> 0; then a
    // missing fourth node (0) becomes -999.
    let mut face_nodes = vec![0i32; 4 * no_faces];
    for k in 0..no_faces {
        for j in 0..max_nodes {
            face_nodes[k * 4 + j] = face_nodes_temp[k * max_nodes + j];
        }
    }
    for v in face_nodes.iter_mut() {
        if *v == -1 {
            *v = 0;
        }
    }
    for k in 0..no_faces {
        if face_nodes[k * 4 + 3] == 0 {
            face_nodes[k * 4 + 3] = -999;
        }
    }

    let sferic = if y.first().map(|v| v.abs() > 90.0).unwrap_or(false) { 0 } else { sferic };

    Ok(Mesh {
        no_nodes,
        no_faces,
        max_nodes,
        sferic,
        x,
        y,
        zb,
        msk: vec![1; no_nodes],
        face_nodes,
    })
}

fn anyhow_missing_dim(path: &Path, dim: &str) -> anyhow::Error {
    anyhow::anyhow!("mesh {}: missing dimension '{}'", path.display(), dim)
}

// ----------------------------------------------------------------------
// Comparison against the Fortran oracle (--compare-mesh)
// ----------------------------------------------------------------------

/// Parse the canonical mesh dump produced by `snapwave_mesh_dump_c` and
/// compare it against `rust`. Returns the number of values compared; errors
/// name every mismatching field.
pub fn check(rust: &Mesh, dump_text: &str) -> Result<usize> {
    let dump = parse_dump(dump_text)?;
    let mut mismatches: Vec<String> = Vec::new();
    let mut compared = 0usize;

    let no_nodes = expect_scalar(&dump, "no_nodes")?;
    let no_faces = expect_scalar(&dump, "no_faces")?;
    let max_nodes = expect_scalar(&dump, "max_nodes")?;
    let sferic = expect_scalar(&dump, "sferic")?;
    compared += 4;

    if rust.no_nodes != no_nodes {
        mismatches.push(format!("no_nodes: rust {} vs fortran {}", rust.no_nodes, no_nodes));
    }
    if rust.no_faces != no_faces {
        mismatches.push(format!("no_faces: rust {} vs fortran {}", rust.no_faces, no_faces));
    }
    if rust.max_nodes != max_nodes {
        mismatches.push(format!("max_nodes: rust {} vs fortran {}", rust.max_nodes, max_nodes));
    }
    if rust.sferic != sferic as i32 {
        mismatches.push(format!("sferic: rust {} vs fortran {}", rust.sferic, sferic));
    }

    // Only compare the deterministic rows 1..max_nodes of face_nodes: for a
    // pure-triangle mesh the Fortran fourth row is never written (uninitialized).
    let face_nodes_subset: Vec<i32> = (0..no_faces)
        .flat_map(|f| (0..max_nodes).map(move |j| rust.face_nodes[f * 4 + j]))
        .collect();

    if rust.x.len() == no_nodes {
        let fortran = dump_f64(&dump, "x")?;
        compared += cmp_f64(&mut mismatches, "x", &rust.x, &fortran);
    }
    if rust.y.len() == no_nodes {
        let fortran = dump_f64(&dump, "y")?;
        compared += cmp_f64(&mut mismatches, "y", &rust.y, &fortran);
    }
    if rust.zb.len() == no_nodes {
        let fortran = dump_f32(&dump, "zb")?;
        compared += cmp_f32(&mut mismatches, "zb", &rust.zb, &fortran);
    }
    if rust.msk.len() == no_nodes {
        let fortran = dump_i32(&dump, "msk")?;
        compared += cmp_i32(&mut mismatches, "msk", &rust.msk, &fortran);
    }
    {
        let fortran = dump_i32(&dump, "face_nodes")?;
        compared += cmp_i32(&mut mismatches, "face_nodes", &face_nodes_subset, &fortran);
    }

    if !mismatches.is_empty() {
        bail!(
            "Rust and Fortran mesh reads disagree ({} of {} values):\n  {}",
            mismatches.len(),
            compared,
            mismatches.join("\n  ")
        );
    }
    Ok(compared)
}

// ----------------------------------------------------------------------
// Dump parsing (same conventions as text_compare.rs, but flat — no sections)
// ----------------------------------------------------------------------

/// Keys whose line carries a count followed by that many value lines.
const ARRAY_KEYS: [&str; 5] = ["x", "y", "zb", "msk", "face_nodes"];

fn parse_dump(dump_text: &str) -> Result<std::collections::BTreeMap<String, Vec<String>>> {
    let mut map: std::collections::BTreeMap<String, Vec<String>> = std::collections::BTreeMap::new();
    let mut lines = dump_text.lines().peekable();
    while let Some(line) = lines.next() {
        let line = line.trim_end_matches('\r');
        if line.trim().is_empty() {
            continue;
        }
        let Some((key, rest)) = line.split_once(' ') else { continue };
        let key = key.trim().to_string();
        let rest = rest.trim().to_string();
        if ARRAY_KEYS.contains(&key.as_str()) {
            let count: usize = rest
                .parse()
                .with_context(|| format!("dump array '{key}' count '{rest}' is not an integer"))?;
            let mut vals = Vec::with_capacity(count);
            for _ in 0..count {
                let Some(v) = lines.next() else {
                    bail!("dump truncated in array '{key}'");
                };
                vals.push(v.trim_end_matches('\r').to_string());
            }
            map.insert(key, vals);
        } else {
            map.insert(key, vec![rest]);
        }
    }
    Ok(map)
}

fn expect_scalar(dump: &std::collections::BTreeMap<String, Vec<String>>, key: &str) -> Result<usize> {
    let vals = dump.get(key).ok_or_else(|| anyhow::anyhow!("dump is missing scalar '{key}'"))?;
    vals[0]
        .parse::<usize>()
        .with_context(|| format!("dump scalar '{key}' = '{}' is not an integer", vals[0]))
}

fn dump_values<'a>(
    dump: &'a std::collections::BTreeMap<String, Vec<String>>,
    key: &str,
) -> Result<&'a [String]> {
    dump.get(key)
        .map(|v| v.as_slice())
        .ok_or_else(|| anyhow::anyhow!("dump is missing array '{key}'"))
}

fn dump_f64(dump: &std::collections::BTreeMap<String, Vec<String>>, key: &str) -> Result<Vec<f64>> {
    dump_values(dump, key)?
        .iter()
        .map(|v| {
            u64::from_str_radix(v.trim(), 16)
                .map(f64::from_bits)
                .with_context(|| format!("dump array '{key}' value '{v}' is not a real*8 bit pattern"))
        })
        .collect()
}

fn dump_f32(dump: &std::collections::BTreeMap<String, Vec<String>>, key: &str) -> Result<Vec<f32>> {
    dump_values(dump, key)?
        .iter()
        .map(|v| {
            u32::from_str_radix(v.trim(), 16)
                .map(f32::from_bits)
                .with_context(|| format!("dump array '{key}' value '{v}' is not a real*4 bit pattern"))
        })
        .collect()
}

fn dump_i32(dump: &std::collections::BTreeMap<String, Vec<String>>, key: &str) -> Result<Vec<i32>> {
    dump_values(dump, key)?
        .iter()
        .map(|v| {
            v.trim()
                .parse::<i32>()
                .with_context(|| format!("dump array '{key}' value '{v}' is not an integer"))
        })
        .collect()
}

fn cmp_f64(mismatches: &mut Vec<String>, label: &str, rust: &[f64], fortran: &[f64]) -> usize {
    let n = rust.len().min(fortran.len());
    let mut compared = 0;
    for i in 0..n {
        compared += 1;
        if rust[i].to_bits() != fortran[i].to_bits() {
            mismatches.push(format!(
                "{label}[{i}]: rust {} (0x{:016x}) vs fortran {} (0x{:016x})",
                rust[i],
                rust[i].to_bits(),
                fortran[i],
                fortran[i].to_bits()
            ));
        }
    }
    if rust.len() != fortran.len() {
        mismatches.push(format!("{label}: length {} vs {}", rust.len(), fortran.len()));
    }
    compared
}

fn cmp_f32(mismatches: &mut Vec<String>, label: &str, rust: &[f32], fortran: &[f32]) -> usize {
    let n = rust.len().min(fortran.len());
    let mut compared = 0;
    for i in 0..n {
        compared += 1;
        if rust[i].to_bits() != fortran[i].to_bits() {
            mismatches.push(format!(
                "{label}[{i}]: rust {} (0x{:08x}) vs fortran {} (0x{:08x})",
                rust[i],
                rust[i].to_bits(),
                fortran[i],
                fortran[i].to_bits()
            ));
        }
    }
    if rust.len() != fortran.len() {
        mismatches.push(format!("{label}: length {} vs {}", rust.len(), fortran.len()));
    }
    compared
}

fn cmp_i32(mismatches: &mut Vec<String>, label: &str, rust: &[i32], fortran: &[i32]) -> usize {
    let n = rust.len().min(fortran.len());
    let mut compared = 0;
    for i in 0..n {
        compared += 1;
        if rust[i] != fortran[i] {
            mismatches.push(format!("{label}[{i}]: rust {} vs fortran {}", rust[i], fortran[i]));
        }
    }
    if rust.len() != fortran.len() {
        mismatches.push(format!("{label}: length {} vs {}", rust.len(), fortran.len()));
    }
    compared
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::netcdf::{f32_payload, i32_payload, Attr, NcType, Var, VarData, Writer};

    /// Build a tiny triangle UGRID mesh in memory and read it back, checking
    /// the dimension detection, coordinate widening, `zb` sign flip and the
    /// `-1 -> 0 -> -999` fourth-node chain.
    #[test]
    fn reads_a_triangle_mesh() {
        let dir = std::env::temp_dir().join(format!("snapwave_mesh_unit_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("create scratch dir");
        let path = dir.join("mesh.nc");

        let mut w = Writer::new();
        w.dim("nmesh2d_node", Some(4));
        w.dim("nmesh2d_face", Some(2));
        w.dim("max_nmesh2d_face_nodes", Some(3));
        w.var(Var {
            name: "mesh2d_node_x".into(),
            dims: vec![0],
            typ: NcType::Float,
            attrs: vec![Attr::text("standard_name", "longitude")],
            data: VarData::Fixed(f32_payload(&[0.0, 10.0, 0.0, 10.0])),
        });
        w.var(Var {
            name: "mesh2d_node_y".into(),
            dims: vec![0],
            typ: NcType::Float,
            attrs: vec![],
            data: VarData::Fixed(f32_payload(&[0.0, 0.0, 10.0, 10.0])),
        });
        w.var(Var {
            name: "mesh2d_node_z".into(),
            dims: vec![0],
            typ: NcType::Float,
            attrs: vec![],
            data: VarData::Fixed(f32_payload(&[-5.0, -4.0, -3.0, -2.0])),
        });
        w.var(Var {
            name: "mesh2d_face_nodes".into(),
            dims: vec![1, 2], // C order: (nmesh2d_face, max_nmesh2d_face_nodes)
            typ: NcType::Int,
            attrs: vec![Attr::int("_FillValue", -999)],
            data: VarData::Fixed(i32_payload(&[1, 2, 3, 1, 3, 4])),
        });
        w.write_to(&path).expect("write mesh");

        // posdwn = 1 flips zb sign (zb = -posdwn*zb), exercising the flip
        // logic with a non-trivial multiplier. sferic stays 0 because
        // |y[0]| = 0 (not > 90).
        let m = read_ugrid_netcdf(&path, 1.0, 0).expect("read mesh");
        assert_eq!(m.no_nodes, 4);
        assert_eq!(m.no_faces, 2);
        assert_eq!(m.max_nodes, 3);
        assert_eq!(m.sferic, 0);
        assert_eq!(m.x, vec![0.0, 10.0, 0.0, 10.0]);
        assert_eq!(m.y, vec![0.0, 0.0, 10.0, 10.0]);
        assert_eq!(m.zb, vec![5.0, 4.0, 3.0, 2.0]);
        assert_eq!(m.msk, vec![1, 1, 1, 1]);
        // Face 0: nodes 1,2,3 + missing fourth -> -999; face 1: 1,3,4 + -999.
        assert_eq!(m.face_nodes, vec![1, 2, 3, -999, 1, 3, 4, -999]);

        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn sferic_reset_for_geographic_coordinates() {
        let dir = std::env::temp_dir().join(format!("snapwave_mesh_sferic_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("create scratch dir");

        // One node at |y| > 90 must reset sferic to 0.
        let path = dir.join("mesh.nc");
        let mut w = Writer::new();
        w.dim("nmesh2d_node", Some(1));
        w.dim("nmesh2d_face", Some(0));
        w.dim("max_nmesh2d_face_nodes", Some(3));
        w.var(Var {
            name: "mesh2d_node_x".into(),
            dims: vec![0],
            typ: NcType::Float,
            attrs: vec![],
            data: VarData::Fixed(f32_payload(&[120.0])),
        });
        w.var(Var {
            name: "mesh2d_node_y".into(),
            dims: vec![0],
            typ: NcType::Float,
            attrs: vec![],
            data: VarData::Fixed(f32_payload(&[95.0])),
        });
        w.var(Var {
            name: "mesh2d_node_z".into(),
            dims: vec![0],
            typ: NcType::Float,
            attrs: vec![],
            data: VarData::Fixed(f32_payload(&[-1.0])),
        });
        w.var(Var {
            name: "mesh2d_face_nodes".into(),
            dims: vec![1, 2],
            typ: NcType::Int,
            attrs: vec![],
            data: VarData::Fixed(vec![]),
        });
        w.write_to(&path).expect("write mesh");

        let m = read_ugrid_netcdf(&path, -1.0, 1).expect("read mesh");
        assert_eq!(m.sferic, 0, "|y| > 90 must reset sferic to 0");

        // A node at |y| <= 90 keeps the configured sferic.
        let path2 = dir.join("mesh2.nc");
        let mut w = Writer::new();
        w.dim("nmesh2d_node", Some(1));
        w.dim("nmesh2d_face", Some(0));
        w.dim("max_nmesh2d_face_nodes", Some(3));
        w.var(Var {
            name: "mesh2d_node_x".into(),
            dims: vec![0],
            typ: NcType::Float,
            attrs: vec![],
            data: VarData::Fixed(f32_payload(&[120.0])),
        });
        w.var(Var {
            name: "mesh2d_node_y".into(),
            dims: vec![0],
            typ: NcType::Float,
            attrs: vec![],
            data: VarData::Fixed(f32_payload(&[20.0])),
        });
        w.var(Var {
            name: "mesh2d_node_z".into(),
            dims: vec![0],
            typ: NcType::Float,
            attrs: vec![],
            data: VarData::Fixed(f32_payload(&[-1.0])),
        });
        w.var(Var {
            name: "mesh2d_face_nodes".into(),
            dims: vec![1, 2],
            typ: NcType::Int,
            attrs: vec![],
            data: VarData::Fixed(vec![]),
        });
        w.write_to(&path2).expect("write mesh");

        let m = read_ugrid_netcdf(&path2, -1.0, 1).expect("read mesh");
        assert_eq!(m.sferic, 1, "|y| <= 90 keeps the configured sferic");

        let _ = std::fs::remove_dir_all(&dir);
    }
}
