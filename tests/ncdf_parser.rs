//! Self-tests for the bundled classic-NetCDF reader (`tests/support/ncdf.rs`)
//! against committed testcase fixtures.
//!
//! These run without executing the model: they validate that the Phase-1
//! reader parses real SnapWave files correctly (which the numeric regression
//! in `regression.rs` relies on) and they pin the schema of the committed
//! baseline outputs themselves.

mod support;

use std::path::PathBuf;

use support::ncdf::{NcFile, NcType};

fn fixture(rel: &str) -> PathBuf {
    let p = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(rel);
    assert!(p.is_file(), "committed fixture missing: {}", p.display());
    p
}

/// Treat values at or below ~-1e6 as FILL_VALUE (-999999) writes.
fn is_fill(v: f32) -> bool {
    v <= -9.99e5
}

#[test]
fn parses_committed_31_his_fixture() {
    let f = NcFile::open(&fixture(
        "testcases/31_linear_shoaling_refraction/output/shoalref_coarse_neu_his.nc",
    ))
    .expect("parse committed 31 his fixture");

    assert_eq!(f.dim("stations").expect("stations dim").len, 201);
    assert_eq!(f.dim("pointnamelength").expect("pointnamelength dim").len, 256);
    assert_eq!(f.dim("runtime").expect("runtime dim").len, 1);
    let time = f.dim("time").expect("time dim");
    assert!(time.unlimited, "time must be the record dimension");
    assert!(time.len >= 1, "expected at least one history frame");

    // Fortran declares (stations, time); on disk (C order) this is
    // (time, stations).
    assert_eq!(f.var_dim_names("point_hm0").unwrap(), vec!["time", "stations"]);
    assert_eq!(f.var("point_hm0").unwrap().typ, NcType::Float);
    assert_eq!(f.var_dim_names("point_wavdir").unwrap(), vec!["time", "stations"]);
    assert_eq!(f.var_dim_names("time").unwrap(), vec!["time"]);

    // Station naming is a char variable pinned by schema only.
    let station_name = f.var("station_name").expect("station_name var");
    assert_eq!(station_name.typ, NcType::Char);
    assert_eq!(f.var_dim_names("station_name").unwrap(), vec!["stations", "pointnamelength"]);

    // Fill value attribute from snapwave_data.f90.
    let hm0 = f.var("point_hm0").unwrap();
    let fill = hm0
        .attrs
        .iter()
        .find(|a| a.name == "_FillValue")
        .expect("point_hm0 _FillValue attr");
    match &fill.value {
        support::ncdf::NcAttrValue::Floats(v) => assert_eq!(v, &vec![-999999.0f32]),
        other => panic!("point_hm0 _FillValue has unexpected payload: {other:?}"),
    }

    // Numeric access: full time series is readable and mostly non-fill.
    let vals = f.read_f32("point_hm0").expect("read point_hm0");
    assert_eq!(vals.len() as u64, time.len * 201);
    assert!(
        vals.iter().copied().any(|v| !is_fill(v)),
        "point_hm0 should contain real values"
    );
    let times = f.read_f32("time").expect("read time");
    assert_eq!(times.len() as u64, f.record_count());
}

#[test]
fn parses_committed_32_map_fixture() {
    let f = NcFile::open(&fixture("testcases/32_curvi_island/output/snapwave_map.nc"))
        .expect("parse committed 32 map fixture");

    let nodes = f.dim("nmesh2d_node").expect("nmesh2d_node dim").len;
    let faces = f.dim("nmesh2d_face").expect("nmesh2d_face dim").len;
    let max_face_nodes = f.dim("max_nmesh2d_face_nodes").expect("max_nmesh2d_face_nodes dim").len;
    let ntheta = f.dim("ntheta").expect("ntheta dim").len;
    let time = f.dim("time").expect("time dim");
    assert!(nodes > 0 && faces > 0 && ntheta > 0);
    assert!(time.unlimited);
    assert!(f.record_count() >= 1);

    // Directional energy density: Fortran (ntheta, nmesh2d_node, time)
    // -> C order (time, nmesh2d_node, ntheta).
    assert_eq!(f.var_dim_names("ee").unwrap(), vec!["time", "nmesh2d_node", "ntheta"]);
    assert_eq!(f.var("ee").unwrap().typ, NcType::Float);
    assert_eq!(f.var_dim_names("hm0").unwrap(), vec!["time", "nmesh2d_node"]);

    // Face-node connectivity is int with 1-based indices.
    let face_nodes = f.var("mesh2d_face_nodes").expect("mesh2d_face_nodes var");
    assert_eq!(face_nodes.typ, NcType::Int);
    assert_eq!(
        f.var_dim_names("mesh2d_face_nodes").unwrap(),
        vec!["nmesh2d_face", "max_nmesh2d_face_nodes"]
    );

    // Static mesh geometry is a non-record float variable.
    assert_eq!(f.var_dim_names("mesh2d_node_x").unwrap(), vec!["nmesh2d_node"]);
    let x = f.read_f32("mesh2d_node_x").expect("read mesh2d_node_x");
    assert_eq!(x.len() as u64, nodes);
    let fnvals = f.read_i32("mesh2d_face_nodes").expect("read mesh2d_face_nodes");
    assert_eq!(fnvals.len() as u64, faces * max_face_nodes);

    // Record series layout is consistent with dims and record count.
    let ee = f.read_f32("ee").expect("read ee");
    assert_eq!(ee.len() as u64, f.record_count() * nodes * ntheta);
    let times = f.read_f32("time").expect("read time");
    assert_eq!(times.len() as u64, f.record_count());

    // Global metadata written by ncoutput_map_init.
    let title = f
        .global_attrs
        .iter()
        .find(|a| a.name == "title")
        .expect("global title attr");
    match &title.value {
        support::ncdf::NcAttrValue::Text(t) => assert_eq!(t.trim_end(), "SnapWave map netcdf output"),
        other => panic!("title attr has unexpected payload: {other:?}"),
    }
    assert!(f.global_attrs.iter().any(|a| a.name == "Conventions"));
}
