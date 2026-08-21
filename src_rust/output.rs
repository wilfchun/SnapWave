//! Rust NetCDF map/history writers (plan.md, Phase 7, step 3).
//!
//! These reproduce — dimension names, variable names, attribute strings,
//! fill values, ordering and time indexing — the schema that
//! `src/snapwave_ncoutput.F90` writes. The *data* comes from the capture
//! stream ([`crate::capture`]): the Fortran solver computes the same
//! buffers it would have handed to `nf90_put_var`, so the Rust files are
//! byte-identical in every variable that SnapWave actually writes.
//!
//! The writers are pure "write what Fortran computed": every value,
//! including the `where (depth < hmin)` fill masking and the
//! `modulo(270 - …*rad2deg, 360.)` direction wrapping, is captured verbatim
//! rather than recomputed here (keeping Fortran the numerical authority,
//! AGENTS.md rule 1). Values the Fortran writer *never* writes are emitted
//! exactly as netCDF would leave them: variables with no `_FillValue`
//! attribute (`mesh2d`, `crs`, `station_id`, `total_runtime`, `average_dt`)
//! carry the library default fill value, while `point_zb` — which Fortran
//! defines with `_FillValue = FILL_VALUE` — carries SnapWave's `FILL_VALUE`
//! (-999999.0).

use std::path::Path;

use anyhow::{anyhow, Result};

use crate::capture::{HisRecord, MapRecord, StaticHis, StaticMap};
use crate::input::SnapWaveInput;
use crate::netcdf::{
    f32_payload, f32_single, f64_as_f32_payload, i32_payload, i32_single, Attr, NcType, NC_FILL_FLOAT,
    NC_FILL_INT, Var, VarData, Writer,
};

/// `FILL_VALUE` of `snapwave_data.f90` (`-999999.0`, real*4).
const FILL_VALUE: f32 = -999_999.0;

/// The `Build-Revision-Date-Netcdf-library` attribute value. The Fortran
/// writer embeds `nf90_inq_libvers()`; the pure-Rust writer has no NetCDF
/// library, so it records the hand-rolled classic-format writer instead.
/// The regression harness excludes this attribute from schema comparison
/// (it legitimately differs between machines).
pub const LIBVERS: &str = "Rust classic-format writer (no NetCDF library)";

// ----------------------------------------------------------------------
// Shared helpers
// ----------------------------------------------------------------------

fn fill() -> Attr {
    Attr::float("_FillValue", FILL_VALUE)
}

/// Extract one optional per-record field into record slabs; errors if the
/// field is absent when the schema says it must be present.
fn map_slabs(
    recs: &[MapRecord],
    name: &str,
    get: impl Fn(&MapRecord) -> Option<&[f32]>,
) -> Result<Vec<Vec<u8>>> {
    recs.iter()
        .map(|r| {
            let v = get(r).ok_or_else(|| anyhow!("map record missing field '{name}'"))?;
            Ok(f32_payload(v))
        })
        .collect()
}

/// Build a record float variable whose slabs come from a required per-record
/// field, and push it onto `w`.
fn push_map_record_var(
    w: &mut Writer,
    name: &str,
    dims: Vec<usize>,
    attrs: Vec<Attr>,
    recs: &[MapRecord],
    get: impl Fn(&MapRecord) -> Option<&[f32]>,
) -> Result<()> {
    let slabs = map_slabs(recs, name, get)?;
    w.var(Var { name: name.into(), dims, typ: NcType::Float, attrs, data: VarData::Record(slabs) });
    Ok(())
}

/// Build record slabs for a required per-record `Vec<f32>` field.
fn his_slabs(recs: &[HisRecord], get: impl Fn(&HisRecord) -> &Vec<f32>) -> Vec<Vec<u8>> {
    recs.iter().map(|r| f32_payload(get(r))).collect()
}

// ----------------------------------------------------------------------
// Map output
// ----------------------------------------------------------------------

/// Write the map NetCDF file from the captured static mesh and map records.
pub fn write_map(path: &Path, cfg: &SnapWaveInput, sm: &StaticMap, recs: &[MapRecord]) -> Result<()> {
    let o = &cfg.output;
    let wind = cfg.wind.enabled;
    let ig = cfg.physics.ig == 1;
    let veg = cfg.vegetation.ja_vegetation == 1;

    let mut w = Writer::new();
    let node = w.dim("nmesh2d_node", Some(sm.no_nodes as u64));
    let face = w.dim("nmesh2d_face", Some(sm.no_faces as u64));
    let maxface = w.dim("max_nmesh2d_face_nodes", Some(4)); // hard-coded 4, not max_nodes
    let ntheta = w.dim("ntheta", Some(sm.ntheta as u64));
    let time = w.dim("time", None);

    // The Fortran writer's string has an opening quote but no closing quote
    // (a literal typo in ncoutput_map_init); reproduced verbatim so the
    // Rust-written file matches the oracle byte-for-byte.
    w.global_attr(Attr::text("Conventions", "Conventions = 'CF-1.6, SGRID-0.3"));
    w.global_attr(Attr::text("Build-Revision-Date-Netcdf-library", &sm.libvers));
    w.global_attr(Attr::text("Producer", "SnapWave"));
    w.global_attr(Attr::text("title", "SnapWave map netcdf output"));

    // mesh2d (topology dummy, never written)
    w.var(Var {
        name: "mesh2d".into(),
        dims: vec![],
        typ: NcType::Int,
        attrs: vec![
            Attr::text("cf_role", "mesh_topology"),
            Attr::text("long_name", "Topology data of 2D network"),
            Attr::int("topology_dimension", 2),
            Attr::text("node_coordinates", "mesh2d_node_x mesh2d_node_y"),
            Attr::text("node_dimension", "nmesh2d_node"),
            Attr::text("face_node_connectivity", "mesh2d_face_nodes"),
            Attr::text("face_dimension", "nmesh2d_face"),
            Attr::text("max_face_nodes_dimension", "max_nmesh2d_face_nodes"),
        ],
        data: VarData::Fixed(i32_single(NC_FILL_INT)),
    });

    // mesh2d_face_nodes: rows 1..max_nodes per face, padded with -999 to 4.
    let mut face_nodes_data = Vec::with_capacity(4 * sm.no_faces);
    for f in 0..sm.no_faces {
        for j in 0..4 {
            face_nodes_data.push(if j < sm.max_nodes { sm.face_nodes[f * sm.max_nodes + j] } else { -999 });
        }
    }
    w.var(Var {
        name: "mesh2d_face_nodes".into(),
        dims: vec![face, maxface],
        typ: NcType::Int,
        attrs: vec![
            Attr::text("cf_role", "face_node_connectivity"),
            Attr::text("mesh", "mesh2d"),
            Attr::text("location", "face"),
            Attr::text("long_name", "Mapping from every face to its corner nodes (counterclockwise)"),
            Attr::int("start_index", 1),
            Attr::int("_FillValue", -999),
        ],
        data: VarData::Fixed(i32_payload(&face_nodes_data)),
    });

    // mesh2d_node_x / mesh2d_node_y: units/standard_name switch on sferic.
    let (x_units, x_std) = if sm.sferic == 1 { ("degrees_east", "longitude") } else { ("m", "projection_x_coordinate") };
    let (y_units, y_std) = if sm.sferic == 1 { ("degrees_north", "latitude") } else { ("m", "projection_y_coordinate") };
    w.var(Var {
        name: "mesh2d_node_x".into(),
        dims: vec![node],
        typ: NcType::Float,
        attrs: vec![
            fill(),
            Attr::text("units", x_units),
            Attr::text("standard_name", x_std),
            Attr::text("long_name", "x-coordinate of mesh nodes"),
            Attr::text("location", "node"),
            Attr::text("mesh", "mesh2d"),
        ],
        data: VarData::Fixed(f64_as_f32_payload(&sm.x)),
    });
    w.var(Var {
        name: "mesh2d_node_y".into(),
        dims: vec![node],
        typ: NcType::Float,
        attrs: vec![
            fill(),
            Attr::text("units", y_units),
            Attr::text("standard_name", y_std),
            Attr::text("long_name", "y-coordinate of mesh nodes"),
            Attr::text("location", "node"),
            Attr::text("mesh", "mesh2d"),
        ],
        data: VarData::Fixed(f64_as_f32_payload(&sm.y)),
    });
    w.var(Var {
        name: "mesh2d_node_z".into(),
        dims: vec![node],
        typ: NcType::Float,
        attrs: vec![
            fill(),
            Attr::text("units", "m"),
            Attr::text("standard_name", "projection_z_coordinate"),
            Attr::text("long_name", "z-coordinate of mesh nodes"),
            Attr::text("location", "node"),
            Attr::text("mesh", "mesh2d"),
        ],
        data: VarData::Fixed(f32_payload(&sm.zb)),
    });

    // crs (EPSG code, never written)
    w.var(Var {
        name: "crs".into(),
        dims: vec![],
        typ: NcType::Int,
        attrs: vec![Attr::text("EPSG", "-")],
        data: VarData::Fixed(i32_single(NC_FILL_INT)),
    });

    // time
    let units = format!("seconds since {}", sm.tref_iso8601);
    let long_name = format!("time_in_seconds_since_{}", sm.tref_iso8601);
    w.var(Var {
        name: "time".into(),
        dims: vec![time],
        typ: NcType::Float,
        attrs: vec![
            Attr::text("units", &units),
            Attr::text("standard_name", "time"),
            Attr::text("long_name", &long_name),
        ],
        data: VarData::Record(recs.iter().map(|r| f32_single(r.t as f32)).collect()),
    });

    // Time-varying map variables, in ncoutput_map_init definition order.
    if o.map_depth == 1 {
        push_map_record_var(&mut w, "depth", vec![time, node], depth_attrs(), recs, |r| r.depth.as_deref())?;
    }
    if o.map_Hm0 == 1 {
        push_map_record_var(
            &mut w,
            "hm0",
            vec![time, node],
            map_attrs("m", "sea_surface_wind_wave_significant_height", "Wave height Hm0"),
            recs,
            |r| r.hm0.as_deref(),
        )?;
    }
    if ig && o.map_Hig == 1 {
        push_map_record_var(
            &mut w,
            "hm0_ig",
            vec![time, node],
            map_attrs("m", "sea_surface_infragravity_wave_significant_height", "Infragravity wave height Hm0ig"),
            recs,
            |r| r.hm0_ig.as_deref(),
        )?;
    }
    if o.map_Tp == 1 {
        push_map_record_var(
            &mut w,
            "tp",
            vec![time, node],
            map_attrs("s", "sea_surface_wave_period_at_variance_spectral_density_maximum", "Peak period Tp"),
            recs,
            |r| r.tp.as_deref(),
        )?;
    }
    if o.map_dir == 1 {
        push_map_record_var(
            &mut w,
            "wd",
            vec![time, node],
            map_attrs("degree", "sea_surface_wave_from_direction", "Mean wave from direction"),
            recs,
            |r| r.wd.as_deref(),
        )?;
    }
    if o.map_dirspr == 1 {
        push_map_record_var(
            &mut w,
            "wdspr",
            vec![time, node],
            map_attrs("degree", "sea_surface_wave_directional_spread", "Wave mean directional spread"),
            recs,
            |r| r.wdspr.as_deref(),
        )?;
    }
    if o.map_cg == 1 {
        push_map_record_var(
            &mut w,
            "cg",
            vec![time, node],
            map_attrs("m/s", "", "Wave group velocity"),
            recs,
            |r| r.cg.as_deref(),
        )?;
    }
    if o.map_Dw == 1 {
        push_map_record_var(
            &mut w,
            "dw",
            vec![time, node],
            map_attrs("W m-2", "", "Depth-induced wave breaking dissipation"),
            recs,
            |r| r.dw.as_deref(),
        )?;
    }
    if o.map_Df == 1 {
        push_map_record_var(
            &mut w,
            "df",
            vec![time, node],
            map_attrs("W m-2", "", "Bottom friction dissipation"),
            recs,
            |r| r.df.as_deref(),
        )?;
    }
    if wind && o.map_SwE == 1 {
        push_map_record_var(
            &mut w,
            "SwE",
            vec![time, node],
            map_attrs("W m-2", "", "Wind input short wave energy"),
            recs,
            |r| r.sw.as_deref(),
        )?;
    }
    if wind && o.map_SwA == 1 {
        push_map_record_var(
            &mut w,
            "SwA",
            vec![time, node],
            map_attrs("W m-2", "", "Wind input short wave action"),
            recs,
            |r| r.st.as_deref(),
        )?;
    }
    if wind && o.map_sig == 1 {
        push_map_record_var(
            &mut w,
            "sig",
            vec![time, node],
            map_attrs("Hz", "wave_frequency", "Relative wave frequency"),
            recs,
            |r| r.sig.as_deref(),
        )?;
    }
    if wind && o.map_u10 == 1 {
        push_map_record_var(
            &mut w,
            "u10",
            vec![time, node],
            map_attrs("m/s", "wind_speed", "Wind speed"),
            recs,
            |r| r.u10.as_deref(),
        )?;
        push_map_record_var(
            &mut w,
            "u10dir",
            vec![time, node],
            map_attrs("degree", "wind_from_direction", "Wind from direction"),
            recs,
            |r| r.u10dir.as_deref(),
        )?;
    }
    if veg && o.map_Dveg == 1 {
        push_map_record_var(
            &mut w,
            "mesh2d_veg_Dveg",
            vec![time, node],
            map_attrs("J/m2", "vegetation_dissipation", "Short wave dissipation by vegetation"),
            recs,
            |r| r.dveg.as_deref(),
        )?;
    }
    if o.map_ee == 1 {
        push_map_record_var(
            &mut w,
            "ee",
            vec![time, node, ntheta],
            map_attrs("J/m2/rad", "", "Wave energy density"),
            recs,
            |r| r.ee.as_deref(),
        )?;
    }
    if o.map_ee == 1 || o.map_ctheta == 1 {
        push_map_record_var(
            &mut w,
            "theta",
            vec![time, ntheta],
            vec![
                Attr::text("long_name", "Wave directional grid "),
                Attr::int("start_index", 1),
                fill(),
            ],
            recs,
            |r| r.theta_deg.as_deref(),
        )?;
    }
    if o.map_ctheta == 1 {
        push_map_record_var(
            &mut w,
            "ctheta",
            vec![time, node, ntheta],
            map_attrs("rad/s", "", "Wave refraction speed"),
            recs,
            |r| r.ctheta.as_deref(),
        )?;
    }

    // Static friction factors.
    w.var(Var {
        name: "fw".into(),
        dims: vec![node],
        typ: NcType::Float,
        attrs: vec![fill(), Attr::text("units", "-"), Attr::text("standard_name", ""), Attr::text("long_name", "Short wave friction factor")],
        data: VarData::Fixed(f32_payload(&sm.fw)),
    });
    w.var(Var {
        name: "fw_ig".into(),
        dims: vec![node],
        typ: NcType::Float,
        attrs: vec![fill(), Attr::text("units", "-"), Attr::text("standard_name", ""), Attr::text("long_name", "IG wave friction factor")],
        data: VarData::Fixed(f32_payload(&sm.fw_ig)),
    });

    // Static vegetation parameters.
    if veg {
        if let Some(v) = &sm.veg {
            w.var(Var {
                name: "mesh2d_veg_ah".into(),
                dims: vec![node],
                typ: NcType::Float,
                attrs: veg_attrs("m", "vegetation height", "Height of vegetation at mesh nodes"),
                data: VarData::Fixed(f32_payload(&v.ah)),
            });
            w.var(Var {
                name: "mesh2d_veg_bstems".into(),
                dims: vec![node],
                typ: NcType::Float,
                attrs: veg_attrs("m", "vegetation width", "Width of vegetation at mesh nodes"),
                data: VarData::Fixed(f32_payload(&v.bstems)),
            });
            w.var(Var {
                name: "mesh2d_veg_Nstems".into(),
                dims: vec![node],
                typ: NcType::Float,
                attrs: veg_attrs("plants/m2", "vegetation density", "Density of vegetation at mesh nodes"),
                data: VarData::Fixed(f32_payload(&v.nstems)),
            });
        }
    }

    w.write_to(path)
}

fn depth_attrs() -> Vec<Attr> {
    vec![
        fill(),
        Attr::text("units", "m"),
        Attr::text("standard_name", "sea_floor_depth_below_sea_surface"),
        Attr::text("long_name", "Water depth"),
    ]
}

fn map_attrs(units: &str, standard_name: &str, long_name: &str) -> Vec<Attr> {
    vec![fill(), Attr::text("units", units), Attr::text("standard_name", standard_name), Attr::text("long_name", long_name)]
}

fn veg_attrs(units: &str, standard_name: &str, long_name: &str) -> Vec<Attr> {
    vec![
        fill(),
        Attr::text("units", units),
        Attr::text("standard_name", standard_name),
        Attr::text("long_name", long_name),
        Attr::text("location", "node"),
        Attr::text("mesh", "mesh2d"),
    ]
}

// ----------------------------------------------------------------------
// History output
// ----------------------------------------------------------------------

/// Write the history NetCDF file from the captured observation points and
/// history records.
pub fn write_his(path: &Path, cfg: &SnapWaveInput, sh: &StaticHis, recs: &[HisRecord]) -> Result<()> {
    let wind = cfg.wind.enabled;
    let ig = cfg.physics.ig == 1;

    let mut w = Writer::new();
    let time = w.dim("time", None);
    let stations = w.dim("stations", Some(sh.nobs as u64));
    let namelen = w.dim("pointnamelength", Some(256));
    let runtime = w.dim("runtime", Some(1));

    // Same typo as the map writer: opening quote, no closing quote.
    w.global_attr(Attr::text("Conventions", "Conventions = 'CF-1.6, SGRID-0.3"));
    w.global_attr(Attr::text("Build-Revision-Date-Netcdf-library", &sh.libvers));
    w.global_attr(Attr::text("Producer", "SnapWave"));
    w.global_attr(Attr::text("title", "Snapwave his point netcdf output"));

    // station_id (never written)
    w.var(Var {
        name: "station_id".into(),
        dims: vec![stations],
        typ: NcType::Float,
        attrs: vec![],
        data: VarData::Fixed(f32_payload(&vec![NC_FILL_FLOAT; sh.nobs])),
    });

    // station_name: character*32 names padded to pointnamelength (256).
    let mut names = Vec::with_capacity(sh.nobs * 256);
    for name in &sh.names {
        let mut buf = name.as_bytes().to_vec();
        buf.resize(256, b' ');
        names.extend_from_slice(&buf);
    }
    w.var(Var { name: "station_name".into(), dims: vec![stations, namelen], typ: NcType::Char, attrs: vec![], data: VarData::Fixed(names) });

    w.var(Var {
        name: "station_x".into(),
        dims: vec![stations],
        typ: NcType::Float,
        attrs: vec![
            Attr::text("units", "m"),
            Attr::text("standard_name", "projection_x_coordinate"),
            Attr::text("long_name", "original_x_coordinate_of_station"),
            Attr::text("grid_mapping", "crs"),
        ],
        data: VarData::Fixed(f64_as_f32_payload(&sh.xobs)),
    });
    w.var(Var {
        name: "station_y".into(),
        dims: vec![stations],
        typ: NcType::Float,
        attrs: vec![
            Attr::text("units", "m"),
            Attr::text("standard_name", "projection_y_coordinate"),
            Attr::text("long_name", "original_y_coordinate_of_station"),
            Attr::text("grid_mapping", "crs"),
        ],
        data: VarData::Fixed(f64_as_f32_payload(&sh.yobs)),
    });

    w.var(Var {
        name: "crs".into(),
        dims: vec![],
        typ: NcType::Int,
        attrs: vec![Attr::text("EPSG", "-")],
        data: VarData::Fixed(i32_single(NC_FILL_INT)),
    });

    // point_zb has _FillValue = FILL_VALUE and is never written by Fortran,
    // so netCDF pre-fills it with FILL_VALUE (-999999.0) — NOT the default
    // float fill value.
    w.var(Var {
        name: "point_zb".into(),
        dims: vec![stations],
        typ: NcType::Float,
        attrs: vec![
            fill(),
            Attr::text("units", "m"),
            Attr::text("standard_name", "altitude"),
            Attr::text("long_name", "bed_level_above_reference_level"),
        ],
        data: VarData::Fixed(f32_payload(&vec![FILL_VALUE; sh.nobs])),
    });

    let units = format!("seconds since {}", sh.tref_iso8601);
    let long_name = format!("time_in_seconds_since_{}", sh.tref_iso8601);
    w.var(Var {
        name: "time".into(),
        dims: vec![time],
        typ: NcType::Float,
        attrs: vec![
            Attr::text("units", &units),
            Attr::text("standard_name", "time"),
            Attr::text("long_name", &long_name),
        ],
        data: VarData::Record(recs.iter().map(|r| f32_single(r.t as f32)).collect()),
    });

    let his_attrs = |units: &str, std: &str, long: &str| {
        vec![fill(), Attr::text("units", units), Attr::text("standard_name", std), Attr::text("long_name", long)]
    };
    let his_var = |w: &mut Writer, name: &str, attrs: Vec<Attr>, slabs: Vec<Vec<u8>>| {
        w.var(Var { name: name.into(), dims: vec![time, stations], typ: NcType::Float, attrs, data: VarData::Record(slabs) });
    };

    his_var(&mut w, "point_zs", his_attrs("m", "sea_surface_height_above_reference_level", "Water level zs"), his_slabs(recs, |r| &r.zs));
    his_var(&mut w, "point_hm0", his_attrs("m", "sea_surface_wave_significant_height", "Significant wave height Hm0"), his_slabs(recs, |r| &r.hm0));
    his_var(&mut w, "point_tp", his_attrs("s", "sea_surface_wave_period_at_variance_spectral_density_maximum", "Peak wave period Tp"), his_slabs(recs, |r| &r.tp));
    his_var(&mut w, "point_wavdir", his_attrs("degree", "sea_surface_wave_from_direction_at_variance_spectral_density_maximum", "Peak wave direction"), his_slabs(recs, |r| &r.wavdir));
    his_var(&mut w, "point_dirspr", his_attrs("degree", "sea_surface_wave_directional_spread", "Wave directional spread"), his_slabs(recs, |r| &r.dirspr));
    if ig {
        let slabs = recs.iter().map(|r| f32_payload(r.hm0ig.as_deref().unwrap_or(&[]))).collect();
        his_var(&mut w, "point_hm0ig", his_attrs("m", "sea_surface_infragravity_wave_significant_height", "Significant infragravity wave height Hm0"), slabs);
    }
    his_var(&mut w, "point_dw", his_attrs("W m-2", "", "Depth induced wave breaking"), his_slabs(recs, |r| &r.dw));
    his_var(&mut w, "point_df", his_attrs("W m-2", "", "Bottom friction"), his_slabs(recs, |r| &r.df));
    if wind {
        let sw_slabs = recs.iter().map(|r| f32_payload(r.sw.as_deref().unwrap_or(&[]))).collect();
        his_var(&mut w, "point_Sw", his_attrs("W m-2", "", "Wind input short wave energy"), sw_slabs);
        let st_slabs = recs.iter().map(|r| f32_payload(r.st.as_deref().unwrap_or(&[]))).collect();
        his_var(&mut w, "point_St", his_attrs("W m-2", "", "Wind input short wave action"), st_slabs);
    }

    // total_runtime / average_dt (never written)
    w.var(Var {
        name: "total_runtime".into(),
        dims: vec![runtime],
        typ: NcType::Float,
        attrs: vec![Attr::text("units", "s"), Attr::text("long_name", "total_model_runtime_in_seconds")],
        data: VarData::Fixed(f32_single(NC_FILL_FLOAT)),
    });
    w.var(Var {
        name: "average_dt".into(),
        dims: vec![runtime],
        typ: NcType::Float,
        attrs: vec![Attr::text("units", "s"), Attr::text("long_name", "model_average_timestep_in_seconds")],
        data: VarData::Fixed(f32_single(NC_FILL_FLOAT)),
    });

    w.write_to(path)
}
