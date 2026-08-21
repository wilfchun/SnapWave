//! Rust port of the mesh-preprocessing and boundary/observation-mapping
//! geometry in `src/snapwave_domain.f90`, `src/snapwave_boundaries.f90`
//! and `src/snapwave_obspoints.f90` (plan.md, Phase 9: "Surrounding-point
//! and upwind-neighbour preprocessing" and "Boundary and observation point
//! interpolation").
//!
//! # What is ported
//!
//! The derived geometry the Fortran core computes during
//! `initialize_snapwave_domain` / `read_obs_points` / `read_boundary_data`,
//! reproduced here as pure functions over Rust-owned buffers:
//!
//! * [`fm_surrounding_points`] — the sorted ring of surrounding nodes per
//!   node (`kp`) and the bed slopes (`dhdx`/`dhdy`) via a `real*4` plane
//!   fit ([`plane_fit`] → [`solve_linear_system`]).
//! * [`find_upwind_neighbours`] / [`intersect_angle`] — the two upwind
//!   neighbours and interpolation weight per node and direction.
//! * [`make_map_fm`] (in `crate::interp`) — observation-point interpolation
//!   weights and references.
//! * [`find_boundary_indices`] — the two nearest boundary support points
//!   and interpolation factor per grid boundary node.
//! * [`neuboundaries`] — the Neumann-boundary connection table.
//!
//! These are validated against the unchanged Fortran routines through the
//! temporary `snapwave_geometry_dump_c` hook (`crate::geometry_compare`,
//! driven by the wrapper's `--compare-geometry` mode). The run path still
//! lets Fortran compute them (the Fortran remains the numerical authority;
//! plan.md Phase 9 step 5 only removes the Fortran utilities once parity is
//! proven and a later phase wires the Rust results into the state handoff).
//!
//! # K-d tree / Triangle decision (plan.md Phase 9, step 4)
//!
//! The *sample-point* interpolation path — `read_interpolate_map_input`
//! feeding `triintfast`/`findtri_kdtree`/`dlaun`, which in turn drive the
//! bundled C Triangle triangulation and the Fortran `kdtree2` wrapper — is
//! deliberately **not** migrated in this phase:
//!
//! 1. It is reachable only through the `fw`/`fwig`/`u10`/`u10dir`
//!    value-or-file inputs when those strings name a *file*, and no
//!    checked-in testcase exercises that branch (every mesh is NetCDF and
//!    the friction/wind inputs are uniform values).
//! 2. Triangle and `kdtree2` are read-only third-party code (AGENTS.md
//!    rule 6); replacing them with Rust-native crates (a Delaunay crate +
//!    a k-d tree crate) would add runtime dependencies against the repo's
//!    "keep dependencies minimal" rule, with no oracle testcase to justify
//!    the licensing/parity review the plan requires.
//!
//! So Triangle + `kdtree2` + `triintfast` stay in Fortran for now, and the
//! decision is recorded here and in `plan.md`. When a sample-backed
//! testcase exists (or the value-or-file readers migrate), the parity tests
//! and licensing review the plan demands can precede any replacement.

// ----------------------------------------------------------------------
// Small helpers shared with the Fortran routines
// ----------------------------------------------------------------------

/// `cosd` (gfortran extension used by `fm_surrounding_points` /
/// `find_upwind_neighbours` for spherical grids): cosine of an angle in
/// degrees.
pub fn cosd(deg: f64) -> f64 {
    (deg * std::f64::consts::PI / 180.0).cos()
}

/// `nint(360./dtheta)` of `initialize_snapwave_domain`: `nint` rounds half
/// away from zero, which is exactly Rust's `f32::round`. (The sector-based
/// `ntheta = nint(sector/dtheta)` uses the same `nint`, so it is not given
/// its own helper.)
pub fn ntheta360(dtheta: f32) -> usize {
    (360.0f32 / dtheta).round() as usize
}

/// The `real*8` directional grid `theta360d0` of `initialize_snapwave_domain`
/// (degrees → radians). `0.0 + 0.5*dtheta + (itheta-1)*dtheta` is evaluated
/// in `real*4` and widened, then scaled by the `real*4` `deg2rad` parameter.
pub fn theta360d0(dtheta: f32, ntheta360: usize) -> Vec<f64> {
    let deg2rad = crate::text_input::deg2rad_f32() as f64;
    (0..ntheta360)
        .map(|i| (0.5f32 * dtheta + (i as f32) * dtheta) as f64 * deg2rad)
        .collect()
}

/// `findloc` of `snapwave_domain.f90`: first index (zero-based) of `b` in
/// `a[..n]`, or `None` (Fortran returns `-1`).
fn findloc(a: &[i32], b: i32) -> Option<usize> {
    (0..a.len()).find(|&i| a[i] == b)
}

/// `findloc` over a length-prefixed slice (the Fortran `findloc(a, n, b, i)`).
fn findloc_range(a: &[i32], n: usize, b: i32) -> Option<usize> {
    (0..n).find(|&i| a[i] == b)
}

// ----------------------------------------------------------------------
// Plane fit (real*4) — exact Gauss-Jordan, matching the Fortran
// ----------------------------------------------------------------------

/// `solve_linear_system` of `snapwave_domain.f90` (real*4 Gauss-Jordan).
/// `a` is the 3x3 system, `b` the right-hand side; the division of the
/// pivot row is skipped when the pivot's magnitude is `< 1e-10`, exactly
/// like the Fortran. All arithmetic is `f32`, in the same operation order.
pub fn solve_linear_system(a: [[f32; 3]; 3], b: [f32; 3]) -> [f32; 3] {
    let mut m = a;
    let mut x = b;
    for i in 1..=3usize {
        let mut factor = m[i - 1][i - 1];
        // `1.d-10` is a real*8 literal; `abs(factor)` (real*4) is widened
        // before the comparison.
        if (factor.abs() as f64) >= 1.0e-10 {
            for c in 0..3 {
                m[i - 1][c] /= factor;
            }
            x[i - 1] /= factor;
        }
        for j in 1..=3usize {
            if i != j {
                factor = m[j - 1][i - 1];
                for c in 0..3 {
                    m[j - 1][c] -= factor * m[i - 1][c];
                }
                x[j - 1] -= factor * x[i - 1];
            }
        }
    }
    x
}

/// `plane_fit` of `snapwave_domain.f90`: least-squares plane `z = a*x +
/// b*y + c` over `n` points; returns `(dzdx, dzdy)` = `(a, b)`. The
/// design matrix `XX(k, 1..3) = [x, y, 1.0]` (the `1.0d0` widens to real*4)
/// and the `XT_X`/`XT_z` accumulation are `real*4`, in the Fortran order.
pub fn plane_fit(x: &[f32], y: &[f32], z: &[f32], n: usize) -> (f32, f32) {
    let col = |k: usize| [x[k], y[k], 1.0f32];
    let mut xtx = [[0.0f32; 3]; 3];
    let mut xtz = [0.0f32; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..n {
                xtx[i][j] += col(k)[i] * col(k)[j];
            }
        }
    }
    for i in 0..3 {
        for k in 0..n {
            xtz[i] += col(k)[i] * z[k];
        }
    }
    let coeffs = solve_linear_system(xtx, xtz);
    (coeffs[0], coeffs[1])
}

// ----------------------------------------------------------------------
// Surrounding points (fm_surrounding_points)
// ----------------------------------------------------------------------

/// Result of [`fm_surrounding_points`]: `kp(np, no_nodes)` flattened
/// `[k*np + ip]`, plus the bed slopes `dhdx`/`dhdy` (`real*4`).
#[derive(Debug, Clone, PartialEq)]
pub struct SurroundingPoints {
    pub np: usize,
    /// `kp(np, no_nodes)`, `[k*np + ip]`; `0` for "no more surrounding node".
    pub kp: Vec<i32>,
    pub dhdx: Vec<f32>,
    pub dhdy: Vec<f32>,
}

/// `fm_surrounding_points` of `snapwave_domain.f90`. `face_nodes` is
/// `(4, no_faces)` node-major with the `-999` "no fourth node" sentinel;
/// `zn` is the bed level (`real*4`). The edge-connection logic and the
/// `real*4` plane fit reproduce the Fortran exactly (including the
/// `xp(j) = xn(surr)-xn(kn)` double→single truncation).
pub fn fm_surrounding_points(
    xn: &[f64],
    yn: &[f64],
    zn: &[f32],
    sferic: i32,
    face_nodes: &[i32],
    no_faces: usize,
    np: usize,
) -> SurroundingPoints {
    let no_nodes = xn.len();
    // `40075017.` / `40007863.` are single-precision literals widened to
    // real*8 (only used on the sferic==1 path).
    let circumf_eq = 40075017.0f32 as f64;
    let circumf_pole = 40007863.0f32 as f64;

    // Build node -> connected-cell table (Fortran `connected_cells(12,:)`).
    let mut connected: Vec<Vec<i32>> = vec![Vec::new(); no_nodes];
    for k in 0..no_faces {
        for inode in 0..4 {
            let knode = face_nodes[k * 4 + inode];
            if knode != -999 {
                connected[(knode - 1) as usize].push((k + 1) as i32);
            }
        }
    }

    let mut kp = vec![0i32; np * no_nodes];
    let mut dhdx = vec![0.0f32; no_nodes];
    let mut dhdy = vec![0.0f32; no_nodes];

    for kn in 0..no_nodes {
        let kn1 = (kn + 1) as i32;
        let mut surr_points = vec![0i32; np];
        let mut surr_pts = vec![0i32; np];
        let mut isp = 0usize; // count of edge endpoints collected

        for &cell1 in &connected[kn] {
            let k = (cell1 - 1) as usize;
            let kpts = [
                face_nodes[k * 4],
                face_nodes[k * 4 + 1],
                face_nodes[k * 4 + 2],
                face_nodes[k * 4 + 3],
            ];
            if kpts[3] == -999 {
                // Triangle: the two other nodes (in face order).
                let mut jj = 0usize;
                for j in 0..3 {
                    if kpts[j] != kn1 {
                        surr_points[isp + jj] = kpts[j];
                        jj += 1;
                    }
                }
                isp += 2;
            } else {
                // Quad: the two edges adjacent to `kn`, in wrapped order.
                let ip = (0..4).find(|&j| kpts[j] == kn1).expect("kn in its own face");
                let edge = match ip {
                    0 => [kpts[1], kpts[2], kpts[2], kpts[3]],
                    1 => [kpts[2], kpts[3], kpts[3], kpts[0]],
                    2 => [kpts[3], kpts[0], kpts[0], kpts[1]],
                    _ => [kpts[0], kpts[1], kpts[1], kpts[2]],
                };
                surr_points[isp..isp + 4].copy_from_slice(&edge);
                isp += 4;
            }
        }

        let mut count_pts = 0usize; // isp2 (number of ring entries)
        if isp >= 2 {
            surr_pts[0] = surr_points[0];
            surr_pts[1] = surr_points[1];
            // surr_points = [surr_points(3:), 0, 0]
            shift_left2(&mut surr_points);
            count_pts = 2;
            let mut remaining = isp - 2;
            while remaining >= 2 {
                let target = surr_pts[count_pts - 1]; // surr_pts(isp2)
                match findloc_range(&surr_points, remaining, target) {
                    Some(n) => {
                        if (n + 1) % 2 == 1 {
                            surr_pts[count_pts] = surr_points[n + 1];
                            count_pts += 1;
                            remove_two(&mut surr_points, n);
                        } else {
                            surr_pts[count_pts] = surr_points[n - 1];
                            count_pts += 1;
                            remove_two(&mut surr_points, n - 1);
                        }
                    }
                    None => {
                        let first = surr_pts[0];
                        match findloc_range(&surr_points, remaining, first) {
                            Some(n) => {
                                if (n + 1) % 2 == 1 {
                                    prepend(&mut surr_pts, count_pts, surr_points[n + 1]);
                                    count_pts += 1;
                                    remove_two(&mut surr_points, n);
                                } else {
                                    prepend(&mut surr_pts, count_pts, surr_points[n - 1]);
                                    count_pts += 1;
                                    remove_two(&mut surr_points, n - 1);
                                }
                            }
                            None => {
                                surr_pts.iter_mut().for_each(|v| *v = 0);
                            }
                        }
                    }
                }
                remaining -= 2;
            }
        }
        // else: surr_pts stays all-zero (isp < 2).

        for ip in 0..np {
            kp[kn * np + ip] = surr_pts[ip];
        }

        if surr_pts.iter().copied().sum::<i32>() > 0 {
            let nfit = count_pts - 1;
            let mut xp = vec![0.0f32; np];
            let mut yp = vec![0.0f32; np];
            let mut zp = vec![0.0f32; np];
            for j in 0..nfit {
                let sj = (surr_pts[j] - 1) as usize;
                xp[j] = (xn[sj] - xn[kn]) as f32;
                yp[j] = (yn[sj] - yn[kn]) as f32;
                zp[j] = zn[sj];
                if sferic == 1 {
                    xp[j] = (xp[j] as f64 * circumf_eq / 360.0 * cosd(yn[kn])) as f32;
                    yp[j] = (yp[j] as f64 * circumf_pole / 360.0) as f32;
                }
            }
            let (dzdx, dzdy) = plane_fit(&xp[..nfit], &yp[..nfit], &zp[..nfit], nfit);
            dhdx[kn] = -dzdx;
            dhdy[kn] = -dzdy;
        } else {
            dhdx[kn] = 0.0;
            dhdy[kn] = 0.0;
        }
    }

    SurroundingPoints { np, kp, dhdx, dhdy }
}

/// `surr_points = [surr_points(3:), 0, 0]`: shift the whole fixed-length
/// array left by two and zero the last two slots.
fn shift_left2(arr: &mut [i32]) {
    let n = arr.len();
    for idx in 0..n {
        arr[idx] = if idx + 2 < n { arr[idx + 2] } else { 0 };
    }
}

/// `surr_points = [surr_points(1:at-1), surr_points(at+2:), 0, 0]`: remove
/// the two elements at (zero-based) `at` and `at+1`, shift left, zero the
/// last two slots.
fn remove_two(arr: &mut [i32], at: usize) {
    let n = arr.len();
    for idx in at..n {
        arr[idx] = if idx + 2 < n { arr[idx + 2] } else { 0 };
    }
}

/// `surr_pts(1:count+1) = [new, surr_pts(1:count)]`: prepend one element.
fn prepend(arr: &mut [i32], count: usize, val: i32) {
    for idx in (0..count).rev() {
        arr[idx + 1] = arr[idx];
    }
    arr[0] = val;
}

// ----------------------------------------------------------------------
// Upwind neighbours (find_upwind_neighbours / intersect_angle)
// ----------------------------------------------------------------------

/// `intersect_angle` of `snapwave_domain.f90`: intersection of the ray from
/// `(x0, y0)` at angle `phi` with the segment `(x, y)`. Returns the two
/// interpolation weights `w`, the upwind distance `ds` (0.0 when the ray
/// misses), and the intersection point `(xi, yi)`.
pub fn intersect_angle(x0: f64, y0: f64, phi: f64, x: [f64; 2], y: [f64; 2]) -> ([f64; 2], f64, f64, f64) {
    // `eps = 1.0e-2` is a single-precision literal widened to real*8.
    let eps = 1.0e-2f32 as f64;
    let (xi, yi);
    if (x[1] - x[0]).abs() > eps {
        let m = (y[1] - y[0]) / (x[1] - x[0]);
        let a = y[0] - m * x[0];
        let n = phi.tan();
        let b = y0 - n * x0;
        xi = (b - a) / (m - n);
        yi = a + m * xi;
    } else {
        yi = (x[0] - x0) * phi.tan() + y0;
        xi = x[0];
    }
    let l = (x[1] - x[0]).hypot(y[1] - y[0]);
    let d1 = (xi - x[0]).hypot(yi - y[0]);
    let d2 = (xi - x[1]).hypot(yi - y[1]);
    let ds = (xi - x0).hypot(yi - y0);
    let err = ((x0 - xi) + ds * phi.cos()).hypot((y0 - yi) + ds * phi.sin());
    if (l - d1 - d2).abs() < eps && err < eps {
        ([d2 / l, d1 / l], ds, xi, yi)
    } else {
        ([0.0, 0.0], 0.0, xi, yi)
    }
}

/// Result of [`find_upwind_neighbours`]. All arrays are Fortran
/// column-major over `(2, ntheta, no_nodes)` / `(ntheta, no_nodes)`,
/// flattened `[k*(2*ntheta) + itheta*2 + i]` and `[k*ntheta + itheta]`.
#[derive(Debug, Clone, PartialEq)]
pub struct UpwindNeighbours {
    pub ntheta: usize,
    pub w: Vec<f64>,
    pub prev: Vec<i32>,
    pub ds: Vec<f64>,
}

/// `find_upwind_neighbours` of `snapwave_domain.f90`: the two upwind
/// neighbours, their weights and the upwind distance for each node and
/// direction. `kp`/`np` come from [`fm_surrounding_points`]; `theta` is the
/// `real*8` directional grid (radians). Nodes with no surrounding points
/// keep the zero initialisation; directions with no intersection get the
/// `prev = 1, w = 0, ds = 0` "closed boundary" sentinel.
pub fn find_upwind_neighbours(
    x: &[f64],
    y: &[f64],
    sferic: i32,
    theta: &[f64],
    kp: &[i32],
    np: usize,
) -> UpwindNeighbours {
    let no_nodes = x.len();
    let ntheta = theta.len();
    // `pi = 4.0*atan(1.0)` in the Fortran is real*4 arithmetic widened to
    // real*8; `40075017.`/`40007863.`/`0.001` are single-precision literals
    // widened to real*8. Reproduce the exact values, not their f64 twins.
    let pi = (4.0f32 * 1.0f32.atan()) as f64;
    let circumf_eq = 40075017.0f32 as f64;
    let circumf_pole = 40007863.0f32 as f64;
    let epsdist = 0.001f32 as f64;

    let mut w = vec![0.0f64; 2 * ntheta * no_nodes];
    let mut prev = vec![0i32; 2 * ntheta * no_nodes];
    let mut ds = vec![0.0f64; ntheta * no_nodes];

    for k in 0..no_nodes {
        let col = &kp[k * np..(k + 1) * np];
        // nploc = findloc(kp(:,k), 0) - 1; the first zero bounds the list.
        let nploc = findloc(col, 0).unwrap_or(np);
        if col[0] == 0 {
            continue;
        }
        for itheta in 0..ntheta {
            let mut dss = 0.0f64;
            for ip in 0..nploc.saturating_sub(1) {
                let ind1 = col[ip] as usize;
                let ind2 = col[ip + 1] as usize;
                let (ww, d, _xi, _yi) = if sferic == 0 {
                    intersect_angle(
                        x[k],
                        y[k],
                        theta[itheta] + pi,
                        [x[ind1 - 1], x[ind2 - 1]],
                        [y[ind1 - 1], y[ind2 - 1]],
                    )
                } else {
                    let xsect = [
                        (x[ind1 - 1] - x[k]) * circumf_eq / 360.0 * cosd(y[k]),
                        (x[ind2 - 1] - x[k]) * circumf_eq / 360.0 * cosd(y[k]),
                    ];
                    let ysect = [
                        (y[ind1 - 1] - y[k]) * circumf_pole / 360.0,
                        (y[ind2 - 1] - y[k]) * circumf_pole / 360.0,
                    ];
                    intersect_angle(0.0, 0.0, theta[itheta] + pi, xsect, ysect)
                };
                dss = d;
                if dss > epsdist {
                    w[k * (2 * ntheta) + itheta * 2] = ww[0];
                    w[k * (2 * ntheta) + itheta * 2 + 1] = ww[1];
                    ds[k * ntheta + itheta] = dss;
                    prev[k * (2 * ntheta) + itheta * 2] = col[ip];
                    prev[k * (2 * ntheta) + itheta * 2 + 1] = col[ip + 1];
                    break;
                }
            }
            if dss <= epsdist {
                prev[k * (2 * ntheta) + itheta * 2] = 1;
                prev[k * (2 * ntheta) + itheta * 2 + 1] = 1;
                w[k * (2 * ntheta) + itheta * 2] = 0.0;
                w[k * (2 * ntheta) + itheta * 2 + 1] = 0.0;
                // ds stays 0 (its initialisation).
            }
        }
    }

    UpwindNeighbours { ntheta, w, prev, ds }
}

// ----------------------------------------------------------------------
// Neumann boundary connection (neuboundaries)
// ----------------------------------------------------------------------

/// `neuboundaries` of `snapwave_domain.f90`: for every node on a Neumann
/// polyline segment, the index of the node on the opposite side of the
/// line (0 when none). `neumannconnected(kmin) = k` with 1-based indices;
/// the returned slice is `neumannconnected(k)` for `k = 1..no_nodes`.
pub fn neuboundaries(x: &[f64], y: &[f64], xneu: &[f64], yneu: &[f64], tol: f32) -> Vec<i32> {
    let no_nodes = x.len();
    let n_neu = xneu.len();
    let tol_f64 = tol as f64;
    let mut conn = vec![0i32; no_nodes];

    for ib in 0..n_neu.saturating_sub(1) {
        if xneu[ib] != -999.0 && xneu[ib + 1] != -999.0 {
            let alpha = (yneu[ib + 1] - yneu[ib]).atan2(xneu[ib + 1] - xneu[ib]);
            let cosa = alpha.cos();
            let sina = alpha.sin();
            let xend = (xneu[ib + 1] - xneu[ib]) * cosa + (yneu[ib + 1] - yneu[ib]) * sina;
            for k in 0..no_nodes {
                let x1 = (x[k] - xneu[ib]) * cosa + (y[k] - yneu[ib]) * sina;
                let y1 = -(x[k] - xneu[ib]) * sina + (y[k] - yneu[ib]) * cosa;
                if x1 >= 0.0 && x1 <= xend && y1.abs() < tol_f64 {
                    let mut distmin = 1.0e10f64;
                    let mut kmin = 0i32;
                    for k2 in 0..no_nodes {
                        let x2 = (x[k2] - xneu[ib]) * cosa + (y[k2] - yneu[ib]) * sina;
                        let y2 = -(x[k2] - xneu[ib]) * sina + (y[k2] - yneu[ib]) * cosa;
                        if (x2 - x1).abs() < tol_f64 && k2 != k && (y2 - y1).abs() < distmin {
                            kmin = (k2 + 1) as i32;
                            distmin = (y2 - y1).abs();
                        }
                    }
                    if kmin > 0 {
                        conn[(kmin - 1) as usize] = (k + 1) as i32;
                    }
                }
            }
        }
    }
    conn
}

// ----------------------------------------------------------------------
// Boundary support-point interpolation (find_boundary_indices)
// ----------------------------------------------------------------------

/// Result of [`find_boundary_indices`].
#[derive(Debug, Clone, PartialEq)]
pub struct BoundaryIndices {
    /// Number of grid boundary nodes (`nb`).
    pub nb: usize,
    /// `nmindbnd(nb)`: the grid node index of each boundary node (1-based).
    pub nmindbnd: Vec<i32>,
    /// `ind1_bwv_cst(nb)`: nearest support point (1-based) or `1`.
    pub ind1: Vec<i32>,
    /// `ind2_bwv_cst(nb)`: second-nearest support point or `1`.
    pub ind2: Vec<i32>,
    /// `fac_bwv_cst(nb)`: weight of the nearest point (`real*4`).
    pub fac: Vec<f32>,
}

/// `find_boundary_indices` of `snapwave_boundaries.f90`: for each grid
/// boundary node (`msk == 2`) the two nearest wave-boundary support points
/// and the interpolation factor. With a single support point (`nwbnd == 1`)
/// the trivial `ind1 = ind2 = 1, fac = 1.0` branch is taken, exactly as in
/// the Fortran. `xgb`/`ygb` are `real*4` truncations of `x`/`y`, and the
/// distance is a `real*8` `sqrt` truncated back to `real*4`, matching the
/// mixed-precision comparisons.
pub fn find_boundary_indices(
    x: &[f64],
    y: &[f64],
    msk: &[i32],
    x_bwv: &[f64],
    y_bwv: &[f64],
    nwbnd: usize,
) -> BoundaryIndices {
    let mut temp: Vec<i32> = Vec::new();
    for (k, &m) in msk.iter().enumerate() {
        if m == 2 {
            temp.push((k + 1) as i32);
        }
    }
    let nb = temp.len();
    let mut nmindbnd = vec![0i32; nb];
    let mut ind1 = vec![0i32; nb];
    let mut ind2 = vec![0i32; nb];
    let mut fac = vec![0.0f32; nb];

    for i in 0..nb {
        let k = temp[i] as usize - 1;
        nmindbnd[i] = temp[i];
        if nwbnd > 1 {
            let xgb = x[k] as f32;
            let ygb = y[k] as f32;
            let mut dst1 = 1.0e10f32;
            let mut dst2 = 1.0e10f32;
            let mut ib1 = 0i32;
            let mut ib2 = 0i32;
            for ic in 0..nwbnd {
                let dx = x_bwv[ic] - xgb as f64;
                let dy = y_bwv[ic] - ygb as f64;
                let dst = (dx * dx + dy * dy).sqrt() as f32;
                if dst < dst1 {
                    dst2 = dst1;
                    ib2 = ib1;
                    dst1 = dst;
                    ib1 = (ic + 1) as i32;
                } else if dst < dst2 {
                    dst2 = dst;
                    ib2 = (ic + 1) as i32;
                }
            }
            ind1[i] = ib1;
            ind2[i] = ib2;
            fac[i] = dst2 / (dst1 + dst2);
        } else if nwbnd == 1 {
            ind1[i] = 1;
            ind2[i] = 1;
            fac[i] = 1.0;
        }
        // nwbnd == 0: ind1/ind2/fac stay 0 — the Fortran leaves those
        // globals unallocated in that case (find_boundary_indices is guarded
        // by `if (nwbnd > 0)`), so the dump omits them and the comparison
        // skips them.
    }

    BoundaryIndices { nb, nmindbnd, ind1, ind2, fac }
}

// ----------------------------------------------------------------------
// Full domain-geometry orchestration (mirrors initialize_snapwave_domain
// + read_obs_points + read_boundary_data, geometry only)
// ----------------------------------------------------------------------

/// Everything the geometry comparison needs from one mesh + its derived
/// quantities, matching the Fortran globals the dump hook emits.
#[derive(Debug, Clone, PartialEq)]
pub struct DomainGeometry {
    pub no_nodes: usize,
    pub np: usize,
    pub ntheta360: usize,
    pub kp: Vec<i32>,
    pub dhdx: Vec<f32>,
    pub dhdy: Vec<f32>,
    /// `w360 = w360d0*1.0` (`real*4`), column-major `(2, ntheta360, no_nodes)`.
    pub w360: Vec<f32>,
    pub prev360: Vec<i32>,
    /// `ds360 = ds360d0*1.0` (`real*4`), column-major `(ntheta360, no_nodes)`.
    pub ds360: Vec<f32>,
    /// `msk` after enclosure + Neumann refinement (0/1/2/3).
    pub msk: Vec<i32>,
    pub neumannconnected: Vec<i32>,
    pub nb: usize,
    pub nnmb: usize,
    pub nmindbnd: Vec<i32>,
    pub neubnd: Vec<i32>,
}

/// Compute the domain-derived geometry from a mesh and the polylines, exactly
/// mirroring `initialize_snapwave_domain` (surrounding points, upwind
/// neighbours, enclosure→mask and Neumann→mask refinement, boundary-node
/// lists). `msk` is the as-read mask (all `1` for a NetCDF mesh).
///
/// `enclosure`/`neumann` are the polylines (`None` = absent). The enclosure
/// refinement uses the `real*8` upwind distances (`ds360d0 == 0` means "no
/// intersection in this direction"), so [`UpwindNeighbours`]' `ds` (f64) is
/// used here *before* the `real*4` truncation.
pub fn compute_domain_geometry(
    x: &[f64],
    y: &[f64],
    zb: &[f32],
    sferic: i32,
    face_nodes: &[i32],
    no_faces: usize,
    dtheta: f32,
    tol: f32,
    enclosure: Option<(&[f64], &[f64])>,
    neumann: Option<(&[f64], &[f64])>,
    msk: &[i32],
) -> DomainGeometry {
    let no_nodes = x.len();
    let np = 22;
    let ntheta360 = ntheta360(dtheta);

    let surround = fm_surrounding_points(x, y, zb, sferic, face_nodes, no_faces, np);
    let theta = theta360d0(dtheta, ntheta360);
    let upwind = find_upwind_neighbours(x, y, sferic, &theta, &surround.kp, np);

    // w360/ds360 single-precision conversions.
    let w360: Vec<f32> = upwind.w.iter().map(|&v| v as f32).collect();
    let ds360: Vec<f32> = upwind.ds.iter().map(|&v| v as f32).collect();

    // Enclosure refinement (uses ds360d0 == 0, the f64 upwind distances).
    let mut msk = msk.to_vec();
    if let Some((xb, yb)) = enclosure {
        if !xb.is_empty() {
            for k in 0..no_nodes {
                for itheta in 0..ntheta360 {
                    if upwind.ds[k * ntheta360 + itheta] == 0.0 {
                        if crate::interp::ipon(xb, yb, x[k], y[k]) > 0 {
                            msk[k] = 2;
                        }
                    }
                }
            }
        }
    }

    // Neumann refinement.
    let mut neumannconnected = match neumann {
        Some((xn, yn)) => neuboundaries(x, y, xn, yn, tol),
        None => vec![0i32; no_nodes],
    };
    for k in 0..no_nodes {
        if neumannconnected[k] > 0 {
            if msk[k] == 1 {
                msk[(neumannconnected[k] - 1) as usize] = 3;
            } else {
                neumannconnected[k] = 0;
            }
        }
    }

    let mut nmindbnd = Vec::new();
    let mut neubnd = Vec::new();
    for (k, &m) in msk.iter().enumerate() {
        if m == 2 {
            nmindbnd.push((k + 1) as i32);
        } else if m == 3 {
            neubnd.push((k + 1) as i32);
        }
    }
    let nb = nmindbnd.len();
    let nnmb = neubnd.len();

    DomainGeometry {
        no_nodes,
        np,
        ntheta360,
        kp: surround.kp,
        dhdx: surround.dhdx,
        dhdy: surround.dhdy,
        w360,
        prev360: upwind.prev,
        ds360,
        msk,
        neumannconnected,
        nb,
        nnmb,
        nmindbnd,
        neubnd,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ntheta_ranges_match_fortran_nint() {
        assert_eq!(ntheta360(10.0), 36);
        assert_eq!(ntheta360(5.0), 72);
    }

    #[test]
    fn solve_linear_system_solves_a_known_system() {
        // [[2,0,0],[0,2,0],[0,0,2]] x = [4,6,8] -> x = [2,3,4].
        let a = [[2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [0.0, 0.0, 2.0]];
        let x = solve_linear_system(a, [4.0, 6.0, 8.0]);
        assert_eq!(x, [2.0, 3.0, 4.0]);
    }

    #[test]
    fn plane_fit_recovers_a_flat_plane() {
        // z = 2x - 3y + 5 over four points.
        let pts = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0)];
        let x: Vec<f32> = pts.iter().map(|p| p.0).collect();
        let y: Vec<f32> = pts.iter().map(|p| p.1).collect();
        let z: Vec<f32> = pts.iter().map(|&(px, py)| 2.0 * px - 3.0 * py + 5.0).collect();
        let (dzdx, dzdy) = plane_fit(&x, &y, &z, 4);
        assert!((dzdx - 2.0).abs() < 1e-4, "dzdx = {dzdx}");
        assert!((dzdy - -3.0).abs() < 1e-4, "dzdy = {dzdy}");
    }

    #[test]
    fn intersect_angle_hits_and_misses() {
        // Segment from (1,0) to (1,2); ray from (0,1) eastward (phi=0).
        let (w, ds, xi, yi) = intersect_angle(0.0, 1.0, 0.0, [1.0, 1.0], [0.0, 2.0]);
        assert!((ds - 1.0).abs() < 1e-9, "ds = {ds}");
        assert!((xi - 1.0).abs() < 1e-9 && (yi - 1.0).abs() < 1e-9);
        // Midpoint of the segment: equal weights.
        assert!((w[0] - 0.5).abs() < 1e-9 && (w[1] - 0.5).abs() < 1e-9, "w = {w:?}");
        // Ray pointing away (phi = pi) misses.
        let (_, ds2, _, _) = intersect_angle(0.0, 1.0, std::f64::consts::PI, [1.0, 1.0], [0.0, 2.0]);
        assert_eq!(ds2, 0.0);
    }

    #[test]
    fn find_boundary_indices_two_point_case() {
        // Two support points at x=0 and x=10; one boundary node at x=2.
        let x = [2.0f64];
        let y = [0.0f64];
        let msk = [2i32];
        let x_bwv = [0.0f64, 10.0];
        let y_bwv = [0.0f64, 0.0];
        let r = find_boundary_indices(&x, &y, &msk, &x_bwv, &y_bwv, 2);
        assert_eq!(r.nb, 1);
        assert_eq!(r.nmindbnd, vec![1]);
        assert_eq!(r.ind1, vec![1]); // nearest = point 1 (x=0)
        assert_eq!(r.ind2, vec![2]);
        assert!((r.fac[0] - (8.0 / (2.0 + 8.0))).abs() < 1e-6, "fac = {}", r.fac[0]);
    }

    #[test]
    fn find_boundary_indices_single_point_is_trivial() {
        let r = find_boundary_indices(&[1.0, 2.0], &[0.0, 0.0], &[2, 2], &[5.0], &[5.0], 1);
        assert_eq!(r.nb, 2);
        assert_eq!(r.ind1, vec![1, 1]);
        assert_eq!(r.ind2, vec![1, 1]);
        assert_eq!(r.fac, vec![1.0, 1.0]);
    }

    #[test]
    fn neuboundaries_connects_opposite_nodes() {
        // A horizontal Neumann segment at y=0 from x=0 to x=10; nodes at
        // (5, -1) and (5, +1) are a mirrored pair.
        let x = [5.0, 5.0];
        let y = [-1.0, 1.0];
        let xneu = [0.0, 10.0];
        let yneu = [0.0, 0.0];
        let conn = neuboundaries(&x, &y, &xneu, &yneu, 10.0);
        // Node 1 (kmin=1) is on the segment; its opposite is node 2.
        assert_eq!(conn[1], 1, "conn = {conn:?}"); // neumannconnected(node2) = node1
    }
}
