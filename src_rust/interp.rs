//! Rust port of the interpolation helpers in `src/interp.F90`
//! (plan.md, Phase 9: "Generic interpolation helpers from interp.F90").
//!
//! # What is ported
//!
//! Every routine of `interp.F90` that is either (a) on the SnapWave runtime
//! path — [`make_map_fm`] and its helpers [`sort`]/[`indexx`], [`hunt`],
//! [`ipon`], [`bilin5`], [`triangle_intp`] — or (b) a small, deterministic
//! generic helper the plan lists in this phase's scope: [`binary_search`],
//! [`linear_interp`], [`linear_interp_2d`], the trapezoidal family
//! ([`trapezoidal`], [`trapezoidal_cyclic`],
//! [`interp_using_trapez_rule`], [`interp_in_cyclic_function`]), and the
//! curvilinear-grid mapping routines ([`make_map`], [`mkmap_step`],
//! [`grmap`], [`grmap2`], [`grmap_sg`]).
//!
//! Not ported (documented decision, plan.md Phase 9 step 4): `triintfast`,
//! `findtri_kdtree`, `dlaun`, `linear`, `getdx`/`getdy`/`getdxdy`, `cross`
//! — the sample-point interpolation path backed by the bundled C Triangle
//! and the Fortran `kdtree2` wrapper. See the module docs of
//! `crate::geometry` for the decision.
//!
//! # Floating-point fidelity
//!
//! These are ports, not clean reimplementations: every `real*8`/`real*4`
//! width, every operation order and every early-return (`ier`) is preserved
//! so the outputs match the Fortran oracle bit-for-bit on the same toolchain
//! (both runtimes call the same libm for `sqrt`/`hypot`/`tan`/`atan`). Where
//! Fortran stores an out-parameter only on the success path, the Rust
//! equivalent returns `Option`/`0` so the caller leaves the slot untouched,
//! matching the Fortran "leave it as initialised" behaviour.
//!
//! # Indexing conventions
//!
//! Public slices are zero-based. Where a routine's Fortran analogue uses a
//! sentinel (`j = 0` meaning "below the array", `j = n` "at/above the top")
//! the Rust return value keeps that convention and documents it.

// Only [`make_map_fm`] and its helpers ([`sort`]/[`indexx`], [`hunt`],
// [`ipon`], [`bilin5`], [`triangle_intp`]) are on the SnapWave runtime path
// today (via `crate::geometry_compare`). The remaining generic helpers are
// ported for the Phase 9 scope ("generic interpolation helpers from
// interp.F90") and exercised by unit tests, but have no production caller
// until the curvilinear/waveparams code migrates — the same position
// `crate::ffi_layout` documents for its own one-based helpers.
#![allow(dead_code)]

// ----------------------------------------------------------------------
// Search and 1-D interpolation
// ----------------------------------------------------------------------

/// `binary_search(xx, n, x, j)` of `interp.F90` (real*8). Returns the
/// Fortran `j` **verbatim** (1-based, in `0..=n`): the largest index with
/// `xx(j) < x` — i.e. the lower bracket of the interval `xx(j) ..= xx(j+1)`
/// containing `x` — with the out-of-range sentinels `0` (when `x <= xx(1)`)
/// and `n` (when `x > xx(n)`). Works for ascending or descending input (the
/// `xx(n) > xx(1)` test). Integer division truncates toward zero exactly
/// like Fortran. Callers mirror the Fortran `linear_interp` guard logic
/// (`j <= 0` / `j >= n`); see [`linear_interp`].
pub fn binary_search(xx: &[f64], x: f64) -> usize {
    let n = xx.len();
    if n == 0 {
        return 0;
    }
    let mut jl = 0usize;
    let mut ju = n + 1;
    while ju - jl > 1 {
        let jm = (ju + jl) / 2;
        let l1 = xx[n - 1] > xx[0];
        let l2 = x > xx[jm - 1];
        if (l1 && l2) || !(l1 || l2) {
            jl = jm;
        } else {
            ju = jm;
        }
    }
    jl
}

/// `linear_interp(x, y, n, xx, yy, indint)` of `interp.F90` (real*8).
/// Returns `(yy, j)` where `j` is the Fortran `indint` **verbatim** (1-based,
/// in `0..=n`): the lower-bracket index `x(j) ..= x(j+1)` (0-indexed slots
/// `j-1` and `j`), with `0`/`n` when `xx` is outside the array. `n == 0`
/// yields `(0.0, 0)` and `n == 1` yields `(y[0], 0)` (both `indint = 0`),
/// matching the Fortran early returns.
pub fn linear_interp(x: &[f64], y: &[f64], xx: f64) -> (f64, usize) {
    let n = x.len();
    if n == 0 {
        return (0.0, 0);
    }
    if n == 1 {
        return (y[0], 0);
    }
    let j = binary_search(x, xx); // Fortran `indint`: 1-based, 0..=n
    let yy = if j == 0 {
        y[0]
    } else if j >= n {
        y[n - 1]
    } else {
        // Fortran interpolates between x(j) (lower, slot j-1) and x(j+1)
        // (upper, slot j): a = x(j+1), b = x(j).
        let a = x[j];
        let b = x[j - 1];
        let dyy = if a == b { 0.0 } else { (y[j] - y[j - 1]) / (a - b) };
        y[j - 1] + (xx - x[j - 1]) * dyy
    };
    (yy, j)
}

/// First index `i` (zero-based) of the minimum of `arr` subject to `mask`,
/// mirroring `minloc(arr, 1, mask)`; returns `None` when the mask selects
/// nothing (Fortran's result would be 0, treated as an error upstream).
fn minloc_masked(arr: &[f64], mask: impl Fn(f64) -> bool) -> Option<usize> {
    let mut best: Option<(f64, usize)> = None;
    for (i, &v) in arr.iter().enumerate() {
        if mask(v) {
            match best {
                Some((b, _)) if b <= v => {}
                _ => best = Some((v, i)),
            }
        }
    }
    best.map(|(_, i)| i)
}

/// First index `i` (zero-based) of the maximum of `arr` subject to `mask`,
/// mirroring `maxloc(arr, 1, mask)`.
fn maxloc_masked(arr: &[f64], mask: impl Fn(f64) -> bool) -> Option<usize> {
    let mut best: Option<(f64, usize)> = None;
    for (i, &v) in arr.iter().enumerate() {
        if mask(v) {
            match best {
                Some((b, _)) if b >= v => {}
                _ => best = Some((v, i)),
            }
        }
    }
    best.map(|(_, i)| i)
}

/// `linear_interp_2d` of `interp.F90`: bilinear interpolation into a
/// curvilinear field `z(nx, ny)` (column-major: element `(ix, iy)` at
/// `ix + nx*iy`). `method` is `"interp"` (out-of-range → `exception`) or
/// `"extendclosest"` (out-of-range → nearest point); any other value falls
/// through to `exception`, matching the Fortran `select case`.
pub fn linear_interp_2d(
    x: &[f64],
    y: &[f64],
    z: &[f64],
    xx: f64,
    yy: f64,
    method: &str,
    exception: f64,
) -> f64 {
    let nx = x.len();
    let at = |ix: usize, iy: usize| z[ix + nx * iy];

    let xmin = x.iter().copied().fold(f64::INFINITY, f64::min);
    let xmax = x.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let ymin = y.iter().copied().fold(f64::INFINITY, f64::min);
    let ymax = y.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let interp_x = xx >= xmin && xx <= xmax;
    let interp_y = yy >= ymin && yy <= ymax;

    if interp_x && interp_y {
        let i1 = minloc_masked(x, |v| v >= xx).expect("xx within range");
        let i2 = maxloc_masked(x, |v| v <= xx).expect("xx within range");
        let i3 = minloc_masked(y, |v| v >= yy).expect("yy within range");
        let i4 = maxloc_masked(y, |v| v <= yy).expect("yy within range");
        let disx = x[i1] - x[i2];
        let disy = y[i3] - y[i4];
        let modx = if disx > 0.0 { (xx - x[i2]) / disx } else { 0.0 };
        let mody = if disy > 0.0 { (yy - y[i4]) / disy } else { 0.0 };
        if disy > 0.0 {
            let yint1 = (1.0 - modx) * at(i2, i3) + modx * at(i1, i3);
            let yint2 = (1.0 - modx) * at(i2, i4) + modx * at(i1, i4);
            (1.0 - mody) * yint2 + mody * yint1
        } else {
            (1.0 - modx) * at(i2, i3) + modx * at(i1, i3)
        }
    } else {
        match method {
            "extendclosest" => {
                // `minloc(abs(X-xx))`: nearest x, nearest y.
                let (i1, _) = x.iter().enumerate().fold((0usize, f64::INFINITY), |a, (i, &v)| {
                    let d = (v - xx).abs();
                    if d < a.1 { (i, d) } else { a }
                });
                let (i3, _) = y.iter().enumerate().fold((0usize, f64::INFINITY), |a, (i, &v)| {
                    let d = (v - yy).abs();
                    if d < a.1 { (i, d) } else { a }
                });
                at(i1, i3)
            }
            _ => exception,
        }
    }
}

/// `hunt` of `interp.F90` (Numerical Recipes): given `jlo` (the previous
/// bracket index, 1-based, `0..=n`) locate the bracket of `x` in the sorted
/// array `xx` and store it back into `jlo`. The cached-`jlo` start is
/// preserved for the same reason Fortran keeps it (`make_map_fm` reuses the
/// table position across cells); the *result* is independent of the cache.
pub fn hunt(xx: &[f64], x: f64, jlo: &mut usize) {
    let n = xx.len();
    let ascnd = xx[n - 1] >= xx[0];
    let mut jhi: usize;
    if *jlo == 0 || *jlo > n {
        *jlo = 0;
        jhi = n + 1;
    } else {
        let mut inc = 1usize;
        if (x >= xx[*jlo - 1]) == ascnd {
            loop {
                jhi = *jlo + inc;
                if jhi > n {
                    jhi = n + 1;
                    break;
                } else if (x >= xx[jhi - 1]) == ascnd {
                    *jlo = jhi;
                    inc += inc;
                } else {
                    break;
                }
            }
        } else {
            jhi = *jlo;
            loop {
                *jlo = jhi.saturating_sub(inc);
                if *jlo < 1 {
                    *jlo = 0;
                    break;
                } else if (x < xx[*jlo - 1]) == ascnd {
                    jhi = *jlo;
                    inc += inc;
                } else {
                    break;
                }
            }
        }
    }
    loop {
        if jhi - *jlo == 1 {
            return;
        }
        let jm = (jhi + *jlo) / 2;
        if (x > xx[jm - 1]) == ascnd {
            *jlo = jm;
        } else {
            jhi = jm;
        }
    }
}

// ----------------------------------------------------------------------
// Sorting (Numerical Recipes indexx / sort)
// ----------------------------------------------------------------------

/// `indexx` of `interp.F90`: the permutation that sorts `arr` ascending.
/// Returns zero-based indices (`sorted` in `sort` reads `arr[perm[j]]`);
/// equal values keep their input order only as far as the heap-sort does
/// (which is *not* stable, exactly like the Fortran).
pub fn indexx(arr: &[f64]) -> Vec<usize> {
    let n = arr.len();
    let mut indx: Vec<usize> = (0..n).collect();
    if n <= 1 {
        return indx;
    }
    let mut l = n / 2 + 1; // Fortran `l = n/2 + 1`
    let mut ir = n; // Fortran `ir = n`
    loop {
        let indxt: usize;
        let q: f64;
        if l > 1 {
            l -= 1;
            indxt = indx[l - 1];
            q = arr[indxt];
        } else {
            indxt = indx[ir - 1];
            q = arr[indxt];
            indx[ir - 1] = indx[0];
            ir -= 1;
            if ir == 1 {
                indx[0] = indxt;
                return indx;
            }
        }
        let mut i = l;
        let mut j = l + l;
        loop {
            if j <= ir {
                if j < ir && arr[indx[j - 1]] < arr[indx[j]] {
                    j += 1;
                }
                if q < arr[indx[j - 1]] {
                    indx[i - 1] = indx[j - 1];
                    i = j;
                    j += j;
                } else {
                    j = ir + 1;
                }
                continue;
            }
            break;
        }
        indx[i - 1] = indxt;
    }
}

/// `sort` of `interp.F90`: `(wksp, iwksp)` where `wksp[j] = ra[iwksp(j)]`
/// (ascending `wksp`, `iwksp` the zero-based permutation).
pub fn sort(arr: &[f64]) -> (Vec<f64>, Vec<usize>) {
    let perm = indexx(arr);
    let sorted = perm.iter().map(|&i| arr[i]).collect();
    (sorted, perm)
}

// ----------------------------------------------------------------------
// Point-in-polygon and cell interpolation weights
// ----------------------------------------------------------------------

/// `ipon` of `interp.F90`: point-in-polygon. Returns `-1` (outside),
/// `0` (on an edge) or `1` (inside), with the exact same `real*4`
/// coordinate truncation and `1.0E-8` edge tolerance as the Fortran.
pub fn ipon(poly_x: &[f64], poly_y: &[f64], qx: f64, qy: f64) -> i32 {
    let n = poly_x.len();
    // x(i) = real(xq(i)-xp, 4): double -> single truncation of the
    // polygon vertex minus the query point.
    let x: Vec<f32> = poly_x.iter().map(|&v| (v - qx) as f32).collect();
    let y: Vec<f32> = poly_y.iter().map(|&v| (v - qy) as f32).collect();
    let mut nunder = 0i32;
    for i in 0..n {
        let (xi, xi1) = (x[i], x[(i + 1) % n]);
        let (yi, yi1) = (y[i], y[(i + 1) % n]);
        if (xi < 0.0 && xi1 >= 0.0) || (xi1 < 0.0 && xi >= 0.0) {
            if yi < 0.0 && yi1 < 0.0 {
                nunder += 1;
            } else if (yi <= 0.0 && yi1 >= 0.0) || (yi1 <= 0.0 && yi >= 0.0) {
                let ysn = (yi * xi1 - xi * yi1) / (xi1 - xi);
                if ysn < 0.0 {
                    nunder += 1;
                } else if ysn <= 0.0 {
                    return 0; // edge
                }
            }
        } else if xi.abs() < 1.0e-8 && xi1.abs() < 1.0e-8 {
            if (yi <= 0.0 && yi1 >= 0.0) || (yi1 <= 0.0 && yi >= 0.0) {
                return 0; // edge
            }
        }
    }
    if nunder % 2 == 0 {
        -1
    } else {
        1
    }
}

/// `bilin5` of `interp.F90`: bilinear interpolation weights for a
/// (possibly distorted) quadrangle. Returns `None` on the `ier = 1` error
/// paths, leaving the caller's weights untouched — exactly the Fortran
/// behaviour (`goto 99999` skips the `w` assignment).
pub fn bilin5(xa: &[f64; 4], ya: &[f64; 4], x0: f64, y0: f64) -> Option<[f64; 4]> {
    // The tolerance literals in the Fortran (`1.0e-20`, `1.0e-7`, `1.0e-6`)
    // are single-precision and are widened to real*8 for the comparisons;
    // `1.0d0` is real*8. Reproduce the widened values exactly.
    let eps20 = 1.0e-20f32 as f64;
    let eps7 = 1.0e-7f32 as f64;
    let eps6 = 1.0e-6f32 as f64;
    let (x1, x2, x3, x4) = (xa[0], xa[1], xa[2], xa[3]);
    let (y1, y2, y3, y4) = (ya[0], ya[1], ya[2], ya[3]);
    let (x, y) = (x0, y0);

    let a21 = x2 - x1;
    let a22 = y2 - y1;
    let a31 = x3 - x1;
    let a32 = y3 - y1;
    let a41 = x4 - x1;
    let a42 = y4 - y1;
    let det = a21 * a42 - a22 * a41;
    if det.abs() < eps20 {
        return None;
    }
    let x3t = (a42 * a31 - a41 * a32) / det;
    let y3t = (-a22 * a31 + a21 * a32) / det;
    let xt = (a42 * (x - x1) - a41 * (y - y1)) / det;
    let yt = (-a22 * (x - x1) + a21 * (y - y1)) / det;
    if x3t < 0.0 || y3t < 0.0 {
        return None;
    }
    let (xi, eta);
    if (x3t - 1.0).abs() < eps7 {
        xi = xt;
        if (y3t - 1.0).abs() < eps7 {
            eta = yt;
        } else if (1.0 + (y3t - 1.0) * xt).abs() < eps6 {
            return None;
        } else {
            eta = yt / (1.0 + (y3t - 1.0) * xt);
        }
    } else if (y3t - 1.0).abs() < eps6 {
        eta = yt;
        if (1.0 + (x3t - 1.0) * yt).abs() < eps6 {
            return None;
        } else {
            xi = xt / (1.0 + (x3t - 1.0) * yt);
        }
    } else {
        let a = y3t - 1.0;
        let b = 1.0 + (x3t - 1.0) * yt - (y3t - 1.0) * xt;
        let c = -xt;
        let discr = b * b - 4.0 * a * c;
        if discr < eps6 {
            return None;
        }
        xi = (-b + discr.sqrt()) / (2.0 * a);
        eta = ((y3t - 1.0) * (xi - xt) + (x3t - 1.0) * yt) / (x3t - 1.0);
    }
    let mut w = [0.0f64; 4];
    w[0] = (1.0 - xi) * (1.0 - eta);
    w[1] = xi * (1.0 - eta);
    w[2] = xi * eta;
    w[3] = eta * (1.0 - xi);
    Some(w)
}

/// `triangle_intp` of `interp.F90`: barycentric weights of a point in a
/// triangle (no degenerate-denominator guard, exactly like the Fortran).
pub fn triangle_intp(x: &[f64; 3], y: &[f64; 3], xt: f64, yt: f64) -> [f64; 3] {
    let den = (y[1] - y[2]) * (x[0] - x[2]) + (x[2] - x[1]) * (y[0] - y[2]);
    let mut w = [0.0f64; 3];
    w[0] = ((y[1] - y[2]) * (xt - x[2]) + (x[2] - x[1]) * (yt - y[2])) / den;
    w[1] = ((y[2] - y[0]) * (xt - x[2]) + (x[0] - x[2]) * (yt - y[2])) / den;
    w[2] = 1.0 - w[0] - w[1];
    w
}

// ----------------------------------------------------------------------
// make_map_fm (unstructured mesh -> arbitrary points)
// ----------------------------------------------------------------------

/// Result of [`make_map_fm`]: `w(4, n2)` and `iref(4, n2)` flattened
/// point-major (`[ip + 4*i2]`), plus `nref(no_nodes)`.
#[derive(Debug, Clone, PartialEq)]
pub struct MapFmResult {
    /// `w(4, n2)`, `[ip + 4*i2]`; `real*8`.
    pub w: Vec<f64>,
    /// `iref(4, n2)`, `[ip + 4*i2]`; node indices (1-based) or `0`.
    pub iref: Vec<i32>,
    /// `nref(no_nodes)`; per-node reference count.
    pub nref: Vec<i32>,
}

/// `make_map_fm` of `interp.F90`: interpolation weights and node references
/// from an unstructured mesh (triangles and/or quads, `face_nodes(4,
/// no_faces)` node-major) to arbitrary points. `x1`/`y1` are node
/// coordinates (`real*8`); `face_nodes` carries the `-999` "no fourth node"
/// sentinel. `iref`/`w`/`nref` all start zeroed, so a point never covered by
/// a valid face keeps zero weights and zero references.
pub fn make_map_fm(
    x1: &[f64],
    y1: &[f64],
    face_nodes: &[i32],
    no_faces: usize,
    x2: &[f64],
    y2: &[f64],
) -> MapFmResult {
    let no_nodes1 = x1.len();
    let n2 = x2.len();

    let (xs, nrx) = sort(x2);
    let (ys, nry) = sort(y2);

    let mut iflag = vec![0i32; n2];
    let mut nrin = vec![0i32; n2];
    let mut w = vec![0.0f64; 4 * n2];
    let mut iref = vec![0i32; 4 * n2];
    let mut nref = vec![0i32; no_nodes1];

    let mut lomnx = 1usize; // 1-based `jlo`
    let mut lomny = 1usize;

    for i1 in 0..no_faces {
        let node_at = |ip: usize| face_nodes[i1 * 4 + ip] as usize; // 1-based
        let np = if face_nodes[i1 * 4 + 3] == -999 { 3 } else { 4 };

        let mut xp = [0.0f64; 4];
        let mut yp = [0.0f64; 4];
        for ip in 0..np {
            let k1 = node_at(ip) - 1; // 0-based node
            xp[ip] = x1[k1];
            yp[ip] = y1[k1];
        }

        let xpmin = xp[..np].iter().copied().fold(f64::INFINITY, f64::min);
        let xpmax = xp[..np].iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let ypmin = yp[..np].iter().copied().fold(f64::INFINITY, f64::min);
        let ypmax = yp[..np].iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let xpmean = 0.5 * (xpmin + xpmax);
        let ypmean = 0.5 * (ypmin + ypmax);

        hunt(&xs, xpmean, &mut lomnx);
        hunt(&ys, ypmean, &mut lomny);

        // For points with X between xpmin and xpmax, set iflag.
        let mut lominx = lomnx;
        let mut lomaxx = lomnx;
        // do i = lomnx, 1, -1
        let mut i = lomnx as isize;
        loop {
            if i < 1 {
                break;
            }
            let idx = i as usize - 1;
            if xs[idx] >= xpmin {
                lominx = idx + 1;
                iflag[nrx[idx]] = 1;
            } else {
                break;
            }
            i -= 1;
        }
        // do i = lomnx + 1, n2
        for i in (lomnx + 1)..=n2 {
            if xs[i - 1] <= xpmax {
                lomaxx = i;
                iflag[nrx[i - 1]] = 1;
            } else {
                break;
            }
        }

        // For points with Y between ypmin and ypmax (and iflag set), store.
        let mut iin = 1usize;
        // do i = lomny, 1, -1
        let mut i = lomny as isize;
        loop {
            if i < 1 {
                break;
            }
            let idx = i as usize - 1;
            if ys[idx] >= ypmin {
                let p = nry[idx]; // 0-based point index
                nrin[iin - 1] = (p as i32 + 1) * iflag[p];
                iin += iflag[p] as usize;
            } else {
                break;
            }
            i -= 1;
        }
        // do i = lomny + 1, n2
        for i in (lomny + 1)..=n2 {
            if ys[i - 1] <= ypmax {
                let p = nry[i - 1];
                nrin[iin - 1] = (p as i32 + 1) * iflag[p];
                iin += iflag[p] as usize;
            } else {
                break;
            }
        }
        let nin = iin - 1;

        // Reset iflag.
        for i in lominx..=lomaxx {
            if i != 0 {
                iflag[nrx[i - 1]] = 0;
            }
        }

        for iin_idx in 1..=nin {
            let i2 = nrin[iin_idx - 1] as usize; // 1-based point index
            if i2 == 0 {
                continue;
            }
            let i2_0 = i2 - 1;
            let inout = ipon(&xp[..np], &yp[..np], x2[i2_0], y2[i2_0]);
            if inout < 0 {
                continue;
            }
            if np == 4 {
                // `ier` is ignored: on error the weights stay as-is (0).
                if let Some(ww) = bilin5(&xp, &yp, x2[i2_0], y2[i2_0]) {
                    for ip in 0..4 {
                        w[ip + 4 * i2_0] = ww[ip];
                    }
                }
                for ip in 0..4 {
                    let node = node_at(ip) as i32;
                    iref[ip + 4 * i2_0] = node;
                    nref[node as usize - 1] += 1;
                }
            } else {
                let ww = triangle_intp(&[xp[0], xp[1], xp[2]], &[yp[0], yp[1], yp[2]], x2[i2_0], y2[i2_0]);
                for ip in 0..3 {
                    w[ip + 4 * i2_0] = ww[ip];
                    let node = node_at(ip) as i32;
                    iref[ip + 4 * i2_0] = node;
                    nref[node as usize - 1] += 1;
                }
                w[3 + 4 * i2_0] = 0.0;
                iref[3 + 4 * i2_0] = 1;
            }
        }
    }

    MapFmResult { w, iref, nref }
}

// ----------------------------------------------------------------------
// Trapezoidal / cyclic integration (Delft3D waveparams, dead in SnapWave)
// ----------------------------------------------------------------------

/// `trapezoidal` of `interp.F90`: integral of `y(x)` from `x1` to `x2` on
/// an equidistant `x` grid; `0.0` when the interval is fully outside.
pub fn trapezoidal(x: &[f64], y: &[f64], x1_in: f64, x2_in: f64) -> f64 {
    let n = x.len();
    if n == 0 {
        return 0.0;
    }
    if x1_in > x[n - 1] || x2_in < x[0] {
        return 0.0;
    }
    let x1 = x1_in.max(x[0] + 1.0e-60);
    let x2 = x2_in.min(x[n - 1] - 1.0e-60);
    let dx = (x[n - 1] - x[0]) / (n as f64 - 1.0);
    // floor((x1-x(1))/dx)+1
    let i1 = ((x1 - x[0]) / dx).floor() as isize + 1;
    let i2 = ((x2 - x[0]) / dx).floor() as isize + 1;
    let i1p1 = (i1 + 1).min(n as isize);
    let i2p1 = (i2 + 1).min(n as isize);
    let at = |i: isize| y[(i - 1) as usize];
    let xat = |i: isize| x[(i - 1) as usize];

    // first partial trapezoid
    let y1 = at(i1) + (x1 - xat(i1)) / dx * (at(i1p1) - at(i1));
    let ifirst = 0.5 * (xat(i1p1) - x1) * (at(i1p1) + y1);
    // middle part
    let mut imid = 0.5 * at(i1p1);
    let mut i = i1 + 2;
    while i <= i2 - 1 {
        imid += at(i);
        i += 1;
    }
    imid += 0.5 * at(i2);
    imid *= dx;
    // last partial trapezoid
    let y2 = at(i2) + (x2 - xat(i2)) / dx * (at(i2p1) - at(i2));
    let ilast = 0.5 * (x2 - xat(i2)) * (y2 + at(i2));
    ifirst + imid + ilast
}

/// `interp_in_cyclic_function` of `interp.F90`: interpolate a cyclic
/// function (`x` equidistant, period `xcycle`) at points `xp`.
pub fn interp_in_cyclic_function(x: &[f64], y: &[f64], xcycle: f64, xp: &[f64]) -> Vec<f64> {
    let n = x.len();
    let dx = x[1] - x[0];
    let icycle = (xcycle / dx).round() as usize;
    let mut xc = vec![0.0f64; icycle + 1];
    let mut yc = vec![0.0f64; icycle + 1];
    if n > icycle + 1 {
        // nonsense; error (Fortran leaves xc/yc empty — replicate with
        // no valid data, but the branch is unreachable for sane input).
    } else if n == icycle + 1 {
        xc[..n].copy_from_slice(x);
        yc[..n].copy_from_slice(y);
    } else if n == icycle {
        xc[..n].copy_from_slice(x);
        yc[..n].copy_from_slice(y);
        xc[n] = xc[n - 1] + dx;
        yc[n] = yc[0];
    } else {
        xc[..n].copy_from_slice(x);
        yc[..n].copy_from_slice(y);
        for i in n + 1..=icycle {
            xc[i - 1] = xc[i - 2] + dx;
            yc[i - 1] = 0.0;
        }
        xc[icycle] = xc[icycle - 1] + dx;
        yc[icycle] = yc[0];
    }
    let nc = icycle + 1;

    xp.iter()
        .map(|&xpp| {
            let mut ileft = ((xpp - xc[0]) / dx).floor() as isize + 1;
            while ileft < 1 {
                ileft += icycle as isize;
            }
            while ileft > nc as isize {
                ileft -= icycle as isize;
            }
            let yleft = if ileft > nc as isize || ileft < 1 { 0.0 } else { yc[(ileft - 1) as usize] };
            let mut iright = ileft + 1;
            while iright < 1 {
                iright += icycle as isize;
            }
            while iright > nc as isize {
                iright -= icycle as isize;
            }
            let yright = if iright > nc as isize || iright < 1 { 0.0 } else { yc[(iright - 1) as usize] };
            let facright = (xpp % dx) / dx;
            let facleft = 1.0 - facright;
            facleft * yleft + facright * yright
        })
        .collect()
}

/// `trapezoidal_cyclic` of `interp.F90`: integral over `[x1, x2]` of a
/// cyclic/linear function.
pub fn trapezoidal_cyclic(x: &[f64], y: &[f64], xcycle: f64, x1: f64, x2: f64) -> f64 {
    let dx = x[1] - x[0];
    let np = (x2 / dx).floor() as isize - ((x1 / dx).floor() as isize + 1) + 2;
    let np = np as usize;
    let mut xp = vec![0.0f64; np];
    let yp: Vec<f64>;
    xp[0] = x1;
    for ip in 2..np {
        xp[ip - 1] = ((x1 / dx).floor() as isize + (ip as isize) - 1) as f64 * dx;
    }
    xp[np - 1] = x2;
    if xcycle > 0.0 {
        yp = interp_in_cyclic_function(x, y, xcycle, &xp);
    } else {
        yp = xp.iter().map(|&v| linear_interp(x, y, v).0).collect();
    }
    let mut integ = 0.0;
    for ip in 1..np {
        integ += 0.5 * (xp[ip] - xp[ip - 1]) * (yp[ip] + yp[ip - 1]);
    }
    integ
}

/// `interp_using_trapez_rule` of `interp.F90`: cell averages of a cyclic
/// function on the target grid `xtarg`.
pub fn interp_using_trapez_rule(x: &[f64], y: &[f64], xcycle: f64, xtarg: &[f64]) -> Vec<f64> {
    let ntarg = xtarg.len();
    if ntarg < 2 {
        return vec![0.0; ntarg];
    }
    let dx = xtarg[1] - xtarg[0];
    (0..ntarg)
        .map(|itarg| {
            let x1 = xtarg[itarg] - 0.5 * dx;
            let x2 = xtarg[itarg] + 0.5 * dx;
            let integ = trapezoidal_cyclic(x, y, xcycle, x1, x2);
            integ / dx
        })
        .collect()
}

// ----------------------------------------------------------------------
// Curvilinear-grid mapping (Delft3D, dead in SnapWave)
// ----------------------------------------------------------------------

/// `grmap` of `interp.F90`: interpolate `f1` (grid 1) onto `f2` (grid 2)
/// using the reference table `iref(np, n2)` and weights `w(np, n2)`; points
/// with `iref(1, i2) <= 0` keep their previous `f2` value.
pub fn grmap(f1: &[f32], f2: &mut [f32], iref: &[i32], w: &[f64], np: usize) {
    let n2 = f2.len();
    for i2 in 0..n2 {
        let i = iref[i2];
        if i > 0 {
            // Fortran `f2(i2) = f2(i2) + w*f1` truncates to real*4 on every
            // step, so accumulate in f32 (not f64 then a single cast).
            let mut acc = 0.0f32;
            for ip in 0..np {
                let ii = iref[ip * n2 + i2];
                let i1 = ii.max(1) as usize;
                acc = (acc as f64 + w[ip * n2 + i2] * f1[i1 - 1] as f64) as f32;
            }
            f2[i2] = acc;
        }
    }
}

/// `grmap2` of `interp.F90`: the inverse mapping — accumulate grid-2 values
/// `f2` back onto grid 1 `f1` using the same weights (area-weighted).
pub fn grmap2(f1: &mut [f64], cellsz1i: &[f64], f2: &[f64], cellsz2: f64, iref: &[i32], w: &[f64], np: usize) {
    let n2 = f2.len();
    for ip in 0..np {
        for i2 in 0..n2 {
            let i1 = iref[ip * n2 + i2];
            if i1 > 0 {
                f1[i1 as usize - 1] += w[ip * n2 + i2] * f2[i2] * cellsz2 * cellsz1i[i1 as usize - 1];
            }
        }
    }
}

/// `grmap_sg` of `interp.F90`: like [`grmap2`] but normalised by the node
/// reference count `nref` instead of a cell size.
pub fn grmap_sg(f1: &mut [f64], nref: &[i32], f2: &[f64], iref: &[i32], w: &[f64], np: usize) {
    let n2 = f2.len();
    for v in f1.iter_mut() {
        *v = 0.0;
    }
    for i2 in 0..n2 {
        for ip in 0..np {
            let i1 = iref[ip * n2 + i2];
            if i1 > 0 {
                f1[i1 as usize - 1] += f2[i2] / nref[i1 as usize - 1] as f64;
            }
        }
    }
    // `w` is accepted but unused, exactly like the Fortran (whose weighted
    // branch is commented out in favour of the 1/nref normalisation).
    let _ = w;
}

/// `mkmap_step` of `interp.F90`: update weights/references from a previous
/// cell using the local grid directions `alfaz`/`dsu`/`dnv` (ordered grid).
pub fn mkmap_step(
    x1: &[f64],
    y1: &[f64],
    m1: usize,
    alfaz: &[f64],
    dsu: &[f64],
    dnv: &[f64],
    x2: &[f64],
    y2: &[f64],
    iref: &mut [i32],
    w: &mut [f64],
) {
    let n2 = x2.len();
    let at = |i: isize, j: isize| (j - 1) as usize * m1 + (i - 1) as usize;
    for i2 in 0..n2 {
        if iref[i2] > 0 {
            let i1 = ((iref[i2] - 1) % (m1 as i32)) as isize + 1;
            let j1 = (iref[i2] as isize - i1) / (m1 as isize) + 1;
            let dx = x2[i2] - x1[at(i1, j1)];
            let dy = y2[i2] - y1[at(i1, j1)];
            let ds = dx * alfaz[at(i1, j1)].cos() + dy * alfaz[at(i1, j1)].sin();
            let dn = -dx * alfaz[at(i1, j1)].sin() + dy * alfaz[at(i1, j1)].cos();
            let dsrel = ds / dsu[at(i1, j1)];
            let dnrel = dn / dnv[at(i1, j1)];
            let i1 = i1 + dsrel.floor() as isize;
            let j1 = j1 + dnrel.floor() as isize;
            let dsrel = dsrel - dsrel.floor();
            let dnrel = dnrel - dnrel.floor();
            let i1c = i1 as i32;
            let j1c = j1 as i32;
            let m1c = m1 as i32;
            iref[0 * n2 + i2] = i1c + (j1c - 1) * m1c;
            iref[1 * n2 + i2] = i1c + 1 + (j1c - 1) * m1c;
            iref[2 * n2 + i2] = i1c + 1 + j1c * m1c;
            iref[3 * n2 + i2] = i1c + j1c * m1c;
            w[0 * n2 + i2] = (1.0 - dsrel) * (1.0 - dnrel);
            w[1 * n2 + i2] = dsrel * (1.0 - dnrel);
            w[2 * n2 + i2] = dsrel * dnrel;
            w[3 * n2 + i2] = (1.0 - dsrel) * dnrel;
        }
    }
}

/// Result of [`make_map`] (curvilinear ordered grid -> arbitrary points).
#[derive(Debug, Clone, PartialEq)]
pub struct MakeMapResult {
    /// `w(4, n2)`, `[ip + 4*i2]`; single precision (Fortran `real(hbuff,4)`).
    pub w: Vec<f32>,
    /// `iref(4, n2)`, `[ip + 4*i2]`; 1-based linear grid-1 indices or `0`.
    pub iref: Vec<i32>,
    /// `covered(n2)`: 1 valid, 0 boundary-adjacent, -1 invalid.
    pub covered: Vec<i32>,
}

/// `make_map` of `interp.F90`: bilinear interpolation weights from an
/// ordered curvilinear grid `x1(m1, n1)` (column-major, index `(i-1) +
/// m1*(j-1)`) to arbitrary points, with a `code` mask distinguishing valid
/// (1/3), boundary (>0) and invalid (<=0) cells.
pub fn make_map(
    code: &[i32],
    x1: &[f64],
    y1: &[f64],
    m1: usize,
    n1: usize,
    x2: &[f64],
    y2: &[f64],
    xymiss: f32,
) -> MakeMapResult {
    let n2 = x2.len();
    let at = |i: usize, j: usize| (j - 1) * m1 + (i - 1);
    let code_at = |i: usize, j: usize| code[(j - 1) * m1 + (i - 1)];

    let (xs, nrx) = sort(x2);
    let (ys, nry) = sort(y2);

    let mut iflag = vec![0i32; n2];
    let mut nrin = vec![0i32; n2];
    let mut w = vec![0.0f32; 4 * n2];
    let mut iref = vec![0i32; 4 * n2];
    let mut covered = vec![0i32; n2];

    let mut lomnx = 1usize;
    let mut lomny = 1usize;

    let xymiss = xymiss as f64; // real*4 missing value widened for the checks
    let eps = 0.00001f32 as f64; // `eps = 0.00001` (single literal -> real*8)

    for j1 in 1..n1 {
        for i1 in 1..m1 {
            let xp = [x1[at(i1, j1)], x1[at(i1 + 1, j1)], x1[at(i1 + 1, j1 + 1)], x1[at(i1, j1 + 1)]];
            let yp = [y1[at(i1, j1)], y1[at(i1 + 1, j1)], y1[at(i1 + 1, j1 + 1)], y1[at(i1, j1 + 1)]];
            let has_miss = xp.iter().chain(yp.iter()).any(|&v| v > xymiss - eps && v < xymiss + eps);
            if has_miss {
                continue;
            }
            let xpmin = xp.iter().copied().fold(f64::INFINITY, f64::min);
            let xpmax = xp.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let ypmin = yp.iter().copied().fold(f64::INFINITY, f64::min);
            let ypmax = yp.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let xpmean = 0.5 * (xpmin + xpmax);
            let ypmean = 0.5 * (ypmin + ypmax);

            hunt(&xs, xpmean, &mut lomnx);
            hunt(&ys, ypmean, &mut lomny);

            let mut lominx = lomnx;
            let mut lomaxx = lomnx;
            let mut i = lomnx as isize;
            loop {
                if i < 1 {
                    break;
                }
                let idx = i as usize - 1;
                if xs[idx] >= xpmin {
                    lominx = idx + 1;
                    iflag[nrx[idx]] = 1;
                } else {
                    break;
                }
                i -= 1;
            }
            for i in (lomnx + 1)..=n2 {
                if xs[i - 1] <= xpmax {
                    lomaxx = i;
                    iflag[nrx[i - 1]] = 1;
                } else {
                    break;
                }
            }

            let mut iin = 1usize;
            let mut i = lomny as isize;
            loop {
                if i < 1 {
                    break;
                }
                let idx = i as usize - 1;
                if ys[idx] >= ypmin {
                    let p = nry[idx];
                    nrin[iin - 1] = (p as i32 + 1) * iflag[p];
                    iin += iflag[p] as usize;
                } else {
                    break;
                }
                i -= 1;
            }
            for i in (lomny + 1)..=n2 {
                if ys[i - 1] <= ypmax {
                    let p = nry[i - 1];
                    nrin[iin - 1] = (p as i32 + 1) * iflag[p];
                    iin += iflag[p] as usize;
                } else {
                    break;
                }
            }
            let nin = iin - 1;

            for i in lominx..=lomaxx {
                if i != 0 {
                    iflag[nrx[i - 1]] = 0;
                }
            }

            for iin_idx in 1..=nin {
                let i2 = nrin[iin_idx - 1] as usize;
                if i2 == 0 {
                    continue;
                }
                let i2_0 = i2 - 1;
                let inout = ipon(&xp, &yp, x2[i2_0], y2[i2_0]);
                if inout < 0 {
                    continue;
                }
                let c = [code_at(i1, j1), code_at(i1 + 1, j1), code_at(i1 + 1, j1 + 1), code_at(i1, j1 + 1)];
                let all_valid = c.iter().all(|&v| v == 1 || v == 3);
                let all_positive = c.iter().all(|&v| v > 0);
                if all_valid {
                    covered[i2_0] = 1;
                    if let Some(h) = bilin5(&xp, &yp, x2[i2_0], y2[i2_0]) {
                        for ip in 0..4 {
                            w[ip + 4 * i2_0] = h[ip] as f32;
                        }
                    }
                    let m1c = m1 as i32;
                    let (i1c, j1c) = (i1 as i32, j1 as i32);
                    iref[0 + 4 * i2_0] = i1c + (j1c - 1) * m1c;
                    iref[1 + 4 * i2_0] = i1c + 1 + (j1c - 1) * m1c;
                    iref[2 + 4 * i2_0] = i1c + 1 + j1c * m1c;
                    iref[3 + 4 * i2_0] = i1c + j1c * m1c;
                } else if all_positive {
                    covered[i2_0] = 0;
                } else {
                    covered[i2_0] = -1;
                }
            }
        }
    }

    MakeMapResult { w, iref, covered }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn binary_search_brackets_ascending_and_descending() {
        // Returns the Fortran 1-based `j`: the lower bracket of the interval
        // `xx(j) ..= xx(j+1)` holding `x`, with 0/n sentinels at the edges.
        let x = [0.0, 1.0, 2.0, 3.0, 4.0];
        assert_eq!(binary_search(&x, 0.5), 1); // 0.0 < 0.5 <= 1.0
        assert_eq!(binary_search(&x, 1.0), 1); // 0.0 < 1.0 <= 1.0
        assert_eq!(binary_search(&x, 2.5), 3); // 2.0 < 2.5 <= 3.0
        assert_eq!(binary_search(&x, 3.0), 3); // 2.0 < 3.0 <= 3.0
        assert_eq!(binary_search(&x, 99.0), 5); // above the last element
        assert_eq!(binary_search(&x, -1.0), 0); // at/below the first element
        // Descending: `j` is still the lower bracket in the descending sense.
        let xd = [4.0, 3.0, 2.0, 1.0, 0.0];
        assert_eq!(binary_search(&xd, 3.5), 1); // 4.0 >= 3.5 > 3.0
        assert_eq!(binary_search(&xd, 2.5), 2); // 3.0 >= 2.5 > 2.0
        assert_eq!(binary_search(&xd, 0.5), 4); // 1.0 >= 0.5 > 0.0
    }

    #[test]
    fn linear_interp_matches_hand_computed_values() {
        let x = [0.0, 1.0, 2.0];
        let y = [0.0, 10.0, 20.0];
        let (v, j) = linear_interp(&x, &y, 1.5);
        assert_eq!(v, 15.0);
        assert_eq!(j, 2); // Fortran `indint`: bracket x(2)=1 ..= x(3)=2
        // Edge clamping.
        assert_eq!(linear_interp(&x, &y, -5.0).0, 0.0);
        assert_eq!(linear_interp(&x, &y, 99.0).0, 20.0);
        // Single element.
        assert_eq!(linear_interp(&[1.0], &[7.0], 0.0).0, 7.0);
        // Empty.
        assert_eq!(linear_interp(&[], &[], 0.0).0, 0.0);
    }

    #[test]
    fn indexx_sorts_ascending() {
        let arr = [3.0, 1.0, 2.0];
        let perm = indexx(&arr);
        let sorted: Vec<f64> = perm.iter().map(|&i| arr[i]).collect();
        assert_eq!(sorted, vec![1.0, 2.0, 3.0]);
        let (s, _) = sort(&arr);
        assert_eq!(s, vec![1.0, 2.0, 3.0]);
    }

    #[test]
    fn hunt_matches_hand_computed_brackets() {
        // `hunt` is the Numerical Recipes "locate with a guess"; its result
        // is the same 1-based lower bracket as `binary_search` for values
        // strictly inside the array, but it keeps the initial guess as a
        // bracket when `x` falls exactly on the first array element (0.0
        // below → jlo = 1, whereas binary_search reports 0). Pin the exact
        // behaviour make_map relies on, starting from the make_map guess
        // `jlo = 1`.
        let x = [0.0, 2.0, 4.0, 6.0, 8.0, 10.0];
        let cases = [
            (-1.0, 0usize),
            (0.0, 1usize), // exactly at xx(1): keeps the guess bracket
            (1.0, 1usize),
            (5.0, 3usize),
            (10.0, 5usize),
            (11.0, 6usize),
        ];
        for &(v, expect) in &cases {
            let mut jlo = 1usize;
            hunt(&x, v, &mut jlo);
            assert_eq!(jlo, expect, "hunt value {v}");
        }
    }

    #[test]
    fn ipon_classifies_points() {
        // Unit square.
        let px = [0.0, 1.0, 1.0, 0.0];
        let py = [0.0, 0.0, 1.0, 1.0];
        assert_eq!(ipon(&px, &py, 0.5, 0.5), 1); // inside
        assert_eq!(ipon(&px, &py, 2.0, 0.5), -1); // outside
        assert_eq!(ipon(&px, &py, 0.5, 0.0), 0); // on edge
    }

    #[test]
    fn bilin5_recovers_unit_square_weights() {
        let xa = [0.0, 1.0, 1.0, 0.0];
        let ya = [0.0, 0.0, 1.0, 1.0];
        // Centre: each corner gets 1/4.
        let w = bilin5(&xa, &ya, 0.5, 0.5).unwrap();
        assert!(w.iter().all(|&v| (v - 0.25).abs() < 1e-12), "weights: {w:?}");
        // At corner 1: w(1)=1.
        let w = bilin5(&xa, &ya, 0.0, 0.0).unwrap();
        assert!((w[0] - 1.0).abs() < 1e-12 && w[1..].iter().all(|&v| v.abs() < 1e-12));
    }

    #[test]
    fn triangle_intp_recovers_barycentric_weights() {
        let x = [0.0, 1.0, 0.0];
        let y = [0.0, 0.0, 1.0];
        let w = triangle_intp(&x, &y, 0.25, 0.25);
        assert_eq!(w[2], 0.25); // node 3 weight
        assert_eq!(w[1], 0.25); // node 2 weight
        assert_eq!(w[0], 0.5); // node 1 weight
        // Sum to 1.
        assert!((w.iter().sum::<f64>() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn make_map_fm_interpolates_a_single_triangle() {
        // One triangle: nodes 1(0,0) 2(1,0) 3(0,1).
        let x1 = [0.0, 1.0, 0.0];
        let y1 = [0.0, 0.0, 1.0];
        let face_nodes = [1, 2, 3, -999];
        let x2 = [0.25, 10.0];
        let y2 = [0.25, 0.0];
        let r = make_map_fm(&x1, &y1, &face_nodes, 1, &x2, &y2);
        // Point 0 inside: iref = [1,2,3,1], weights [0.5,0.25,0.25,0].
        assert_eq!(&r.iref[..4], &[1, 2, 3, 1]);
        assert!((r.w[0] - 0.5).abs() < 1e-12);
        assert!((r.w[1] - 0.25).abs() < 1e-12);
        assert!((r.w[2] - 0.25).abs() < 1e-12);
        assert_eq!(r.w[3], 0.0);
        // nref: nodes 1,2,3 each referenced once.
        assert_eq!(r.nref, vec![1, 1, 1]);
        // Point 1 outside: no reference.
        assert_eq!(&r.iref[4..8], &[0, 0, 0, 0]);
    }

    #[test]
    fn trapezoidal_computes_a_known_integral() {
        // y = x on [0,1,2,3]: integral 0..2 = 2.0.
        let x = [0.0, 1.0, 2.0, 3.0];
        let y = [0.0, 1.0, 2.0, 3.0];
        assert!((trapezoidal(&x, &y, 0.0, 2.0) - 2.0).abs() < 1e-9);
        assert_eq!(trapezoidal(&x, &y, 10.0, 20.0), 0.0);
    }

    #[test]
    fn linear_interp_2d_recovers_grid_values() {
        // Z = x + y on a 2x2 grid.
        let x = [0.0, 1.0];
        let y = [0.0, 1.0];
        let z = [0.0, 1.0, 1.0, 2.0]; // column-major: (0,0)=0,(1,0)=1,(0,1)=1,(1,1)=2
        assert_eq!(linear_interp_2d(&x, &y, &z, 0.5, 0.5, "interp", -999.0), 1.0);
        assert_eq!(linear_interp_2d(&x, &y, &z, 0.0, 0.0, "interp", -999.0), 0.0);
        assert_eq!(linear_interp_2d(&x, &y, &z, 1.0, 1.0, "interp", -999.0), 2.0);
        // Out of range: 'interp' -> exception, 'extendclosest' -> nearest.
        assert_eq!(linear_interp_2d(&x, &y, &z, 5.0, 0.0, "interp", -999.0), -999.0);
        assert_eq!(linear_interp_2d(&x, &y, &z, 5.0, 0.0, "extendclosest", -999.0), 1.0);
    }
}
