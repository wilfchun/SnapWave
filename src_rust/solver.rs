//! Rust port of the SnapWave numerical solver internals
//! (plan.md, Phase 11: "Solver Internals In Rust").
//!
//! # What is ported
//!
//! Every routine of `src/snapwave_solver.f90` and
//! `src/snapwave_windsource.f90` that is on the runtime path:
//!
//! * [`solve_tridiag`] — Thomas algorithm for tridiagonal linear systems.
//! * [`baldock`] — Baldock wave breaking dissipation model.
//! * [`hpsort_eps_epw`] — Heapsort with epsilon tolerance for sweep ordering.
//! * [`disper_nr`] — Newton-Raphson solver for the linear dispersion relation.
//! * [`compute_celerities`] — Wave celerity, group velocity, refraction speed.
//! * [`numerical_limiter`] — Depth-limits energy/action, bounds frequency.
//! * [`windinput`] — Wind growth source terms (Kahma & Calkoen / Breugem &
//!   Holthuijsen).
//! * [`vegatt`] / [`swvegatt`] / [`bulkdragcoeff`] — Vegetation dissipation
//!   (Suzuki et al. 2012, Mendez & Losada 2004 / Ozeren et al. 2013).
//! * [`compute_wave_field`] — Top-level solver orchestrator called each
//!   timestep.
//! * [`solve_energy_balance2Dstat`] — The core implicit 4-sweep solver on
//!   unstructured grids.
//!
//! # Floating-point fidelity
//!
//! These are ports, not clean reimplementations: every `real*4`/`real*8`
//! width, every operation order and every literal is preserved so the
//! outputs match the Fortran oracle. Where Fortran uses `real*4` (the
//! default real kind), Rust uses `f32`; where Fortran uses `real*8`, Rust
//! uses `f64`. The Fortran `parameter` constants (`pi`, `g`, `rho`) are
//! reproduced exactly.
//!
//! # Indexing conventions
//!
//! All public slices are zero-based. Internal loops mirror the Fortran
//! one-based indexing where it aids readability of the port, but the
//! public API is zero-based throughout.

use std::f32::consts::PI as PI_F32;

// ---------------------------------------------------------------------------
// Physical constants (mirror snapwave_data parameters)
// ---------------------------------------------------------------------------

/// Water density [kg/m³] — `snapwave_data::rho`.
pub const RHO: f32 = 1025.0;
/// Gravitational acceleration [m/s²] — `snapwave_data::g` (used by
/// `compute_wave_field` for celerities/refraction speed).
pub const G: f32 = 9.813;
/// Gravitational acceleration used *inside* `solve_energy_balance2Dstat`:
/// the Fortran routine declares its own local `g = 9.81` parameter (distinct
/// from the module `g = 9.813`), so the two must not be collapsed.
pub const G_SOLVER: f32 = 9.81;
/// π in single precision — `snapwave_data::pi`.
pub const PI: f32 = PI_F32;
/// Degrees to radians — `snapwave_data::deg2rad`.
pub const DEG2RAD: f32 = PI_F32 / 180.0;
/// Radians to degrees — `snapwave_data::rad2deg`.
pub const RAD2DEG: f32 = 180.0 / PI_F32;

/// A snapshot of the solver state at one per-iteration map output point
/// (`ja_save_each_iter == 1`). Enough for the model to build one
/// [`crate::capture::MapRecord`]; the static fields (depth, `theta`, wind)
/// are already owned by the model.
#[derive(Clone, Debug, Default)]
pub struct IterSnapshot {
    pub h: Vec<f32>,
    pub thetam: Vec<f32>,
    pub df: Vec<f32>,
    pub dw: Vec<f32>,
    pub f: Vec<f32>,
    pub h_ig: Vec<f32>,
    pub tp: Vec<f32>,
    pub sig: Vec<f32>,
    pub swe: Vec<f32>,
    pub swa: Vec<f32>,
    pub ee: Vec<f32>,
    pub ctheta: Vec<f32>,
    pub cg: Vec<f32>,
}

/// One per-iteration map output event: the record index (`ntmapout`, the
/// Fortran `iter`) and the output time (`time + iter`).
#[derive(Clone, Debug)]
pub struct IterOutput {
    pub ntmapout: i32,
    pub time: f64,
    pub snapshot: IterSnapshot,
}

// ---------------------------------------------------------------------------
// 1. solve_tridiag — Thomas algorithm
// ---------------------------------------------------------------------------

/// Thomas algorithm for solving a tridiagonal linear system.
///
/// Port of `solve_tridiag` in `src/snapwave_solver.f90` (lines 1014–1050).
///
/// Solves `A x = d` where:
/// * `a` — sub-diagonal (below the main diagonal), length `n`
/// * `b` — main diagonal, length `n`
/// * `c` — super-diagonal (above the main diagonal), length `n`
/// * `d` — right-hand side, length `n`
///
/// Returns the solution `x` of length `n`.
///
/// # Panics
///
/// Panics if `n == 0` or if any of the input slices have length < `n`.
pub fn solve_tridiag(a: &[f32], b: &[f32], c: &[f32], d: &[f32], n: usize) -> Vec<f32> {
    assert!(n > 0, "solve_tridiag: n must be > 0");
    let mut cp = vec![0.0f32; n];
    let mut dp = vec![0.0f32; n];
    let mut x = vec![0.0f32; n];

    // initialize c-prime and d-prime
    cp[0] = c[0] / b[0];
    dp[0] = d[0] / b[0];

    // solve for vectors c-prime and d-prime
    for i in 1..n {
        let m = b[i] - cp[i - 1] * a[i];
        cp[i] = c[i] / m;
        dp[i] = (d[i] - dp[i - 1] * a[i]) / m;
    }

    // initialize x
    x[n - 1] = dp[n - 1];

    // solve for x from the vectors c-prime and d-prime
    for i in (0..n - 1).rev() {
        x[i] = dp[i] - cp[i] * x[i + 1];
    }

    x
}

// ---------------------------------------------------------------------------
// 2. baldock — Baldock wave breaking dissipation
// ---------------------------------------------------------------------------

/// Baldock wave breaking dissipation model.
///
/// Port of `baldock` in `src/snapwave_solver.f90` (lines 1052–1075).
///
/// Computes the dissipation due to depth-induced wave breaking.
///
/// * `opt == 1`: `Dw = 0.28 * alfa * rho * g / T * exp(-(Hmax/Hloc)²) * (Hmax² + Hloc²)`
/// * `opt == 2`: `Dw = 0.28 * alfa * rho * g / T * exp(-(Hmax/Hloc)²) * (Hmax³ + Hloc³) / gamma / depth`
pub fn baldock(rho: f32, g: f32, alfa: f32, gamma: f32, depth: f32, h: f32, t: f32, opt: i32, hmax: f32) -> f32 {
    let hloc = h.max(1.0e-6);
    let ratio = hmax / hloc;
    let exp_term = (-ratio * ratio).exp();

    if opt == 1 {
        0.28 * alfa * rho * g / t * exp_term * (hmax * hmax + hloc * hloc)
    } else {
        0.28 * alfa * rho * g / t * exp_term * (hmax.powi(3) + hloc.powi(3)) / gamma / depth
    }
}

// ---------------------------------------------------------------------------
// 3. hpsort_eps_epw — Heapsort with epsilon tolerance
// ---------------------------------------------------------------------------

/// Heapsort with epsilon tolerance for sweep-direction ordering.
///
/// Port of `hpsort_eps_epw` in `src/snapwave_solver.f90` (lines 1077–1216).
///
/// Sorts `ra` in-place into ascending order, considering two elements equal
/// if their values differ by less than `eps`. The index array `ind` tracks
/// the original positions.
///
/// If `ind[0] == 0` on input, indices are initialized to `1..n`; otherwise
/// the existing indices are carried through the sort.
pub fn hpsort_eps_epw(n: usize, ra: &mut [f32], ind: &mut [i32], eps: f32) {
    // initialize index array
    if ind[0] == 0 {
        for i in 0..n {
            ind[i] = (i + 1) as i32;
        }
    }

    if n < 2 {
        return;
    }

    let mut l = n / 2;
    let mut ir = n - 1; // zero-based: ir = n - 1

    loop {
        let (rra, iind): (f32, i32);
        if l > 0 {
            // still in hiring phase
            l -= 1;
            rra = ra[l];
            iind = ind[l];
        } else {
            // in retirement-promotion phase
            rra = ra[ir];
            iind = ind[ir];
            // retire the top of the heap into it
            ra[ir] = ra[0];
            ind[ir] = ind[0];
            // decrease the size of the corporation
            if ir == 0 {
                // the least competent worker at all!
                ra[0] = rra;
                ind[0] = iind;
                break;
            }
            ir -= 1;
        }

        // set up to place rra in its proper level
        let mut i = l;
        let mut j = 2 * l + 1; // zero-based: j = l + l + 1

        while j <= ir {
            if j < ir {
                // compare to better underling
                if hslt(ra[j], ra[j + 1], eps) {
                    j += 1;
                }
            }
            // demote rra
            if hslt(rra, ra[j], eps) {
                ra[i] = ra[j];
                ind[i] = ind[j];
                i = j;
                j = 2 * j + 1; // zero-based: j = j + j + 1
            } else {
                // set j to terminate do-while loop
                j = ir + 1;
            }
        }
        ra[i] = rra;
        ind[i] = iind;
    }
}

/// Internal comparison function: returns true if `a < b` (beyond epsilon).
fn hslt(a: f32, b: f32, eps: f32) -> bool {
    if (a - b).abs() < eps {
        false
    } else {
        a < b
    }
}

// ---------------------------------------------------------------------------
// 4. disper_nr — Newton-Raphson dispersion relation solver
// ---------------------------------------------------------------------------

/// Newton-Raphson solver for the linear dispersion relation ω² = gk tanh(kh).
///
/// Port of `disper_nr` in `src/snapwave_windsource.f90` (lines 214–252).
///
/// Returns `(k, cg)` where `k` is the wavenumber [rad/m] and `cg` is the
/// group velocity [m/s].
pub fn disper_nr(h: f32, t: f32) -> (f32, f32) {
    let tol: f32 = 1.0e-5;
    let max_iter = 20;

    // Calculate angular frequency from period
    let omega = 2.0 * PI_F32 / t;

    // Initial guess using deep water approximation: k = ω²/g
    let mut k = omega * omega / G;

    // Newton-Raphson iteration
    for _ in 0..max_iter {
        let kh = k * h;
        let tanh_kh = kh.tanh();
        let f = G * k * tanh_kh - omega * omega; // residual
        let df = G * (tanh_kh + kh * (1.0 - tanh_kh * tanh_kh)); // derivative
        let k_old = k;
        k = k - f / df;

        // Check convergence
        if (k - k_old).abs() < tol * k.abs() {
            break;
        }
    }

    let kh = k * h;
    let c = omega / k;
    let arg = (2.0 * kh).min(50.0);
    let n = 0.5 + kh / arg.sinh();
    let cg = n * c;

    (k, cg)
}

// ---------------------------------------------------------------------------
// 5. compute_celerities — Wave celerity and refraction speed
// ---------------------------------------------------------------------------

/// Compute wave celerities and refraction speed for a single node.
///
/// Port of `compute_celerities` in `src/snapwave_windsource.f90` (lines 67–94).
///
/// Called per-node during the iteration loop when wind is active.
///
/// Returns `(sinhkh, hmx, kwav, cg, ctheta)` where `ctheta` has length `ntheta`.
pub fn compute_celerities(
    hh: f32,
    sig: f32,
    sinth: &[f32],
    costh: &[f32],
    ntheta: usize,
    gamma: f32,
    dhdx: f32,
    dhdy: f32,
) -> (f32, f32, f32, f32, Vec<f32>) {
    // Compute celerities and refraction speed
    let (kwav, cg) = disper_nr(hh, 2.0 * PI_F32 / sig);

    let sinhkh = (kwav * hh).min(50.0).sinh();
    let hmx = gamma * hh;

    let mut ctheta = vec![0.0f32; ntheta];
    for itheta in 0..ntheta {
        ctheta[itheta] = sig / (2.0 * kwav * hh).min(50.0).sinh()
            * (dhdx * sinth[itheta] - dhdy * costh[itheta]);
        // Limit unrealistic refraction speed to 1/2 pi per wave period
        ctheta[itheta] = ctheta[itheta].signum() * ctheta[itheta].abs().min(sig / 4.0);
    }

    (sinhkh, hmx, kwav, cg, ctheta)
}

// ---------------------------------------------------------------------------
// 6. numerical_limiter — Depth-limits energy/action, bounds frequency
// ---------------------------------------------------------------------------

/// Depth-limits energy and action density, bounds frequency to `[sigmin, sigmax]`.
///
/// Port of `numerical_limiter` in `src/snapwave_windsource.f90` (lines 7–65).
///
/// Returns `(H, E, A, sig)` and mutates `ee` and `aa` in place.
pub fn numerical_limiter(
    ee: &mut [f32],
    aa: &mut [f32],
    waveps: f32,
    hh: f32,
    dtheta: f32,
    rho: f32,
    g: f32,
    gamma: f32,
    sigmin: f32,
    sigmax: f32,
) -> (f32, f32, f32, f32) {
    let ntheta = ee.len();

    // ee(:) = max(ee(:), waveps)
    for v in ee.iter_mut() {
        *v = (*v).max(waveps);
    }
    // aa(:) = max(ee(:), waveps) / sigmax
    for i in 0..ntheta {
        aa[i] = ee[i].max(waveps) / sigmax;
    }

    // sigt = ee / aa  (not stored, used only for the initial aa assignment above)

    let e: f32 = ee.iter().sum::<f32>() * dtheta;
    let a: f32 = aa.iter().sum::<f32>() * dtheta;

    // depth limitation of energy and corresponding wave period
    let h = (8.0 * e / rho / g).sqrt();
    let depthlimfac = (1.0f32).max((h / (gamma * hh)).powi(2));
    let h = h.min(gamma * hh);

    let e = e / depthlimfac;
    let a = a / depthlimfac;

    for i in 0..ntheta {
        ee[i] /= depthlimfac;
        aa[i] /= depthlimfac;
    }

    // limit frequency to range
    let mut sig = e / a;
    sig = sig.max(sigmin).min(sigmax);
    for i in 0..ntheta {
        aa[i] = ee[i] / sig;
    }
    let a = aa.iter().sum::<f32>() * dtheta;
    let sig = e / a;

    (h, e, a, sig)
}

// ---------------------------------------------------------------------------
// 7. windinput — Wind growth source terms
// ---------------------------------------------------------------------------

/// Wind growth source terms based on Kahma & Calkoen (1992) and
/// Breugem & Holthuijsen (2007).
///
/// Port of `windinput` in `src/snapwave_windsource.f90` (lines 96–212).
///
/// Returns `(wsor_e, wsor_a)` vectors of length `ntheta`.
pub fn windinput(
    u10: f32,
    rho: f32,
    g: f32,
    hh: f32,
    ntheta: usize,
    windspreadfac: &[f32],
    e: f32,
    a: f32,
    cg: f32,
    eeprev: &[f32],
    aaprev: &[f32],
    ds: &[f32],
    jadcgdx: i32,
) -> (Vec<f32>, Vec<f32>) {
    // Fully developed dimensionless wave energy (Pierson Moskowitz 1964)
    let eful: f32 = 0.0036;
    // Fully developed dimensionless peak period (Pierson Moskowitz 1964)
    let tful: f32 = 7.69;
    // Shape parameters (Kahma Calkoen 1992)
    let aa1: f32 = 0.00288;
    let bb1: f32 = 0.45;
    let aa2: f32 = 0.459;
    let bb2: f32 = 0.27;
    // Shape parameters (Breugem Holthuijzen 2007)
    let aa3: f32 = 0.13;
    let bb3: f32 = 0.65;
    let aa4: f32 = 5.0;
    let bb4: f32 = 0.375;

    let pi = PI_F32;

    // compute dimensionless wave state, maximized by fully developed sea states
    let ddmlss = g * hh / (u10 * u10);
    let emaxddmlss = (aa3.powi(2) / 16.0 * ddmlss.powf(2.0 * bb3)).min(eful);
    let tmaxddmlss = (aa4 * ddmlss.powf(bb4)).min(tful);

    let cgdmlss = cg.abs() / u10;
    let edmlss = (e * g / rho / u10.powi(4)).min(emaxddmlss);
    let t = 2.0 * pi * a / e;
    let tdmlss = (g * t / u10).min(tmaxddmlss);

    // dimensionless magnitude of source terms, based on Kahma and Calkoen
    let fe = 16.0 / (2.0 * bb1 * aa1 * aa1) * (16.0 * edmlss / aa1 / aa1).powf(0.5 / bb1 - 1.0);
    let de = cgdmlss / fe;

    let ft = 1.0 / aa2 / bb2 * (tdmlss / aa2).powf(1.0 / bb2 - 1.0);
    let mut dt = cgdmlss / ft;

    let mut wsor_e = vec![0.0f32; ntheta];
    let mut wsor_a = vec![0.0f32; ntheta];

    for itheta in 0..ntheta {
        // gradT component computed from dimensional parameters
        let tprev = 2.0 * pi * aaprev[itheta] / eeprev[itheta].max(0.001);
        let tprev = tprev.max(1.0);
        let deltat = t - tprev;
        let (_kdum, cg1) = disper_nr(hh, t);
        let (_kdum, cg2) = disper_nr(hh, tprev);
        let dcgdt = if deltat.abs() > 1e-6 {
            (cg1 - cg2) / deltat
        } else {
            0.0
        };
        let dtdx = (jadcgdx as f32) * ((t - tprev) / ds[itheta]).max(0.0);
        let mut dcgdx = e * dcgdt * dtdx / u10.powi(3) / rho; // dimensionless
        dcgdx = de.min(dcgdx.abs());

        dt = dt * (0.0f32).max((1000.0 * (tful - tdmlss) / tful).tanh());

        let detmp = (de + dcgdx) * (0.0f32).max((2.0 * (eful - edmlss) / eful).tanh());
        let wsor_adlss = windspreadfac[itheta] * 0.5 / pi * (tdmlss * detmp + edmlss * dt);
        let wsor_edlss = windspreadfac[itheta] * detmp;

        // make dimensional growth rates
        wsor_e[itheta] = (u10.powi(3) * rho * wsor_edlss).max(0.0);
        wsor_a[itheta] = (u10.powi(4) * rho / g * wsor_adlss).max(0.0);
    }

    (wsor_e, wsor_a)
}

// ---------------------------------------------------------------------------
// 8. bulkdragcoeff — Bulk drag coefficient for vegetation
// ---------------------------------------------------------------------------

/// Compute bulk drag coefficient for short wave energy dissipation based on
/// the Keulegan-Carpenter number.
///
/// Port of `bulkdragcoeff` in `src/snapwave_solver.f90` (lines 1354–1424).
///
/// Uses Mendez and Losada (2004) formulation (`myflag = 2`).
pub fn bulkdragcoeff(
    ahh: f32,
    _m: i32,
    _no_nodes: i32,
    _no_secveg: i32,
    depth: f32,
    h: f32,
    kwav: f32,
    veg_bstems: f32,
    sigm: f32,
) -> f32 {
    let myflag: i32 = 2;
    let pi = PI_F32;

    // Representative wave period
    let tp = 2.0 * pi / sigm;

    // Coefficient alfa
    let alfav = if ahh >= depth { 1.0 } else { ahh / depth };

    // Representative orbital velocity
    let um = 0.5 * h * sigm * (kwav * alfav * depth).cosh() / (kwav * depth).sinh();

    // Keulegan-Carpenter number
    let kc = um * tp / veg_bstems;

    // Bulk drag coefficient
    if myflag == 1 {
        // Approach from Ozeren et al. (2013)
        if kc >= 10.0 {
            0.036 + 50.0 / kc.powf(0.926)
        } else {
            0.036 + 50.0 / 10.0f32.powf(0.926)
        }
    } else {
        // Approach from Mendez and Losada (2004), eq. 40
        let q = kc / alfav.powf(0.76);
        if q >= 7.0 {
            (-0.0138f32 * q).exp() / q.powf(0.3)
        } else {
            (-0.0138f32 * 7.0f32).exp() / 7.0f32.powf(0.3)
        }
    }
}

// ---------------------------------------------------------------------------
// 9. swvegatt — Short wave dissipation by vegetation
// ---------------------------------------------------------------------------

/// Short wave dissipation by vegetation (Suzuki et al. 2012).
///
/// Port of `swvegatt` in `src/snapwave_solver.f90` (lines 1290–1352).
pub fn swvegatt(
    sigm: f32,
    _no_nodes: i32,
    kwav: f32,
    no_secveg: usize,
    veg_ah: &[f32],
    veg_bstems: &[f32],
    veg_nstems: &[f32],
    veg_cd: &[f32],
    depth: f32,
    rho: f32,
    g: f32,
    h: f32,
) -> f32 {
    let pi = PI_F32;
    let kmr = kwav.max(0.01).min(100.0);

    let mut dvg: f32 = 0.0;
    let mut htermold: f32 = 0.0;
    let mut ahtold: f32 = 0.0;

    if no_secveg > 0 {
        for m in 0..no_secveg {
            // Determine height of vegetation section
            let aht = veg_ah[m] + ahtold;

            // restrict vegetation height to local water depth
            let aht = aht.min(depth);

            // compute hterm based on ah
            let sinh_kmr_aht = (kmr * aht).sinh();
            let cosh_kmr_depth = (kmr * depth).cosh();
            let hterm = (sinh_kmr_aht.powi(3) + 3.0 * sinh_kmr_aht)
                / (3.0 * kmr * cosh_kmr_depth.powi(3));

            // compute dissipation based on aht (Suzuki et al. 2012)
            let dvgt = 0.5 / pi.sqrt()
                * rho
                * veg_cd[m]
                * veg_bstems[m]
                * veg_nstems[m]
                * (0.5 * kmr * g / sigm).powi(3)
                * (hterm - htermold)
                * h.powi(3);

            // save hterm to htermold for next vegetation section
            htermold = hterm;
            ahtold = aht;

            // add dissipation current vegetation section
            dvg += dvgt;
        }
    }

    dvg
}

// ---------------------------------------------------------------------------
// 10. vegatt — Vegetation attenuation (top-level dispatcher)
// ---------------------------------------------------------------------------

/// Vegetation attenuation: computes bulk drag coefficient (if needed) and
/// then short wave dissipation by vegetation.
///
/// Port of `vegatt` in `src/snapwave_solver.f90` (lines 1240–1288).
pub fn vegatt(
    sigm: f32,
    no_nodes: i32,
    kwav: f32,
    no_secveg: usize,
    veg_ah: &[f32],
    veg_bstems: &[f32],
    veg_nstems: &[f32],
    veg_cd: &mut [f32],
    depth: f32,
    rho: f32,
    g: f32,
    h: f32,
) -> f32 {
    // First compute drag coefficient (if not user-defined)
    if no_secveg > 0 {
        for m in 0..no_secveg {
            if veg_cd[m] < 0.0 {
                // Cd is not user specified: compute via bulkdragcoeff
                veg_cd[m] = bulkdragcoeff(
                    veg_ah[m],
                    m as i32,
                    no_nodes,
                    no_secveg as i32,
                    depth,
                    h,
                    kwav,
                    veg_bstems[m],
                    sigm,
                );
            }
        }
    }

    // Short wave dissipation by vegetation
    swvegatt(
        sigm,
        no_nodes,
        kwav,
        no_secveg,
        veg_ah,
        veg_bstems,
        veg_nstems,
        veg_cd,
        depth,
        rho,
        g,
        h,
    )
}

// ---------------------------------------------------------------------------
// 11. compute_wave_field — Top-level solver orchestrator
// ---------------------------------------------------------------------------

/// Top-level solver orchestrator called each timestep.
///
/// Port of `compute_wave_field` in `src/snapwave_solver.f90` (lines 5–170).
///
/// This is the entry point that:
/// 1. Initializes energies on first call
/// 2. Computes celerities and refraction speed
/// 3. Calls `solve_energy_balance2Dstat`
/// 4. Computes wave forces
///
/// Returns the updated solver state arrays. The caller is responsible for
/// providing the pre-allocated mutable arrays.
#[allow(clippy::too_many_arguments)]
pub fn compute_wave_field(
    // Inputs
    time: f64,
    restart: bool,
    ig: i32,
    wind: i32,
    ja_vegetation: i32,
    ja_save_each_iter: i32,
    ntheta: usize,
    no_nodes: usize,
    no_secveg: usize,
    // Mesh / geometry (immutable inputs)
    x: &[f64],
    y: &[f64],
    dhdx: &[f32],
    dhdy: &[f32],
    msk: &[i8],
    neumannconnected: &[i32],
    theta: &[f32],
    thetamean: f32,
    // Mean boundary wave period (`tpmean_bwv`), used for the Fortran
    // `Tp = max(tpmean_bwv, Tpini)` broadcast to every node.
    tpmean_bwv: f32,
    depth: &[f32],
    fw: &[f32],
    fw_ig: &[f32],
    // Upwind geometry
    w: &[f32],    // shape (2, ntheta, no_nodes) column-major
    ds: &[f32],   // shape (ntheta, no_nodes) column-major
    prev: &[i32], // shape (2, ntheta, no_nodes) column-major
    // Solver parameters
    dt: f32,
    rho: f32,
    alfa: f32,
    gamma: f32,
    gammax: f32,
    u10: &[f32],
    niter: i32,
    crit: f32,
    upwindref: i32,
    tpini: f32,
    // Wind parameters
    windspreadfac: &[f32], // shape (ntheta, no_nodes) column-major
    jadcgdx: i32,
    sigmin: f32,
    sigmax: f32,
    c_dispt: f32,
    // IG parameters
    kwav_ig_in: &[f32],
    cg_ig_in: &[f32],
    ctheta_ig_in: &[f32], // shape (ntheta, no_nodes) column-major
    hmx_ig_in: &[f32],
    // Vegetation
    veg_ah: &[f32],     // shape (no_nodes, no_secveg) column-major
    veg_bstems: &[f32], // shape (no_nodes, no_secveg) column-major
    veg_nstems: &[f32], // shape (no_nodes, no_secveg) column-major
    veg_cd: &[f32],     // shape (no_nodes, no_secveg) column-major
    // Mutable outputs / inouts
    kwav: &mut [f32],
    cg: &mut [f32],
    ctheta: &mut [f32], // shape (ntheta, no_nodes) column-major
    ee: &mut [f32],     // shape (ntheta, no_nodes) column-major
    ee_ig: &mut [f32],  // shape (ntheta, no_nodes) column-major
    sinhkh: &mut [f32],
    hmx: &mut [f32],
    tp: &mut [f32],
    sig: &mut [f32],
    aa: &mut [f32], // shape (ntheta, no_nodes) column-major
    wsor_e: &mut [f32], // shape (ntheta, no_nodes) column-major
    wsor_a: &mut [f32], // shape (ntheta, no_nodes) column-major
    swe: &mut [f32],
    swa: &mut [f32],
    // Outputs
    h_out: &mut [f32],
    h_ig_out: &mut [f32],
    dw_out: &mut [f32],
    df_out: &mut [f32],
    f_out: &mut [f32],
    thetam_out: &mut [f32],
    dveg_out: &mut [f32],
    // Forces
    fx: &mut [f32],
    fy: &mut [f32],
    // Per-iteration map output sink (ja_save_each_iter == 1)
    iter_outputs: &mut Vec<IterOutput>,
) {
    let waveps: f32 = 1e-5;
    let g = G;

    // Allocate local IG arrays
    let mut sigm_ig = vec![0.0f32; no_nodes];
    let mut tp_ig = vec![0.0f32; no_nodes];

    if !restart {
        // Set energies to waveps for inner points
        for k in 0..no_nodes {
            if msk[k] == 1 {
                for itheta in 0..ntheta {
                    ee[itheta + k * ntheta] = waveps;
                }
            }
        }
        for v in ee_ig.iter_mut() {
            *v = waveps;
        }
    }

    // Initialize wave period
    for k in 0..no_nodes {
        if msk[k] == 1 {
            tp[k] = tpini;
        }
        if neumannconnected[k] > 0 {
            let kn = neumannconnected[k] as usize - 1; // one-based -> zero-based
            tp[kn] = tpini;
        }
    }

    // Fortran: `Tp = max(tpmean_bwv, Tpini)` broadcast to every node
    // (tpmean_bwv is the mean boundary period from update_boundary_conditions).
    for k in 0..no_nodes {
        tp[k] = tpmean_bwv.max(tpini);
    }

    for k in 0..no_nodes {
        sig[k] = 2.0 * PI_F32 / tp[k];
    }
    for k in 0..no_nodes {
        tp_ig[k] = 9.0 * tp[k];
        sigm_ig[k] = 2.0 * PI_F32 / tp_ig[k];
    }

    // kwav = sig²/g * (1 - exp(-(sig*sqrt(depth/g))^2.5))^(-0.4)
    for k in 0..no_nodes {
        let expon = -(sig[k] * (depth[k] / g).sqrt()).powf(2.5);
        kwav[k] = sig[k].powi(2) / g * (1.0 - expon.exp()).powf(-0.4);
        // C = sig / kwav (not stored globally in the port)
        let c = sig[k] / kwav[k];
        let nwav = 0.5 + kwav[k] * depth[k] / (2.0 * kwav[k] * depth[k]).min(50.0).sinh();
        cg[k] = nwav * c;
    }

    if ig == 1 {
        for k in 0..no_nodes {
            cg_ig_out(k, cg[k]); // placeholder — see below
        }
        // Actually cg_ig is an input in this port; the Fortran sets cg_ig = Cg
        // We'll handle this in the caller
    }

    // sinhkh and Hmx
    for k in 0..no_nodes {
        sinhkh[k] = (kwav[k] * depth[k]).min(50.0).sinh();
        hmx[k] = gamma * depth[k];
    }

    // ctheta
    for itheta in 0..ntheta {
        for k in 0..no_nodes {
            let idx = itheta + k * ntheta;
            ctheta[idx] = sig[k] / (2.0 * kwav[k] * depth[k]).min(50.0).sinh()
                * (dhdx[k] * theta[itheta].sin() - dhdy[k] * theta[itheta].cos());
        }
    }

    // Limit unrealistic refraction speed
    for k in 0..no_nodes {
        for itheta in 0..ntheta {
            let idx = itheta + k * ntheta;
            ctheta[idx] = ctheta[idx].signum() * ctheta[idx].abs().min(sig[k] / 4.0);
        }
    }

    // Solve the directional wave energy balance
    solve_energy_balance2Dstat(
        x, y, dhdx, dhdy, no_nodes, msk, w, ds, prev, neumannconnected, theta, ntheta,
        thetamean, depth, kwav, cg, ctheta, fw, tp, dt, rho, alfa, gamma, gammax, wind,
        h_out, dw_out, f_out, df_out, thetam_out, sinhkh, hmx, ee, windspreadfac, u10,
        niter, crit, upwindref, aa, sig, jadcgdx, sigmin, sigmax, c_dispt, wsor_e, wsor_a,
        swe, swa, tpini, ig, kwav_ig_in, cg_ig_in, h_ig_out, ctheta_ig_in, hmx_ig_in,
        ee_ig, fw_ig, time, ja_save_each_iter, ja_vegetation, no_secveg, veg_ah,
        veg_bstems, veg_nstems, veg_cd, dveg_out, iter_outputs,
    );

    // Wave forces
    for k in 0..no_nodes {
        fx[k] = f_out[k] * thetam_out[k].cos();
        fy[k] = f_out[k] * thetam_out[k].sin();
    }
}

/// Helper: write to cg_ig_out (not used in the port since cg_ig is an input).
fn cg_ig_out(_k: usize, _val: f32) {}

// ---------------------------------------------------------------------------
// 12. solve_energy_balance2Dstat — Core implicit 4-sweep solver
// ---------------------------------------------------------------------------

/// The core implicit 4-sweep solver on unstructured grids.
///
/// Port of `solve_energy_balance2Dstat` in `src/snapwave_solver.f90`
/// (lines 172–1012).
///
/// This is the heart of SnapWave: it solves the directional wave energy
/// balance equation using a 4-sweep implicit scheme with tridiagonal
/// solves per grid point.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn solve_energy_balance2Dstat(
    x: &[f64],
    y: &[f64],
    dhdx: &[f32],
    dhdy: &[f32],
    no_nodes: usize,
    msk: &[i8],
    w: &[f32],
    ds: &[f32],
    prev: &[i32],
    neumannconnected: &[i32],
    theta: &[f32],
    ntheta: usize,
    thetamean: f32,
    depth: &[f32],
    kwav: &mut [f32],
    cg: &mut [f32],
    ctheta: &mut [f32],
    fw: &[f32],
    tp: &mut [f32],
    dt: f32,
    rho: f32,
    alfa: f32,
    gamma: f32,
    gammax: f32,
    wind: i32,
    h_out: &mut [f32],
    dw_out: &mut [f32],
    f_out: &mut [f32],
    df_out: &mut [f32],
    thetam_out: &mut [f32],
    sinhkh: &mut [f32],
    hmx: &mut [f32],
    ee: &mut [f32],
    windspreadfac: &[f32],
    u10: &[f32],
    niter: i32,
    crit: f32,
    upwindref: i32,
    aa: &mut [f32],
    sig: &mut [f32],
    jadcgdx: i32,
    sigmin: f32,
    sigmax: f32,
    c_dispt: f32,
    wsor_e: &mut [f32],
    wsor_a: &mut [f32],
    swe: &mut [f32],
    swa: &mut [f32],
    tpini: f32,
    ig: i32,
    // kwav_ig is not read in the port (Hmx_ig is passed in by the caller);
    // keep the parameter for call-site symmetry with the Fortran signature.
    _kwav_ig: &[f32],
    cg_ig: &[f32],
    h_ig_out: &mut [f32],
    ctheta_ig: &[f32],
    hmx_ig: &[f32],
    ee_ig: &mut [f32],
    fw_ig: &[f32],
    time: f64,
    ja_save_each_iter: i32,
    ja_vegetation: i32,
    no_secveg: usize,
    veg_ah: &[f32],
    veg_bstems: &[f32],
    veg_nstems: &[f32],
    veg_cd: &[f32],
    dveg_out: &mut [f32],
    iter_outputs: &mut Vec<IterOutput>,
) {
    // Note: the Fortran `solve_energy_balance2Dstat` declares a *local*
    // `g = 9.81` (not the module `g = 9.813`), and a local `pi`.
    let g = G_SOLVER;
    let pi = PI_F32;
    let hmin: f32 = 0.1;
    let fac: f32 = 1.0; // underrelaxation factor for DoverA
    let epsdist: f64 = 0.001;
    let waveps: f32 = 0.0001;

    // Allocate local arrays
    let mut ok = vec![0i32; no_nodes];
    let mut indx = vec![0i32; no_nodes * 4]; // (no_nodes, 4) column-major
    let mut eeold = vec![0.0f32; ntheta * no_nodes];
    let mut dee = vec![0.0f32; ntheta];
    let mut eeprev = vec![0.0f32; ntheta];
    let mut cgprev = vec![0.0f32; ntheta];
    let mut a_arr = vec![0.0f32; ntheta];
    let mut b_arr = vec![0.0f32; ntheta];
    let mut c_arr = vec![0.0f32; ntheta];
    let mut r_arr = vec![0.0f32; ntheta];
    let mut dover_e = vec![0.0f32; no_nodes];
    let mut e_arr = vec![waveps; no_nodes];
    let mut eold = vec![0.0f32; no_nodes];

    // IG arrays
    let mut a_ig: Vec<f32> = Vec::new();
    let mut b_ig: Vec<f32> = Vec::new();
    let mut c_ig: Vec<f32> = Vec::new();
    let mut r_ig: Vec<f32> = Vec::new();
    let mut eeprev_ig: Vec<f32> = Vec::new();
    let mut cgprev_ig: Vec<f32> = Vec::new();
    let mut dover_e_ig: Vec<f32> = Vec::new();
    let mut e_ig_arr: Vec<f32> = Vec::new();
    let mut t_ig: Vec<f32> = Vec::new();
    let mut sigm_ig: Vec<f32> = Vec::new();

    if ig == 1 {
        a_ig = vec![0.0f32; ntheta];
        b_ig = vec![0.0f32; ntheta];
        c_ig = vec![0.0f32; ntheta];
        r_ig = vec![0.0f32; ntheta];
        eeprev_ig = vec![0.0f32; ntheta];
        cgprev_ig = vec![0.0f32; ntheta];
        dover_e_ig = vec![0.0f32; no_nodes];
        e_ig_arr = vec![waveps; no_nodes];
        t_ig = vec![0.0f32; no_nodes];
        sigm_ig = vec![0.0f32; no_nodes];
    }

    // Wind arrays
    let mut b_aa: Vec<f32> = Vec::new();
    let mut r_aa: Vec<f32> = Vec::new();
    let mut dover_a: Vec<f32> = Vec::new();
    let mut aaprev: Vec<f32> = Vec::new();

    if wind == 1 {
        b_aa = vec![0.0f32; ntheta];
        r_aa = vec![0.0f32; ntheta];
        dover_a = vec![0.0f32; no_nodes];
        aaprev = vec![0.0f32; ntheta];
    }

    let mut diff = vec![0.0f32; no_nodes];
    let mut ra = vec![0.0f32; no_nodes];
    let mut srcsh = vec![0.0f32; ntheta * no_nodes];

    // Precompute sin/cos of theta
    let mut sinth = vec![0.0f32; ntheta];
    let mut costh = vec![0.0f32; ntheta];
    for itheta in 0..ntheta {
        sinth[itheta] = theta[itheta].sin();
        costh[itheta] = theta[itheta].cos();
    }

    // Initialize outputs
    for v in df_out.iter_mut() { *v = 0.0; }
    for v in dw_out.iter_mut() { *v = 0.0; }
    for v in f_out.iter_mut() { *v = 0.0; }

    let mut eemax = ee.iter().cloned().fold(0.0f32, f32::max);
    let mut dtheta = theta[1] - theta[0];
    if dtheta < 0.0 {
        dtheta += 2.0 * pi;
    }

    if wind == 1 {
        for k in 0..no_nodes {
            sig[k] = 2.0 * pi / tpini;
        }
    } else {
        for k in 0..no_nodes {
            sig[k] = 2.0 * pi / tp[k];
        }
    }

    let oneoverdt = 1.0 / dt;
    let oneoverdtheta = 1.0 / dtheta;
    let oneover2dtheta = 1.0 / 2.0 / dtheta;
    let rhog8 = 0.125 * rho * g;

    for v in thetam_out.iter_mut() { *v = 0.0; }
    for v in h_out.iter_mut() { *v = 0.0; }
    for v in dveg_out.iter_mut() { *v = 0.0; }

    let (_finc2ig, shinc2ig) = if ig == 1 {
        for k in 0..no_nodes {
            t_ig[k] = 7.0 * tp[k];
            sigm_ig[k] = 2.0 * pi / t_ig[k];
        }
        (0.20f32, 0.8f32)
    } else {
        (0.0, 0.0)
    };

    let (ndissip, mut ak_val) = if wind == 1 {
        (3.0f32, waveps / sigmax)
    } else {
        (0.0, 0.0)
    };

    // ---- Sort coordinates in sweep directions ----
    // Fortran `ra(k) = x(k)*cos(thetamean + 0.5*pi*shift) + y(k)*sin(...)`:
    // x,y are real*8, the angle/cos/sin are real*4, so the product is
    // real*8 and only the final `ra` is truncated to real*4. Truncating
    // x,y to f32 first (as an earlier port did) shifts the sort key by a
    // few ulp, which flips the order of near-tied points and changes the
    // per-iteration intermediate states saved by `ja_save_each_iter`.
    let shift = [0i32, 1, -1, 2];
    for sweep in 0..4usize {
        let sweep_idx = sweep;
        let angle = thetamean + 0.5 * pi * shift[sweep_idx] as f32;
        let cosang = angle.cos();
        let sinang = angle.sin();
        for k in 0..no_nodes {
            ra[k] = (x[k] * cosang as f64 + y[k] * sinang as f64) as f32;
        }
        let col_start = sweep_idx * no_nodes;
        hpsort_eps_epw(no_nodes, &mut ra, &mut indx[col_start..col_start + no_nodes], 1.0e-6);
    }

    // ---- Boundary condition at sea side + inner-point initialization ----
    for k in 0..no_nodes {
        if msk[k] == 2 {
            // Boundary node
            for itheta in 0..ntheta {
                let idx = itheta + k * ntheta;
                ee[idx] = ee[idx].max(waveps);
            }
            e_arr[k] = ee.iter().skip(k * ntheta).take(ntheta).sum::<f32>() * dtheta;
            h_out[k] = (8.0 * e_arr[k] / rho / g).sqrt();
            let sum_sin: f32 = (0..ntheta).map(|i| ee[i + k * ntheta] * sinth[i]).sum();
            let sum_cos: f32 = (0..ntheta).map(|i| ee[i + k * ntheta] * costh[i]).sum();
            thetam_out[k] = sum_sin.atan2(sum_cos);

            if ig == 1 {
                for itheta in 0..ntheta {
                    ee_ig[itheta + k * ntheta] = 0.01 * ee[itheta + k * ntheta];
                }
                e_ig_arr[k] = ee_ig.iter().skip(k * ntheta).take(ntheta).sum::<f32>() * dtheta;
                h_ig_out[k] = (8.0 * e_ig_arr[k] / rho / g).sqrt();
            }

            if wind == 1 {
                sig[k] = 2.0 * pi / tp[k];
                for itheta in 0..ntheta {
                    let idx = itheta + k * ntheta;
                    aa[idx] = ee[idx].max(waveps) / sig[k];
                }
                ak_val = e_arr[k] / sig[k];
            }
        }

        if msk[k] == 1 {
            // Inner point
            if ig == 1 {
                for itheta in 0..ntheta {
                    let k1 = prev[0 + itheta * 2 + k * (2 * ntheta)] as usize;
                    let k2 = prev[1 + itheta * 2 + k * (2 * ntheta)] as usize;
                    // Note: prev is 1-based; k1=0 or k2=0 means no upwind point
                    if k1 > 0 && k2 > 0 {
                        let k1z = k1 - 1;
                        let k2z = k2 - 1;
                        let beta = ((w[0 + itheta * 2 + k * (2 * ntheta)] * (depth[k1z] - depth[k])
                            + w[1 + itheta * 2 + k * (2 * ntheta)] * (depth[k2z] - depth[k]))
                            / ds[itheta + k * ntheta].max(epsdist as f32))
                            .max(0.0);
                        let betan = (beta / sigm_ig[k]) * (9.81f32 / depth[k].max(0.1)).sqrt();
                        let fbr = 1.0;
                        let fsh = fbr * (-4.0 * betan.sqrt()).exp();

                        cgprev_ig[itheta] = w[0 + itheta * 2 + k * (2 * ntheta)] * cg_ig[k1z]
                            + w[1 + itheta * 2 + k * (2 * ntheta)] * cg_ig[k2z];

                        let srcsh_val = -shinc2ig * fsh
                            * ((cg[k] - cgprev_ig[itheta]) / ds[itheta + k * ntheta].max(epsdist as f32));
                        srcsh[itheta + k * ntheta] = srcsh_val.max(0.0);
                    } else {
                        srcsh[itheta + k * ntheta] = 0.0;
                    }
                }
            } else {
                for itheta in 0..ntheta {
                    srcsh[itheta + k * ntheta] = 0.0;
                }
            }

            if wind == 1 {
                for itheta in 0..ntheta {
                    ee[itheta + k * ntheta] = waveps;
                }
                sig[k] = 2.0 * pi / tpini;
                for itheta in 0..ntheta {
                    aa[itheta + k * ntheta] = ee[itheta + k * ntheta] / sig[k];
                }
            } else {
                for itheta in 0..ntheta {
                    ee[itheta + k * ntheta] = waveps;
                }
            }

            // Make sure DoverE is filled based on previous ee
            let ek = ee.iter().skip(k * ntheta).take(ntheta).sum::<f32>() * dtheta;
            let hk = (ek / rhog8).sqrt().min(gamma * depth[k]);
            let ek = rhog8 * hk * hk;

            if wind != 1 {
                let uorbi = 0.5 * sig[k] * hk / sinhkh[k];
                let dfk = 0.28 * rho * fw[k] * uorbi.powi(3);
                let dwk = baldock(rho, g, alfa, gamma, depth[k], hk, tp[k], 1, hmx[k]);
                dover_e[k] = (dwk + dfk) / ek.max(1.0e-6);
            }

            if wind == 1 {
                let (sh, hmx_val, kwav_val, cg_val, ctheta_val) = compute_celerities(
                    depth[k], sig[k], &sinth, &costh, ntheta, gamma, dhdx[k], dhdy[k],
                );
                sinhkh[k] = sh;
                hmx[k] = hmx_val;
                kwav[k] = kwav_val;
                cg[k] = cg_val;
                for itheta in 0..ntheta {
                    ctheta[itheta + k * ntheta] = ctheta_val[itheta];
                }

                let uorbi = 0.5 * sig[k] * hk / sinhkh[k];
                let dfk = 0.28 * rho * fw[k] * uorbi.powi(3);
                let dwk = baldock(rho, g, alfa, gamma, depth[k], hk, 2.0 * pi / sig[k], 1, hmx[k]);
                dover_e[k] = (dwk + dfk) / ek.max(1.0e-6);

                // initial conditions are not equal to bc conditions
                let dwt = -c_dispt / (1.0 - ndissip) * (2.0 * pi) / sig[k].powi(2)
                    * cg[k] * kwav[k] * dover_e[k];
                let dwak = 0.5 / pi * (e_arr[k] * dwt + 2.0 * pi * ak_val * dover_e[k]);
                dover_a[k] = dwak / ak_val.max(1e-6);
            }
        }
    }

    // ---- Start iteration ----
    for iter in 1..=niter {
        let mut sweep = iter % 4;
        if sweep == 0 {
            sweep = 4;
        }
        let sweep_idx = (sweep - 1) as usize;

        if sweep == 1 {
            eeold.copy_from_slice(ee);
            for k in 0..no_nodes {
                eold[k] = eeold.iter().skip(k * ntheta).take(ntheta).sum::<f32>();
            }
        }

        // Loop over all points depending on sweep direction
        for count in 0..no_nodes {
            let k = indx[count + sweep_idx * no_nodes] as usize - 1; // one-based -> zero-based

            if msk[k] == 1 {
                if depth[k] > 1.1 * hmin {
                    if ok[k] != 1 {
                        // Only perform computations on wet inner points not yet converged

                        // Interpolate upwind values
                        for itheta in 0..ntheta {
                            let k1 = prev[0 + itheta * 2 + k * (2 * ntheta)] as usize;
                            let k2 = prev[1 + itheta * 2 + k * (2 * ntheta)] as usize;
                            let k1z = if k1 > 0 { k1 - 1 } else { 0 };
                            let k2z = if k2 > 0 { k2 - 1 } else { 0 };

                            eeprev[itheta] = w[0 + itheta * 2 + k * (2 * ntheta)] * ee[itheta + k1z * ntheta]
                                + w[1 + itheta * 2 + k * (2 * ntheta)] * ee[itheta + k2z * ntheta];
                            cgprev[itheta] = w[0 + itheta * 2 + k * (2 * ntheta)] * cg[k1z]
                                + w[1 + itheta * 2 + k * (2 * ntheta)] * cg[k2z];

                            if ig == 1 {
                                eeprev_ig[itheta] = w[0 + itheta * 2 + k * (2 * ntheta)] * ee_ig[itheta + k1z * ntheta]
                                    + w[1 + itheta * 2 + k * (2 * ntheta)] * ee_ig[itheta + k2z * ntheta];
                                cgprev_ig[itheta] = w[0 + itheta * 2 + k * (2 * ntheta)] * cg_ig[k1z]
                                    + w[1 + itheta * 2 + k * (2 * ntheta)] * cg_ig[k2z];
                            }

                            if wind == 1 {
                                aaprev[itheta] = w[0 + itheta * 2 + k * (2 * ntheta)] * aa[itheta + k1z * ntheta]
                                    + w[1 + itheta * 2 + k * (2 * ntheta)] * aa[itheta + k2z * ntheta];
                            }
                        }

                        let ek = eeprev.iter().sum::<f32>() * dtheta;

                        let depthlimfac = (1.0f32).max(((ek / rhog8).sqrt() / (gammax * depth[k])).powi(2));
                        let hk = (ek / rhog8).sqrt().min(gamma * depth[k]);
                        let ek = ek / depthlimfac;

                        let mut ak: f32 = 0.0;
                        if wind == 1 {
                            ak = aaprev.iter().sum::<f32>() * dtheta;
                            ak /= depthlimfac;
                            for itheta in 0..ntheta {
                                ee[itheta + k * ntheta] /= depthlimfac;
                                aa[itheta + k * ntheta] /= depthlimfac;
                            }
                            sig[k] = ek / ak;
                            sig[k] = sig[k].max(sigmin).min(sigmax);
                            ak = ek / sig[k]; // to avoid small T in windinput

                            if wind == 1 {
                                for itheta in 0..ntheta {
                                    aaprev[itheta] = aaprev[itheta].min(eeprev[itheta] / sigmin);
                                    aaprev[itheta] = aaprev[itheta].max(eeprev[itheta] / sigmax);
                                }
                            }

                            let (sh, hmx_val, kwav_val, cg_val, ctheta_val) = compute_celerities(
                                depth[k], sig[k], &sinth, &costh, ntheta, gamma, dhdx[k], dhdy[k],
                            );
                            sinhkh[k] = sh;
                            hmx[k] = hmx_val;
                            kwav[k] = kwav_val;
                            cg[k] = cg_val;
                            for itheta in 0..ntheta {
                                ctheta[itheta + k * ntheta] = ctheta_val[itheta];
                            }
                        }

                        // Fill DoverE
                        let uorbi = 0.5 * sig[k] * hk / sinhkh[k];
                        let dfk = 0.28 * rho * fw[k] * uorbi.powi(3);
                        let dwk = if hk > 0.0 {
                            baldock(rho, g, alfa, gamma, depth[k], hk, 2.0 * pi / sig[k], 1, hmx[k])
                        } else {
                            0.0
                        };

                        let dvegk = if ja_vegetation == 1 {
                            let veg_ah_k = &veg_ah[k * no_secveg..(k + 1) * no_secveg];
                            let veg_bstems_k = &veg_bstems[k * no_secveg..(k + 1) * no_secveg];
                            let veg_nstems_k = &veg_nstems[k * no_secveg..(k + 1) * no_secveg];
                            let veg_cd_k = &veg_cd[k * no_secveg..(k + 1) * no_secveg];
                            let mut veg_cd_mut: Vec<f32> = veg_cd_k.to_vec();
                            vegatt(
                                sig[k], no_nodes as i32, kwav[k], no_secveg,
                                veg_ah_k, veg_bstems_k, veg_nstems_k, &mut veg_cd_mut,
                                depth[k], rho, g, hk,
                            )
                        } else {
                            0.0
                        };

                        dover_e[k] = (dwk + dfk + dvegk) / ek.max(1.0e-6);

                        if wind == 1 {
                            if iter == 1 {
                                let (ws_e, ws_a) = windinput(
                                    u10[k], rho, g, depth[k], ntheta,
                                    &windspreadfac[k * ntheta..(k + 1) * ntheta],
                                    ek, ak, cg[k], &eeprev, &aaprev,
                                    &ds[k * ntheta..(k + 1) * ntheta], jadcgdx,
                                );
                                for itheta in 0..ntheta {
                                    wsor_e[itheta + k * ntheta] = ws_e[itheta];
                                    wsor_a[itheta + k * ntheta] = ws_a[itheta];
                                }
                            } else {
                                let (ws_e, ws_a) = windinput(
                                    u10[k], rho, g, depth[k], ntheta,
                                    &windspreadfac[k * ntheta..(k + 1) * ntheta],
                                    ek, ak, cg[k],
                                    &ee[k * ntheta..(k + 1) * ntheta],
                                    &aa[k * ntheta..(k + 1) * ntheta],
                                    &ds[k * ntheta..(k + 1) * ntheta], jadcgdx,
                                );
                                for itheta in 0..ntheta {
                                    wsor_e[itheta + k * ntheta] = ws_e[itheta];
                                    wsor_a[itheta + k * ntheta] = ws_a[itheta];
                                }
                            }

                            let dwt = -c_dispt / (1.0 - ndissip) * (2.0 * pi) / sig[k].powi(2)
                                * cg[k] * kwav[k] * dover_e[k];
                            let dwak = 1.0 / 2.0 / pi * (e_arr[k] * dwt + 2.0 * pi * ak_val * dover_e[k]);

                            if iter == 1 {
                                dover_a[k] = dwak / ak_val.max(1e-6);
                            } else {
                                dover_a[k] = (1.0 - fac) * dover_a[k] + fac * dwak / ak_val.max(1e-6);
                            }

                            let (sh, hmx_val, kwav_val, cg_val, ctheta_val) = compute_celerities(
                                depth[k], sig[k], &sinth, &costh, ntheta, gamma, dhdx[k], dhdy[k],
                            );
                            sinhkh[k] = sh;
                            hmx[k] = hmx_val;
                            kwav[k] = kwav_val;
                            cg[k] = cg_val;
                            for itheta in 0..ntheta {
                                ctheta[itheta + k * ntheta] = ctheta_val[itheta];
                            }
                        }

                        // Build RHS
                        for itheta in 0..ntheta {
                            r_arr[itheta] = oneoverdt * ee[itheta + k * ntheta]
                                + cgprev[itheta] * eeprev[itheta]
                                    / ds[itheta + k * ntheta].max(epsdist as f32);
                        }

                        // Build tridiagonal system
                        if upwindref == 0 {
                            // central scheme
                            for itheta in 1..ntheta - 1 {
                                a_arr[itheta] = -ctheta[itheta - 1 + k * ntheta] * oneover2dtheta;
                                b_arr[itheta] = oneoverdt
                                    + cg[k] / ds[itheta + k * ntheta].max(epsdist as f32)
                                    + dover_e[k]
                                    + srcsh[itheta + k * ntheta];
                                c_arr[itheta] = ctheta[itheta + 1 + k * ntheta] * oneover2dtheta;
                            }

                            a_arr[0] = -ctheta[ntheta - 1 + k * ntheta] * oneover2dtheta;
                            b_arr[0] = oneoverdt
                                + cg[k] / ds[0 + k * ntheta].max(epsdist as f32)
                                + dover_e[k]
                                + srcsh[0 + k * ntheta];
                            c_arr[0] = ctheta[1 + k * ntheta] * oneover2dtheta;

                            a_arr[ntheta - 1] = -ctheta[ntheta - 2 + k * ntheta] * oneover2dtheta;
                            b_arr[ntheta - 1] = oneoverdt
                                + cg[k] / ds[ntheta - 1 + k * ntheta].max(epsdist as f32)
                                + dover_e[k]
                                + srcsh[ntheta - 1 + k * ntheta];
                            c_arr[ntheta - 1] = ctheta[0 + k * ntheta] * oneover2dtheta;
                        } else {
                            // upwind scheme
                            for itheta in 1..ntheta - 1 {
                                if ctheta[itheta + k * ntheta] < 0.0 {
                                    a_arr[itheta] = 0.0;
                                    b_arr[itheta] = oneoverdt
                                        - ctheta[itheta + k * ntheta] * oneoverdtheta
                                        + cg[k] / ds[itheta + k * ntheta].max(epsdist as f32)
                                        + dover_e[k]
                                        + srcsh[itheta + k * ntheta];
                                    c_arr[itheta] = ctheta[itheta + 1 + k * ntheta] * oneoverdtheta;
                                } else {
                                    a_arr[itheta] = -ctheta[itheta - 1 + k * ntheta] / dtheta;
                                    b_arr[itheta] = oneoverdt
                                        + ctheta[itheta + k * ntheta] * oneoverdtheta
                                        + cg[k] / ds[itheta + k * ntheta].max(epsdist as f32)
                                        + dover_e[k]
                                        + srcsh[itheta + k * ntheta];
                                    c_arr[itheta] = 0.0;
                                }
                            }

                            if ctheta[0 + k * ntheta] < 0.0 {
                                a_arr[0] = 0.0;
                                b_arr[0] = oneoverdt
                                    - ctheta[0 + k * ntheta] * oneoverdtheta
                                    + cg[k] / ds[0 + k * ntheta].max(epsdist as f32)
                                    + dover_e[k]
                                    + srcsh[0 + k * ntheta];
                                c_arr[0] = ctheta[1 + k * ntheta] * oneoverdtheta;
                            } else {
                                a_arr[0] = 0.0;
                                b_arr[0] = oneoverdt;
                                c_arr[0] = 0.0;
                                r_arr[0] = 0.0;
                            }

                            if ctheta[ntheta - 1 + k * ntheta] > 0.0 {
                                a_arr[ntheta - 1] = -ctheta[ntheta - 2 + k * ntheta] / dtheta;
                                b_arr[ntheta - 1] = oneoverdt
                                    + ctheta[ntheta - 1 + k * ntheta] * oneoverdtheta
                                    + cg[k] / ds[ntheta - 1 + k * ntheta].max(epsdist as f32)
                                    + dover_e[k]
                                    + srcsh[ntheta - 1 + k * ntheta];
                                c_arr[ntheta - 1] = 0.0;
                            } else {
                                a_arr[ntheta - 1] = 0.0;
                                b_arr[ntheta - 1] = oneoverdt;
                                c_arr[ntheta - 1] = 0.0;
                                r_arr[ntheta - 1] = 0.0;
                            }
                        }

                        // Solve tridiagonal system per point
                        if wind == 1 {
                            for itheta in 0..ntheta {
                                r_aa[itheta] = oneoverdt * aa[itheta + k * ntheta]
                                    + cgprev[itheta] * aaprev[itheta]
                                        / ds[itheta + k * ntheta].max(epsdist as f32);
                            }

                            if upwindref == 0 {
                                for itheta in 0..ntheta {
                                    b_aa[itheta] = oneoverdt
                                        + cg[k] / ds[itheta + k * ntheta].max(epsdist as f32)
                                        + dover_a[k];
                                }
                            } else {
                                for itheta in 1..ntheta - 1 {
                                    if ctheta[itheta + k * ntheta] < 0.0 {
                                        b_aa[itheta] = oneoverdt
                                            - ctheta[itheta + k * ntheta] * oneoverdtheta
                                            + cg[k] / ds[itheta + k * ntheta].max(epsdist as f32)
                                            + dover_a[k];
                                    } else {
                                        b_aa[itheta] = oneoverdt
                                            + ctheta[itheta + k * ntheta] * oneoverdtheta
                                            + cg[k] / ds[itheta + k * ntheta].max(epsdist as f32)
                                            + dover_a[k];
                                    }
                                }

                                if ctheta[0 + k * ntheta] < 0.0 {
                                    b_aa[0] = oneoverdt
                                        - ctheta[0 + k * ntheta] / dtheta
                                        + cg[k] / ds[0 + k * ntheta].max(epsdist as f32)
                                        + dover_a[k];
                                } else {
                                    b_aa[0] = oneoverdt;
                                    r_aa[0] = 0.0;
                                }

                                if ctheta[ntheta - 1 + k * ntheta] > 0.0 {
                                    b_aa[ntheta - 1] = oneoverdt
                                        + ctheta[ntheta - 1 + k * ntheta] / dtheta
                                        + cg[k] / ds[ntheta - 1 + k * ntheta].max(epsdist as f32)
                                        + dover_a[k];
                                } else {
                                    b_aa[ntheta - 1] = oneoverdt;
                                    r_aa[ntheta - 1] = 0.0;
                                }
                            }

                            // Add wind source terms to RHS
                            for itheta in 0..ntheta {
                                r_arr[itheta] += wsor_e[itheta + k * ntheta];
                                r_aa[itheta] += wsor_a[itheta + k * ntheta];
                            }

                            let ee_sol = solve_tridiag(&a_arr, &b_arr, &c_arr, &r_arr, ntheta);
                            let aa_sol = solve_tridiag(&a_arr, &b_aa, &c_arr, &r_aa, ntheta);

                            for itheta in 0..ntheta {
                                ee[itheta + k * ntheta] = ee_sol[itheta].max(waveps);
                                aa[itheta + k * ntheta] = aa_sol[itheta].max(waveps / sigmax);
                                aa[itheta + k * ntheta] = aa[itheta + k * ntheta].max(waveps / sig[k]);
                            }

                            let ek = ee.iter().skip(k * ntheta).take(ntheta).sum::<f32>() * dtheta;
                            let ak = aa.iter().skip(k * ntheta).take(ntheta).sum::<f32>() * dtheta;

                            let depthlimfac = (1.0f32).max(((ek / rhog8).sqrt() / (gammax * depth[k])).powi(2));
                            let _hk = (ek / rhog8 / depthlimfac).sqrt();

                            let ek = ek / depthlimfac;
                            let ak = ak / depthlimfac;
                            for itheta in 0..ntheta {
                                ee[itheta + k * ntheta] /= depthlimfac;
                                aa[itheta + k * ntheta] /= depthlimfac;
                            }

                            sig[k] = ek / ak;
                            sig[k] = sig[k].max(sigmin).min(sigmax);

                            let (sh, hmx_val, kwav_val, cg_val, ctheta_val) = compute_celerities(
                                depth[k], sig[k], &sinth, &costh, ntheta, gamma, dhdx[k], dhdy[k],
                            );
                            sinhkh[k] = sh;
                            hmx[k] = hmx_val;
                            kwav[k] = kwav_val;
                            cg[k] = cg_val;
                            for itheta in 0..ntheta {
                                ctheta[itheta + k * ntheta] = ctheta_val[itheta];
                            }
                        } else {
                            // No wind: solve tridiagonal for ee only
                            let ee_sol = solve_tridiag(&a_arr, &b_arr, &c_arr, &r_arr, ntheta);
                            for itheta in 0..ntheta {
                                ee[itheta + k * ntheta] = ee_sol[itheta].max(waveps);
                            }
                        }

                        // IG
                        if ig == 1 {
                            let ek_ig = eeprev_ig.iter().sum::<f32>() * dtheta;
                            let hk_ig = (ek_ig / rhog8).sqrt();
                            let ek_ig = rhog8 * hk_ig * hk_ig;

                            // Bottom friction Henderson and Bowen (2002)
                            let dfk_ig = fw_ig[k] * 0.0361 * (9.81 / depth[k]).powf(1.5) * hk * ek_ig;

                            // Dissipation of infragravity waves
                            let dwk_ig = if hk_ig > 0.2 * hmx_ig[k] {
                                baldock(rho, g, alfa, gamma, depth[k], hk_ig, t_ig[k], 1, hmx_ig[k])
                            } else {
                                0.0
                            };

                            dover_e_ig[k] = (dwk_ig + dfk_ig) / ek_ig.max(1.0e-6);

                            for itheta in 0..ntheta {
                                r_ig[itheta] = oneoverdt * ee_ig[itheta + k * ntheta]
                                    + cgprev_ig[itheta] * eeprev_ig[itheta]
                                        / ds[itheta + k * ntheta].max(epsdist as f32)
                                    + srcsh[itheta + k * ntheta] * ee[itheta + k * ntheta];
                            }

                            // IG always uses central scheme for interior, upwind for boundaries
                            for itheta in 1..ntheta - 1 {
                                a_ig[itheta] = -ctheta_ig[itheta - 1 + k * ntheta] * oneover2dtheta;
                                b_ig[itheta] = oneoverdt
                                    + cg_ig[k] / ds[itheta + k * ntheta].max(epsdist as f32)
                                    + dover_e_ig[k];
                                c_ig[itheta] = ctheta_ig[itheta + 1 + k * ntheta] * oneover2dtheta;
                            }

                            if ctheta_ig[0 + k * ntheta] < 0.0 {
                                a_ig[0] = 0.0;
                                b_ig[0] = oneoverdt
                                    - ctheta_ig[0 + k * ntheta] / dtheta
                                    + cg_ig[k] / ds[0 + k * ntheta].max(epsdist as f32)
                                    + dover_e_ig[k];
                                c_ig[0] = ctheta_ig[1 + k * ntheta] / dtheta;
                            } else {
                                a_ig[0] = 0.0;
                                b_ig[0] = 1.0 / dt
                                    + cg_ig[k] / ds[0 + k * ntheta].max(epsdist as f32)
                                    + dover_e_ig[k];
                                c_ig[0] = 0.0;
                            }

                            if ctheta_ig[ntheta - 1 + k * ntheta] > 0.0 {
                                a_ig[ntheta - 1] = -ctheta_ig[ntheta - 2 + k * ntheta] / dtheta;
                                b_ig[ntheta - 1] = oneoverdt
                                    + ctheta_ig[ntheta - 1 + k * ntheta] / dtheta
                                    + cg_ig[k] / ds[ntheta - 1 + k * ntheta].max(epsdist as f32)
                                    + dover_e_ig[k];
                                c_ig[ntheta - 1] = 0.0;
                            } else {
                                a_ig[ntheta - 1] = 0.0;
                                b_ig[ntheta - 1] = oneoverdt
                                    + cg_ig[k] / ds[ntheta - 1 + k * ntheta].max(epsdist as f32)
                                    + dover_e_ig[k];
                                c_ig[ntheta - 1] = 0.0;
                            }

                            let ee_ig_sol = solve_tridiag(&a_ig, &b_ig, &c_ig, &r_ig, ntheta);
                            for itheta in 0..ntheta {
                                ee_ig[itheta + k * ntheta] = ee_ig_sol[itheta].max(0.0);
                            }
                        } else {
                            for itheta in 0..ntheta {
                                ee_ig[itheta + k * ntheta] = 0.0;
                            }
                        }
                    }
                } else {
                    // depth <= 1.1 * hmin: dry point
                    for itheta in 0..ntheta {
                        ee[itheta + k * ntheta] = 0.0;
                    }
                    if wind == 1 {
                        for itheta in 0..ntheta {
                            aa[itheta + k * ntheta] = 0.0;
                        }
                    }
                    for itheta in 0..ntheta {
                        ee_ig[itheta + k * ntheta] = 0.0;
                    }
                }
            }

            // Neumann boundary: copy values from connected inner point
            if neumannconnected[k] != 0 {
                let kn = neumannconnected[k] as usize - 1;
                sinhkh[kn] = sinhkh[k];
                kwav[kn] = kwav[k];
                hmx[kn] = hmx[k];
                for itheta in 0..ntheta {
                    ee[itheta + kn * ntheta] = ee[itheta + k * ntheta];
                    ctheta[itheta + kn * ntheta] = ctheta[itheta + k * ntheta];
                }
                cg[kn] = cg[k];
                if wind == 1 {
                    sig[kn] = sig[k];
                    tp[kn] = 2.0 * pi / sig[k];
                    for itheta in 0..ntheta {
                        wsor_e[itheta + kn * ntheta] = wsor_e[itheta + k * ntheta];
                        wsor_a[itheta + kn * ntheta] = wsor_a[itheta + k * ntheta];
                        aa[itheta + kn * ntheta] = aa[itheta + k * ntheta];
                    }
                }
                df_out[kn] = df_out[k];
                dw_out[kn] = dw_out[k];
            }
        }

        // Check convergence after all 4 sweeps
        if sweep == 4 {
            for k in 0..no_nodes {
                for itheta in 0..ntheta {
                    dee[itheta] = ee[itheta + k * ntheta] - eeold[itheta + k * ntheta];
                }
                diff[k] = dee.iter().map(|v| v.abs()).fold(0.0f32, f32::max);

                if diff[k] / eemax < crit {
                    ok[k] = 1;
                }
            }

            let _percok = ok.iter().sum::<i32>() as f64 / no_nodes as f64 * 100.0;
            eemax = ee.iter().cloned().fold(0.0f32, f32::max);

            let error = diff.iter().cloned().fold(0.0f32, f32::max) / eemax;

            if error < crit {
                break;
            }
        }

        // ---- Per-iteration map output (ja_save_each_iter == 1) ----
        // Mirrors the Fortran block at the end of each sweep iteration:
        // recompute H/thetam/Df/Dw/F/H_ig/Tp/SwE/SwA from the current ee,
        // then emit a map record at time+iter with index iter.
        if ja_save_each_iter == 1 {
            let mut h_snap = vec![0.0f32; no_nodes];
            let mut thetam_snap = vec![0.0f32; no_nodes];
            let mut df_snap = vec![0.0f32; no_nodes];
            let mut dw_snap = vec![0.0f32; no_nodes];
            let mut f_snap = vec![0.0f32; no_nodes];
            let mut h_ig_snap = vec![0.0f32; no_nodes];
            let mut tp_snap = tp.to_vec();
            let mut swe_snap = vec![0.0f32; no_nodes];
            let mut swa_snap = vec![0.0f32; no_nodes];

            for k in 0..no_nodes {
                let ek = ee.iter().skip(k * ntheta).take(ntheta).sum::<f32>() * dtheta;
                let _ak = aa.iter().skip(k * ntheta).take(ntheta).sum::<f32>() * dtheta;
                h_snap[k] = (8.0 * ek / rho / g).sqrt();
                let sum_sin: f32 = (0..ntheta).map(|i| ee[i + k * ntheta] * sinth[i]).sum();
                let sum_cos: f32 = (0..ntheta).map(|i| ee[i + k * ntheta] * costh[i]).sum();
                thetam_snap[k] = sum_sin.atan2(sum_cos);
                let uorbi = 0.5 * sig[k] * h_snap[k] / sinhkh[k];
                df_snap[k] = 0.28 * rho * fw[k] * uorbi.powi(3);
                dw_snap[k] = baldock(rho, g, alfa, gamma, depth[k], h_snap[k], tp_snap[k], 1, hmx[k]);
                f_snap[k] = dw_snap[k] * kwav[k] / sig[k] / rho / depth[k];

                if ig == 1 {
                    let ek_ig = ee_ig.iter().skip(k * ntheta).take(ntheta).sum::<f32>() * dtheta;
                    h_ig_snap[k] = (8.0 * ek_ig / rho / g).sqrt();
                }
                if wind == 1 {
                    tp_snap[k] = 2.0 * pi / sig[k];
                    swe_snap[k] = wsor_e.iter().skip(k * ntheta).take(ntheta).sum::<f32>() * dtheta;
                    swa_snap[k] = wsor_a.iter().skip(k * ntheta).take(ntheta).sum::<f32>() * dtheta;
                }
            }

            iter_outputs.push(IterOutput {
                ntmapout: iter,
                time: time + iter as f64,
                snapshot: IterSnapshot {
                    h: h_snap,
                    thetam: thetam_snap,
                    df: df_snap,
                    dw: dw_snap,
                    f: f_snap,
                    h_ig: h_ig_snap,
                    tp: tp_snap,
                    sig: sig.to_vec(),
                    swe: swe_snap,
                    swa: swa_snap,
                    ee: ee.to_vec(),
                    ctheta: ctheta.to_vec(),
                    cg: cg.to_vec(),
                },
            });
        }
    }

    // ---- Final output computation ----
    for k in 0..no_nodes {
        if depth[k] > 1.1 * hmin {
            e_arr[k] = ee.iter().skip(k * ntheta).take(ntheta).sum::<f32>() * dtheta;
            h_out[k] = (8.0 * e_arr[k] / rho / g).sqrt();
            let sum_sin: f32 = (0..ntheta).map(|i| ee[i + k * ntheta] * sinth[i]).sum();
            let sum_cos: f32 = (0..ntheta).map(|i| ee[i + k * ntheta] * costh[i]).sum();
            thetam_out[k] = sum_sin.atan2(sum_cos);

            let uorbi = if wind == 1 {
                0.5 * sig[k] * h_out[k] / sinhkh[k] // Note: uses hk from last iteration, approximate
            } else {
                PI_F32 * h_out[k] / tp[k] / sinhkh[k]
            };
            df_out[k] = 0.28 * rho * fw[k] * uorbi.powi(3);

            if wind == 1 {
                dw_out[k] = baldock(rho, g, alfa, gamma, depth[k], h_out[k], 2.0 * pi / sig[k], 1, hmx[k]);
                f_out[k] = dw_out[k] * kwav[k] / sig[k] / rho / depth[k];
            } else {
                dw_out[k] = baldock(rho, g, alfa, gamma, depth[k], h_out[k], tp[k], 1, hmx[k]);
                f_out[k] = dw_out[k] * kwav[k] / sig[k] / rho / depth[k];
            }

            if ja_vegetation == 1 {
                let veg_ah_k = &veg_ah[k * no_secveg..(k + 1) * no_secveg];
                let veg_bstems_k = &veg_bstems[k * no_secveg..(k + 1) * no_secveg];
                let veg_nstems_k = &veg_nstems[k * no_secveg..(k + 1) * no_secveg];
                let veg_cd_k = &veg_cd[k * no_secveg..(k + 1) * no_secveg];
                let mut veg_cd_mut: Vec<f32> = veg_cd_k.to_vec();
                dveg_out[k] = vegatt(
                    sig[k], no_nodes as i32, kwav[k], no_secveg,
                    veg_ah_k, veg_bstems_k, veg_nstems_k, &mut veg_cd_mut,
                    depth[k], rho, g, h_out[k],
                );
            } else {
                dveg_out[k] = 0.0;
            }

            if ig == 1 {
                e_ig_arr[k] = ee_ig.iter().skip(k * ntheta).take(ntheta).sum::<f32>() * dtheta;
                h_ig_out[k] = (8.0 * e_ig_arr[k] / rho / g).sqrt();
            }

            if wind == 1 {
                tp[k] = 2.0 * pi / sig[k];
                swe[k] = wsor_e.iter().skip(k * ntheta).take(ntheta).sum::<f32>() * dtheta;
                swa[k] = wsor_a.iter().skip(k * ntheta).take(ntheta).sum::<f32>() * dtheta;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_solve_tridiag_simple() {
        // Solve a simple 3x3 tridiagonal system
        // 2x₁ + x₂         = 4
        // x₁  + 2x₂ + x₃   = 8
        //        x₂ + 2x₃  = 6
        // Solution: x = [0.5, 3.0, 1.5]
        let a = [0.0f32, 1.0, 1.0];
        let b = [2.0f32, 2.0, 2.0];
        let c = [1.0f32, 1.0, 0.0];
        let d = [4.0f32, 8.0, 6.0];
        let x = solve_tridiag(&a, &b, &c, &d, 3);
        assert!((x[0] - 0.5).abs() < 1e-5, "x[0] = {}", x[0]);
        assert!((x[1] - 3.0).abs() < 1e-5, "x[1] = {}", x[1]);
        assert!((x[2] - 1.5).abs() < 1e-5, "x[2] = {}", x[2]);
    }

    #[test]
    fn test_solve_tridiag_identity() {
        // Identity matrix: x = d
        let a = [0.0f32, 0.0, 0.0];
        let b = [1.0f32, 1.0, 1.0];
        let c = [0.0f32, 0.0, 0.0];
        let d = [5.0f32, 7.0, 9.0];
        let x = solve_tridiag(&a, &b, &c, &d, 3);
        assert!((x[0] - 5.0).abs() < 1e-5);
        assert!((x[1] - 7.0).abs() < 1e-5);
        assert!((x[2] - 9.0).abs() < 1e-5);
    }

    #[test]
    fn test_baldock_opt1() {
        let dw = baldock(1025.0, 9.813, 1.0, 0.6, 10.0, 2.0, 8.0, 1, 6.0);
        // Just verify it's positive and finite
        assert!(dw > 0.0);
        assert!(dw.is_finite());
    }

    #[test]
    fn test_baldock_opt2() {
        let dw = baldock(1025.0, 9.813, 1.0, 0.6, 10.0, 2.0, 8.0, 2, 6.0);
        assert!(dw > 0.0);
        assert!(dw.is_finite());
    }

    #[test]
    fn test_baldock_zero_height() {
        // H=0 should use Hloc=1e-6
        let dw = baldock(1025.0, 9.813, 1.0, 0.6, 10.0, 0.0, 8.0, 1, 6.0);
        assert!(dw.is_finite());
    }

    #[test]
    fn test_hpsort_eps_epw_basic() {
        let mut ra = vec![3.0f32, 1.0, 2.0, 5.0, 4.0];
        let mut ind = vec![0i32; 5];
        hpsort_eps_epw(5, &mut ra, &mut ind, 1e-6);
        assert_eq!(ra, vec![1.0, 2.0, 3.0, 4.0, 5.0]);
        // Check that indices track original positions (1-based)
        let sorted_indices: Vec<i32> = ind.iter().copied().collect();
        // The sorted values came from original positions: 1.0@2, 2.0@3, 3.0@1, 4.0@5, 5.0@4
        assert_eq!(sorted_indices, vec![2, 3, 1, 5, 4]);
    }

    #[test]
    fn test_hpsort_eps_epw_with_duplicates() {
        let mut ra = vec![1.0f32, 1.0, 2.0];
        let mut ind = vec![0i32; 3];
        hpsort_eps_epw(3, &mut ra, &mut ind, 1e-6);
        assert_eq!(ra, vec![1.0, 1.0, 2.0]);
    }

    #[test]
    fn test_hpsort_eps_epw_single_element() {
        let mut ra = vec![42.0f32];
        let mut ind = vec![0i32; 1];
        hpsort_eps_epw(1, &mut ra, &mut ind, 1e-6);
        assert_eq!(ra, vec![42.0]);
        assert_eq!(ind, vec![1]);
    }

    #[test]
    fn test_disper_nr_deep_water() {
        // Deep water: h=100m, T=10s -> k ≈ ω²/g = (2π/10)²/9.813 ≈ 0.040
        let (k, cg) = disper_nr(100.0, 10.0);
        let omega = 2.0 * PI_F32 / 10.0;
        let k_deep = omega * omega / G;
        assert!((k - k_deep).abs() < 0.01 * k_deep, "k={}, k_deep={}", k, k_deep);
        // In deep water, cg ≈ C/2
        let c = omega / k;
        assert!((cg - c * 0.5).abs() < 0.1 * c, "cg={}, c/2={}", cg, c * 0.5);
    }

    #[test]
    fn test_disper_nr_shallow_water() {
        // Shallow water: h=1m, T=10s
        let (k, cg) = disper_nr(1.0, 10.0);
        assert!(k > 0.0);
        assert!(cg > 0.0);
        // In shallow water, cg ≈ C ≈ sqrt(gh)
        let c_shallow = (G * 1.0).sqrt();
        assert!((cg - c_shallow).abs() < 0.5, "cg={}, c_shallow={}", cg, c_shallow);
    }

    #[test]
    fn test_compute_celerities() {
        let sinth = [0.0f32, 1.0, 0.0, -1.0];
        let costh = [1.0f32, 0.0, -1.0, 0.0];
        let (sinhkh, hmx, kwav, cg, ctheta) = compute_celerities(
            10.0, 0.5, &sinth, &costh, 4, 0.6, 0.01, 0.0,
        );
        assert!(sinhkh > 0.0);
        assert!(hmx > 0.0);
        assert!(kwav > 0.0);
        assert!(cg > 0.0);
        assert_eq!(ctheta.len(), 4);
        // Refraction speed should be bounded by sig/4
        for &ct in &ctheta {
            assert!(ct.abs() <= 0.5 / 4.0 + 1e-5);
        }
    }

    #[test]
    fn test_numerical_limiter() {
        let mut ee = vec![10.0f32, 20.0, 15.0];
        let mut aa = vec![1.0f32, 2.0, 1.5];
        let (h, e, a, sig) = numerical_limiter(
            &mut ee, &mut aa, 1e-5, 5.0, 0.1745, 1025.0, 9.813, 0.6, 0.2, 2.0,
        );
        assert!(h > 0.0);
        assert!(e > 0.0);
        assert!(a > 0.0);
        assert!(sig >= 0.2 && sig <= 2.0);
    }

    #[test]
    fn test_windinput_no_wind() {
        // u10=0 should produce zero source terms
        let windspreadfac = vec![0.25f32; 4];
        let eeprev = vec![1.0f32; 4];
        let aaprev = vec![0.5f32; 4];
        let ds = vec![100.0f32; 4];
        let (ws_e, ws_a) = windinput(
            0.0, 1025.0, 9.813, 10.0, 4, &windspreadfac,
            10.0, 5.0, 5.0, &eeprev, &aaprev, &ds, 0,
        );
        // With u10=0, ddmlss is infinite, so source terms should be zero
        for &v in &ws_e {
            assert!(v.is_finite());
        }
        for &v in &ws_a {
            assert!(v.is_finite());
        }
    }

    #[test]
    fn test_bulkdragcoeff_mendez_losada() {
        let cd = bulkdragcoeff(1.0, 0, 100, 3, 5.0, 1.0, 0.5, 0.01, 0.6);
        assert!(cd > 0.0);
        assert!(cd.is_finite());
    }

    #[test]
    fn test_swvegatt_no_vegetation() {
        let dveg = swvegatt(0.5, 100, 0.5, 0, &[], &[], &[], &[], 5.0, 1025.0, 9.813, 1.0);
        assert_eq!(dveg, 0.0);
    }

    #[test]
    fn test_swvegatt_with_vegetation() {
        let veg_ah = [1.0f32];
        let veg_bstems = [0.01f32];
        let veg_nstems = [100.0f32];
        let veg_cd = [1.0f32];
        let dveg = swvegatt(0.5, 100, 0.5, 1, &veg_ah, &veg_bstems, &veg_nstems, &veg_cd, 5.0, 1025.0, 9.813, 1.0);
        assert!(dveg >= 0.0);
        assert!(dveg.is_finite());
    }
}