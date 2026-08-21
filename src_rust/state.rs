//! Rust-owned time loop and output scheduling (plan.md Phase 10), now the
//! single scheduling authority of the pure-Rust model run (Phase 12).
//!
//! Since Phase 12 the model core is Rust-owned end to end (`crate::model`),
//! so the `snapwave_data` globals and the FFI state handoff
//! (`SnapWaveStateC`/`FfiState`) are gone. What remains here is the
//! scheduling state [`ModelState`] — current time, iteration counter,
//! next-output times, output counters — whose predicates mirror the legacy
//! Fortran `run_time_loop` logic exactly, plus [`rust_owns_gridfile`], the
//! grid-file extension dispatch used to decide whether the Rust mesh reader
//! owns a case.

use crate::input::SnapWaveInput;

/// The time-loop scheduling state. Mirrors the Fortran `run_time_loop` in
/// `src/snapwave_c_api.f90` (now retired), which the Rust `execute()` path
/// replaced in Phase 10 and `crate::model` now drives entirely in Rust.
///
/// # Scheduling invariants (mirror `run_time_loop`)
///
/// 1. Both output schedules start at `tstart` (first iteration outputs).
/// 2. Output fires when `t >= next_*_output - output_tol`.
/// 3. After output, `next_*_output` advances by the interval in a
///    `do while` loop (handles timestep > interval).
/// 4. History output additionally requires `nobs > 0` and a non-empty
///    `his_filename`.
/// 5. Map output is suppressed when `ja_save_each_iter != 0`.
/// 6. Time advances by `timestep` after output checks.
/// 7. The loop runs while `t <= tstop`.
#[derive(Debug, Clone)]
pub struct ModelState {
    /// Current model time (seconds since reference).
    pub t: f64,
    /// Iteration counter (1-based, incremented at the top of each iteration).
    pub it: i32,
    /// End time; the loop stops when `t > tstop`.
    pub tstop: f64,
    /// Time step (seconds).
    pub timestep: f32,
    /// Next scheduled map output time.
    pub next_map_output: f64,
    /// Next scheduled history output time.
    pub next_his_output: f64,
    /// Map output record counter (1-based, passed to the map writer).
    pub map_output_count: i32,
    /// History output record counter (1-based, passed to the his writer).
    pub his_output_count: i32,
    /// `max(1.0d-6, abs(dble(timestep))*1.0d-6)` — the Fortran loop's
    /// floating-point tolerance for output-time comparisons.
    pub output_tol: f64,
}

impl ModelState {
    /// Initialise the scheduling state exactly the way `run_time_loop` did.
    pub fn new(cfg: &SnapWaveInput) -> Self {
        ModelState {
            t: cfg.time.tstart,
            it: 0,
            tstop: cfg.time.tstop,
            timestep: cfg.time.timestep,
            next_map_output: cfg.time.tstart,
            next_his_output: cfg.time.tstart,
            map_output_count: 0,
            his_output_count: 0,
            output_tol: 1.0e-6_f64.max((cfg.time.timestep as f64).abs() * 1.0e-6),
        }
    }

    /// Whether the time loop should continue (`t <= tstop`).
    pub fn is_running(&self) -> bool {
        self.t <= self.tstop
    }

    /// Advance the iteration counter (called at the top of each iteration,
    /// before the solver step).
    pub fn advance_iteration(&mut self) {
        self.it += 1;
    }

    /// Advance model time by one timestep (called at the bottom of each
    /// iteration, after output checks).
    pub fn advance_time(&mut self) {
        self.t += self.timestep as f64;
    }

    /// Whether map output should fire at the current time. Mirrors the
    /// Fortran condition:
    ///   `ja_save_each_iter == 0 .and. map_filename /= '' .and.
    ///    t >= next_map_output - output_tol`
    pub fn should_output_map(&self, ja_save_each_iter: i32, map_file_nonempty: bool) -> bool {
        ja_save_each_iter == 0 && map_file_nonempty && self.t >= self.next_map_output - self.output_tol
    }

    /// Whether history output should fire at the current time. Mirrors the
    /// Fortran condition:
    ///   `his_filename /= '' .and. nobs > 0 .and.
    ///    t >= next_his_output - output_tol`
    pub fn should_output_his(&self, his_file_nonempty: bool, nobs: i32) -> bool {
        his_file_nonempty && nobs > 0 && self.t >= self.next_his_output - self.output_tol
    }

    /// Record a map output event: increment the counter and advance
    /// `next_map_output` by `map_interval` in a `do while` loop (handles
    /// the case where the timestep is larger than the output interval).
    pub fn record_map_output(&mut self, map_interval: f32) {
        self.map_output_count += 1;
        let interval = map_interval as f64;
        while self.next_map_output <= self.t + self.output_tol {
            self.next_map_output += interval;
        }
    }

    /// Record a history output event: increment the counter and advance
    /// `next_his_output` by `his_interval` in a `do while` loop.
    pub fn record_his_output(&mut self, his_interval: f32) {
        self.his_output_count += 1;
        let interval = his_interval as f64;
        while self.next_his_output <= self.t + self.output_tol {
            self.next_his_output += interval;
        }
    }
}

/// Does the Rust mesh reader own this gridfile? Mirrors the extension
/// dispatch of `initialize_snapwave_domain` *verbatim*, including the
/// quirks: the two characters after the **first** `.` decide (`a.ncd`
/// takes the NetCDF branch), and a name without any `.` compares its
/// first two characters (Fortran `index` returns 0, so `gridfile(1:2)`
/// is inspected). Since Phase 12 the Rust build only supports NetCDF
/// meshes; other formats are rejected with a clear error.
pub fn rust_owns_gridfile(gridfile: &str) -> bool {
    let after_dot = match gridfile.find('.') {
        Some(j) => &gridfile[j + 1..],
        None => gridfile,
    };
    after_dot.get(..2) == Some("nc")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::input::parse_str;

    #[test]
    fn model_state_initialises_like_the_fortran_loop() {
        let cfg = parse_str("timestep = 60\nmap_interval = 3600\nhis_interval = 1800\n").unwrap();
        let m = ModelState::new(&cfg);
        assert_eq!(m.t, cfg.time.tstart);
        assert_eq!(m.it, 0);
        assert_eq!(m.tstop, cfg.time.tstop);
        assert_eq!(m.timestep, 60.0);
        assert_eq!(m.next_map_output, cfg.time.tstart);
        assert_eq!(m.next_his_output, cfg.time.tstart);
        assert_eq!(m.map_output_count, 0);
        assert_eq!(m.his_output_count, 0);
        assert_eq!(m.output_tol, 60.0f64 * 1.0e-6);
        assert!(m.is_running());
    }

    #[test]
    fn model_state_output_tol_floor() {
        // For very small timesteps the absolute floor (1e-6) wins.
        let cfg = parse_str("timestep = 0.001\n").unwrap();
        let m = ModelState::new(&cfg);
        assert_eq!(m.output_tol, 1.0e-6);
    }

    #[test]
    fn model_state_iteration_and_time_advance() {
        let cfg = parse_str("timestep = 10\n").unwrap();
        let mut m = ModelState::new(&cfg);
        assert_eq!(m.it, 0);
        m.advance_iteration();
        assert_eq!(m.it, 1);
        assert_eq!(m.t, cfg.time.tstart);
        m.advance_time();
        assert_eq!(m.t, cfg.time.tstart + 10.0);
    }

    #[test]
    fn model_state_output_predicates() {
        let cfg = parse_str("timestep = 60\nmap_interval = 3600\nhis_interval = 1800\n").unwrap();
        let mut m = ModelState::new(&cfg);

        assert!(m.should_output_map(0, true), "map output at t=tstart");
        assert!(m.should_output_his(true, 5), "his output at t=tstart with nobs>0");

        assert!(!m.should_output_map(1, true), "ja_save_each_iter suppresses map");
        assert!(!m.should_output_map(0, false), "empty map_file suppresses map");
        assert!(!m.should_output_his(false, 5), "empty his_file suppresses his");
        assert!(!m.should_output_his(true, 0), "nobs==0 suppresses his");

        m.record_map_output(3600.0);
        assert_eq!(m.map_output_count, 1);
        assert!(m.next_map_output > m.t, "next_map_output advanced past current t");
        assert!(!m.should_output_map(0, true), "map output not due again yet");

        m.record_his_output(1800.0);
        assert_eq!(m.his_output_count, 1);
        assert!(m.next_his_output > m.t);
        assert!(!m.should_output_his(true, 5));
    }

    #[test]
    fn model_state_handles_timestep_larger_than_interval() {
        let cfg = parse_str("timestep = 7200\nmap_interval = 3600\n").unwrap();
        let mut m = ModelState::new(&cfg);

        assert!(m.should_output_map(0, true));
        m.record_map_output(3600.0);
        assert_eq!(m.map_output_count, 1);
        assert!(!m.should_output_map(0, true));

        m.advance_time();
        assert!(m.should_output_map(0, true));
    }

    #[test]
    fn model_state_stops_at_tstop() {
        let cfg = parse_str("timestep = 10\n").unwrap();
        let mut m = ModelState::new(&cfg);
        m.t = m.tstop + 1.0;
        assert!(!m.is_running());
    }

    #[test]
    fn gridfile_dispatch_mirrors_the_fortran_extension_check() {
        assert!(rust_owns_gridfile("mesh.nc"), "'nc' extension is Rust-owned");
        assert!(rust_owns_gridfile("run/mesh.nc"));
        assert!(rust_owns_gridfile(".nc"), "bare '.nc' works");
        assert!(rust_owns_gridfile("a.ncd"), "quirk: only two chars after the first dot are compared");
        assert!(rust_owns_gridfile("ncmesh"), "quirk: no dot inspects the first two characters");
        assert!(!rust_owns_gridfile("mesh.txt"));
        assert!(!rust_owns_gridfile(".txt"));
        assert!(!rust_owns_gridfile("mesh.NET"));
        assert!(!rust_owns_gridfile("n"));
        assert!(!rust_owns_gridfile("snapshot.vnc"));
        assert!(!rust_owns_gridfile(""));
    }
}
