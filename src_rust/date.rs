//! Date/time conversion utilities, ported from `src/snapwave_date.f90`
//! (plan.md, Phase 9: "Date/time conversion utilities").
//!
//! # What this module covers
//!
//! `snapwave_date` holds four routines; three are ported here and one is
//! deliberately not:
//!
//! * [`julian_date`] — the Fliegel & Van Flandern Julian-day formula
//!   (`julian_date` in `snapwave_date.f90`, also inlined by the input
//!   parser since Phase 3).
//! * [`parse_date15`] — the fixed-position `'(I4,2I2,1X,3I2)'` read used by
//!   every date consumer (`time_difference`, `date_to_iso8601`).
//! * [`seconds_between`] — `time_difference(datespw, datesim, dtsec)`, the
//!   seconds between two `yyyymmdd hhmmss` strings.
//! * [`date_to_iso8601`] — `date_to_iso8601`, used by `snapwave_ncoutput`
//!   to build the NetCDF time reference string (`20240417 000000` →
//!   `2024-04-17 00:00:00`).
//!
//! `convert_fewsdate` is **not** ported: it has no callers anywhere in the
//! tree and reads its `trefstr` from an uninitialised local `character*41`
//! (it is dead code with a latent bug), so there is no oracle behaviour to
//! preserve. It is documented here rather than silently reproduced.
//!
//! # Why the port exists now
//!
//! Until Phase 7 the Rust side only needed `julian_date` + `seconds_between`
//! (for `tstart`/`tstop`), which lived inline in `src_rust/input.rs`. Phase 9
//! moves the *date utilities* as a unit (plan.md Phase 9, step 1: "split
//! pure utility routines from file-reading and global-state routines"), so
//! the canonical implementations now live here and `input` delegates to them.
//! `date_to_iso8601` is ported and unit-tested for parity even though the run
//! path still receives the formatted string from Fortran's capture stream —
//! the NetCDF header generation is the Phase 7 design, not this phase's.

use anyhow::{bail, Result};

/// Width of the Fortran `character*15` date globals (`trefstr`,
/// `tstartstr`, `tstopstr`).
const WIDTH_DATE: usize = 15;

/// Calendar date and time fields (`yyyy, mm, dd, hh, mn, ss`).
pub type DateFields = (i32, i32, i32, i32, i32, i32);

/// Parse the fixed-position date read by `snapwave_date` with the format
/// `'(I4,2I2,1X,3I2)'` from a `character*15` value: `yyyymmdd hhmmss`.
///
/// Position 9 is skipped entirely (`1X`); blanks inside numeric fields are
/// ignored (Fortran `BLANK='NULL'` default), blank fields read as zero, and
/// any other non-digit character is an error (a Fortran formatted-read
/// runtime error). A value shorter than 15 characters is blank-padded first,
/// matching the fixed `character*15` receiver.
pub fn parse_date15(s: &str) -> Result<DateFields> {
    let padded = format!("{s:<width$}", width = WIDTH_DATE);
    let b = padded.as_bytes();
    let field = |range: std::ops::Range<usize>| -> Result<i32> {
        let mut digits = String::new();
        for &c in &b[range] {
            match c {
                b' ' => continue,
                b'0'..=b'9' => digits.push(c as char),
                _ => bail!("invalid character '{}' in date '{s}' (expected yyyymmdd hhmmss)", c as char),
            }
        }
        if digits.is_empty() {
            Ok(0)
        } else {
            Ok(digits.parse::<i32>()?)
        }
    };
    Ok((
        field(0..4)?,   // yyyy
        field(4..6)?,   // mm
        field(6..8)?,   // dd
        field(9..11)?,  // hh (position 8 is the skipped separator)
        field(11..13)?, // mn
        field(13..15)?, // ss
    ))
}

/// Fliegel & Van Flandern Julian day number, identical to `julian_date`
/// in `src/snapwave_date.f90`. Both Fortran and Rust integer division
/// truncate toward zero, which this formula relies on for months < 3.
/// (Computed in i64 so pathological date ranges cannot overflow; Fortran
/// 32-bit integers would wrap only beyond ~68-year spans.)
pub fn julian_date(yyyy: i32, mm: i32, dd: i32) -> i64 {
    let (yyyy, mm, dd) = (yyyy as i64, mm as i64, dd as i64);
    dd - 32075 + 1461 * (yyyy + 4800 + (mm - 14) / 12) / 4
        + 367 * (mm - 2 - ((mm - 14) / 12) * 12) / 12
        - 3 * ((yyyy + 4900 + (mm - 14) / 12) / 100) / 4
}

/// Seconds between two `yyyymmdd hhmmss` strings (date2 - date1), as
/// `time_difference` in `src/snapwave_date.f90` computes for the globals
/// `tstart`/`tstop` (seconds relative to `tref`).
pub fn seconds_between(date1: &str, date2: &str) -> Result<f64> {
    let (y1, m1, d1, h1, n1, s1) = parse_date15(date1)?;
    let (y2, m2, d2, h2, n2, s2) = parse_date15(date2)?;
    let jul1 = julian_date(y1, m1, d1);
    let jul2 = julian_date(y2, m2, d2);
    let sec1 = (h1 as i64) * 3600 + (n1 as i64) * 60 + s1 as i64;
    let sec2 = (h2 as i64) * 3600 + (n2 as i64) * 60 + s2 as i64;
    Ok(((jul2 - jul1) * 86400 + sec2 - sec1) as f64)
}

/// `date_to_iso8601` of `src/snapwave_date.f90`: render a `yyyymmdd
/// hhmmss` value as `yyyy-mm-dd HH:MM:SS` (space-separated, not `T` — the
/// DFM/FEWS variant, see the commented-out alternative in the Fortran).
///
/// The Fortran write format is `'(I4,A1,I0.2,A1,I0.2,A1,I0.2,A1,I0.2,A1,I0.2)'`:
/// the year is right-justified in width 4 (space-padded, not zero-padded)
/// and the five other fields are zero-padded to at least two digits.
// `allow(dead_code)`: the run path still receives the formatted string from
// Fortran's capture stream (the Phase 7 NetCDF-header design), so this port
// is pinned by unit tests only until the header generation moves to Rust.
#[allow(dead_code)]
pub fn date_to_iso8601(date_string: &str) -> Result<String> {
    let (yyyy, mm, dd, hh, mn, ss) = parse_date15(date_string)?;
    Ok(format!("{yyyy:>4}-{mm:02}-{dd:02} {hh:02}:{mn:02}:{ss:02}"))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn julian_day_matches_reference_values() {
        // Reference values from the snapwave_date.f90 header and the
        // Fliegel & Van Flandern paper: 1970-01-01 -> 2440588,
        // 2000-01-01 -> 2451545.
        assert_eq!(julian_date(1970, 1, 1), 2440588);
        assert_eq!(julian_date(2000, 1, 1), 2451545);
        // Consistency across the March boundary where the formula's
        // truncating division matters (Jan/Feb belong to the previous
        // "Roman" year).
        assert_eq!(julian_date(2000, 2, 28) - julian_date(2000, 1, 31), 28);
        assert_eq!(julian_date(2000, 3, 1) - julian_date(2000, 2, 28), 2); // leap day
        assert_eq!(julian_date(2001, 3, 1) - julian_date(2001, 2, 28), 1);
    }

    #[test]
    fn date_parsing_mirrors_the_fortran_format() {
        // Position 9 is skipped entirely (the '1X' in '(I4,2I2,1X,3I2)').
        assert_eq!(parse_date15("20240417 000000").unwrap(), (2024, 4, 17, 0, 0, 0));
        assert_eq!(parse_date15("20240417T010203").unwrap(), (2024, 4, 17, 1, 2, 3));
        // Blanks are ignored inside numeric fields; blank fields read 0,
        // so a date-only value means midnight.
        assert_eq!(parse_date15("20240417").unwrap(), (2024, 4, 17, 0, 0, 0));
        assert_eq!(parse_date15("  240417 000000").unwrap(), (24, 4, 17, 0, 0, 0));
        // Non-digits are an error (a Fortran formatted-read abort).
        assert!(parse_date15("notadate1234567").is_err());
        assert!(parse_date15("2024041X 000000").is_err());
    }

    #[test]
    fn seconds_between_computes_signed_differences() {
        assert_eq!(seconds_between("20240417 000000", "20240417 000000").unwrap(), 0.0);
        assert_eq!(seconds_between("20240417 000000", "20240418 000000").unwrap(), 86400.0);
        assert_eq!(seconds_between("20240417 000000", "20240416 235959").unwrap(), -1.0);
        // Cross month and year boundaries.
        assert_eq!(seconds_between("20240131 235959", "20240201 000000").unwrap(), 1.0);
        assert_eq!(seconds_between("20231231 235959", "20240101 000000").unwrap(), 1.0);
    }

    #[test]
    fn iso8601_renders_the_fews_variant() {
        assert_eq!(date_to_iso8601("20240417 000000").unwrap(), "2024-04-17 00:00:00");
        assert_eq!(date_to_iso8601("20240417 010203").unwrap(), "2024-04-17 01:02:03");
        // A date-only value is midnight.
        assert_eq!(date_to_iso8601("20240417").unwrap(), "2024-04-17 00:00:00");
        // Single-digit fields are zero-padded to width 2.
        assert_eq!(date_to_iso8601("20240101 010101").unwrap(), "2024-01-01 01:01:01");
    }
}
