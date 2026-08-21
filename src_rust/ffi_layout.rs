//! Localized Fortran array-layout conversions (plan.md Phase 8, step 3).
//!
//! SnapWave's Fortran core stores every array column-major with one-based
//! indices; Rust is row-major and zero-based. Rather than sprinkling
//! `+1`s and stride arithmetic through the migration code, the layout
//! facts live here once:
//!
//! * [`fortran_index`] / [`rust_index`] convert one array index between
//!   the zero-based Rust and one-based Fortran conventions. Fortran index
//!   `0` (and negatives) are legal *stored values* in SnapWave (the
//!   `-999` "no fourth node" sentinel of `face_nodes`), but they are
//!   never valid *subscripts*, so [`rust_index`] rejects them.
//! * [`ColMajor`] maps a multi-dimensional index onto the flat offset a
//!   Fortran array of that shape occupies in memory: element `(i1, i2,
//! * ...)` of `A(d1, d2, ...)` sits at `i1 + d1*(i2-1) + ...` — the first
//!   dimension is contiguous. Both zero-based and one-based entry points
//!   are provided; the one-based one is the greppable mirror of the
//!   Fortran declaration.
//!
//! A pleasant consequence pinned by the tests below: two of the layouts
//! this crate already uses are *bit-compatible* with the Fortran memory
//! order and need no data movement —
//!
//! * `mesh.face_nodes` (`[face*4 + node]`, node-major within a face) is
//!   exactly the column-major flattening of `face_nodes(4, no_faces)`;
//! * the time-major series layout of `text_input::BoundarySeries`
//!   (`[itb*nwbnd + ib]`) is exactly the column-major flattening of
//!   `hs_bwv(nwbnd, ntwbnd)`.
//!
//! The helpers therefore document and verify those equivalences instead
//! of reshuffling buffers, and the FFI handoff (`state.rs`) can hand the
//! raw buffers to Fortran directly.

// The one-based index helpers and the general `from_row_major` reshuffle
// have no production consumer yet: the Phase 8 handoff buffers are
// already Fortran-compatible, so only the tests exercise them. Phases
// 9-10 (geometry, interpolation, solver state) consume them heavily;
// like `paths::RunPaths` before Phase 6, the API exists so those phases
// do not re-derive the layout rules.
#![allow(dead_code)]

/// Convert a zero-based Rust subscript to its one-based Fortran value.
pub fn fortran_index(i0: usize) -> i64 {
    i0 as i64 + 1
}

/// Convert a one-based Fortran subscript (or stored sentinel value) to a
/// zero-based Rust subscript. Subscripts below 1 are not valid array
/// positions and yield `None` (callers decide whether that is an error
/// or a sentinel).
pub fn rust_index(i1: i64) -> Option<usize> {
    if i1 >= 1 {
        Some((i1 - 1) as usize)
    } else {
        None
    }
}

/// Column-major layout of a Fortran array `A(d1, d2, ...)`: the first
/// dimension is contiguous, exactly as Fortran stores it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ColMajor {
    dims: Vec<usize>,
}

impl ColMajor {
    /// Layout of an array with the given dimensions (in Fortran
    /// declaration order). Panics only in the blatant contract violation
    /// of a zero rank.
    pub fn new(dims: &[usize]) -> Self {
        assert!(!dims.is_empty(), "Fortran arrays have at least one dimension");
        ColMajor { dims: dims.to_vec() }
    }

    /// The Fortran declaration-order dimensions.
    pub fn dims(&self) -> &[usize] {
        &self.dims
    }

    /// Total number of elements.
    pub fn len(&self) -> usize {
        self.dims.iter().product()
    }

    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Flat offset of a **zero-based** multi-index (Rust order is
    /// irrelevant here — the tuple follows Fortran declaration order).
    /// `None` when the index is out of bounds or has the wrong rank.
    pub fn offset(&self, idx0: &[usize]) -> Option<usize> {
        if idx0.len() != self.dims.len() {
            return None;
        }
        let mut offset = 0usize;
        let mut stride = 1usize;
        for (&i, &d) in idx0.iter().zip(self.dims.iter()) {
            if i >= d {
                return None;
            }
            offset += i * stride;
            stride *= d;
        }
        Some(offset)
    }

    /// Flat offset of a **one-based** Fortran subscript tuple, mirroring
    /// `A(i1, i2, ...)`. `None` when any subscript is outside `1..=d`.
    pub fn offset_fortran(&self, idx1: &[i64]) -> Option<usize> {
        if idx1.len() != self.dims.len() {
            return None;
        }
        let idx0: Vec<usize> = idx1
            .iter()
            .zip(self.dims.iter())
            .map(|(&i, &d)| if i >= 1 && (i as usize) <= d { (i - 1) as usize } else { usize::MAX })
            .collect();
        // usize::MAX marks out-of-range subscripts; offset() re-validates.
        self.offset(&idx0)
    }

    /// Inverse of [`ColMajor::offset`]: the zero-based multi-index stored
    /// at a flat offset. `None` when the offset is out of range.
    pub fn indices(&self, offset: usize) -> Option<Vec<usize>> {
        if offset >= self.len() {
            return None;
        }
        let mut idx = vec![0usize; self.dims.len()];
        let mut rem = offset;
        for (k, &d) in self.dims.iter().enumerate() {
            idx[k] = rem % d;
            rem /= d;
        }
        Some(idx)
    }

    /// Re-lay a flat buffer stored in **C row-major order** (first
    /// dimension slowest) with the *same logical dimensions, in the same
    /// order*, as this column-major layout: element `(i1, i2, ...)` moves
    /// from `i1*d2 + i2 + ...` to `i1 + d1*i2 + ...`. This is the
    /// general reshuffle for arrays whose Rust layout is NOT already
    /// Fortran-compatible; the SnapWave handoffs that are compatible
    /// skip it (see the module docs). `None` when the shapes disagree.
    pub fn from_row_major<T: Clone>(&self, data: &[T], dims: &[usize]) -> Option<Vec<T>> {
        if dims.len() != self.dims.len() || dims != self.dims || data.len() != self.len() {
            return None;
        }
        let rank = self.dims.len();
        let mut out = vec![None; data.len()];
        for (flat, value) in data.iter().enumerate() {
            // Row-major decode: the first dimension is the slowest.
            let mut idx = vec![0usize; rank];
            let mut rem = flat;
            for k in (0..rank).rev() {
                idx[k] = rem % self.dims[k];
                rem /= self.dims[k];
            }
            let foff = self.offset(&idx)?;
            out[foff] = Some(value.clone());
        }
        out.into_iter().collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn index_conversions_round_trip() {
        for i0 in [0usize, 1, 2, 41, 999999] {
            let i1 = fortran_index(i0);
            assert_eq!(i1, i0 as i64 + 1);
            assert_eq!(rust_index(i1), Some(i0));
        }
    }

    #[test]
    fn fortran_subscripts_below_one_are_rejected() {
        // 0 and negatives are stored sentinels (e.g. face_nodes -999),
        // never valid subscripts.
        assert_eq!(rust_index(0), None);
        assert_eq!(rust_index(-999), None);
        assert_eq!(rust_index(-1), None);
    }

    #[test]
    fn column_major_first_dimension_is_contiguous() {
        // face_nodes(4, no_faces): (j, k) at (j-1) + 4*(k-1), the exact
        // formula a Fortran compiler uses.
        let cm = ColMajor::new(&[4, 7]);
        assert_eq!(cm.len(), 28);
        assert_eq!(cm.offset_fortran(&[1, 1]), Some(0));
        assert_eq!(cm.offset_fortran(&[4, 1]), Some(3));
        assert_eq!(cm.offset_fortran(&[1, 2]), Some(4));
        assert_eq!(cm.offset_fortran(&[4, 7]), Some(27));
        // The Rust node-major-within-face layout [face*4 + node].
        for k in 1..=7i64 {
            for j in 1..=4i64 {
                let rust_flat = ((k - 1) * 4 + (j - 1)) as usize;
                assert_eq!(cm.offset_fortran(&[j, k]), Some(rust_flat));
            }
        }
    }

    #[test]
    fn column_major_time_major_equivalence() {
        // hs_bwv(nwbnd, ntwbnd): (ib, itb) at (ib-1) + nwbnd*(itb-1),
        // which equals the time-major Rust layout [itb*nwbnd + ib].
        let (nwbnd, ntwbnd) = (3usize, 5usize);
        let cm = ColMajor::new(&[nwbnd, ntwbnd]);
        for itb in 0..ntwbnd {
            for ib in 0..nwbnd {
                let rust_flat = itb * nwbnd + ib;
                assert_eq!(
                    cm.offset_fortran(&[fortran_index(ib), fortran_index(itb)]),
                    Some(rust_flat)
                );
            }
        }
    }

    #[test]
    fn offset_rejects_out_of_bounds_and_wrong_rank() {
        let cm = ColMajor::new(&[4, 7]);
        assert_eq!(cm.offset_fortran(&[0, 1]), None, "subscript 0 is not a position");
        assert_eq!(cm.offset_fortran(&[5, 1]), None, "first dim too large");
        assert_eq!(cm.offset_fortran(&[1, 8]), None, "second dim too large");
        assert_eq!(cm.offset_fortran(&[-999, 1]), None, "sentinel is not a position");
        assert_eq!(cm.offset_fortran(&[1]), None, "wrong rank");
        assert_eq!(cm.offset(&[4, 7]), None, "zero-based out of bounds");
        assert_eq!(cm.offset(&[0, 6]), Some(24));
    }

    #[test]
    fn indices_invert_offset_in_3d() {
        // A 3-D case mirrors the solver's prev(2, ntheta, no_nodes).
        let cm = ColMajor::new(&[2, 3, 4]);
        assert_eq!(cm.len(), 24);
        for off in 0..24 {
            let idx = cm.indices(off).expect("every offset decodes");
            assert_eq!(cm.offset(&idx), Some(off), "round trip at {off}");
        }
        assert_eq!(cm.indices(24), None);
        // Fortran formula: (i-1) + 2*(j-1) + 6*(k-1).
        assert_eq!(cm.offset_fortran(&[2, 3, 4]), Some(1 + 2 * 2 + 6 * 3));
    }

    #[test]
    fn zero_sized_dimensions_are_handled() {
        // Empty inputs associate as zero-extent Fortran arrays (the
        // absent-data convention of the Phase 8 state handoff).
        let cm = ColMajor::new(&[0, 5]);
        assert!(cm.is_empty());
        assert_eq!(cm.offset_fortran(&[1, 1]), None);
        assert_eq!(cm.indices(0), None);
    }

    #[test]
    fn row_major_sources_are_permuted() {
        // A logical 2x3 matrix stored row-major (C order, first
        // dimension slowest) becomes the same matrix stored
        // column-major: M = [[1,2,3],[4,5,6]] has A(1,1)=1, A(2,1)=4,
        // A(1,2)=2, ... so the Fortran flat buffer is [1,4,2,5,3,6].
        let cm = ColMajor::new(&[2, 3]);
        let src = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0];
        let out = cm.from_row_major(&src, &[2, 3]).expect("shapes must match");
        assert_eq!(out, vec![1.0, 4.0, 2.0, 5.0, 3.0, 6.0]);
        // Round trip: every logical element keeps its value.
        for i1 in 0..2usize {
            for i2 in 0..3usize {
                let foff = cm.offset(&[i1, i2]).unwrap();
                assert_eq!(out[foff], src[i1 * 3 + i2], "element ({i1},{i2})");
            }
        }
        // Mismatched shapes are refused rather than mis-permuted.
        assert!(cm.from_row_major(&src, &[3, 2]).is_none());
        assert!(cm.from_row_major(&src, &[2, 3, 1]).is_none());
        assert!(cm.from_row_major(&src[..5], &[2, 3]).is_none());
    }
}
