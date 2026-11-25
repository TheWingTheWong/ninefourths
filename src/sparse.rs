//! Sparse matrices with entries in F_p[ε, ε⁻¹]: the encoders and decoders of
//! a packing.

use crate::laurent::Laurent;
use std::collections::HashMap;

#[derive(Clone, Debug)]
pub struct Sparse {
    pub nrows: usize,
    pub ncols: usize,
    /// Row-major: for each row, (column, entry) with distinct columns, nonzero entries.
    pub rows: Vec<Vec<(usize, Laurent)>>,
}

impl Sparse {
    pub fn zeros(nrows: usize, ncols: usize) -> Self {
        Sparse { nrows, ncols, rows: vec![vec![]; nrows] }
    }
    pub fn identity(n: usize) -> Self {
        let mut m = Self::zeros(n, n);
        for i in 0..n {
            m.rows[i].push((i, Laurent::constant(1)));
        }
        m
    }
    /// Add `x` to entry (i, j).
    pub fn add_entry(&mut self, i: usize, j: usize, x: &Laurent) {
        assert!(i < self.nrows && j < self.ncols, "entry ({i},{j}) outside {}x{}", self.nrows, self.ncols);
        if x.is_zero() {
            return;
        }
        let row = &mut self.rows[i];
        if let Some(pos) = row.iter().position(|(c, _)| *c == j) {
            row[pos].1.add_assign(x);
            if row[pos].1.is_zero() {
                row.swap_remove(pos);
            }
        } else {
            row.push((j, x.clone()));
        }
    }
    pub fn nnz(&self) -> usize {
        self.rows.iter().map(|r| r.len()).sum()
    }
    pub fn apply(&self, x: &[Laurent]) -> Vec<Laurent> {
        assert_eq!(x.len(), self.ncols, "apply: dimension mismatch");
        self.rows
            .iter()
            .map(|row| {
                let mut acc = Laurent::zero();
                for (j, a) in row {
                    if !x[*j].is_zero() {
                        acc.add_assign(&a.mul(&x[*j]));
                    }
                }
                acc
            })
            .collect()
    }
    /// self ∘ other (apply `other` first).
    pub fn compose(&self, other: &Sparse) -> Sparse {
        assert_eq!(self.ncols, other.nrows, "compose: dimension mismatch");
        let mut out = Sparse::zeros(self.nrows, other.ncols);
        for (i, row) in self.rows.iter().enumerate() {
            let mut acc: HashMap<usize, Laurent> = HashMap::new();
            for (k, a) in row {
                for (j, b) in &other.rows[*k] {
                    acc.entry(*j).or_default().add_assign(&a.mul(b));
                }
            }
            let mut r: Vec<_> = acc.into_iter().filter(|(_, v)| !v.is_zero()).collect();
            r.sort_by_key(|e| e.0);
            out.rows[i] = r;
        }
        out
    }
    pub fn transpose(&self) -> Sparse {
        let mut out = Sparse::zeros(self.ncols, self.nrows);
        for (i, row) in self.rows.iter().enumerate() {
            for (j, a) in row {
                out.rows[*j].push((i, a.clone()));
            }
        }
        out
    }
    /// Kronecker product: entry ((i,i'),(j,j')) = a_ij · b_i'j', index i·b.nrows + i'.
    pub fn kron(&self, b: &Sparse) -> Sparse {
        let mut out = Sparse::zeros(self.nrows * b.nrows, self.ncols * b.ncols);
        for (i, row) in self.rows.iter().enumerate() {
            for ib in 0..b.nrows {
                let r = &mut out.rows[i * b.nrows + ib];
                for (j, x) in row {
                    for (jb, y) in &b.rows[ib] {
                        r.push((j * b.ncols + jb, x.mul(y)));
                    }
                }
            }
        }
        out
    }
    /// Block-diagonal sum.
    pub fn block_diag(blocks: &[Sparse]) -> Sparse {
        let nr = blocks.iter().map(|b| b.nrows).sum();
        let nc = blocks.iter().map(|b| b.ncols).sum();
        let mut out = Sparse::zeros(nr, nc);
        let (mut r0, mut c0) = (0, 0);
        for b in blocks {
            for (i, row) in b.rows.iter().enumerate() {
                out.rows[r0 + i] = row.iter().map(|(j, x)| (j + c0, x.clone())).collect();
            }
            r0 += b.nrows;
            c0 += b.ncols;
        }
        out
    }
    /// Permutation/selection matrix: output i = input src[i] (None = zero row).
    pub fn select(src: &[Option<usize>], ncols: usize) -> Sparse {
        let mut m = Sparse::zeros(src.len(), ncols);
        for (i, s) in src.iter().enumerate() {
            if let Some(j) = s {
                m.rows[i].push((*j, Laurent::constant(1)));
            }
        }
        m
    }
    /// Substitute ε → ε^n in every entry.
    pub fn subst(&self, n: i64) -> Sparse {
        Sparse {
            nrows: self.nrows,
            ncols: self.ncols,
            rows: self.rows.iter().map(|r| r.iter().map(|(j, x)| (*j, x.subst(n))).collect()).collect(),
        }
    }
    /// (min, max) ε exponent over all entries; (0, 0) for an ε-free map.
    pub fn eps_range(&self) -> (i64, i64) {
        let mut lo = 0;
        let mut hi = 0;
        for row in &self.rows {
            for (_, x) in row {
                lo = lo.min(x.valuation().unwrap());
                hi = hi.max(x.degree().unwrap());
            }
        }
        (lo, hi)
    }
}

impl Sparse {
    /// Evaluate every entry at ε = e (nonzero), giving an ε-free map.
    pub fn at_eps(&self, e: u64) -> Sparse {
        Sparse {
            nrows: self.nrows,
            ncols: self.ncols,
            rows: self
                .rows
                .iter()
                .map(|r| r.iter().map(|(j, x)| (*j, Laurent::constant(x.eval(e)))).filter(|(_, x)| !x.is_zero()).collect())
                .collect(),
        }
    }
    pub fn scale(&self, x: &Laurent) -> Sparse {
        Sparse {
            nrows: self.nrows,
            ncols: self.ncols,
            rows: self.rows.iter().map(|r| r.iter().map(|(j, y)| (*j, y.mul(x))).filter(|(_, y)| !y.is_zero()).collect()).collect(),
        }
    }
}
