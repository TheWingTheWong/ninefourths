//! Laurent polynomials in the formal parameter ε with coefficients in F_p.
//!
//! Packings "modulo ε" have encoder and decoder entries in F_p[ε, ε⁻¹]; the
//! verifier carries ε symbolically, so "the error is divisible by ε" is
//! checked exactly rather than numerically.

use crate::field as f;

/// `c[i]` is the coefficient of ε^(lo + i). Always trimmed: no zero
/// coefficients at either end; the zero polynomial has empty `c`.
#[derive(Clone, Debug, PartialEq, Eq, Default)]
pub struct Laurent {
    pub lo: i64,
    pub c: Vec<u64>,
}

impl Laurent {
    pub fn zero() -> Self {
        Laurent { lo: 0, c: vec![] }
    }
    pub fn constant(x: u64) -> Self {
        Self::mono(x, 0)
    }
    pub fn mono(x: u64, e: i64) -> Self {
        if x == 0 { Self::zero() } else { Laurent { lo: e, c: vec![x] } }
    }
    pub fn int(x: i64) -> Self {
        Self::constant(f::from_i64(x))
    }
    pub fn is_zero(&self) -> bool {
        self.c.is_empty()
    }
    /// Lowest exponent present (None for zero).
    pub fn valuation(&self) -> Option<i64> {
        if self.is_zero() { None } else { Some(self.lo) }
    }
    /// Highest exponent present (None for zero).
    pub fn degree(&self) -> Option<i64> {
        if self.is_zero() { None } else { Some(self.lo + self.c.len() as i64 - 1) }
    }
    pub fn coeff(&self, e: i64) -> u64 {
        let i = e - self.lo;
        if i < 0 || i >= self.c.len() as i64 { 0 } else { self.c[i as usize] }
    }
    fn trim(mut self) -> Self {
        while let Some(&0) = self.c.last() {
            self.c.pop();
        }
        let lead = self.c.iter().take_while(|&&x| x == 0).count();
        if lead == self.c.len() {
            return Self::zero();
        }
        if lead > 0 {
            self.c.drain(..lead);
            self.lo += lead as i64;
        }
        self
    }
    pub fn add(&self, o: &Laurent) -> Laurent {
        if self.is_zero() {
            return o.clone();
        }
        if o.is_zero() {
            return self.clone();
        }
        let lo = self.lo.min(o.lo);
        let hi = self.degree().unwrap().max(o.degree().unwrap());
        let mut c = vec![0u64; (hi - lo + 1) as usize];
        for (i, &x) in self.c.iter().enumerate() {
            let k = (self.lo - lo) as usize + i;
            c[k] = f::add(c[k], x);
        }
        for (i, &x) in o.c.iter().enumerate() {
            let k = (o.lo - lo) as usize + i;
            c[k] = f::add(c[k], x);
        }
        Laurent { lo, c }.trim()
    }
    pub fn add_assign(&mut self, o: &Laurent) {
        if o.is_zero() {
            return;
        }
        if self.is_zero() {
            *self = o.clone();
            return;
        }
        // fast path: o inside self's range
        let (slo, shi) = (self.lo, self.degree().unwrap());
        let (olo, ohi) = (o.lo, o.degree().unwrap());
        if olo >= slo && ohi <= shi {
            for (i, &x) in o.c.iter().enumerate() {
                let k = (olo - slo) as usize + i;
                self.c[k] = f::add(self.c[k], x);
            }
            let t = std::mem::take(self);
            *self = t.trim();
        } else {
            *self = self.add(o);
        }
    }
    pub fn neg(&self) -> Laurent {
        Laurent { lo: self.lo, c: self.c.iter().map(|&x| f::neg(x)).collect() }
    }
    pub fn sub(&self, o: &Laurent) -> Laurent {
        self.add(&o.neg())
    }
    pub fn scale(&self, x: u64) -> Laurent {
        if x == 0 {
            return Self::zero();
        }
        Laurent { lo: self.lo, c: self.c.iter().map(|&y| f::mul(x, y)).collect() }
    }
    pub fn shift(&self, e: i64) -> Laurent {
        if self.is_zero() { Self::zero() } else { Laurent { lo: self.lo + e, c: self.c.clone() } }
    }
    pub fn mul(&self, o: &Laurent) -> Laurent {
        if self.is_zero() || o.is_zero() {
            return Self::zero();
        }
        if self.c.len() == 1 {
            return Laurent { lo: self.lo + o.lo, c: o.c.iter().map(|&y| f::mul(self.c[0], y)).collect() };
        }
        if o.c.len() == 1 {
            return Laurent { lo: self.lo + o.lo, c: self.c.iter().map(|&y| f::mul(o.c[0], y)).collect() };
        }
        let mut c = vec![0u64; self.c.len() + o.c.len() - 1];
        for (i, &x) in self.c.iter().enumerate() {
            if x == 0 {
                continue;
            }
            for (j, &y) in o.c.iter().enumerate() {
                c[i + j] = f::add(c[i + j], f::mul(x, y));
            }
        }
        Laurent { lo: self.lo + o.lo, c }.trim()
    }
    /// Substitute ε → ε^n (n ≥ 1).
    pub fn subst(&self, n: i64) -> Laurent {
        assert!(n >= 1);
        if self.is_zero() || n == 1 {
            return self.clone();
        }
        let mut c = vec![0u64; (self.c.len() - 1) * n as usize + 1];
        for (i, &x) in self.c.iter().enumerate() {
            c[i * n as usize] = x;
        }
        Laurent { lo: self.lo * n, c }
    }
    /// Evaluate at a nonzero field element.
    pub fn eval(&self, e: u64) -> u64 {
        if self.is_zero() {
            return 0;
        }
        let mut acc = 0u64;
        for &x in self.c.iter().rev() {
            acc = f::add(f::mul(acc, e), x);
        }
        let s = if self.lo >= 0 { f::pow(e, self.lo as u64) } else { f::inv(f::pow(e, (-self.lo) as u64)) };
        f::mul(acc, s)
    }
}

/// Inverse of a Laurent *monomial-times-unit* is not generally Laurent; we only
/// ever divide by monomials times scalars, which this helper supports.
pub fn inv_monomial(x: &Laurent) -> Laurent {
    assert!(x.c.len() == 1, "only monomials are invertible in F_p[ε, ε⁻¹]");
    Laurent::mono(f::inv(x.c[0]), -x.lo)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::field::Rng;

    fn rand(r: &mut Rng) -> Laurent {
        let lo = r.below(7) as i64 - 3;
        let n = r.below(5) as usize;
        Laurent { lo, c: (0..n).map(|_| r.elem()).collect() }.trim()
    }

    #[test]
    fn ring_laws_and_eval() {
        let mut r = Rng(7);
        for _ in 0..500 {
            let (a, b, c) = (rand(&mut r), rand(&mut r), rand(&mut r));
            assert_eq!(a.mul(&b.add(&c)), a.mul(&b).add(&a.mul(&c)));
            assert_eq!(a.mul(&b), b.mul(&a));
            let e = r.elem().max(1);
            assert_eq!(a.mul(&b).eval(e), f::mul(a.eval(e), b.eval(e)));
            assert_eq!(a.subst(3).eval(e), a.eval(f::pow(e, 3)));
            let mut s = a.clone();
            s.add_assign(&b);
            assert_eq!(s, a.add(&b));
        }
    }
}
