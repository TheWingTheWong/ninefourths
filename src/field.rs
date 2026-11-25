//! Arithmetic in F_p with p = 2^61 − 1 (a Mersenne prime).
//!
//! All identities are rational, so they hold over F_p for every p that does
//! not divide a denominator; checking them at random points of F_p is the
//! Schwartz–Zippel test (a wrong identity of degree D survives one random
//! point with probability at most D/p ≈ D·4.3e-19).

pub const P: u64 = (1u64 << 61) - 1;

#[inline]
pub fn reduce128(x: u128) -> u64 {
    let lo = (x as u64) & P;
    let hi = (x >> 61) as u64;
    let mut s = lo + (hi & P) + ((x >> 122) as u64);
    while s >= P {
        s -= P;
    }
    s
}

#[inline]
pub fn add(a: u64, b: u64) -> u64 {
    let s = a + b;
    if s >= P { s - P } else { s }
}

#[inline]
pub fn sub(a: u64, b: u64) -> u64 {
    if a >= b { a - b } else { a + P - b }
}

#[inline]
pub fn neg(a: u64) -> u64 {
    if a == 0 { 0 } else { P - a }
}

#[inline]
pub fn mul(a: u64, b: u64) -> u64 {
    reduce128(a as u128 * b as u128)
}

pub fn pow(mut a: u64, mut e: u64) -> u64 {
    let mut r = 1;
    while e > 0 {
        if e & 1 == 1 {
            r = mul(r, a);
        }
        a = mul(a, a);
        e >>= 1;
    }
    r
}

pub fn inv(a: u64) -> u64 {
    assert!(a != 0, "inverse of zero in F_p");
    pow(a, P - 2)
}

pub fn from_i64(x: i64) -> u64 {
    if x >= 0 { (x as u64) % P } else { neg(((-x) as u64) % P) }
}

/// Map a field element back to the symmetric integer range (−p/2, p/2).
pub fn to_i64(x: u64) -> i64 {
    if x > P / 2 { -((P - x) as i64) } else { x as i64 }
}

/// Small deterministic PRNG (SplitMix64) so every test run is reproducible.
#[derive(Clone)]
pub struct Rng(pub u64);

impl Rng {
    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    pub fn elem(&mut self) -> u64 {
        self.next_u64() % P
    }
    pub fn below(&mut self, n: u64) -> u64 {
        self.next_u64() % n
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn field_laws() {
        let mut r = Rng(1);
        for _ in 0..1000 {
            let (a, b, c) = (r.elem(), r.elem(), r.elem());
            assert_eq!(mul(a, add(b, c)), add(mul(a, b), mul(a, c)));
            if a != 0 {
                assert_eq!(mul(a, inv(a)), 1);
            }
            assert_eq!(add(sub(a, b), b), a);
        }
        assert_eq!(to_i64(from_i64(-5)), -5);
        assert_eq!(reduce128((P as u128) * (P as u128)), 0);
    }
}
