//! Selected matrix-polynomial products (Definition 2.1).
//!
//! A shape is a triple of monomial supports (I, J, K) in `nvars` variables
//! with matrix dimensions [p, q, r]: left coefficients are p×q, right q×r,
//! output p×r, and output coefficient w is Σ_{u+v=w} A_u B_v (exponent
//! vectors add). Coordinates of the three coefficient lists:
//!
//! * left  (u, a, b) ↦ u·p·q + a·q + b
//! * right (v, b, c) ↦ v·q·r + b·r + c
//! * out   (w, a, c) ↦ w·p·r + a·r + c

use crate::laurent::Laurent;
use std::collections::HashMap;

pub type Mono = Vec<i32>;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Shape {
    pub nvars: usize,
    pub i: Vec<Mono>,
    pub j: Vec<Mono>,
    pub k: Vec<Mono>,
    pub p: usize,
    pub q: usize,
    pub r: usize,
}

/// The three roles of a selected product.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Role {
    Left,
    Right,
    Out,
}

impl Shape {
    /// Full_{a,b}[p,q,r]: all a+b−1 coefficients of a length-a times length-b product.
    pub fn full(a: usize, b: usize, p: usize, q: usize, r: usize) -> Shape {
        let s = a + b - 1;
        Shape {
            nvars: 1,
            i: (0..a as i32).map(|e| vec![e]).collect(),
            j: (0..b as i32).map(|e| vec![e]).collect(),
            k: (0..s as i32).map(|e| vec![e]).collect(),
            p,
            q,
            r,
        }
    }
    /// Mid_{a,b} := dual of Full_{a,b} (equation (2)); supports from the table in §2.1.
    pub fn mid(a: usize, b: usize, p: usize, q: usize, r: usize) -> Shape {
        let s = (a + b - 1) as i32;
        Shape {
            nvars: 1,
            i: (0..a as i32).map(|e| vec![e]).collect(),
            j: (1 - s..=0).map(|e| vec![e]).collect(),
            k: (1 - b as i32..=0).map(|e| vec![e]).collect(),
            p,
            q,
            r,
        }
    }
    /// Ordinary matrix multiplication [p,q,r] in zero variables.
    pub fn matmul(p: usize, q: usize, r: usize) -> Shape {
        Shape { nvars: 0, i: vec![vec![]], j: vec![vec![]], k: vec![vec![]], p, q, r }
    }
    pub fn with_dims(&self, p: usize, q: usize, r: usize) -> Shape {
        Shape { p, q, r, ..self.clone() }
    }
    pub fn size(&self, role: Role) -> usize {
        match role {
            Role::Left => self.i.len() * self.p * self.q,
            Role::Right => self.j.len() * self.q * self.r,
            Role::Out => self.k.len() * self.p * self.r,
        }
    }

    /// Evaluate the selected product.
    pub fn eval(&self, a: &[Laurent], b: &[Laurent]) -> Vec<Laurent> {
        assert_eq!(a.len(), self.size(Role::Left));
        assert_eq!(b.len(), self.size(Role::Right));
        let kidx: HashMap<&Mono, usize> = self.k.iter().enumerate().map(|(i, m)| (m, i)).collect();
        let (p, q, r) = (self.p, self.q, self.r);
        let mut out = vec![Laurent::zero(); self.size(Role::Out)];
        let mut sum = vec![0i32; self.nvars];
        for (ui, u) in self.i.iter().enumerate() {
            for (vi, v) in self.j.iter().enumerate() {
                for t in 0..self.nvars {
                    sum[t] = u[t] + v[t];
                }
                let Some(&wi) = kidx.get(&sum) else { continue };
                for x in 0..p {
                    for y in 0..q {
                        let av = &a[ui * p * q + x * q + y];
                        if av.is_zero() {
                            continue;
                        }
                        for z in 0..r {
                            let bv = &b[vi * q * r + y * r + z];
                            if bv.is_zero() {
                                continue;
                            }
                            out[wi * p * r + x * r + z].add_assign(&av.mul(bv));
                        }
                    }
                }
            }
        }
        out
    }

    /// Role exchange (12): Prod_{I,J,K}[p,q,r]^{12} = Prod_{J,I,K}[r,q,p] (BᵀAᵀ = (AB)ᵀ).
    /// Returns the new shape and, for each new role, the map new coordinate →
    /// old coordinate together with the old role it comes from.
    pub fn ex12(&self) -> (Shape, [(Role, Vec<usize>); 3]) {
        let (p, q, r) = (self.p, self.q, self.r);
        let n = Shape { nvars: self.nvars, i: self.j.clone(), j: self.i.clone(), k: self.k.clone(), p: r, q, r: p };
        // new left (v, c, b) [r×q] ← old right (v, b, c)
        let mut left = vec![0; n.size(Role::Left)];
        for v in 0..self.j.len() {
            for c in 0..r {
                for b in 0..q {
                    left[v * r * q + c * q + b] = v * q * r + b * r + c;
                }
            }
        }
        // new right (u, b, a) [q×p] ← old left (u, a, b)
        let mut right = vec![0; n.size(Role::Right)];
        for u in 0..self.i.len() {
            for b in 0..q {
                for a in 0..p {
                    right[u * q * p + b * p + a] = u * p * q + a * q + b;
                }
            }
        }
        // new out (w, c, a) [r×p] ← old out (w, a, c)
        let mut out = vec![0; n.size(Role::Out)];
        for w in 0..self.k.len() {
            for c in 0..r {
                for a in 0..p {
                    out[w * r * p + c * p + a] = w * p * r + a * r + c;
                }
            }
        }
        (n, [(Role::Right, left), (Role::Left, right), (Role::Out, out)])
    }

    /// Role exchange (23): Prod_{I,J,K}[p,q,r]^{23} = Prod_{I,K⁻¹,J⁻¹}[q,p,r] (AᵀZ).
    pub fn ex23(&self) -> (Shape, [(Role, Vec<usize>); 3]) {
        let (p, q, r) = (self.p, self.q, self.r);
        let inv = |s: &Vec<Mono>| s.iter().map(|m| m.iter().map(|e| -e).collect()).collect::<Vec<Mono>>();
        let n = Shape { nvars: self.nvars, i: self.i.clone(), j: inv(&self.k), k: inv(&self.j), p: q, q: p, r };
        // new left (u, b, a) [q×p] ← old left (u, a, b)
        let mut left = vec![0; n.size(Role::Left)];
        for u in 0..self.i.len() {
            for b in 0..q {
                for a in 0..p {
                    left[u * q * p + b * p + a] = u * p * q + a * q + b;
                }
            }
        }
        // new right (w⁻¹, a, c) [p×r] ← old out (w, a, c): identical layout
        let right = (0..n.size(Role::Right)).collect();
        // new out (v⁻¹, b, c) [q×r] ← old right (v, b, c): identical layout
        let out = (0..n.size(Role::Out)).collect();
        (n, [(Role::Left, left), (Role::Out, right), (Role::Right, out)])
    }

    /// Tensor product (Definition 2.6), variables of `o` appended after ours.
    /// Coordinates: support index (x, x') ↦ x·|S'| + x', matrix indices likewise.
    pub fn tensor(&self, o: &Shape) -> Shape {
        let prod = |a: &Vec<Mono>, b: &Vec<Mono>| {
            let mut v = Vec::with_capacity(a.len() * b.len());
            for x in a {
                for y in b {
                    let mut m = x.clone();
                    m.extend_from_slice(y);
                    v.push(m);
                }
            }
            v
        };
        Shape {
            nvars: self.nvars + o.nvars,
            i: prod(&self.i, &o.i),
            j: prod(&self.j, &o.j),
            k: prod(&self.k, &o.k),
            p: self.p * o.p,
            q: self.q * o.q,
            r: self.r * o.r,
        }
    }

    /// For the tensor shape `self ⊗ o`, map each coordinate of `role` to its
    /// index in the Kronecker layout (x·size' + x').
    pub fn tensor_to_kron(&self, o: &Shape, role: Role) -> Vec<usize> {
        let (n1, d1a, d1b, n2, d2a, d2b) = match role {
            Role::Left => (self.i.len(), self.p, self.q, o.i.len(), o.p, o.q),
            Role::Right => (self.j.len(), self.q, self.r, o.j.len(), o.q, o.r),
            Role::Out => (self.k.len(), self.p, self.r, o.k.len(), o.p, o.r),
        };
        let size2 = n2 * d2a * d2b;
        let (da, db) = (d1a * d2a, d1b * d2b);
        let mut m = vec![0; n1 * n2 * da * db];
        for s1 in 0..n1 {
            for s2 in 0..n2 {
                for a1 in 0..d1a {
                    for a2 in 0..d2a {
                        for b1 in 0..d1b {
                            for b2 in 0..d2b {
                                let t = (s1 * n2 + s2) * da * db + (a1 * d2a + a2) * db + (b1 * d2b + b2);
                                let x1 = s1 * d1a * d1b + a1 * d1b + b1;
                                let x2 = s2 * d2a * d2b + a2 * d2b + b2;
                                m[t] = x1 * size2 + x2;
                            }
                        }
                    }
                }
            }
        }
        m
    }

    /// Effort (4): Π_k min(a_k+b_k−1, a_k+h_k−1, b_k+h_k−1) · pqr.
    pub fn effort(&self) -> u128 {
        if self.i.is_empty() || self.j.is_empty() || self.k.is_empty() {
            return 0;
        }
        let range = |s: &Vec<Mono>, t: usize| {
            let lo = s.iter().map(|m| m[t]).min().unwrap();
            let hi = s.iter().map(|m| m[t]).max().unwrap();
            (hi - lo + 1) as u128
        };
        let mut e: u128 = (self.p * self.q * self.r) as u128;
        for t in 0..self.nvars {
            let (a, b, h) = (range(&self.i, t), range(&self.j, t), range(&self.k, t));
            e *= (a + b - 1).min(a + h - 1).min(b + h - 1);
        }
        e
    }

    /// Rename variables: new variable t is old variable perm[t].
    pub fn rename(&self, perm: &[usize]) -> Shape {
        let f = |s: &Vec<Mono>| s.iter().map(|m| perm.iter().map(|&t| m[t]).collect()).collect();
        Shape { nvars: self.nvars, i: f(&self.i), j: f(&self.j), k: f(&self.k), p: self.p, q: self.q, r: self.r }
    }

    /// Same supports as sets (orderings may differ) and same dims.
    pub fn same_as_sets(&self, o: &Shape) -> bool {
        let s = |v: &Vec<Mono>| {
            let mut v = v.clone();
            v.sort();
            v
        };
        self.nvars == o.nvars
            && (self.p, self.q, self.r) == (o.p, o.q, o.r)
            && s(&self.i) == s(&o.i)
            && s(&self.j) == s(&o.j)
            && s(&self.k) == s(&o.k)
    }

    /// Find a variable renaming π with self.rename(π) equal to `o` as sets
    /// (Lemma 2.10 isomorphism). Brute force over permutations; fine for ≤ 8 vars.
    pub fn find_isomorphism(&self, o: &Shape) -> Option<Vec<usize>> {
        if self.nvars != o.nvars {
            return None;
        }
        let mut perm: Vec<usize> = (0..self.nvars).collect();
        loop {
            if self.rename(&perm).same_as_sets(o) {
                return Some(perm);
            }
            if !next_permutation(&mut perm) {
                return None;
            }
        }
    }

    /// Coordinate map `role` of self → role of `o`, for shapes equal as sets:
    /// result[x_in_o] = matching coordinate in self.
    pub fn reorder_to(&self, o: &Shape, role: Role) -> Vec<usize> {
        let (s1, s2, da, db) = match role {
            Role::Left => (&self.i, &o.i, self.p, self.q),
            Role::Right => (&self.j, &o.j, self.q, self.r),
            Role::Out => (&self.k, &o.k, self.p, self.r),
        };
        let idx: HashMap<&Mono, usize> = s1.iter().enumerate().map(|(i, m)| (m, i)).collect();
        let mut m = vec![0; s2.len() * da * db];
        for (x, mono) in s2.iter().enumerate() {
            let y = idx[mono];
            for e in 0..da * db {
                m[x * da * db + e] = y * da * db + e;
            }
        }
        m
    }
}

pub fn next_permutation(v: &mut [usize]) -> bool {
    let n = v.len();
    if n < 2 {
        return false;
    }
    let mut i = n - 1;
    while i > 0 && v[i - 1] >= v[i] {
        i -= 1;
    }
    if i == 0 {
        return false;
    }
    let mut j = n - 1;
    while v[j] <= v[i - 1] {
        j -= 1;
    }
    v.swap(i - 1, j);
    v[i..].reverse();
    true
}
