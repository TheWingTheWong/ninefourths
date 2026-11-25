//! The constructions of Sections 3 and 4, built as explicit packings.

use crate::field as f;
use crate::laurent::Laurent;
use crate::packing::{lagrange_at, offsets, Packing};
use crate::shape::{Role, Shape};
use crate::sparse::Sparse;

/// Theorem 3.1 (primal packing): an SL-PP modulo ε of c+1 copies of
/// Full_{a,b}[p,q,r] and one Mid_{a,b}[p,q,r] using c copies of
/// Full_{a,b+t}[p,q,r], valid when ct ≥ a + 2b − 1.
///
/// Targets are ordered Full_1..Full_{c+1}, Mid. The construction is built
/// whether or not the condition holds, so tests can show it fails without it.
pub fn primal(a: usize, b: usize, c: usize, t: usize, p: usize, q: usize, r: usize) -> Packing {
    let s = a + b - 1;
    let full = Shape::full(a, b, p, q, r);
    let mid = Shape::mid(a, b, p, q, r);
    let src = Shape::full(a, b + t, p, q, r);
    let mut targets = vec![full.clone(); c + 1];
    targets.push(mid);
    let sources = vec![src; c];

    // Scalar (1×1) maps on polynomial coefficient indices, later tensored with
    // the identity on matrix entries.
    // Target right space (scalar): P_1..P_{c+1} (b coeffs each), then Q (s coeffs,
    // exponents 1−s..0 in that order, so Q̃ = z^{s−1}Q has coefficient e at index e).
    let tr_n = (c + 1) * b + s;
    let p_idx = |j: usize, e: usize| (j - 1) * b + e; // j in 1..=c+1
    let q_idx = |e: usize| (c + 1) * b + e; // Q̃ coefficient e
    let sr_n = c * (b + t);
    let mut er = Sparse::zeros(sr_n, tr_n);
    let one = Laurent::constant(1);
    let minus = Laurent::int(-1);
    let eps = Laurent::mono(1, 1);
    for j in 1..=c {
        let base = (j - 1) * (b + t);
        for e in 0..b {
            er.add_entry(base + e, p_idx(j, e), &one); // P_j
            er.add_entry(base + t + e, p_idx(j + 1, e), &minus); // − z^t P_{j+1}
        }
        for v in 0..t {
            let qe = t * (j - 1) + v; // [z^{t(j−1)+v}] Q̃
            if qe < s {
                er.add_entry(base + b + v, q_idx(qe), &eps); // + ε z^b Q^{(j)}
            }
        }
    }
    // Decoder: source outputs W_j (s+t coeffs each) → targets.
    let so_n = c * (s + t);
    let w_idx = |j: usize, e: usize| (j - 1) * (s + t) + e;
    // W = Σ_j z^{t(j−1)} W_j as linear forms over source outputs.
    let wlen = t * (c - 1) + s + t;
    let mut wform: Vec<Vec<(usize, Laurent)>> = vec![vec![]; wlen];
    for j in 1..=c {
        for e in 0..s + t {
            wform[t * (j - 1) + e].push((w_idx(j, e), one.clone()));
        }
    }
    let to_n = (c + 1) * s + b;
    let mut d = Sparse::zeros(to_n, so_n);
    // Z_1 = first s coefficients of W; Z_{j+1}[k] = (Z_j − W_j)[k + t].
    let mut z: Vec<Vec<(usize, Laurent)>> = (0..s).map(|k| wform[k].clone()).collect();
    for j in 1..=c + 1 {
        for (k, form) in z.iter().enumerate() {
            for (col, x) in form {
                d.add_entry((j - 1) * s + k, *col, x);
            }
        }
        if j == c + 1 {
            break;
        }
        let mut next = vec![];
        for k in 0..s {
            let e = k + t;
            let mut form = if e < s { z[e].clone() } else { vec![] };
            if e < s + t {
                form.push((w_idx(j, e), minus.clone()));
            }
            next.push(form);
        }
        z = next;
    }
    // Middle product: ε⁻¹ [z^{s+k}] W for 0 ≤ k < b.
    let einv = Laurent::mono(1, -1);
    for k in 0..b {
        if s + k < wlen {
            for (col, x) in &wform[s + k] {
                d.add_entry((c + 1) * s + k, *col, &x.mul(&einv));
            }
        }
    }
    // Lift to matrix coefficients.
    let el = Sparse::block_diag(&vec![Sparse::identity(a); c]).compose(&copy_map(a, c)).kron(&Sparse::identity(p * q));
    let er = er.kron(&Sparse::identity(q * r));
    let d = d.kron(&Sparse::identity(p * r));
    let pk = Packing { targets, sources, same_left: true, el, er, d };
    pk.check_dims();
    pk
}

/// The SL copying map: one left input of length a copied to `n` sources.
fn copy_map(a: usize, n: usize) -> Sparse {
    let mut m = Sparse::zeros(a * n, a);
    for j in 0..n {
        for e in 0..a {
            m.rows[j * a + e].push((e, Laurent::constant(1)));
        }
    }
    m
}

/// Is `s` (a tensor of univariate Full/Mid factors) balanced: at most `d`
/// definite Mid factors and at most `d` definite Full factors? A factor in
/// variable t is Mid iff its right support has a negative exponent, Full iff its
/// output support has a positive one; when a + b − 1 = 1 the two shapes coincide and the factor counts
/// as either.
pub fn balanced(s: &Shape, d: usize) -> bool {
    let mid = (0..s.nvars).filter(|&t| s.j.iter().any(|m| m[t] < 0)).count();
    let full = (0..s.nvars).filter(|&t| s.k.iter().any(|m| m[t] > 0)).count();
    mid <= d && full <= d
}

/// Lemma 3.2 (same-left amplification): tensor the SL-PP with its dual, take
/// d copies, and retain the targets with exactly d Full and d Mid factors
/// (Claim 3.3 guarantees at least m^{2d}/(d+1) of them).
pub fn amplify(pk: &Packing, d: usize) -> Packing {
    let pair = pk.tensor(&pk.ex23());
    let mut acc = pair.clone();
    for _ in 1..d {
        acc = acc.tensor(&pair);
    }
    let keep: Vec<usize> = (0..acc.targets.len()).filter(|&i| balanced(&acc.targets[i], d)).collect();
    acc.retain(&keep)
}

/// T_{a,b} := (Full_{a,b} ⊗ Mid_{a,b})^{⊗d} (equation (7)) at dims [1,1,1].
pub fn t_shape(a: usize, b: usize, d: usize) -> Shape {
    let pair = Shape::full(a, b, 1, 1, 1).tensor(&Shape::mid(a, b, 1, 1, 1));
    let mut s = Shape::matmul(1, 1, 1);
    for _ in 0..d {
        s = s.tensor(&pair);
    }
    s
}

/// Lemma 3.4 (from SL-PP to PP): an SL-PP of w targets F_i[p,q,r] using u
/// copies of S[p,q,r] becomes a PP of F_i[p,q,wr] (independent left inputs)
/// using 2w − 1 copies of S[p,q,ur].
///
/// Interpolation nodes: A(i) = A_i for i = 1..w, B'_i(i + εj) = B_{i,j},
/// evaluation points k = 0..2w−2 (all in F_p).
pub fn sl_to_pp(sl: &Packing) -> Packing {
    assert!(sl.same_left);
    let w = sl.targets.len();
    let u = sl.sources.len();
    let s0 = &sl.sources[0];
    for s in &sl.sources {
        assert_eq!(s, s0, "Lemma 3.4 needs identical sources");
    }
    let (p, q, r) = (s0.p, s0.q, s0.r);
    for t in &sl.targets {
        assert_eq!((t.p, t.q, t.r), (p, q, r));
    }
    let eta = sl.subst(w as i64); // parameter η = ε^w
    let npts = 2 * w - 1;
    let pts: Vec<u64> = (0..npts as u64).collect();
    let new_targets: Vec<Shape> = sl.targets.iter().map(|t| t.with_dims(p, q, w * r)).collect();
    let big_src = s0.with_dims(p, q, u * r);
    let new_sources = vec![big_src.clone(); npts];

    // --- left encoder: A(k) = Σ_i A_i ℓ_i(k), nodes 1..w
    let nodes_a: Vec<u64> = (1..=w as u64).collect();
    let la = sl.targets[0].size(Role::Left);
    let mut el = Sparse::zeros(npts * la, w * la);
    for (kk, &k) in pts.iter().enumerate() {
        let coef = lagrange_at(&nodes_a, k);
        for (i, &cf) in coef.iter().enumerate() {
            for x in 0..la {
                el.add_entry(kk * la + x, i * la + x, &Laurent::constant(cf));
            }
        }
    }

    // --- right encoder.
    // Target i right input (coords of F_i[p,q,wr]): (v, row, col) with col = j·r + c'.
    // B'_i(k) = Σ_j λ_{ij}(k) B_{i,j}, λ_{ij}(k) = Π_{j'≠j} (k − i − εj') / (ε(j − j')).
    let lambda = |i: usize, j: usize, k: u64| -> Laurent {
        let mut num = Laurent::constant(1);
        let mut den = 1u64;
        for jp in 0..w {
            if jp == j {
                continue;
            }
            // (k − i) − ε j'
            let lin = Laurent::constant(f::sub(k, i as u64)).add(&Laurent::mono(f::neg(jp as u64), 1));
            num = num.mul(&lin);
            den = f::mul(den, f::from_i64(j as i64 - jp as i64));
        }
        num.mul(&Laurent::mono(f::inv(den), -((w - 1) as i64)))
    };
    let t_off = offsets(&new_targets, Role::Right);
    let ntr = *t_off.last().unwrap();
    // inner SL right input space: concat over targets of F_i[p,q,r] right
    let in_off = offsets(&sl.targets, Role::Right);
    let nin = *in_off.last().unwrap();
    let src_r = s0.size(Role::Right); // per inner source (q×r coefficients)
    let big_r = big_src.size(Role::Right);
    // M_k: new target right → inner target right (B'_i(k)).
    // Then eta.er maps inner target right → u inner source rights; then concat
    // columns into big source right.
    // concat map C: inner source rights (u blocks of S[p,q,r] right) → S[p,q,ur] right.
    let mut concat_r = Sparse::zeros(big_r, u * src_r);
    for nu in 0..u {
        for v in 0..s0.j.len() {
            for row in 0..q {
                for col in 0..r {
                    let from = nu * src_r + v * q * r + row * r + col;
                    let to = v * q * (u * r) + row * (u * r) + nu * r + col;
                    concat_r.add_entry(to, from, &Laurent::constant(1));
                }
            }
        }
    }
    let mut er_blocks = vec![];
    for &k in &pts {
        let mut mk = Sparse::zeros(nin, ntr);
        for (i, t) in sl.targets.iter().enumerate() {
            for j in 0..w {
                let lam = lambda(i + 1, j, k);
                for v in 0..t.j.len() {
                    for row in 0..q {
                        for col in 0..r {
                            let from = t_off[i] + v * q * (w * r) + row * (w * r) + j * r + col;
                            let to = in_off[i] + v * q * r + row * r + col;
                            mk.add_entry(to, from, &lam);
                        }
                    }
                }
            }
        }
        er_blocks.push(concat_r.compose(&eta.er).compose(&mk));
    }
    let er = vstack(&er_blocks);

    // --- decoder. Source k output: S[p,q,ur] out (w', a, col) with col = ν r + c'.
    // split into u inner source outputs, apply eta.d → inner target outputs P_i(k).
    let src_o = s0.size(Role::Out);
    let big_o = big_src.size(Role::Out);
    let mut split_o = Sparse::zeros(u * src_o, big_o);
    for nu in 0..u {
        for wv in 0..s0.k.len() {
            for a in 0..p {
                for col in 0..r {
                    let to = nu * src_o + wv * p * r + a * r + col;
                    let from = wv * p * (u * r) + a * (u * r) + nu * r + col;
                    split_o.add_entry(to, from, &Laurent::constant(1));
                }
            }
        }
    }
    let dec_k = eta.d.compose(&split_o); // inner target outs ← one big source out
    let inner_to = offsets(&sl.targets, Role::Out);
    let new_to = offsets(&new_targets, Role::Out);
    let mut d = Sparse::zeros(*new_to.last().unwrap(), npts * big_o);
    // P̃_i(ξ) interpolated from values at pts, evaluated at ξ = i + εj:
    // μ_k(i + εj) = Π_{k'≠k} (i + εj − k') / (k − k').
    for (i, t) in sl.targets.iter().enumerate() {
        for j in 0..w {
            for (kk, &k) in pts.iter().enumerate() {
                let mut mu = Laurent::constant(1);
                let mut den = 1u64;
                for &kp in &pts {
                    if kp == k {
                        continue;
                    }
                    let lin = Laurent::constant(f::sub((i + 1) as u64, kp)).add(&Laurent::mono(j as u64 % f::P, 1));
                    mu = mu.mul(&lin);
                    den = f::mul(den, f::sub(k, kp));
                }
                let mu = mu.scale(f::inv(den));
                for wv in 0..t.k.len() {
                    for a in 0..p {
                        for col in 0..r {
                            let row_new = new_to[i] + wv * p * (w * r) + a * (w * r) + j * r + col;
                            let row_inner = inner_to[i] + wv * p * r + a * r + col;
                            for (c2, x) in &dec_k.rows[row_inner] {
                                d.add_entry(row_new, kk * big_o + c2, &x.mul(&mu));
                            }
                        }
                    }
                }
            }
        }
    }
    let pk = Packing { targets: new_targets, sources: new_sources, same_left: false, el, er, d };
    pk.check_dims();
    pk
}

/// Stack maps vertically (same number of columns).
pub fn vstack(blocks: &[Sparse]) -> Sparse {
    let ncols = blocks[0].ncols;
    let mut out = Sparse::zeros(0, ncols);
    for b in blocks {
        assert_eq!(b.ncols, ncols);
        out.rows.extend(b.rows.iter().cloned());
        out.nrows += b.nrows;
    }
    out
}
