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
///
/// Types follow the paper's labels, not the shapes: in the primal list the
/// last target is Mid and the rest are Full; the dual swaps every type. (When
/// a + b − 1 = 1 the shapes coincide but the labels still matter.)
pub fn amplify(pk: &Packing, d: usize) -> Packing {
    let m = pk.targets.len();
    let is_mid_primal = |i: usize| i == m - 1;
    let pair = pk.tensor(&pk.ex23());
    let mut acc = pair.clone();
    for _ in 1..d {
        acc = acc.tensor(&pair);
    }
    // target index of acc = base-(m²) digits (pair index = i·m + k, i primal, k dual)
    let keep: Vec<usize> = (0..acc.targets.len())
        .filter(|&idx| {
            let mut x = idx;
            let mut mids = 0;
            for _ in 0..d {
                let pi = x % (m * m);
                x /= m * m;
                let (i, k) = (pi / m, pi % m);
                mids += is_mid_primal(i) as usize; // primal letter of type Mid
                mids += (!is_mid_primal(k)) as usize; // dual of a Full is a Mid
            }
            mids == d
        })
        .collect();
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

// ---------------------------------------------------------------------------
// Elementary packings
// ---------------------------------------------------------------------------

/// Exact packing of `t` into `s` by zero padding, valid when every support
/// of `t` is contained in the corresponding support of `s` and the matrix
/// dims of `t` are ≤ those of `s` (inputs embedded, outputs selected).
pub fn embedding(t: &Shape, s: &Shape) -> Packing {
    assert_eq!(t.nvars, s.nvars);
    let emb = |ts: &Vec<crate::shape::Mono>, ss: &Vec<crate::shape::Mono>, ta: usize, tb: usize, sa: usize, sb: usize| {
        // returns for each source coordinate the target coordinate (or None)
        let idx: std::collections::HashMap<&crate::shape::Mono, usize> = ts.iter().enumerate().map(|(i, m)| (m, i)).collect();
        let mut src = vec![None; ss.len() * sa * sb];
        for (si, m) in ss.iter().enumerate() {
            if let Some(&ti) = idx.get(m) {
                for x in 0..ta {
                    for y in 0..tb {
                        src[si * sa * sb + x * sb + y] = Some(ti * ta * tb + x * tb + y);
                    }
                }
            }
        }
        for m in ts {
            assert!(ss.contains(m), "embedding: support not contained");
        }
        src
    };
    assert!(t.p <= s.p && t.q <= s.q && t.r <= s.r);
    let el = Sparse::select(&emb(&t.i, &s.i, t.p, t.q, s.p, s.q), t.size(Role::Left));
    let er = Sparse::select(&emb(&t.j, &s.j, t.q, t.r, s.q, s.r), t.size(Role::Right));
    let dsel = Sparse::select(&emb(&t.k, &s.k, t.p, t.r, s.p, s.r), t.size(Role::Out));
    let pk = Packing { targets: vec![t.clone()], sources: vec![s.clone()], same_left: false, el, er, d: dsel.transpose() };
    pk.check_dims();
    pk
}

/// Identity packing of one shape.
pub fn identity(t: &Shape) -> Packing {
    embedding(t, t)
}

/// Observation 2.5: append unused sources.
pub fn add_unused_sources(pk: &Packing, extra: &[Shape]) -> Packing {
    let mut out = pk.clone();
    for s in extra {
        out.sources.push(s.clone());
        out.el.nrows += s.size(Role::Left);
        out.el.rows.extend(vec![vec![]; s.size(Role::Left)]);
        out.er.nrows += s.size(Role::Right);
        out.er.rows.extend(vec![vec![]; s.size(Role::Right)]);
        out.d.ncols += s.size(Role::Out);
    }
    out.check_dims();
    out
}

/// Stack packings with the same sources' structure: the direct sum, used to
/// replace every source by a packing of it (compose with a block of inner packings).
pub fn direct_sum(parts: &[Packing]) -> Packing {
    let targets = parts.iter().flat_map(|p| p.targets.iter().cloned()).collect();
    let sources = parts.iter().flat_map(|p| p.sources.iter().cloned()).collect();
    for p in parts {
        assert!(!p.same_left);
    }
    let bd = |f: &dyn Fn(&Packing) -> Sparse| Sparse::block_diag(&parts.iter().map(f).collect::<Vec<_>>());
    let pk = Packing { targets, sources, same_left: false, el: bd(&|p| p.el.clone()), er: bd(&|p| p.er.clone()), d: bd(&|p| p.d.clone()) };
    pk.check_dims();
    pk
}

/// Replace each source of `pk` by its padded version `new[j]` (Observation 2.5: zero-pad source matrices).
pub fn pad_sources(pk: &Packing, new: &[Shape]) -> Packing {
    let parts: Vec<Packing> = pk.sources.iter().zip(new).map(|(s, n)| embedding(s, n)).collect();
    pk.compose(&direct_sum(&parts))
}

// ---------------------------------------------------------------------------
// Lemma 2.3: exact evaluation of a selected product by scalar multiplications
// ---------------------------------------------------------------------------

/// Exact packing of the univariate product with supports given as exponent
/// intervals (natural direction: evaluate both inputs at a+b−1 points,
/// interpolate the product, select K), into 1×1 matmul sources.
fn interp_natural(t: &Shape) -> Packing {
    assert_eq!(t.nvars, 1);
    assert_eq!((t.p, t.q, t.r), (1, 1, 1));
    let ilo = t.i.iter().map(|m| m[0]).min().unwrap();
    let ihi = t.i.iter().map(|m| m[0]).max().unwrap();
    let jlo = t.j.iter().map(|m| m[0]).min().unwrap();
    let jhi = t.j.iter().map(|m| m[0]).max().unwrap();
    let a = (ihi - ilo + 1) as usize;
    let b = (jhi - jlo + 1) as usize;
    let n = a + b - 1;
    let pts: Vec<u64> = (1..=n as u64).collect();
    let src = Shape::matmul(1, 1, 1);
    let mut el = Sparse::zeros(n, t.i.len());
    let mut er = Sparse::zeros(n, t.j.len());
    for (k, &x) in pts.iter().enumerate() {
        for (ui, u) in t.i.iter().enumerate() {
            el.add_entry(k, ui, &Laurent::constant(f::pow(x, (u[0] - ilo) as u64)));
        }
        for (vi, v) in t.j.iter().enumerate() {
            er.add_entry(k, vi, &Laurent::constant(f::pow(x, (v[0] - jlo) as u64)));
        }
    }
    // inverse Vandermonde: coefficient e of the product = Σ_k V^{-1}[e][k] y_k
    let vinv = vandermonde_inverse(&pts);
    let mut d = Sparse::zeros(t.k.len(), n);
    for (wi, w) in t.k.iter().enumerate() {
        let e = w[0] - ilo - jlo;
        if e < 0 || e >= n as i32 {
            continue;
        }
        for k in 0..n {
            d.add_entry(wi, k, &Laurent::constant(vinv[e as usize][k]));
        }
    }
    let pk = Packing { targets: vec![t.clone()], sources: vec![src; n], same_left: false, el, er, d };
    pk.check_dims();
    pk
}

/// Inverse of the Vandermonde matrix V[k][e] = x_k^e, over F_p (Gauss–Jordan).
pub fn vandermonde_inverse(xs: &[u64]) -> Vec<Vec<u64>> {
    let n = xs.len();
    let mut m: Vec<Vec<u64>> = (0..n).map(|k| (0..n).map(|e| f::pow(xs[k], e as u64)).collect()).collect();
    let mut inv: Vec<Vec<u64>> = (0..n).map(|i| (0..n).map(|j| (i == j) as u64).collect()).collect();
    for col in 0..n {
        let piv = (col..n).find(|&r| m[r][col] != 0).unwrap();
        m.swap(col, piv);
        inv.swap(col, piv);
        let iv = f::inv(m[col][col]);
        for j in 0..n {
            m[col][j] = f::mul(m[col][j], iv);
            inv[col][j] = f::mul(inv[col][j], iv);
        }
        for r in 0..n {
            if r != col && m[r][col] != 0 {
                let fct = m[r][col];
                for j in 0..n {
                    m[r][j] = f::sub(m[r][j], f::mul(fct, m[col][j]));
                    inv[r][j] = f::sub(inv[r][j], f::mul(fct, inv[col][j]));
                }
            }
        }
    }
    // V c = y  =>  c = V^{-1} y; we built V as rows k, cols e, so inv is (e, k).
    inv
}

/// Lemma 2.3 for one univariate scalar factor: pick the role exchange with
/// the smallest interpolation size, as in the proof.
fn interp_univariate(t: &Shape) -> Packing {
    let len = |s: &Vec<crate::shape::Mono>| {
        (s.iter().map(|m| m[0]).max().unwrap() - s.iter().map(|m| m[0]).min().unwrap() + 1) as usize
    };
    let (a, b, h) = (len(&t.i), len(&t.j), len(&t.k));
    let best = (a + b - 1).min(a + h - 1).min(b + h - 1);
    let pk = if a + b - 1 == best {
        interp_natural(t)
    } else if a + h - 1 == best {
        // (23): P^{23} has input lengths (a, h)
        let (t23, _) = t.ex23();
        interp_natural(&t23).ex23()
    } else {
        // (13): P^{13} has input lengths (h, b)
        let t13 = t.ex12().0.ex23().0.ex12().0;
        interp_natural(&t13).ex13()
    };
    // restore the target's support ordering
    let tt = pk.targets[0].clone();
    let n = pk.rename_target(0, t, &[0]);
    let _ = tt;
    // sources of exchanged scalar packings are 1×1 matmuls (supports {()})
    n
}

/// Schoolbook matmul [p,q,r] as an exact packing into pqr scalar products.
pub fn schoolbook(p: usize, q: usize, r: usize) -> Packing {
    let t = Shape::matmul(p, q, r);
    let n = p * q * r;
    let mut el = Sparse::zeros(n, p * q);
    let mut er = Sparse::zeros(n, q * r);
    let mut d = Sparse::zeros(p * r, n);
    let one = Laurent::constant(1);
    let mut s = 0;
    for a in 0..p {
        for b in 0..q {
            for c in 0..r {
                el.add_entry(s, a * q + b, &one);
                er.add_entry(s, b * r + c, &one);
                d.add_entry(a * r + c, s, &one);
                s += 1;
            }
        }
    }
    Packing { targets: vec![t], sources: vec![Shape::matmul(1, 1, 1); n], same_left: false, el, er, d }
}

/// Lemma 2.3: exact packing of a selected product whose supports are products
/// of univariate intervals, into eff(P) scalar multiplications.
pub fn evaluate(t: &Shape) -> Packing {
    // univariate factors
    let factor = |v: usize| -> Shape {
        let col = |s: &Vec<crate::shape::Mono>| {
            let mut e: Vec<i32> = s.iter().map(|m| m[v]).collect();
            e.sort();
            e.dedup();
            e.into_iter().map(|x| vec![x]).collect::<Vec<_>>()
        };
        Shape { nvars: 1, i: col(&t.i), j: col(&t.j), k: col(&t.k), p: 1, q: 1, r: 1 }
    };
    let mut pk = schoolbook(t.p, t.q, t.r);
    for v in 0..t.nvars {
        pk = interp_univariate(&factor(v)).tensor(&pk);
    }
    // the tensor's target is the product shape with variable order (v_{n-1}, ..., v_0, matrix)
    let tgt = pk.targets[0].clone();
    assert_eq!(tgt.nvars, t.nvars);
    let perm: Vec<usize> = (0..t.nvars).rev().collect();
    assert!(tgt.rename(&perm).same_as_sets(t), "Lemma 2.3 needs product supports");
    let mut out = pk.rename_target(0, t, &perm);
    // sources are 1×1 matmuls in t.nvars "variables" of exponent 0 — normalise to Shape::matmul(1,1,1)
    let m = Shape::matmul(1, 1, 1);
    for j in 0..out.sources.len() {
        out.sources[j] = m.clone();
    }
    out
}

// ---------------------------------------------------------------------------
// Lemma 4.1 / 4.2: the round
// ---------------------------------------------------------------------------

/// Q_{a,b}[p,q,r] := T_{a,b} ⊗ T_{a,b}^{12} ⊗ T_{a,b}^{13} ⊗ (matmul [p,q,r]) (Lemma 4.1).
pub fn q_shape(a: usize, b: usize, d: usize, p: usize, q: usize, r: usize) -> Shape {
    let t = t_shape(a, b, d);
    let t12 = t.ex12().0;
    let t13 = t.ex12().0.ex23().0.ex12().0;
    t.tensor(&t12).tensor(&t13).tensor(&Shape::matmul(p, q, r))
}

/// Rename every target and source of `pk` to the given canonical shapes
/// (each must be isomorphic; Lemma 2.10).
pub fn canonicalize(pk: &Packing, targets: &[Shape], sources: &[Shape]) -> Packing {
    let mut out = pk.clone();
    for (i, t) in targets.iter().enumerate() {
        if &out.targets[i] != t {
            let perm = out.targets[i].find_isomorphism(t).unwrap_or_else(|| panic!("target {i} not isomorphic to canonical shape"));
            out = out.rename_target(i, t, &perm);
        }
    }
    for (j, s) in sources.iter().enumerate() {
        if &out.sources[j] != s {
            let perm = out.sources[j].find_isomorphism(s).unwrap_or_else(|| panic!("source {j} not isomorphic to canonical shape"));
            out = out.rename_source(j, s, &perm);
        }
    }
    out
}

/// Parameters of one round: lengths L → L', amplification depth d, counts
/// w ≥ u, Λ = ⌊w/u⌋², side multiplier x.
#[derive(Clone, Copy, Debug)]
pub struct RoundParams {
    pub l: usize,
    pub l2: usize,
    pub d: usize,
    pub w: usize,
    pub u: usize,
    pub x: usize,
}

/// Lemma 4.2. `p1` packs w copies of T_{L,L}[1,1,w] using ≤ 2w copies of
/// T_{L,L'}[1,1,u]; `p2` the same for (a, b, b') = (L', L, L'). Returns a PP
/// of w³ copies of Q_{L,L}[w²Λx] using 64w³ copies of Q_{L',L'}[w²x].
pub fn round(p1: &Packing, p2: &Packing, rp: RoundParams) -> Packing {
    let RoundParams { l, l2, d, w, u, x } = rp;
    assert!(w >= u);
    let lam = (w / u) * (w / u);
    // tensor P ⊗ P^{12} ⊗ P^{13} ⊗ id(matmul[s,s,s]); canonical shapes; pad sources to 8w³
    let lift = |pk: &Packing, a: usize, b: usize, b2: usize, side: usize| -> Packing {
        assert!(pk.sources.len() <= 2 * w);
        let t = pk.tensor(&pk.ex12()).tensor(&pk.ex13()).tensor(&identity(&Shape::matmul(side, side, side)));
        let tgt = vec![q_shape(a, b, d, w * side, w * side, w * side); t.targets.len()];
        let src = vec![q_shape(a, b2, d, u * side, u * side, u * side); t.sources.len()];
        let t = canonicalize(&t, &tgt, &src);
        let extra = vec![src[0].clone(); 8 * w * w * w - t.sources.len()];
        add_unused_sources(&t, &extra)
    };
    let first = lift(p1, l, l, l2, w * lam * x); // Q_{L,L}[w²Λx] ← 8w³ Q_{L,L'}[wuΛx]
    // rename the sources Q_{L,L'} → Q_{L',L} (Lemma 4.1 symmetry)
    let mid = q_shape(l2, l, d, w * u * lam * x, w * u * lam * x, w * u * lam * x);
    let first = canonicalize(&first, &first.targets.clone(), &vec![mid.clone(); first.sources.len()]);
    let second = lift(p2, l2, l, l2, u * lam * x); // Q_{L',L}[wuΛx] ← 8w³ Q_{L',L'}[u²Λx]
    let second = second.repeat(8);
    let composed = first.compose(&second);
    // pad final sources to side w²x
    let side = w * w * x;
    let padded = vec![q_shape(l2, l2, d, side, side, side); composed.sources.len()];
    pad_sources(&composed, &padded)
}

// ---------------------------------------------------------------------------
// Lemma 2.4: Bini's approximate-to-exact conversion
// ---------------------------------------------------------------------------

/// Turn a packing modulo ε whose output has ε-degree ≤ `deg` into an exact
/// packing: evaluate the maps at deg + 1 nonzero values of ε and interpolate
/// the constant coefficient. Uses (deg + 1)·#sources sources.
pub fn bini(pk: &Packing, deg: usize) -> Packing {
    let pts: Vec<u64> = (1..=(deg + 1) as u64).collect();
    let w0 = lagrange_at(&pts, 0);
    let mut sources = vec![];
    let mut els = vec![];
    let mut ers = vec![];
    let mut ds = vec![];
    for (l, &e) in pts.iter().enumerate() {
        sources.extend(pk.sources.iter().cloned());
        els.push(pk.el.at_eps(e));
        ers.push(pk.er.at_eps(e));
        ds.push(pk.d.at_eps(e).scale(&Laurent::constant(w0[l])));
    }
    // decoder: horizontal concatenation
    let mut d = Sparse::zeros(pk.d.nrows, pk.d.ncols * pts.len());
    for (l, dl) in ds.iter().enumerate() {
        for (i, row) in dl.rows.iter().enumerate() {
            for (j, x) in row {
                d.rows[i].push((l * pk.d.ncols + j, x.clone()));
            }
        }
    }
    let out = Packing { targets: pk.targets.clone(), sources, same_left: pk.same_left, el: vstack(&els), er: vstack(&ers), d };
    out.check_dims();
    out
}

/// Upper bound on the ε-degree of a packing's output: EL.hi + ER.hi + D.hi.
pub fn degree_bound(pk: &Packing) -> usize {
    pk.profile().iter().map(|&(_, hi)| hi.max(0) as usize).sum()
}

// ---------------------------------------------------------------------------
// Theorem 4.3 chain (executable at minimal retention)
// ---------------------------------------------------------------------------

/// One round's input packings at retention w = u = 1: the exact embedding of
/// T_{a,b} into T_{a,b'} (Observation 2.5), which satisfies Lemma 4.2's
/// hypothesis "w copies of T_{a,b}[1,1,w] using ≤ 2w copies of T_{a,b'}[1,1,u]".
pub fn minimal_round_input(a: usize, b: usize, b2: usize, d: usize) -> Packing {
    embedding(&t_shape(a, b, d), &t_shape(a, b2, d))
}

/// L_{j+1} = L_j + ⌈(3L_j − 1)/(n − 1)⌉ (proof of Theorem 4.3).
pub fn next_length(l: usize, n: usize) -> usize {
    l + (3 * l - 1).div_ceil(n - 1)
}

/// The chain of Theorem 4.3 with k rounds at retention w = u = 1 (so Λ = 1,
/// v_k = 1, and x is the matrix side), followed by Lemma 2.3 evaluation of
/// every source: an exact bilinear algorithm for x×x matrix multiplication.
/// `perturb` uses ε-perturbed inputs instead, giving a packing modulo ε.
pub fn minimal_chain(c: usize, d: usize, k: usize, x: usize, input: &dyn Fn(usize, usize, usize, usize) -> Packing) -> Packing {
    let mut ls = vec![1usize];
    for _ in 0..k {
        let l = *ls.last().unwrap();
        ls.push(next_length(l, c));
    }
    let mut chain: Option<Packing> = None;
    for j in 0..k {
        let (l, l2) = (ls[j], ls[j + 1]);
        let r = round(&input(l, l, l2, d), &input(l2, l, l2, d), RoundParams { l, l2, d, w: 1, u: 1, x });
        chain = Some(match chain {
            None => r,
            Some(ch) => {
                let n = ch.sources.len() / r.targets.len();
                ch.compose(&r.repeat(n))
            }
        });
    }
    let ch = chain.unwrap();
    let ev = evaluate(&ch.sources[0]);
    let n = ch.sources.len();
    ch.compose(&ev.repeat(n))
}
