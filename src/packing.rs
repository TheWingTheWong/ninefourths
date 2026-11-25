//! Polynomial packings (Definition 2.2) and the generic operations on them:
//! role exchange (Lemma 2.8), tensoring (Lemma 2.7), composition (Lemma 2.9),
//! elementary operations (Observation 2.5) and renaming (Lemma 2.10), plus an
//! exact verifier.

use crate::field::{self as f, Rng};
use crate::laurent::Laurent;
use crate::shape::{Role, Shape};
use crate::sparse::Sparse;

/// A packing of `targets` using `sources`.
///
/// Coefficient spaces are concatenations over the list, in order. For a
/// same-left packing (SL-PP) all targets share one left input, so the target
/// left space is a single copy of `targets[0]`'s left space.
#[derive(Clone, Debug)]
pub struct Packing {
    pub targets: Vec<Shape>,
    pub sources: Vec<Shape>,
    pub same_left: bool,
    /// Target left → source left.
    pub el: Sparse,
    /// Target right → source right.
    pub er: Sparse,
    /// Source out → target out.
    pub d: Sparse,
}

pub fn offsets(shapes: &[Shape], role: Role) -> Vec<usize> {
    let mut o = vec![0];
    for s in shapes {
        o.push(o.last().unwrap() + s.size(role));
    }
    o
}

/// Result of [`Packing::verify`].
#[derive(Debug, Clone)]
pub struct Verdict {
    /// Highest ε power seen in the decoded output (0 for exact packings).
    pub out_degree: i64,
    /// Lowest ε power seen in the decoded output (must be ≥ 0).
    pub out_valuation: i64,
}

impl Packing {
    pub fn target_space(&self, role: Role) -> usize {
        if role == Role::Left && self.same_left {
            self.targets[0].size(Role::Left)
        } else {
            *offsets(&self.targets, role).last().unwrap()
        }
    }
    pub fn source_space(&self, role: Role) -> usize {
        *offsets(&self.sources, role).last().unwrap()
    }
    /// Dimension sanity checks; panics with a message on mismatch.
    pub fn check_dims(&self) {
        assert_eq!(self.el.ncols, self.target_space(Role::Left), "EL cols");
        assert_eq!(self.el.nrows, self.source_space(Role::Left), "EL rows");
        assert_eq!(self.er.ncols, self.target_space(Role::Right), "ER cols");
        assert_eq!(self.er.nrows, self.source_space(Role::Right), "ER rows");
        assert_eq!(self.d.ncols, self.source_space(Role::Out), "D cols");
        assert_eq!(self.d.nrows, self.target_space(Role::Out), "D rows");
        if self.same_left {
            let t0 = &self.targets[0];
            for t in &self.targets {
                assert_eq!((&t.i, t.p, t.q), (&t0.i, t0.p, t0.q), "SL targets must share the left space");
            }
        }
    }

    /// Run the packing as a program: encode, multiply sources, decode.
    pub fn run(&self, left: &[Laurent], right: &[Laurent]) -> Vec<Laurent> {
        let sl = self.el.apply(left);
        let sr = self.er.apply(right);
        let lo = offsets(&self.sources, Role::Left);
        let ro = offsets(&self.sources, Role::Right);
        let mut so = Vec::with_capacity(self.source_space(Role::Out));
        for (j, s) in self.sources.iter().enumerate() {
            so.extend(s.eval(&sl[lo[j]..lo[j + 1]], &sr[ro[j]..ro[j + 1]]));
        }
        self.d.apply(&so)
    }

    /// The targets' true outputs.
    pub fn expected(&self, left: &[Laurent], right: &[Laurent]) -> Vec<Laurent> {
        let lo = offsets(&self.targets, Role::Left);
        let ro = offsets(&self.targets, Role::Right);
        let mut out = Vec::new();
        for (i, t) in self.targets.iter().enumerate() {
            let a = if self.same_left { left } else { &left[lo[i]..lo[i + 1]] };
            out.extend(t.eval(a, &right[ro[i]..ro[i + 1]]));
        }
        out
    }

    /// Exact check on random inputs from F_p, with ε symbolic: the decoded
    /// output minus the true output must be divisible by ε (and the decoded
    /// output must have no negative powers of ε).
    pub fn verify(&self, rng: &mut Rng, trials: usize) -> Result<Verdict, String> {
        self.check_dims();
        let mut v = Verdict { out_degree: 0, out_valuation: 0 };
        for _ in 0..trials {
            let left: Vec<Laurent> = (0..self.target_space(Role::Left)).map(|_| Laurent::constant(rng.elem())).collect();
            let right: Vec<Laurent> = (0..self.target_space(Role::Right)).map(|_| Laurent::constant(rng.elem())).collect();
            let got = self.run(&left, &right);
            let want = self.expected(&left, &right);
            for (idx, (g, w)) in got.iter().zip(&want).enumerate() {
                let diff = g.sub(w);
                if let Some(val) = diff.valuation() {
                    if val <= 0 {
                        return Err(format!("output coordinate {idx}: error has ε^{val} term"));
                    }
                }
                if let Some(d) = g.degree() {
                    v.out_degree = v.out_degree.max(d);
                }
                if let Some(l) = g.valuation() {
                    v.out_valuation = v.out_valuation.min(l);
                }
            }
        }
        Ok(v)
    }

    /// ε exponent ranges of (EL, ER, D).
    pub fn profile(&self) -> [(i64, i64); 3] {
        [self.el.eps_range(), self.er.eps_range(), self.d.eps_range()]
    }
    pub fn is_exact(&self) -> bool {
        self.profile().iter().all(|&(lo, hi)| lo == 0 && hi == 0)
    }
    /// N = pL + pR + pD + 1 from the proof of Lemma 2.9.
    pub fn clearing_power(&self) -> i64 {
        self.profile().iter().map(|&(lo, _)| (-lo).max(0)).sum::<i64>() + 1
    }
    pub fn subst(&self, n: i64) -> Packing {
        Packing { el: self.el.subst(n), er: self.er.subst(n), d: self.d.subst(n), ..self.clone() }
    }

    /// The target→source map of a role (D is transposed).
    fn forward(&self, role: Role) -> Sparse {
        match role {
            Role::Left => self.el.clone(),
            Role::Right => self.er.clone(),
            Role::Out => self.d.transpose(),
        }
    }

    fn exchange(&self, per_shape: impl Fn(&Shape) -> (Shape, [(Role, Vec<usize>); 3]), keeps_sl: bool) -> Packing {
        assert!(!self.same_left || keeps_sl, "this exchange does not preserve same-left packings");
        let tex: Vec<_> = self.targets.iter().map(&per_shape).collect();
        let sex: Vec<_> = self.sources.iter().map(&per_shape).collect();
        let new_t: Vec<Shape> = tex.iter().map(|x| x.0.clone()).collect();
        let new_s: Vec<Shape> = sex.iter().map(|x| x.0.clone()).collect();
        let roles = [Role::Left, Role::Right, Role::Out];
        let mut maps = vec![];
        for (ri, &new_role) in roles.iter().enumerate() {
            let old_role = tex.first().map(|x| x.1[ri].0).unwrap_or(Role::Left);
            // selection matrices: new coords ← old concatenated coords
            let sel = |ex: &Vec<(Shape, [(Role, Vec<usize>); 3])>, old: &[Shape], shared: bool| {
                let oo = offsets(old, old_role);
                let mut src = vec![];
                if shared {
                    for &x in &ex[0].1[ri].1 {
                        src.push(Some(x));
                    }
                    return Sparse::select(&src, old[0].size(old_role));
                }
                for (n, e) in ex.iter().enumerate() {
                    for &x in &e.1[ri].1 {
                        src.push(Some(oo[n] + x));
                    }
                }
                Sparse::select(&src, *oo.last().unwrap())
            };
            let t_shared = self.same_left && old_role == Role::Left;
            let st = sel(&tex, &self.targets, t_shared);
            let ss = sel(&sex, &self.sources, false);
            let fwd = self.forward(old_role);
            let m = ss.compose(&fwd).compose(&st.transpose());
            let _ = new_role;
            maps.push(m);
        }
        let d = maps.pop().unwrap().transpose();
        let er = maps.pop().unwrap();
        let el = maps.pop().unwrap();
        let p = Packing { targets: new_t, sources: new_s, same_left: self.same_left, el, er, d };
        p.check_dims();
        p
    }

    /// Lemma 2.8 with σ = (12). Not defined for SL-PPs.
    pub fn ex12(&self) -> Packing {
        self.exchange(|s| s.ex12(), false)
    }
    /// Lemma 2.8 with σ = (23): the dual packing; preserves SL-PPs.
    pub fn ex23(&self) -> Packing {
        self.exchange(|s| s.ex23(), true)
    }
    /// σ = (13) = (12)(23)(12).
    pub fn ex13(&self) -> Packing {
        self.ex12().ex23().ex12()
    }

    /// Lemma 2.7: tensor product, targets ordered (i, k), sources (j, l).
    pub fn tensor(&self, o: &Packing) -> Packing {
        assert_eq!(self.same_left, o.same_left, "tensor of SL with non-SL packing");
        let pair = |a: &[Shape], b: &[Shape]| {
            let mut v = vec![];
            for x in a {
                for y in b {
                    v.push(x.tensor(y));
                }
            }
            v
        };
        let new_t = pair(&self.targets, &o.targets);
        let new_s = pair(&self.sources, &o.sources);
        // selection new concat coords -> Kronecker coords of (concat ⊗ concat')
        let sel = |a: &[Shape], b: &[Shape], role: Role, shared: bool| -> Sparse {
            if shared {
                let m = a[0].tensor_to_kron(&b[0], role);
                let n = a[0].size(role) * b[0].size(role);
                return Sparse::select(&m.into_iter().map(Some).collect::<Vec<_>>(), n);
            }
            let oa = offsets(a, role);
            let ob = offsets(b, role);
            let n2 = *ob.last().unwrap();
            let mut src = vec![];
            for (ia, x) in a.iter().enumerate() {
                for (ib, y) in b.iter().enumerate() {
                    let size2 = y.size(role);
                    for k in x.tensor_to_kron(y, role) {
                        let (x1, x2) = (k / size2, k % size2);
                        src.push(Some((oa[ia] + x1) * n2 + ob[ib] + x2));
                    }
                }
            }
            Sparse::select(&src, *oa.last().unwrap() * n2)
        };
        let mut maps = vec![];
        for role in [Role::Left, Role::Right, Role::Out] {
            let shared = self.same_left && role == Role::Left;
            let st = sel(&self.targets, &o.targets, role, shared);
            let ss = sel(&self.sources, &o.sources, role, false);
            let k = self.forward(role).kron(&o.forward(role));
            maps.push(ss.compose(&k).compose(&st.transpose()));
        }
        let d = maps.pop().unwrap().transpose();
        let er = maps.pop().unwrap();
        let el = maps.pop().unwrap();
        let p = Packing { targets: new_t, sources: new_s, same_left: self.same_left, el, er, d };
        p.check_dims();
        p
    }

    /// Lemma 2.9: `self` packs T using S; `inner` packs S using R. The inner
    /// parameter is replaced by ε^N, N = pL + pR + pD + 1 of the outer maps.
    pub fn compose(&self, inner: &Packing) -> Packing {
        assert!(!self.same_left && !inner.same_left, "compose expects general packings");
        assert_eq!(self.sources, inner.targets, "compose: outer sources must equal inner targets");
        let inner = if inner.is_exact() || self.is_exact() { inner.clone() } else { inner.subst(self.clearing_power()) };
        let p = Packing {
            targets: self.targets.clone(),
            sources: inner.sources.clone(),
            same_left: false,
            el: inner.el.compose(&self.el),
            er: inner.er.compose(&self.er),
            d: self.d.compose(&inner.d),
        };
        p.check_dims();
        p
    }

    /// Observation 2.5: keep only the listed targets (in that order).
    pub fn retain(&self, keep: &[usize]) -> Packing {
        let embed = |role: Role| -> Sparse {
            // old concat coords <- new concat coords (columns = new)
            let o = offsets(&self.targets, role);
            let mut src = vec![None; *o.last().unwrap()];
            let mut n = 0;
            for &i in keep {
                for x in o[i]..o[i + 1] {
                    src[x] = Some(n);
                    n += 1;
                }
            }
            Sparse::select(&src, n)
        };
        let el = if self.same_left { self.el.clone() } else { self.el.compose(&embed(Role::Left)) };
        let er = self.er.compose(&embed(Role::Right));
        let d = embed(Role::Out).transpose().compose(&self.d);
        let p = Packing {
            targets: keep.iter().map(|&i| self.targets[i].clone()).collect(),
            sources: self.sources.clone(),
            same_left: self.same_left,
            el,
            er,
            d,
        };
        p.check_dims();
        p
    }

    /// Observation 2.5: k independent copies.
    pub fn repeat(&self, k: usize) -> Packing {
        assert!(!self.same_left);
        let rep = |s: &[Shape]| (0..k).flat_map(|_| s.iter().cloned()).collect::<Vec<_>>();
        let bd = |m: &Sparse| Sparse::block_diag(&vec![m.clone(); k]);
        Packing {
            targets: rep(&self.targets),
            sources: rep(&self.sources),
            same_left: false,
            el: bd(&self.el),
            er: bd(&self.er),
            d: bd(&self.d),
        }
    }

    /// Lemma 2.10: replace target `i` by `new`, which must equal
    /// `targets[i].rename(perm)` as sets.
    pub fn rename_target(&self, i: usize, new: &Shape, perm: &[usize]) -> Packing {
        let old = &self.targets[i];
        let ren = old.rename(perm);
        assert!(ren.same_as_sets(new), "rename_target: not isomorphic via perm");
        let mut out = self.clone();
        out.targets[i] = new.clone();
        // coordinate maps between old (= ren layout) and new
        let reorder = |role: Role| ren.reorder_to(new, role); // new coord -> old local coord
        let o_r = offsets(&self.targets, Role::Right);
        let o_o = offsets(&self.targets, Role::Out);
        let perm_space = |o: &Vec<usize>, total: usize, role: Role| -> Sparse {
            // P: old concat <- new concat (identity outside block i)
            let mut src: Vec<Option<usize>> = (0..total).map(Some).collect();
            let m = reorder(role);
            for (xn, &xo) in m.iter().enumerate() {
                src[o[i] + xo] = Some(o[i] + xn);
            }
            Sparse::select(&src, total)
        };
        if !self.same_left {
            let o_l = offsets(&self.targets, Role::Left);
            out.el = self.el.compose(&perm_space(&o_l, *o_l.last().unwrap(), Role::Left));
        } else {
            assert!(reorder(Role::Left).iter().enumerate().all(|(a, &b)| a == b), "SL rename must keep left layout");
        }
        out.er = self.er.compose(&perm_space(&o_r, *o_r.last().unwrap(), Role::Right));
        out.d = perm_space(&o_o, *o_o.last().unwrap(), Role::Out).transpose().compose(&self.d);
        out.check_dims();
        out
    }

    /// Lemma 2.10: replace source `j` by `new` (= sources[j].rename(perm) as sets).
    pub fn rename_source(&self, j: usize, new: &Shape, perm: &[usize]) -> Packing {
        let old = &self.sources[j];
        let ren = old.rename(perm);
        assert!(ren.same_as_sets(new), "rename_source: not isomorphic via perm");
        let mut out = self.clone();
        out.sources[j] = new.clone();
        let sp = |role: Role| -> Sparse {
            let o = offsets(&self.sources, role);
            let total = *o.last().unwrap();
            let m = ren.reorder_to(new, role);
            let mut src: Vec<Option<usize>> = (0..total).map(Some).collect();
            for (xn, &xo) in m.iter().enumerate() {
                src[o[j] + xn] = Some(o[j] + xo);
            }
            // rows = new coords, value = old coord
            Sparse::select(&src, total)
        };
        out.el = sp(Role::Left).compose(&self.el);
        out.er = sp(Role::Right).compose(&self.er);
        out.d = self.d.compose(&sp(Role::Out).transpose());
        out.check_dims();
        out
    }
}

/// Lemma 2.3 (scalar version): evaluate a univariate Full_{a,b}[1,1,1] by
/// interpolation at a+b−1 points, as an exact packing into scalar products.
pub fn interpolation_points(n: usize) -> Vec<u64> {
    (1..=n as u64).collect()
}

/// Lagrange basis values: ℓ_i(x) for nodes `xs` (all distinct), over F_p.
pub fn lagrange_at(xs: &[u64], x: u64) -> Vec<u64> {
    (0..xs.len())
        .map(|i| {
            let mut num = 1u64;
            let mut den = 1u64;
            for (j, &xj) in xs.iter().enumerate() {
                if j != i {
                    num = f::mul(num, f::sub(x, xj));
                    den = f::mul(den, f::sub(xs[i], xj));
                }
            }
            f::mul(num, f::inv(den))
        })
        .collect()
}
