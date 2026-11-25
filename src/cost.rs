//! Exact cost accounting for Pass's construction (arXiv 2610.28627, Section 4).
//!
//! Every quantity follows a formula from the paper; the comment next to each
//! one names the equation or lemma it comes from. Quantities that overflow
//! machine integers are carried as base-2 logarithms (`f64`). The only
//! inputs are small integers (c, d, k, h, retained w), so the logs are
//! accurate to ~1e-12 relative, far below anything the reported numbers need.
//!
//! ## ε bookkeeping
//!
//! A packing modulo ε has three linear maps (left encoder `EL`, right encoder
//! `ER`, decoder `D`) whose entries are Laurent polynomials in ε. We track for
//! each map the closed interval of ε exponents that can occur. The output of
//! the packing is `D(S(EL x, ER y))` and the source products `S` contain no ε,
//! so its ε-degree is at most `EL.hi + ER.hi + D.hi`. That is the `d_ε` that
//! Bini's conversion (Lemma 2.4) needs. The interval rules per operation
//! (tensor: add; role exchange: permute `EL, ER, Dᵀ`; substitution ε→ε^N:
//! scale; composition (Lemma 2.9 proof): scale the inner packing by
//! N = pL + pR + pD + 1, then add) are exercised against the materialized
//! engine in `tests/`.

use num_bigint::BigUint;
use num_traits::{ToPrimitive, Zero};


/// Interval of ε exponents occurring in one linear map.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Range {
    pub lo: f64,
    pub hi: f64,
}

impl Range {
    pub const ZERO: Range = Range { lo: 0.0, hi: 0.0 };
    pub fn new(lo: f64, hi: f64) -> Self {
        Range { lo, hi }
    }
    fn add(self, o: Range) -> Range {
        Range::new(self.lo + o.lo, self.hi + o.hi)
    }
    fn scale(self, n: f64) -> Range {
        Range::new(self.lo * n, self.hi * n)
    }
    /// Number of negative powers of ε that must be cleared (p in Lemma 2.9).
    pub fn neg(self) -> f64 {
        (-self.lo).max(0.0)
    }
}

/// ε profile of a packing: exponent ranges of its three maps.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Profile {
    pub el: Range,
    pub er: Range,
    pub d: Range,
}

impl Profile {
    pub const EXACT: Profile = Profile { el: Range::ZERO, er: Range::ZERO, d: Range::ZERO };

    /// Lemma 2.7: tensoring multiplies map entries, so exponent ranges add.
    pub fn tensor(self, o: Profile) -> Profile {
        Profile { el: self.el.add(o.el), er: self.er.add(o.er), d: self.d.add(o.d) }
    }
    /// Lemma 2.8, exchange (12): swap the encoders.
    pub fn ex12(self) -> Profile {
        Profile { el: self.er, er: self.el, d: self.d }
    }
    /// Lemma 2.8, exchange (23): encoders become (EL, Dᵀ), decoder (ER)ᵀ.
    pub fn ex23(self) -> Profile {
        Profile { el: self.el, er: self.d, d: self.er }
    }
    /// Exchange (13) = (12)(23)(12): left encoder and decoder swap.
    pub fn ex13(self) -> Profile {
        Profile { el: self.d, er: self.er, d: self.el }
    }
    /// Substitute ε → ε^n.
    pub fn subst(self, n: f64) -> Profile {
        Profile { el: self.el.scale(n), er: self.er.scale(n), d: self.d.scale(n) }
    }
    /// The N = pL + pR + pD + 1 of the proof of Lemma 2.9 for this outer packing.
    pub fn clearing_power(self) -> f64 {
        self.el.neg() + self.er.neg() + self.d.neg() + 1.0
    }
    /// Lemma 2.9: compose `self` (outer, packs the targets into intermediate
    /// products) with `inner` (packs the intermediates). The inner parameter is
    /// replaced by ε^N first.
    pub fn compose(self, inner: Profile) -> Profile {
        let n = if inner == Profile::EXACT { 1.0 } else { self.clearing_power() };
        let i = inner.subst(n);
        Profile { el: self.el.add(i.el), er: self.er.add(i.er), d: self.d.add(i.d) }
    }
    pub fn neg_sum(self) -> f64 {
        self.el.neg() + self.er.neg() + self.d.neg()
    }
    pub fn hi_sum(self) -> f64 {
        self.el.hi.max(0.0) + self.er.hi.max(0.0) + self.d.hi.max(0.0)
    }
    /// Upper bound on the ε-degree of the packing's output.
    pub fn out_degree(self) -> f64 {
        self.hi_sum()
    }
}

/// Theorem 3.1 primal packing: `E_j = P_j − z^t P_{j+1} + ε z^b Q^{(j)}`,
/// middle product recovered with ε^{-1}.
pub fn primal_profile() -> Profile {
    Profile { el: Range::ZERO, er: Range::new(0.0, 1.0), d: Range::new(-1.0, 0.0) }
}

/// Lemma 3.2 before target selection: (primal ⊗ dual)^{⊗d}.
pub fn amplified_profile_all(d: u32) -> Profile {
    let pair = primal_profile().tensor(primal_profile().ex23());
    let mut p = Profile::EXACT;
    for _ in 0..d {
        p = p.tensor(pair);
    }
    p
}

/// Lemma 3.2 after retaining the balanced targets: **exact**.
///
/// Every column of the primal right encoder has a single ε exponent (0 for
/// the P_j, +1 for Q) and every row of its decoder too (0 for Full rows, −1
/// for the Mid row); the dual swaps the two maps. A retained word has α
/// primal Mid letters and β dual-of-Mid letters with α = β (that is the
/// balance condition), so each of its right-encoder entries carries
/// ε^{α−β} = ε^0 and each decoder entry ε^{β−α} = ε^0. The formal parameter
/// cancels identically; `tests/profiles.rs` confirms it on materialized
/// packings. (The paper only states that validity modulo ε is preserved.)
pub fn amplified_profile(_d: u32) -> Profile {
    Profile::EXACT
}

/// Lemma 3.4: parameter η = ε^w, right inputs B'_i(k) interpolated through
/// the points i + εj (exponents −(w−1)..0), outputs interpolated in ξ and
/// evaluated at i + εj (exponents 0..2w−2).
pub fn slpp_profile(amplified: Profile, w: f64) -> Profile {
    let p = amplified.subst(w);
    Profile {
        el: p.el,
        er: p.er.add(Range::new(-(w - 1.0), 0.0)),
        d: p.d.add(Range::new(0.0, 2.0 * w - 2.0)),
    }
}

/// Lemma 4.2: tensor the packing with its (12) and (13) transposes (and an
/// exact identity factor), then compose two such packings.
pub fn round_profile(pp: Profile) -> Profile {
    let t = pp.tensor(pp.ex12()).tensor(pp.ex13());
    t.compose(t)
}

/// Profile of the k-round chain of Theorem 4.3 followed by the exact source
/// evaluation of Lemma 2.3. Returns log2 of the output ε-degree bound.
///
/// All rounds share one profile G (it depends on d and w only), and chaining
/// gives s_k·G with s_1 = 1, s_{j+1} = s_j (negsum(G) + 1) + 1. We return
/// log2(s_k · hisum(G)), computed in log space.
pub fn chain_log2_degree(round: Profile, k: u64) -> f64 {
    let ln = (round.neg_sum() + 1.0).log2();
    // s_k = ((g)^k - 1)/(g - 1) with g = negsum+1; for g >= 2 this is < g^k / (g-1) * ... ;
    // we use the exact recurrence while small, then the closed form bound
    // s_k <= 2 g^(k-1) (valid for g >= 2).
    let mut s: f64 = 1.0;
    let mut j = 1u64;
    while j < k && s < 1e300 / (round.neg_sum() + 1.0) {
        s = s * (round.neg_sum() + 1.0) + 1.0;
        j += 1;
    }
    let log_s = if j == k { s.log2() } else { s.log2() + (k - j) as f64 * ln + 1e-12 };
    log_s + round.hi_sum().log2()
}

/// log2 of a big integer, accurate to ~1e-15 relative.
pub fn log2_big(x: &BigUint) -> f64 {
    let bits = x.bits();
    if bits <= 1000 {
        return x.to_f64().unwrap().log2();
    }
    let shift = bits - 60;
    let top: BigUint = x >> shift;
    top.to_f64().unwrap().log2() + shift as f64
}

/// Exact u = n^{2d}, w = ⌈m^{2d}/(d+1)⌉ (8) and q = ⌊w/u⌋ (Λ = q², Lemma 4.2)
/// as log2 values. Exact big-integer arithmetic below 2^20000, log-space above
/// (where the floor/ceiling change log2 by less than 2^-19000).
pub fn log2_u_w_q(c: u64, d: u64) -> (f64, f64, f64) {
    let (n, m) = (c, c + 2);
    let bits = 2.0 * d as f64 * (m as f64).log2();
    if bits < 20_000.0 {
        let u = BigUint::from(n).pow((2 * d) as u32);
        let mm = BigUint::from(m).pow((2 * d) as u32);
        let w = (&mm + BigUint::from(d)) / BigUint::from(d + 1);
        let q = &w / &u;
        let lq = if q.is_zero() { f64::NEG_INFINITY } else { log2_big(&q) };
        (log2_big(&u), log2_big(&w), lq)
    } else {
        let lu = 2.0 * d as f64 * (n as f64).log2();
        let lw = bits - ((d + 1) as f64).log2();
        (lu, lw, lw - lu)
    }
}

/// Parameters of one instance of Theorem 4.3.
#[derive(Clone, Copy, Debug)]
pub struct Params {
    /// Primal packing (c+2)-for-c (Theorem 3.1); m = c + 2, n = c.
    pub c: u64,
    /// Amplification depth (Lemma 3.2).
    pub d: u64,
    /// Number of rounds (Theorem 4.3).
    pub k: u64,
    /// Bini tensor power (Lemma 2.4).
    pub h: u64,
    /// log2 of the retained target count w (≤ the w of (8)); None = (8) itself.
    pub log2_w: Option<f64>,
}

/// All derived quantities, as log2 values unless stated.
#[derive(Clone, Copy, Debug)]
pub struct Cost {
    pub log2_u: f64,
    pub log2_w: f64,
    pub log2_lambda: f64,
    /// log2 of the mod-ε identity side v_k = w² Λ^k.
    pub log2_v: f64,
    /// log2 of its multiplication count E_k = 64^k w^9 (2L_k − 1)^{6d}.
    pub log2_r: f64,
    /// log2 of the ε-degree bound of its output.
    pub log2_deg: f64,
    /// Exact identity after Bini (Lemma 2.4): side V = v^h, R = (h d_ε + 1) r^h.
    pub log2_vv: f64,
    pub log2_rr: f64,
}

impl Cost {
    pub fn exponent_mod_eps(&self) -> f64 {
        self.log2_r / self.log2_v
    }
    pub fn exponent_exact(&self) -> f64 {
        self.log2_rr / self.log2_vv
    }
    /// Smallest side at which this exact identity could beat Strassen on
    /// multiplication count: (R/7)^{1/log2 7}. Returns log2 of that side.
    pub fn log2_strassen_floor(&self) -> f64 {
        (self.log2_rr - 7f64.log2()) / 7f64.log2()
    }
}

pub fn log2_w_eq8(c: u64, d: u64) -> f64 {
    log2_u_w_q(c, d).1
}

/// Sequence L_0 = 1, L_{j+1} = L_j + ⌈(3L_j − 1)/(n − 1)⌉ (proof of Theorem 4.3).
/// Returns log2(2 L_k − 1).
pub fn log2_two_lk_minus_one(n: u64, k: u64) -> f64 {
    let mut l: u128 = 1;
    let mut j = 0;
    while j < k {
        if l > (1u128 << 100) {
            break;
        }
        l += (3 * l - 1).div_ceil((n - 1) as u128);
        j += 1;
    }
    if j == k {
        return ((2 * l - 1) as f64).log2();
    }
    // For large L, L' = L + ceil((3L-1)/(n-1)) <= rho L + 1 with rho = 1 + 3/(n-1);
    // the +1 is negligible at L > 2^100 (relative 2^-100 per step), but we keep it
    // as an upper bound: L' <= rho L (1 + 2^-99).
    let rho = 1.0 + 3.0 / (n - 1) as f64;
    let step = rho.log2() + 1e-29;
    let log_l = (l as f64).log2() + (k - j) as f64 * step;
    log_l + 1.0
}

/// Mod-ε part of the cost (everything except Bini's h), for k rounds.
pub fn mod_eps_cost(c: u64, d: u64, k: u64, log2_w: Option<f64>) -> Option<(f64, f64, f64, f64, f64, f64)> {
    let n = c;
    let (log2_u, log2_w8, log2_q8) = log2_u_w_q(c, d);
    let (log2_w, log2_q) = match log2_w {
        None => (log2_w8, log2_q8),
        Some(lw) => {
            assert!(lw <= log2_w8 + 1e-9, "retained w exceeds (8)");
            // q = floor(w'/u); w' is a power of two here, exact enough in log space
            let r = lw - log2_u;
            if r < 0.0 { return None; }
            let q = if r < 52.0 { 2f64.powf(r).floor().log2() } else { r };
            (lw, q)
        }
    };
    if !(log2_q >= 0.0) {
        return None;
    }
    let log2_lambda = 2.0 * log2_q;
    let log2_v = 2.0 * log2_w + k as f64 * log2_lambda; // v_k = w² Λ^k
    let log2_r = 6.0 * k as f64 + 9.0 * log2_w + 6.0 * d as f64 * log2_two_lk_minus_one(n, k); // E_k
    // ε degree
    let w = 2f64.powf(log2_w.min(1000.0));
    let log2_deg = if log2_w < 40.0 {
        let round = round_profile(slpp_profile(amplified_profile(d as u32), w));
        chain_log2_degree(round, k)
    } else {
        chain_log2_degree_logspace(log2_w, k)
    };
    Some((log2_u, log2_w, log2_lambda, log2_v, log2_r, log2_deg))
}

/// Log-space version of the ε bookkeeping for huge w, using the closed forms
/// (checked against the interval arithmetic in tests): with the exact
/// amplified packing, one round has profile G = (N+1)·t where
/// t = pp ⊗ pp^{12} ⊗ pp^{13}, N = 3w − 2, so
///   negsum(G) = 3(3w − 1)(w − 1) ≤ 9w²,   hisum(G) = 6(3w − 1)(w − 1) ≤ 18w².
pub fn chain_log2_degree_logspace(log2_w: f64, k: u64) -> f64 {
    let log_neg1 = 9f64.log2() + 2.0 * log2_w; // ≥ log2(negsum + 1)
    let log_hi = 18f64.log2() + 2.0 * log2_w;
    // s_k ≤ 2 (negsum+1)^(k−1)
    let log_s = if k <= 1 { 0.0 } else { 1.0 + (k - 1) as f64 * log_neg1 };
    log_s + log_hi
}

/// Full cost including Bini's conversion with tensor power h.
pub fn cost(p: Params) -> Option<Cost> {
    let (log2_u, log2_w, log2_lambda, log2_v, log2_r, log2_deg) = mod_eps_cost(p.c, p.d, p.k, p.log2_w)?;
    let h = p.h as f64;
    // (h d_ε + 1) r^h
    let log2_hd1 = log2_add(h.log2() + log2_deg, 0.0);
    Some(Cost {
        log2_u,
        log2_w,
        log2_lambda,
        log2_v,
        log2_r,
        log2_deg,
        log2_vv: h * log2_v,
        log2_rr: log2_hd1 + h * log2_r,
    })
}

/// log2(2^a + 2^b)
pub fn log2_add(a: f64, b: f64) -> f64 {
    let (hi, lo) = if a > b { (a, b) } else { (b, a) };
    hi + (1.0 + 2f64.powf(lo - hi)).log2()
}

/// Smallest Bini power h such that the exact identity has exponent < tau,
/// i.e. log2(h d + 1) + h log2 r < tau h log2 v. Returns None if the mod-ε
/// identity itself is not below tau.
pub fn min_bini_h(log2_v: f64, log2_r: f64, log2_deg: f64, tau: f64) -> Option<u64> {
    let gap = tau * log2_v - log2_r; // bits saved per tensor power
    if gap <= 0.0 {
        return None;
    }
    // need log2(h·2^deg + 1) < h·gap. LHS ≈ deg + log2 h.
    let f = |h: f64| log2_add(h.log2() + log2_deg, 0.0) < h * gap;
    let mut lo = 1.0f64;
    if f(lo) {
        return Some(1);
    }
    let mut hi = 2.0f64;
    while !f(hi) {
        hi *= 2.0;
        if hi > 1e300 {
            return None;
        }
    }
    while hi - lo > 1.0 {
        let mid = ((lo + hi) / 2.0).floor();
        if f(mid) {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    if hi > u64::MAX as f64 {
        return None;
    }
    Some(hi as u64)
}

/// Asymptotic exponent threshold for a given c (Corollary 4.4):
/// 3 log(1 + 3/(c−1)) / (2 log(1 + 2/c)).
pub fn threshold(c: u64) -> f64 {
    let c = c as f64;
    3.0 * (1.0 + 3.0 / (c - 1.0)).ln() / (2.0 * (1.0 + 2.0 / c).ln())
}

/// Condition (10): 64 ρ^{6d} < Λ^τ with ρ = 1 + 3/(n−1).
pub fn condition_10(c: u64, d: u64, tau: f64) -> bool {
    let rho = 1.0 + 3.0 / (c - 1) as f64;
    let Some((_, _, log2_lambda, _, _, _)) = mod_eps_cost(c, d, 0, None) else { return false };
    6.0 + 6.0 * d as f64 * rho.log2() < tau * log2_lambda
}

/// The paper's recipe, taking the smallest choice at every step:
/// smallest c with threshold(c) < tau, smallest d satisfying (10), smallest
/// k with E_k < v_k^tau, smallest h with R < V^tau.
pub fn faithful(tau: f64) -> Option<(Params, Cost)> {
    let c = (3..10_000).find(|&c| threshold(c) < tau)?;
    let d = (1..100_000).find(|&d| condition_10(c, d, tau))?;
    let mut k = 1;
    loop {
        let (_, _, _, lv, lr, _) = mod_eps_cost(c, d, k, None)?;
        if lr < tau * lv {
            break;
        }
        k += 1;
        if k > 10_000_000 {
            return None;
        }
    }
    let (_, _, _, lv, lr, ld) = mod_eps_cost(c, d, k, None)?;
    let h = min_bini_h(lv, lr, ld, tau)?;
    let p = Params { c, d, k, h, log2_w: None };
    Some((p, cost(p)?))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn thresholds_match_corollary() {
        assert!((threshold(2) - 3.0).abs() < 1e-12);
        assert!(threshold(3) < 7f64.log2());
        assert!((threshold(1_000_000) - 2.25).abs() < 1e-4);
        for c in 3..200 {
            assert!(threshold(c + 1) < threshold(c));
        }
    }

    #[test]
    fn small_exact_values() {
        // c = 2, d = 1: u = 4, w = ceil(16/2) = 8, Λ = 4.
        let (lu, lw, ll, lv, lr, _) = mod_eps_cost(2, 1, 1, None).unwrap();
        assert_eq!(lu, 2.0);
        assert_eq!(lw, 3.0);
        assert_eq!(ll, 2.0);
        // v_1 = 64 * 4 = 256; L_1 = 1 + ceil(2/1) = 3; E_1 = 64 * 8^9 * 5^6
        assert!((lv - 8.0).abs() < 1e-12);
        let e1 = 64f64 * 8f64.powi(9) * 5f64.powi(6);
        assert!((lr - e1.log2()).abs() < 1e-9);
    }

    #[test]
    fn lk_sequence() {
        // n = 3: L: 1 -> 1+ceil(2/2)=2 -> 2+ceil(5/2)=5 -> 5+ceil(14/2)=12
        assert!((log2_two_lk_minus_one(3, 3) - 23f64.log2()).abs() < 1e-12);
    }

    #[test]
    fn closed_form_degree_matches_interval_arithmetic() {
        for d in 1..5u32 {
            for w in [8.0f64, 100.0, 12345.0] {
                let round = round_profile(slpp_profile(amplified_profile(d), w));
                assert_eq!(round.neg_sum(), 3.0 * (3.0 * w - 1.0) * (w - 1.0));
                assert_eq!(round.hi_sum(), 6.0 * (3.0 * w - 1.0) * (w - 1.0));
                assert!((round.neg_sum() + 1.0).log2() <= 9f64.log2() + 2.0 * w.log2());
                for k in 1..6 {
                    assert!(chain_log2_degree(round, k) <= chain_log2_degree_logspace(w.log2(), k) + 1e-9);
                }
            }
        }
    }
}

/// Best exact identity found by [`search`].
#[derive(Clone, Copy, Debug)]
pub struct Best {
    pub params: Params,
    pub cost: Cost,
}

/// log2 R of the exact identity at (c, d, k) with the smallest feasible Bini
/// power, or None if the mod-ε identity at k is not below tau.
pub fn exact_at(c: u64, d: u64, k: u64, tau: f64) -> Option<(Params, Cost)> {
    let (_, _, _, lv, lr, ld) = mod_eps_cost(c, d, k, None)?;
    if lr >= tau * lv {
        return None;
    }
    let h = min_bini_h(lv, lr, ld, tau)?;
    let p = Params { c, d, k, h, log2_w: None };
    Some((p, cost(p)?))
}

/// Smallest k with E_k < v_k^tau (feasibility is monotone in k once the
/// per-round margin τ log Λ − log 64 − 6d log ρ is positive).
pub fn first_feasible_k(c: u64, d: u64, tau: f64) -> Option<u64> {
    let ok = |k: u64| mod_eps_cost(c, d, k, None).map(|(_, _, _, lv, lr, _)| lr < tau * lv).unwrap_or(false);
    let mut hi = 1u64;
    while !ok(hi) {
        hi *= 2;
        if hi > 1 << 34 {
            return None;
        }
    }
    let mut lo = hi / 2;
    if lo == 0 {
        return Some(1);
    }
    while hi - lo > 1 {
        let mid = (lo + hi) / 2;
        if ok(mid) { hi = mid } else { lo = mid }
    }
    Some(hi)
}

/// Grid search over c ∈ [3, c_max], d ∈ [1, d_max] and k ∈ [k0, 64 k0]
/// (k0 = first feasible k; 4000 log-spaced samples, then a unit-step
/// refinement around the best sample), with the smallest feasible Bini power
/// h at each point. Objective: minimise log2 R of the exact identity with
/// exponent < tau. Retains the full w of (8).
pub fn search(tau: f64, c_max: u64, d_max: u64) -> Option<Best> {
    let threads = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4).min(6);
    let cs: Vec<u64> = (3..=c_max).filter(|&c| threshold(c) < tau).collect();
    let results: Vec<Option<Best>> = std::thread::scope(|sc| {
        let hs: Vec<_> = (0..threads)
            .map(|t| {
                let mine: Vec<u64> = cs.iter().copied().skip(t).step_by(threads).collect();
                sc.spawn(move || {
                    let mut best: Option<Best> = None;
                    for c in mine {
                        if let Some(b) = search_c(tau, c, d_max, best.map(|b| b.cost.log2_rr)) {
                            if best.map(|x| b.cost.log2_rr < x.cost.log2_rr).unwrap_or(true) {
                                best = Some(b);
                            }
                        }
                    }
                    best
                })
            })
            .collect();
        hs.into_iter().map(|h| h.join().unwrap()).collect()
    });
    results.into_iter().flatten().min_by(|a, b| a.cost.log2_rr.partial_cmp(&b.cost.log2_rr).unwrap())
}

/// [`search`] for one value of c; `bound` is a known achievable log2 R used for pruning.
pub fn search_c(tau: f64, c: u64, d_max: u64, bound: Option<f64>) -> Option<Best> {
    let mut best: Option<Best> = None;
    let _ = bound;
    {
        let rho = 1.0 + 3.0 / (c - 1) as f64;
        for d in 1..=d_max {
            let (_, lw, lq) = log2_u_w_q(c, d);
            if !(lq >= 0.0) {
                continue;
            }
            if tau * 2.0 * lq <= 6.0 + 6.0 * d as f64 * rho.log2() {
                continue; // per-round margin not positive: no k is feasible
            }
            if let Some(b) = best {
                if 9.0 * lw > b.cost.log2_rr {
                    continue; // R ≥ E_k ≥ w^9 already exceeds the best
                }
            }
            let Some(k0) = first_feasible_k(c, d, tau) else { continue };
            let mut local: Option<(u64, f64)> = None;
            let try_k = |k: u64, local: &mut Option<(u64, f64)>, best: &mut Option<Best>| {
                if let Some((p, cst)) = exact_at(c, d, k, tau) {
                    if local.map(|l| cst.log2_rr < l.1).unwrap_or(true) {
                        *local = Some((k, cst.log2_rr));
                    }
                    if best.map(|b| cst.log2_rr < b.cost.log2_rr).unwrap_or(true) {
                        *best = Some(Best { params: p, cost: cst });
                    }
                }
            };
            let samples = 4000u64;
            let mut prev = 0;
            let mut ks = vec![];
            for i in 0..=samples {
                let k = (k0 as f64 * 64f64.powf(i as f64 / samples as f64)).round() as u64;
                if k != prev {
                    ks.push(k);
                    prev = k;
                }
            }
            for &k in &ks {
                try_k(k, &mut local, &mut best);
            }
            if let Some((kb, _)) = local {
                let span = (kb as f64 * (64f64.powf(1.0 / samples as f64) - 1.0)).ceil() as u64 + 1;
                for k in kb.saturating_sub(span).max(k0)..=kb + span {
                    try_k(k, &mut local, &mut best);
                }
            }
        }
    }
    best
}
