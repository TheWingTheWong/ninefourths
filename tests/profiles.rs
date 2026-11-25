//! The cost model's ε bookkeeping (src/cost.rs) against ε exponent ranges
//! measured on materialized packings.
use ninefourths::cost::{self, Profile, Range};
use ninefourths::gadgets::*;
use ninefourths::laurent::Laurent;
use ninefourths::packing::Packing;

fn measured(pk: &Packing) -> Profile {
    let p = pk.profile();
    let r = |(lo, hi): (i64, i64)| Range::new(lo as f64, hi as f64);
    Profile { el: r(p[0]), er: r(p[1]), d: r(p[2]) }
}

/// measured ⊆ predicted, interval by interval
fn within(m: Profile, p: Profile) -> bool {
    let ok = |a: Range, b: Range| a.lo >= b.lo.min(0.0) && a.hi <= b.hi.max(0.0);
    ok(m.el, p.el) && ok(m.er, p.er) && ok(m.d, p.d)
}

#[test]
fn primal_profile_exact() {
    for (a, b, c, t) in [(1, 1, 2, 1), (2, 3, 4, 2), (3, 1, 1, 4)] {
        assert_eq!(measured(&primal(a, b, c, t, 1, 1, 1)), cost::primal_profile());
    }
}

#[test]
fn amplified_profile_exact() {
    // full tensor (before selection) has the predicted [−d, d] ranges ...
    let pk = primal(2, 1, 3, 1, 1, 1, 1);
    let pair = pk.tensor(&pk.ex23());
    assert_eq!(measured(&pair), cost::amplified_profile_all(1));
    assert_eq!(measured(&pair.tensor(&pair)), cost::amplified_profile_all(2));
    // ... and the balanced targets of Lemma 3.2 form an exact packing
    for (a, b, c, t, d) in [(1, 1, 2, 1, 1), (1, 1, 2, 1, 2), (2, 1, 3, 1, 1), (1, 2, 2, 2, 2), (1, 1, 3, 1, 2)] {
        let amp = amplify(&primal(a, b, c, t, 1, 1, 1), d);
        assert!(amp.is_exact(), "a={a} b={b} c={c} t={t} d={d}");
        assert_eq!(measured(&amp), cost::amplified_profile(d as u32));
    }
}

#[test]
fn sl_to_pp_profile_within_prediction() {
    let pk = primal(1, 1, 2, 1, 1, 1, 1);
    let amp = amplify(&pk, 1);
    for w in [2usize, 3, 4] {
        let keep: Vec<usize> = (0..w).collect();
        let pp = sl_to_pp(&amp.retain(&keep));
        let m = measured(&pp);
        let p = cost::slpp_profile(cost::amplified_profile(1), w as f64);
        assert!(within(m, p), "w={w}: measured {m:?} predicted {p:?}");
        // the cost model is tight: it uses exactly the measured ranges
        assert_eq!(m, p, "w={w}");
    }
}

fn perturbed(a: usize, b: usize, b2: usize, d: usize) -> Packing {
    let mut pk = minimal_round_input(a, b, b2, d);
    pk.er = pk.er.scale(&Laurent::mono(1, -1));
    pk.d = pk.d.scale(&Laurent::mono(1, 1));
    for i in 0..pk.d.nrows {
        pk.d.add_entry(i, (i * 7) % pk.d.ncols, &Laurent::mono(5, 2));
    }
    pk
}

#[test]
fn round_profile_matches_prediction() {
    let (l, l2, d) = (1, 2, 1);
    let p1 = perturbed(l, l, l2, d);
    let p2 = perturbed(l2, l, l2, d);
    let r = round(&p1, &p2, RoundParams { l, l2, d, w: 1, u: 1, x: 1 });
    let predicted = cost::round_profile(measured(&p1));
    assert_eq!(measured(&p1), measured(&p2));
    assert_eq!(measured(&r), predicted);
}

#[test]
fn chain_count_matches_cost_formula() {
    // E_k = 64^k w^9 (2 L_k − 1)^{6d} at w = 1, k = 1, c = 3 (L_1 = 2), d = 1.
    let lr = 6.0 + 6.0 * cost::log2_two_lk_minus_one(3, 1);
    assert!((2f64.powf(lr) - 64.0 * 729.0).abs() < 1e-6);
    let pk = minimal_chain(3, 1, 1, 1, &|a, b, b2, d| minimal_round_input(a, b, b2, d));
    assert_eq!(pk.sources.len(), 64 * 729);
}
