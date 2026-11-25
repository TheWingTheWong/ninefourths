use ninefourths::field::Rng;
use ninefourths::gadgets::*;
use ninefourths::shape::Shape;

#[test]
fn full_mid_exchange_identities() {
    for a in 1..5 {
        for b in 1..5 {
            // (9): Full_{b,a} = Full^{12}_{a,b}
            let (s, _) = Shape::full(a, b, 2, 3, 4).ex12();
            assert!(s.same_as_sets(&Shape::full(b, a, 4, 3, 2)));
            // (2): Mid_{a,b} = dual of Full_{a,b}
            let (m, _) = Shape::full(a, b, 2, 3, 4).ex23();
            assert!(m.same_as_sets(&Shape::mid(a, b, 3, 2, 4)));
        }
    }
}

#[test]
fn primal_packing_valid_whenever_ct_ge_a_plus_2b_minus_1() {
    let mut rng = Rng(42);
    let mut n = 0;
    for a in 1..=3 {
        for b in 1..=3 {
            for c in 1..=4 {
                for t in 1..=6 {
                    if c * t < a + 2 * b - 1 {
                        continue;
                    }
                    for &(p, q, r) in &[(1, 1, 1), (2, 1, 3), (1, 2, 2)] {
                        let pk = primal(a, b, c, t, p, q, r);
                        let v = pk.verify(&mut rng, 2).unwrap_or_else(|e| panic!("a={a} b={b} c={c} t={t}: {e}"));
                        assert!(v.out_degree <= 1);
                        n += 1;
                    }
                }
            }
        }
    }
    assert!(n > 100);
}

#[test]
fn three_for_one_gadget() {
    // c = 1, t = a + 2b − 1: OpenAI's Lemma 4.2 degeneration.
    let mut rng = Rng(3);
    for a in 1..=4 {
        for b in 1..=4 {
            let pk = primal(a, b, 1, a + 2 * b - 1, 1, 1, 1);
            assert_eq!((pk.targets.len(), pk.sources.len()), (3, 1));
            pk.verify(&mut rng, 3).unwrap();
        }
    }
}

#[test]
fn primal_fails_when_condition_violated() {
    let mut rng = Rng(5);
    let mut failures = 0;
    let mut cases = 0;
    for a in 1..=3 {
        for b in 1..=3 {
            for c in 1..=3 {
                for t in 1..=4 {
                    if c * t >= a + 2 * b - 1 {
                        continue;
                    }
                    cases += 1;
                    if primal(a, b, c, t, 1, 1, 1).verify(&mut rng, 2).is_err() {
                        failures += 1;
                    }
                }
            }
        }
    }
    assert_eq!(failures, cases, "every instance below the shift condition should fail");
}

#[test]
fn corrupted_coefficient_is_caught() {
    let mut rng = Rng(9);
    let mut pk = primal(2, 2, 3, 2, 1, 1, 1);
    // flip one decoder coefficient
    let row = pk.d.rows.iter().position(|r| !r.is_empty()).unwrap();
    pk.d.rows[row][0].1 = pk.d.rows[row][0].1.add(&ninefourths::laurent::Laurent::constant(1));
    assert!(pk.verify(&mut rng, 1).is_err());
}

#[test]
fn exchanges_and_tensor_preserve_validity() {
    let mut rng = Rng(11);
    let pk = primal(2, 1, 2, 2, 1, 2, 1);
    let dual = pk.ex23();
    dual.verify(&mut rng, 2).unwrap();
    // dual of dual returns the original shapes
    for (x, y) in dual.ex23().targets.iter().zip(&pk.targets) {
        assert!(x.same_as_sets(y));
    }
    let t = pk.tensor(&dual);
    t.verify(&mut rng, 1).unwrap();
}

#[test]
fn amplification_d1() {
    let mut rng = Rng(13);
    // a = b = 1, c = 2, t = 1: 4-for-2 primal packing
    let pk = primal(1, 1, 2, 1, 1, 1, 1);
    pk.verify(&mut rng, 2).unwrap();
    let amp = amplify(&pk, 1);
    // m = 4, d = 1: at least m^2/(d+1) = 8 balanced targets
    assert!(amp.targets.len() >= 8, "{}", amp.targets.len());
    assert_eq!(amp.sources.len(), 4);
    amp.verify(&mut rng, 2).unwrap();
    let t = t_shape(1, 1, 1);
    for s in &amp.targets {
        assert!(s.find_isomorphism(&t).is_some());
    }
}

#[test]
fn amplification_distinguishable_types() {
    let mut rng = Rng(19);
    // a = 2, b = 1, c = 3, t = 1: 5-for-3, Full and Mid differ
    let pk = primal(2, 1, 3, 1, 1, 1, 1);
    let amp = amplify(&pk, 1);
    // Claim 3.3: at least m^{2d}/(d+1) = 25/2 balanced words
    assert!(amp.targets.len() * 2 >= 25, "{}", amp.targets.len());
    assert_eq!(amp.sources.len(), 9);
    amp.verify(&mut rng, 1).unwrap();
    let t = t_shape(2, 1, 1);
    for s in &amp.targets {
        assert!(s.find_isomorphism(&t).is_some(), "{s:?}");
    }
}

#[test]
fn sl_to_pp_small() {
    let mut rng = Rng(17);
    let pk = primal(1, 1, 2, 1, 1, 1, 1);
    let keep: Vec<usize> = (0..2).collect();
    let sl = pk.retain(&keep);
    let pp = sl_to_pp(&sl);
    assert_eq!(pp.sources.len(), 3);
    pp.verify(&mut rng, 2).unwrap();
}
