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

#[test]
fn lemma_2_3_evaluation_uses_exactly_effort_products() {
    let mut rng = Rng(23);
    let mut shapes = vec![
        Shape::full(3, 2, 1, 1, 1),
        Shape::mid(3, 2, 1, 1, 1),
        Shape::mid(2, 4, 2, 1, 2),
        Shape::full(1, 1, 2, 3, 2),
        t_shape(2, 2, 1),
        t_shape(2, 1, 2),
    ];
    let t = t_shape(2, 1, 1);
    shapes.push(t.tensor(&t.ex12().0).tensor(&t.ex12().0.ex23().0.ex12().0));
    for s in shapes {
        let pk = evaluate(&s);
        assert_eq!(pk.sources.len() as u128, s.effort(), "{s:?}");
        assert!(pk.is_exact());
        pk.verify(&mut rng, 2).unwrap();
    }
}

fn toy_input(a: usize, b: usize, b2: usize, d: usize, perturb: bool) -> ninefourths::packing::Packing {
    use ninefourths::laurent::Laurent;
    let mut pk = embedding(&t_shape(a, b, d), &t_shape(a, b2, d));
    if perturb {
        // ER·ε⁻¹, D·ε, plus an ε·junk term in D: valid modulo ε, with negative powers
        pk.er = pk.er.scale(&Laurent::mono(1, -1));
        pk.d = pk.d.scale(&Laurent::mono(1, 1));
        for i in 0..pk.d.nrows {
            pk.d.add_entry(i, (i * 7) % pk.d.ncols, &Laurent::mono(5, 2));
        }
    }
    pk
}

#[test]
fn round_lemma_with_minimal_inputs() {
    let mut rng = Rng(29);
    for &perturb in &[false, true] {
        let (l, l2, d) = (1, 2, 1);
        let p1 = toy_input(l, l, l2, d, perturb);
        let p2 = toy_input(l2, l, l2, d, perturb);
        p1.verify(&mut rng, 1).unwrap();
        p2.verify(&mut rng, 1).unwrap();
        let rp = RoundParams { l, l2, d, w: 1, u: 1, x: 1 };
        let r = round(&p1, &p2, rp);
        assert_eq!(r.targets.len(), 1);
        assert_eq!(r.sources.len(), 64);
        assert_eq!(r.targets[0], q_shape(1, 1, 1, 1, 1, 1));
        r.verify(&mut rng, 1).unwrap();
    }
}

#[test]
fn minimal_chain_multiplies_matrices() {
    use ninefourths::field as f;
    use ninefourths::laurent::Laurent;
    let mut rng = Rng(31);
    // c = 3: L_0 = 1, L_1 = 2; one round, d = 1, side x = 2.
    let x = 2;
    let pk = minimal_chain(3, 1, 1, x, &|a, b, b2, d| minimal_round_input(a, b, b2, d));
    assert!(pk.is_exact());
    // 64 sources Q_{2,2}[x], each evaluated with eff = (2·2−1)^{6d} x³ products
    assert_eq!(pk.sources.len(), 64 * 3usize.pow(6) * x * x * x);
    assert_eq!(pk.targets[0], q_shape(1, 1, 1, x, x, x));
    pk.verify(&mut rng, 1).unwrap();
    // run on integer matrices against schoolbook
    let a: Vec<i64> = (0..x * x).map(|i| i as i64 * 3 - 5).collect();
    let b: Vec<i64> = (0..x * x).map(|i| 7 - i as i64 * 2).collect();
    let out = pk.run(&a.iter().map(|&v| Laurent::int(v)).collect::<Vec<_>>(), &b.iter().map(|&v| Laurent::int(v)).collect::<Vec<_>>());
    for i in 0..x {
        for j in 0..x {
            let want: i64 = (0..x).map(|t| a[i * x + t] * b[t * x + j]).sum();
            assert_eq!(f::to_i64(out[i * x + j].coeff(0)), want);
        }
    }
}

#[test]
fn minimal_chain_mod_eps_then_bini() {
    let mut rng = Rng(37);
    let input = |a: usize, b: usize, b2: usize, d: usize| {
        use ninefourths::laurent::Laurent;
        let mut pk = minimal_round_input(a, b, b2, d);
        pk.er = pk.er.scale(&Laurent::mono(1, -1));
        pk.d = pk.d.scale(&Laurent::mono(1, 1));
        for i in 0..pk.d.nrows {
            pk.d.add_entry(i, (i * 7) % pk.d.ncols, &Laurent::mono(5, 2));
        }
        pk
    };
    let pk = minimal_chain(3, 1, 1, 1, &input);
    assert!(!pk.is_exact());
    let v = pk.verify(&mut rng, 1).unwrap();
    let bound = degree_bound(&pk);
    assert!(v.out_degree as usize <= bound);
    let ex = bini(&pk, bound);
    assert!(ex.is_exact());
    assert_eq!(ex.sources.len(), (bound + 1) * pk.sources.len());
    ex.verify(&mut rng, 1).unwrap();
}
