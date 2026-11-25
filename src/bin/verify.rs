//! Builds every construction of Pass's proof as an explicit packing at small
//! parameters and checks it exactly (random points of F_p, ε symbolic).
//! Prints one PASS/FAIL line per lemma. `--quick` skips the end-to-end chain.

use ninefourths::cost;
use ninefourths::field::{self as f, Rng};
use ninefourths::gadgets::*;
use ninefourths::laurent::Laurent;
use ninefourths::packing::Packing;
use ninefourths::shape::Shape;
use std::time::Instant;

struct Line {
    name: &'static str,
    what: String,
    pass: usize,
    total: usize,
    secs: f64,
}

fn check(lines: &mut Vec<Line>, name: &'static str, what: &str, cases: Vec<Box<dyn Fn() -> bool>>) {
    let t = Instant::now();
    let total = cases.len();
    let pass = cases.iter().filter(|c| c()).count();
    let l = Line { name, what: what.to_string(), pass, total, secs: t.elapsed().as_secs_f64() };
    print_line(&l);
    lines.push(l);
}

fn print_line(l: &Line) {
    let status = if l.pass == l.total { "PASS" } else { "FAIL" };
    println!("{status}  {:<22} {:>4}/{:<4} {:>7.2}s  {}", l.name, l.pass, l.total, l.secs, l.what);
}

fn ok(pk: &Packing, seed: u64) -> bool {
    pk.verify(&mut Rng(seed), 2).is_ok()
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

fn main() {
    let quick = std::env::args().any(|a| a == "--quick");
    println!("ninefourths verify: every check is exact over F_p, p = 2^61 - 1, with eps symbolic\n");
    let mut lines = vec![];

    // Shapes: (9) and (2)
    let mut cases: Vec<Box<dyn Fn() -> bool>> = vec![];
    for a in 1..6 {
        for b in 1..6 {
            cases.push(Box::new(move || {
                Shape::full(a, b, 2, 3, 4).ex12().0.same_as_sets(&Shape::full(b, a, 4, 3, 2))
                    && Shape::full(a, b, 2, 3, 4).ex23().0.same_as_sets(&Shape::mid(a, b, 3, 2, 4))
            }));
        }
    }
    check(&mut lines, "shapes (2), (9)", "Full^12_{a,b} = Full_{b,a}; Mid_{a,b} = dual of Full_{a,b}", cases);

    // Theorem 3.1
    let mut cases: Vec<Box<dyn Fn() -> bool>> = vec![];
    let mut neg: Vec<Box<dyn Fn() -> bool>> = vec![];
    for a in 1..=3 {
        for b in 1..=3 {
            for c in 1..=4 {
                for t in 1..=6 {
                    if c * t >= a + 2 * b - 1 {
                        for &(p, q, r) in &[(1, 1, 1), (2, 1, 3)] {
                            cases.push(Box::new(move || ok(&primal(a, b, c, t, p, q, r), (a * 1000 + b * 100 + c * 10 + t) as u64)));
                        }
                    } else {
                        neg.push(Box::new(move || !ok(&primal(a, b, c, t, 1, 1, 1), 7)));
                    }
                }
            }
        }
    }
    check(&mut lines, "Theorem 3.1 primal", "(c+2)-for-c SL-PP mod eps, all a,b<=3, c<=4, t<=6 with ct >= a+2b-1", cases);
    check(&mut lines, "Theorem 3.1 negative", "every instance with ct < a+2b-1 is rejected by the verifier", neg);
    let mut cases: Vec<Box<dyn Fn() -> bool>> = vec![];
    for a in 1..=4 {
        for b in 1..=4 {
            cases.push(Box::new(move || ok(&primal(a, b, 1, a + 2 * b - 1, 1, 1, 1), 11)));
        }
    }
    check(&mut lines, "OpenAI Lemma 4.2", "three-for-one degeneration (c = 1, t = a+2b-1)", cases);
    check(
        &mut lines,
        "negative control",
        "one decoder coefficient changed by 1 is caught",
        vec![Box::new(|| {
            let mut pk = primal(2, 2, 3, 2, 1, 1, 1);
            let row = pk.d.rows.iter().position(|r| !r.is_empty()).unwrap();
            pk.d.rows[row][0].1 = pk.d.rows[row][0].1.add(&Laurent::constant(1));
            !ok(&pk, 3)
        })],
    );

    // Lemma 2.8, 2.7
    let mut cases: Vec<Box<dyn Fn() -> bool>> = vec![];
    for (a, b, c, t) in [(2, 1, 2, 2), (1, 2, 3, 2), (2, 2, 2, 3)] {
        cases.push(Box::new(move || ok(&primal(a, b, c, t, 1, 2, 1).ex23(), 13)));
        cases.push(Box::new(move || {
            let pp = sl_to_pp(&primal(a, b, c, t, 1, 1, 1).retain(&[0, 1]));
            ok(&pp.ex12(), 14) && ok(&pp.ex13(), 15) && ok(&pp.ex23(), 16)
        }));
    }
    check(&mut lines, "Lemma 2.8 exchanges", "(23) on SL-PPs; (12), (13), (23) on PPs", cases);
    let mut cases: Vec<Box<dyn Fn() -> bool>> = vec![];
    for (a, b, c, t) in [(2, 1, 2, 2), (1, 1, 2, 1), (1, 2, 2, 3)] {
        cases.push(Box::new(move || {
            let pk = primal(a, b, c, t, 1, 1, 1);
            ok(&pk.tensor(&pk.ex23()), 17)
        }));
    }
    check(&mut lines, "Lemma 2.7 tensoring", "primal tensor dual", cases);

    // Lemma 3.2 + Claim 3.3 + exactness finding
    let mut cases: Vec<Box<dyn Fn() -> bool>> = vec![];
    for (a, b, c, t, d) in [(1, 1, 2, 1, 1), (1, 1, 2, 1, 2), (2, 1, 3, 1, 1), (1, 2, 2, 2, 2), (2, 2, 3, 2, 1), (1, 1, 3, 1, 2)] {
        cases.push(Box::new(move || {
            let pk = primal(a, b, c, t, 1, 1, 1);
            let m = pk.targets.len();
            let amp = amplify(&pk, d);
            let claim = amp.targets.len() * (d + 1) >= m.pow(2 * d as u32);
            let iso = amp.targets.iter().all(|s| s.find_isomorphism(&t_shape(a, b, d)).is_some());
            claim && iso && ok(&amp, 19)
        }));
    }
    check(&mut lines, "Lemma 3.2 amplify", "d <= 2; >= m^2d/(d+1) targets (Claim 3.3), all ~ T_{a,b}", cases);
    let mut cases: Vec<Box<dyn Fn() -> bool>> = vec![];
    for (a, b, c, t, d) in [(1, 1, 2, 1, 1), (1, 1, 2, 1, 2), (2, 1, 3, 1, 1), (1, 2, 2, 2, 2), (2, 2, 3, 2, 1), (1, 1, 3, 1, 2)] {
        cases.push(Box::new(move || amplify(&primal(a, b, c, t, 1, 1, 1), d).is_exact()));
    }
    check(&mut lines, "finding: 3.2 is exact", "the balanced amplified packing has no eps at all", cases);

    // Lemma 3.4
    let mut cases: Vec<Box<dyn Fn() -> bool>> = vec![];
    for w in 2..=4usize {
        cases.push(Box::new(move || {
            let amp = amplify(&primal(1, 1, 2, 1, 1, 1, 1), 1);
            let keep: Vec<usize> = (0..w).collect();
            let pp = sl_to_pp(&amp.retain(&keep));
            let m = pp.profile();
            let pred = cost::slpp_profile(cost::amplified_profile(1), w as f64);
            pp.sources.len() == 2 * w - 1 && ok(&pp, 23) && m[1] == (pred.er.lo as i64, pred.er.hi as i64) && m[2] == (pred.d.lo as i64, pred.d.hi as i64)
        }));
    }
    for w in 2..=3usize {
        cases.push(Box::new(move || {
            let pk = primal(2, 1, 3, 1, 1, 1, 1);
            let keep: Vec<usize> = (0..w).collect();
            ok(&sl_to_pp(&pk.retain(&keep)), 24)
        }));
    }
    check(&mut lines, "Lemma 3.4 SL -> PP", "w targets via 2w-1 sources; eps ranges match the cost model", cases);

    // Lemma 2.3, 2.9
    let mut cases: Vec<Box<dyn Fn() -> bool>> = vec![];
    for s in [Shape::full(3, 2, 1, 1, 1), Shape::mid(3, 2, 1, 2, 1), t_shape(2, 2, 1), t_shape(2, 1, 2), q_shape(2, 1, 1, 1, 1, 1)] {
        cases.push(Box::new(move || {
            let pk = evaluate(&s);
            pk.sources.len() as u128 == s.effort() && pk.is_exact() && ok(&pk, 29)
        }));
    }
    check(&mut lines, "Lemma 2.3 evaluation", "exact, with exactly eff(P) scalar multiplications", cases);
    check(
        &mut lines,
        "Lemma 2.9 composition",
        "SL->PP (mod eps) composed with exact evaluation of every source",
        vec![Box::new(|| {
            let pp = sl_to_pp(&primal(2, 1, 3, 1, 1, 1, 1).retain(&[0, 1]));
            let ev = evaluate(&pp.sources[0]);
            let all = pp.compose(&ev.repeat(pp.sources.len()));
            ok(&all, 31)
        })],
    );

    // Lemma 4.1
    let mut cases: Vec<Box<dyn Fn() -> bool>> = vec![];
    for (a, b) in [(1, 2), (2, 1), (2, 3), (3, 2)] {
        cases.push(Box::new(move || {
            let q = q_shape(a, b, 1, 1, 1, 1);
            q.effort() == ((a + b - 1) as u128).pow(6) && q.find_isomorphism(&q_shape(b, a, 1, 1, 1, 1)).is_some()
        }));
    }
    cases.push(Box::new(|| {
        let q = q_shape(1, 1, 1, 2, 3, 4);
        q.k.len() == 1 && q.k[0].iter().all(|&e| e == 0) && (q.p, q.q, q.r) == (2, 3, 4)
    }));
    check(&mut lines, "Lemma 4.1 Q_{a,b}", "eff = (a+b-1)^6d, Q_{a,b} ~ Q_{b,a}, Q_{1,1} = matmul", cases);

    // Lemma 4.2
    let mut cases: Vec<Box<dyn Fn() -> bool>> = vec![];
    for &pert in &[false, true] {
        cases.push(Box::new(move || {
            let (l, l2, d) = (1, 2, 1);
            let inp = |a, b, b2, d| if pert { perturbed(a, b, b2, d) } else { minimal_round_input(a, b, b2, d) };
            let (p1, p2) = (inp(l, l, l2, d), inp(l2, l, l2, d));
            let r = round(&p1, &p2, RoundParams { l, l2, d, w: 1, u: 1, x: 1 });
            let m = r.profile();
            let pred = cost::round_profile(cost::Profile {
                el: cost::Range::new(p1.profile()[0].0 as f64, p1.profile()[0].1 as f64),
                er: cost::Range::new(p1.profile()[1].0 as f64, p1.profile()[1].1 as f64),
                d: cost::Range::new(p1.profile()[2].0 as f64, p1.profile()[2].1 as f64),
            });
            r.sources.len() == 64 && ok(&r, 37) && m[0] == (pred.el.lo as i64, pred.el.hi as i64) && m[2] == (pred.d.lo as i64, pred.d.hi as i64)
        }));
    }
    check(&mut lines, "Lemma 4.2 round", "1 Q_{1,1} via 64 Q_{2,2}; exact and eps-perturbed inputs; eps ranges match", cases);

    // Lemma 2.4
    check(
        &mut lines,
        "Lemma 2.4 Bini",
        "mod-eps chain -> exact identity via deg+1 evaluations of eps",
        vec![Box::new(|| {
            let pk = minimal_chain(3, 1, 1, 1, &perturbed);
            let deg = degree_bound(&pk);
            let ex = bini(&pk, deg);
            ex.is_exact() && ex.sources.len() == (deg + 1) * pk.sources.len() && ok(&ex, 41)
        })],
    );

    if !quick {
        // Theorem 4.3 end to end
        check(
            &mut lines,
            "Theorem 4.3 chain",
            "2x2 integer matrices through the whole chain (c=3, d=1, k=1, w=u=1) = schoolbook",
            vec![Box::new(|| {
                let x = 2;
                let pk = minimal_chain(3, 1, 1, x, &|a, b, b2, d| minimal_round_input(a, b, b2, d));
                let mut rng = Rng(43);
                let a: Vec<i64> = (0..x * x).map(|_| rng.below(2001) as i64 - 1000).collect();
                let b: Vec<i64> = (0..x * x).map(|_| rng.below(2001) as i64 - 1000).collect();
                let la: Vec<Laurent> = a.iter().map(|&v| Laurent::int(v)).collect();
                let lb: Vec<Laurent> = b.iter().map(|&v| Laurent::int(v)).collect();
                let out = pk.run(&la, &lb);
                let good = (0..x).all(|i| (0..x).all(|j| f::to_i64(out[i * x + j].coeff(0)) == (0..x).map(|t| a[i * x + t] * b[t * x + j]).sum::<i64>()));
                println!("      A = {a:?}, B = {b:?}\n      A*B via {} scalar multiplications (schoolbook: {}) = {:?}", pk.sources.len(), x * x * x, out.iter().map(|o| f::to_i64(o.coeff(0))).collect::<Vec<_>>());
                good && pk.sources.len() == 64 * 729 * x * x * x
            })],
        );
    }

    let (p, t): (usize, usize) = lines.iter().fold((0, 0), |(p, t), l| (p + l.pass, t + l.total));
    let all = lines.iter().all(|l| l.pass == l.total);
    println!("\n{} checks, {} passed: {}", t, p, if all { "ALL PASS" } else { "FAILURES" });
    if !all {
        std::process::exit(1);
    }
}
