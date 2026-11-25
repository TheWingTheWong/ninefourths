//! Prints the cost of Pass's construction: the faithful parameter choices and
//! a grid search over the same lemmas' free parameters, for several target
//! exponents. `--json` writes one JSON object per target exponent.
use ninefourths::cost::*;

fn l10(x: f64) -> f64 {
    x * std::f64::consts::LOG10_2
}

fn row(kind: &str, tau: f64, p: Params, c: Cost, json: bool) {
    if json {
        println!(
            "{{\"kind\":\"{kind}\",\"tau\":{tau:.6},\"c\":{},\"d\":{},\"k\":{},\"h\":{:.6e},\"log10_v\":{:.6},\"log10_r\":{:.6},\"exp_mod_eps\":{:.8},\"log10_deg\":{:.6},\"log10_V\":{:.6},\"log10_R\":{:.6},\"exp_exact\":{:.8},\"log10_strassen_floor\":{:.6}}}",
            p.c, p.d, p.k, p.h, l10(c.log2_v), l10(c.log2_r), c.exponent_mod_eps(), l10(c.log2_deg), l10(c.log2_vv), l10(c.log2_rr), c.exponent_exact(), l10(c.log2_strassen_floor())
        );
    } else {
        println!("{kind:>9} tau<{tau:.4}  c={:<3} d={:<5} k={:<6} h={:<11.4e} mod-eps side 10^{:<11.0} exp {:.4}  eps-degree 10^{:<11.0} exact side V = 10^{:.4e}  R = 10^{:.4e}",
            p.c, p.d, p.k, p.h, l10(c.log2_v), c.exponent_mod_eps(), l10(c.log2_deg), l10(c.log2_vv), l10(c.log2_rr));
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let json = args.iter().any(|a| a == "--json");
    let env = |k: &str, d: u64| std::env::var(k).ok().and_then(|s| s.parse().ok()).unwrap_or(d);
    let (cm, dm) = (env("CMAX", 8), env("DMAX", 500));
    let taus: Vec<f64> = match std::env::var("TAUS") {
        Ok(s) => s.split(',').map(|x| if x == "strassen" { 7f64.log2() } else { x.parse().unwrap() }).collect(),
        Err(_) => vec![7f64.log2()],
    };
    for tau in taus {
        let mode = std::env::var("MODE").unwrap_or_default();
        if mode != "search" {
            if let Some((p, c)) = faithful(tau) {
                row("faithful", tau, p, c, json);
            }
        }
        if mode == "faithful" {
            continue;
        }
        if let Some(b) = search(tau, cm, dm) {
            row("searched", tau, b.params, b.cost, json);
        }
    }
}
