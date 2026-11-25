use ninefourths::cost::*;

fn main() {
    let tau = 7f64.log2();
    let l10 = |x: f64| x * std::f64::consts::LOG10_2;
    if std::env::var("DBG").is_ok() {
        let c = (3..10_000).find(|&c| threshold(c) < tau).unwrap();
        let d = (1..100_000).find(|&d| condition_10(c, d, tau));
        println!("c={c} d={d:?}");
        let d = d.unwrap();
        let mut k = 1;
        loop { let (_,_,_,lv,lr,ld) = mod_eps_cost(c,d,k,None).unwrap(); if lr < tau*lv { println!("k={k} lv={lv} lr={lr} gap={} ld={ld}", tau*lv-lr); break;} k+=1; }
    }
    match faithful(tau) {
        Some((p, c)) => {
            println!("faithful tau=log2 7: c={} d={} k={} h={}", p.c, p.d, p.k, p.h);
            println!("  mod-eps: log10 v={:.6e} log10 r={:.6e} exp={:.8} log10 deg={:.6e}", l10(c.log2_v), l10(c.log2_r), c.exponent_mod_eps(), l10(c.log2_deg));
            println!("  exact:   log10 V={:.6e} log10 R={:.6e} exp={:.8} floor log10 n={:.6e}", l10(c.log2_vv), l10(c.log2_rr), c.exponent_exact(), l10(c.log2_strassen_floor()));
        }
        None => println!("faithful: none"),
    }
    let cm: u64 = std::env::var("CMAX").ok().and_then(|s| s.parse().ok()).unwrap_or(40);
    let dm: u64 = std::env::var("DMAX").ok().and_then(|s| s.parse().ok()).unwrap_or(400);
    if let Some(b) = search(tau, cm, dm) {
        let (p, c) = (b.params, b.cost);
        println!("searched: c={} d={} k={} h={}", p.c, p.d, p.k, p.h);
        println!("  mod-eps: log10 v={:.6e} log10 r={:.6e} exp={:.8} log10 deg={:.6e}", l10(c.log2_v), l10(c.log2_r), c.exponent_mod_eps(), l10(c.log2_deg));
        println!("  exact:   log10 V={:.6e} log10 R={:.6e} exp={:.8} floor log10 n={:.6e}", l10(c.log2_vv), l10(c.log2_rr), c.exponent_exact(), l10(c.log2_strassen_floor()));
    }
}
#[allow(dead_code)]
fn dbg() {}
