use ninefourths::{field::Rng, gadgets::*};
fn main() {
    let mut rng = Rng(1);
    for (a,b,c,t,d) in [(2,1,3,1,1),(1,1,2,1,2),(1,1,3,1,2),(2,2,3,2,1),(1,2,2,2,2)] {
        let amp = amplify(&primal(a,b,c,t,1,1,1), d);
        println!("a={a} b={b} c={c} t={t} d={d}: targets={} sources={} profile={:?} exact={} verify={:?}", amp.targets.len(), amp.sources.len(), amp.profile(), amp.is_exact(), amp.verify(&mut rng,1).map(|v|v.out_degree));
    }
}
