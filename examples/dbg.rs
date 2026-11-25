use ninefourths::{field::Rng, gadgets::*};
fn main() {
    let mut rng = Rng(1);
    let (a,b)=(t_shape(2,1,1), t_shape(2,2,1));
    println!("{:?}\n{:?}", a, b);
    println!("{:?}", embedding(&a,&b).verify(&mut rng,1));
}
