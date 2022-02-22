use ark_ec::{short_weierstrass_jacobian::GroupAffine, AffineCurve};
use ark_pallas::PallasParameters;
use ipapc::{Init, IpaScheme};
use rand::{thread_rng, Rng};
use std::iter::repeat;

type Fr<P> = <GroupAffine<P> as AffineCurve>::ScalarField;
pub fn commit_iai() {
    const SIZE: u8 = 8;
    let scheme = IpaScheme::<PallasParameters, false>::init(Init::Seed(1), SIZE);
    let mut rng = thread_rng();
    //let poly: [Fr<PallasParameters>; 2_usize.pow(SIZE as u32)] = rng.gen();
    let poly: Vec<Fr<PallasParameters>> = repeat(())
        .map(|_| rng.gen())
        .take(2_usize.pow(SIZE as u32))
        .collect();
    //panic!();
    let coeffs = iai::black_box(poly.to_vec());
    let _ = scheme.commit(coeffs.clone());
}
iai::main!(commit_iai);
