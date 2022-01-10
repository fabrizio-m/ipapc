use crate::{Init, IpaScheme};
use ark_ff::One;
use ark_pallas::{Fr, PallasParameters};
use ark_poly::{Polynomial, UVPolynomial};
use rand::thread_rng;

#[test]
fn test1() {
    let rng = thread_rng();
    let mut scheme = IpaScheme::<PallasParameters, _>::init(Init::Seed(1), 3, rng);
    let poly = [1, 2, 3, 4, 5, 6, 7, 8].map(Fr::from).to_vec();
    let (commit, factor) = scheme.commit(poly.clone());
    let point = Fr::from(5);
    let eval = {
        let poly = ark_poly::univariate::DensePolynomial::<Fr>::from_coefficients_slice(&*poly);
        poly.evaluate(&point)
    };
    println!("commit: {:?}", commit.0);
    println!("eval: {}", eval);
    println!();
    let proof = scheme.open(commit, factor, &poly, point, eval);
    let bad_proof = scheme.open(commit, factor, &poly, point, eval + Fr::one());
    println!("opening: {:?}", proof.a);
    //println!("opening: {:#?}", proof.rounds);
    println!();
    assert_eq!(scheme.verify(commit, proof).unwrap(), eval);
    assert!(scheme.verify(commit, bad_proof).is_none());
}
