use crate::{Init, IpaScheme};
use ark_pallas::{Fr, PallasParameters};
use ark_poly::{Polynomial, UVPolynomial};

#[test]
fn test1() {
    let scheme = IpaScheme::<PallasParameters>::init(Init::Seed(1), 3);
    let poly = [1, 2, 3, 4, 5, 6, 7, 8].map(Fr::from).to_vec();
    let commit = scheme.commit(poly.clone());
    let point = Fr::from(5);
    let eval = {
        let poly = ark_poly::univariate::DensePolynomial::<Fr>::from_coefficients_slice(&*poly);
        poly.evaluate(&point)
    };
    println!("commit: {:?}", commit.0);
    println!("eval: {}", eval);
    println!();
    let proof = scheme.open(commit, &poly, point, eval);
    println!("opening: {:?}", proof.a);
    //println!("opening: {:#?}", proof.rounds);
    println!();
    assert_eq!(scheme.verify(commit, proof).unwrap(), eval);
}
