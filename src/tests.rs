use crate::{Commitment, Fr, HidingOpening, Init, IpaScheme, Opening, UnsafeHidingCommitment};
use ark_ec::SWModelParameters;
use ark_ff::One;
use ark_pallas::{Fr as F, PallasParameters};
use ark_poly::{Polynomial, UVPolynomial};
use rand::{thread_rng, Rng};

#[test]
fn test_hiding() {
    let scheme = IpaScheme::<PallasParameters, _>::init(Init::Seed(1), 3, true, thread_rng());
    let poly = [1, 2, 3, 4, 5, 6, 7, 8].map(F::from).to_vec();
    //let commit = scheme.commit_hiding(poly.clone(), &mut rng);
    let commit: UnsafeHidingCommitment<_> = scheme.commit(poly.clone());
    let point = F::from(5);
    let eval = {
        let poly = ark_poly::univariate::DensePolynomial::<F>::from_coefficients_slice(&*poly);
        poly.evaluate(&point)
    };
    //let proof = scheme.open_hiding(commit.into(), &poly, point, eval, &mut rng);
    let proof: HidingOpening<_> = scheme.open(commit, &poly, point, eval);
    let bad_proof: HidingOpening<_> = scheme.open(commit.into(), &poly, point, eval + F::one());
    assert_eq!(scheme.verify(commit.into(), proof).unwrap(), eval);
    assert!(scheme.verify(commit.into(), bad_proof).is_none());
}

#[test]
fn test_binding() {
    let scheme = IpaScheme::<PallasParameters, _>::init(Init::Seed(1), 3, true, thread_rng());
    let poly = [1, 2, 3, 4, 5, 6, 7, 8].map(F::from).to_vec();
    let commit = scheme.commit(poly.clone());
    let point = F::from(5);
    let eval = {
        let poly = ark_poly::univariate::DensePolynomial::<F>::from_coefficients_slice(&*poly);
        poly.evaluate(&point)
    };
    let proof: Opening<_> = scheme.open(commit, &poly, point, eval);
    let bad_proof: Opening<_> = scheme.open(commit, &poly, point, eval + F::one());
    assert_eq!(scheme.verify(commit, proof).unwrap(), eval);
    assert!(scheme.verify(commit, bad_proof).is_none());
}

pub(crate) fn commit_and_open<P, R>(
    scheme: &IpaScheme<P, R>,
) -> (Commitment<P, false>, Vec<Fr<P>>, Fr<P>, Fr<P>)
where
    P: SWModelParameters,
    Fr<P>: From<i32>,
    R: Rng,
{
    use ark_ff::UniformRand;

    type Poly<P> = ark_poly::univariate::DensePolynomial<Fr<P>>;
    let mut rng = thread_rng();
    let poly = [0; 256].map(|_| Fr::<P>::rand(&mut rng)).to_vec();
    let commit = scheme.commit(poly.clone());
    let point = Fr::<P>::from(5);
    let eval = {
        let poly = Poly::<P>::from_coefficients_slice(&*poly);
        poly.evaluate(&point)
    };
    (commit, poly, point, eval)
}
