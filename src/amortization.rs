use crate::{challenges::ChallengeGenerator, Assert, Commitment, Fr, IpaScheme, IsFalse, Opening};
use ark_ec::{
    short_weierstrass_jacobian::{GroupAffine, GroupProjective},
    AffineCurve, ModelParameters, ProjectiveCurve, SWModelParameters,
};
use ark_ff::{One, Zero};
use ark_poly::{
    univariate::{DenseOrSparsePolynomial, DensePolynomial, SparsePolynomial},
    Polynomial, UVPolynomial,
};
use std::ops::Mul;

pub struct MultiOpening<P: SWModelParameters> {
    openings: Vec<(Opening<P>, GroupAffine<P>)>,
    batch_opening: Opening<P>,
}

impl<P, const HIDING: bool> IpaScheme<P, HIDING>
where
    P: ModelParameters + SWModelParameters,
    Fr<P>: One,
    Commitment<P, HIDING>: Clone + Copy,
{
    pub fn batch_open(
        &self,
        opens: Vec<(Commitment<P, HIDING>, &[Fr<P>], Fr<P>, Fr<P>)>,
    ) -> MultiOpening<P>
    where
        Assert<HIDING>: IsFalse,
    {
        let len = opens.len();
        let (openings, amortization, challenges) = opens
            .into_iter()
            .map(|open| {
                let (commitment, a, point, eval) = open;
                let u = ChallengeGenerator::inner_product_basis(&commitment, &point);
                let basis = &*self.basis;
                let b = self.b(point);
                //let mut rng = &self.rng;
                let first = Self::round(basis, a, &*b, u, Some(vec![]));
                let rounds = vec![(first.lj, first.rj)];
                let open = Self::open_recursive(first, rounds, point, eval, u);
                open
            })
            .fold(
                (
                    Vec::with_capacity(len),
                    Vec::with_capacity(len),
                    ChallengeGenerator::new(),
                ),
                |mut state, open| {
                    let (open, s_challenges, basis) = open;
                    state.0.push(open);
                    state.1.push((s_challenges, basis));
                    state.2.digest_for_amortization(basis);
                    state
                },
            );
        let (combinations, point) = challenges.amortization_elements();
        let mut bs = vec![];
        let (_, s_poly, commitment) = amortization.into_iter().fold(
            (
                Fr::<P>::one(),
                DensePolynomial::from_coefficients_slice(&[]),
                GroupProjective::zero(),
            ),
            |acc, val| {
                let (combination, s_poly, commitment) = acc;
                let (s_challenges, basis) = val;
                let s_poly = Self::challenges_to_poly(s_challenges.unwrap(), combination) + s_poly;
                let commitment = basis.mul(combination) + commitment;
                bs.push(basis);

                (combination * combinations, s_poly, commitment)
            },
        );
        let eval = s_poly.evaluate(&point);
        let commitment = Commitment(commitment.into_affine());
        let batch_opening = self.open(commitment, &*s_poly.coeffs, point, eval);
        debug_assert_eq!(openings.len(), bs.len());
        let openings = openings.into_iter().zip(bs.into_iter()).collect();
        MultiOpening {
            openings,
            batch_opening,
        }
    }
    fn challenges_to_poly(
        challenges: Vec<(Fr<P>, Fr<P>)>,
        combination_element: Fr<P>,
    ) -> DensePolynomial<Fr<P>> {
        let powers = std::iter::successors(Some(1_usize), |a| Some(a * 2));
        let poly = challenges
            .into_iter()
            .rev()
            .zip(powers)
            .map(|((challenge, inverse), exp)| {
                SparsePolynomial::from_coefficients_vec(vec![(0, inverse), (exp, challenge)])
            })
            .reduce(|a, b| a.mul(&b))
            .unwrap();
        let poly: DensePolynomial<_> = DenseOrSparsePolynomial::from(poly).into();
        poly.mul(combination_element)
    }
}
