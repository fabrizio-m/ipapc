use crate::{
    challenges::ChallengeGenerator,
    prove::{Commitment, HidingOpening, Opening},
    Fr, IpaScheme,
};
use ark_ec::{
    short_weierstrass_jacobian::{GroupAffine, GroupProjective},
    AffineCurve, ModelParameters, ProjectiveCurve, SWModelParameters,
};
use ark_ff::{Field, One, PrimeField};
use ark_poly::{
    univariate::{DenseOrSparsePolynomial, DensePolynomial, SparsePolynomial},
    Polynomial,
};
use std::ops::Mul;

impl<P> IpaScheme<P>
where
    P: ModelParameters + SWModelParameters,
    Fr<P>: One,
{
    pub fn verify_hiding(
        &self,
        commitment: Commitment<P, true>,
        open: HidingOpening<P>,
    ) -> Option<Fr<P>> {
        let HidingOpening::<P> {
            point,
            eval,
            a,
            rounds,
            blinding_factor,
        } = open;
        let (final_commit, check) = self.general_verify(commitment, point, eval, a, rounds);
        if final_commit == check + self.blinding_basis.mul(blinding_factor) {
            Some(eval)
        } else {
            None
        }
    }

    pub fn verify(&self, commitment: Commitment<P, false>, open: Opening<P>) -> Option<Fr<P>> {
        let Opening::<P> {
            point,
            eval,
            a,
            rounds,
        } = open;
        let (final_commit, check) = self.general_verify(commitment, point, eval, a, rounds);
        if final_commit == check {
            Some(eval)
        } else {
            None
        }
    }
    pub(crate) fn general_verify<const HIDING: bool>(
        &self,
        commitment: Commitment<P, HIDING>,
        point: Fr<P>,
        eval: Fr<P>,
        a: Fr<P>,
        rounds: Vec<(GroupAffine<P>, GroupAffine<P>)>,
    ) -> (GroupProjective<P>, GroupProjective<P>) {
        let u = ChallengeGenerator::inner_product_basis(&commitment, &point);

        let (final_commit, b_poly) = Self::process_rounds(commitment, point, eval, rounds);
        let b = Self::eval_b_poly(&b_poly, point);
        let s = Self::sparse_to_dense(b_poly).coeffs;
        let basis = Self::s_to_basis(&self, s);

        (final_commit, basis.mul(a) + u.mul(a * b))
    }
    fn sparse_to_dense(polys: Vec<SparsePolynomial<Fr<P>>>) -> DensePolynomial<Fr<P>> {
        let poly = polys.into_iter().reduce(|a, b| a.mul(&b)).unwrap();
        DenseOrSparsePolynomial::from(poly).into()
    }
    pub(crate) fn eval_b_poly(b_poly: &Vec<SparsePolynomial<Fr<P>>>, point: Fr<P>) -> Fr<P> {
        b_poly
            .iter()
            .map(|poly| poly.evaluate(&point))
            .reduce(Mul::mul)
            .unwrap()
    }
    /// compute
    /// final commitment
    /// b_poly
    pub(crate) fn process_rounds<const HIDING: bool>(
        commitment: Commitment<P, HIDING>,
        point: Fr<P>,
        eval: Fr<P>,
        rounds: Vec<(GroupAffine<P>, GroupAffine<P>)>,
    ) -> (GroupProjective<P>, Vec<SparsePolynomial<Fr<P>>>) {
        let u = ChallengeGenerator::inner_product_basis(&commitment, &point);
        let p = commitment.0.into_projective() + u.mul(eval);
        let mut exp = 2_u64.pow(rounds.len() as u32);

        let b_poly = Vec::with_capacity(rounds.len());
        let (final_commit, b_poly) = rounds.iter().fold((p, b_poly), |state, (lj, rj)| {
            let (p, mut b_poly) = state;
            let challenge = <ChallengeGenerator<P>>::round_challenge(lj, rj);
            let inverse = challenge.inverse().unwrap();
            let new_commit = p + (lj.mul(challenge.square()) + rj.mul(inverse.square()));

            exp = exp / 2;
            let term = SparsePolynomial::from_coefficients_vec(vec![
                (0, inverse),
                (exp as usize, challenge),
            ]);
            b_poly.push(term);

            (new_commit, b_poly)
        });
        (final_commit, b_poly)
    }
    fn s_to_basis(&self, s: Vec<Fr<P>>) -> GroupAffine<P> {
        debug_assert_eq!(s.len(), self.max_degree);
        let basis = &*self.basis;
        let coeffs = s.into_iter().map(|e| e.into_repr()).collect::<Vec<_>>();
        let scalars = &*coeffs;
        let result = ark_ec::msm::VariableBaseMSM::multi_scalar_mul(basis, scalars);
        result.into_affine()
    }
}
