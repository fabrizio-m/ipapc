use crate::{
    challenges::ChallengeGenerator,
    prove::{Commitment, HidingOpening, Opening},
    utils::s_vec,
    Assert, Fr, IpaScheme, IsFalse,
};
use ark_ec::{
    short_weierstrass_jacobian::{GroupAffine, GroupProjective},
    AffineCurve, ModelParameters, ProjectiveCurve, SWModelParameters,
};
use ark_ff::{Field, One, PrimeField};
use ark_poly::{univariate::SparsePolynomial, Polynomial};

impl<P, const HIDING: bool> IpaScheme<P, HIDING>
where
    P: ModelParameters + SWModelParameters,
    Fr<P>: One,
    Commitment<P, HIDING>: Clone + Copy,
{
    pub fn verify_hiding(
        &self,
        commitment: Commitment<P, HIDING>,
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

    pub fn verify(&self, commitment: Commitment<P, HIDING>, open: Opening<P>) -> Option<Fr<P>>
    where
        Assert<HIDING>: IsFalse,
    {
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
    fn general_verify(
        &self,
        commitment: Commitment<P, HIDING>,
        point: Fr<P>,
        eval: Fr<P>,
        a: Fr<P>,
        rounds: Vec<(GroupAffine<P>, GroupAffine<P>)>,
    ) -> (GroupProjective<P>, GroupProjective<P>) {
        let u = ChallengeGenerator::inner_product_basis(&commitment, &point);
        let p = commitment.0.into_projective() + u.mul(eval);
        let mut exp = 2_u64.pow(rounds.len() as u32);

        let (final_commit, challenges, b) = rounds.iter().fold(
            (
                p,
                Vec::with_capacity(rounds.len()),
                //SPoly::<P>::new(rounds.len()),
                Fr::<P>::one(),
            ),
            |state, (lj, rj)| {
                let (p, mut challenges, s_eval) = state;
                let challenge = <ChallengeGenerator<P>>::round_challenge(lj, rj);
                let inverse = challenge.inverse().unwrap();
                let new_commit = p + (lj.mul(challenge.square()) + rj.mul(inverse.square()));
                challenges.push((challenge, inverse));

                exp = exp / 2;
                let s_eval = self.eval_term(challenge, inverse, exp, &point) * s_eval;

                (new_commit, challenges, s_eval)
            },
        );
        let s = s_vec::<P>(challenges);
        let basis = self.basis_from_s(s);
        (final_commit, basis.mul(a) + u.mul(a * b))
    }
    fn eval_term(&self, challenge: Fr<P>, inverse: Fr<P>, exp: u64, point: &Fr<P>) -> Fr<P> {
        let term =
            SparsePolynomial::from_coefficients_vec(vec![(0, inverse), (exp as usize, challenge)]);
        term.evaluate(point)
    }
    fn basis_from_s(&self, s: Vec<Fr<P>>) -> GroupAffine<P> {
        debug_assert_eq!(s.len(), self.max_degree);
        let bases = &*self.basis;
        let coeffs = s.into_iter().map(|e| e.into_repr()).collect::<Vec<_>>();
        let scalars = &*coeffs;
        let result = ark_ec::msm::VariableBaseMSM::multi_scalar_mul(bases, scalars);
        result.into_affine()
    }
}
