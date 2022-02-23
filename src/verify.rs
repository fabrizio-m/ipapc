use crate::challenges::ChallengeGenerator;
use crate::prove::{Commitment, HidingOpening, Opening};
use crate::utils::{compress_basis, split};
use crate::{Assert, Fr, IpaScheme, IsFalse};
use ark_ec::{AffineCurve, ModelParameters, SWModelParameters};
use ark_ff::{Field, One};
use std::ops::Mul;

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
        let u = ChallengeGenerator::inner_product_basis(&commitment, &point);
        let p = commitment.0.into_projective() + u.mul(eval);
        let basis = self.basis.clone();
        let mut exp = 2_u64.pow(rounds.len() as u32);

        let (final_commit, basis, b) =
            rounds
                .iter()
                .fold((p, basis, Fr::<P>::one()), |state, (lj, rj)| {
                    let (p, basis, b) = state;
                    let challenge = <ChallengeGenerator<P>>::round_challenge(lj, rj);
                    let inverse = challenge.inverse().unwrap();
                    let new_commit = p + (lj.mul(challenge.square()) + rj.mul(inverse.square()));

                    let (g_l, g_r) = split(&*basis);
                    let basis = compress_basis(g_l, g_r, challenge);

                    exp = exp / 2;
                    let new_b = inverse + challenge * point.pow([exp]);
                    let b = new_b * b;

                    (new_commit, basis, b)
                });
        let final_check = final_commit
            == basis[0].mul(a) + u.mul(a * b) + self.blinding_basis.mul(blinding_factor);
        if final_check {
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
        let u = ChallengeGenerator::inner_product_basis(&commitment, &point);
        let p = commitment.0.into_projective() + u.mul(eval);
        let basis = self.basis.clone();
        let mut exp = 2_u64.pow(rounds.len() as u32);
        struct PolySegment<P: SWModelParameters> {
            inverse: Fr<P>,
            challenge: Fr<P>,
            exp: u64,
        }

        let (final_commit, basis, b) = rounds.iter().fold(
            (p, basis, Vec::<PolySegment<P>>::with_capacity(rounds.len())),
            |state, (lj, rj)| {
                let (p, basis, mut b) = state;
                let challenge = <ChallengeGenerator<P>>::round_challenge(lj, rj);
                let inverse = challenge.inverse().unwrap();
                let new_commit = p + (lj.mul(challenge.square()) + rj.mul(inverse.square()));

                let (g_l, g_r) = split(&*basis);
                let basis = compress_basis(g_l, g_r, challenge);

                exp = exp / 2;
                b.push(PolySegment {
                    inverse,
                    challenge,
                    exp,
                });

                (new_commit, basis, b)
            },
        );
        let b = b
            .into_iter()
            .map(|segment| {
                let PolySegment {
                    inverse,
                    challenge,
                    exp,
                } = segment;
                inverse + challenge * point.pow([exp])
            })
            .reduce(Mul::mul)
            .unwrap();
        let final_check = final_commit == basis[0].mul(a) + u.mul(a * b);
        if final_check {
            Some(eval)
        } else {
            None
        }
    }
}
