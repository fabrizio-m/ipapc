use crate::challenges::ChallengeGenerator;
use crate::prove::{Commitment, HidingOpening, Opening};
use crate::{Assert, Fr, IpaScheme, IsFalse};
use ark_ec::short_weierstrass_jacobian::{GroupAffine, GroupProjective};
use ark_ec::{AffineCurve, ModelParameters, ProjectiveCurve, SWModelParameters};
use ark_ff::{Field, One, PrimeField};
use std::iter::repeat;
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
        struct PolySegment<P: SWModelParameters> {
            inverse: Fr<P>,
            challenge: Fr<P>,
            exp: u64,
        }

        let (final_commit, challenges, b) = rounds.iter().fold(
            (
                p,
                Vec::with_capacity(rounds.len()),
                Vec::<PolySegment<P>>::with_capacity(rounds.len()),
            ),
            |state, (lj, rj)| {
                let (p, mut challenges, mut b) = state;
                let challenge = <ChallengeGenerator<P>>::round_challenge(lj, rj);
                let inverse = challenge.inverse().unwrap();
                let new_commit = p + (lj.mul(challenge.square()) + rj.mul(inverse.square()));
                challenges.push((challenge, inverse));

                exp = exp / 2;
                b.push(PolySegment {
                    inverse,
                    challenge,
                    exp,
                });

                (new_commit, challenges, b)
            },
        );
        let s = s_vec::<P>(challenges);
        let basis = self.basis_from_s(s);
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
        (final_commit, basis.mul(a) + u.mul(a * b))
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

fn s_vec<P: SWModelParameters>(challenges: Vec<(Fr<P>, Fr<P>)>) -> Vec<Fr<P>> {
    let size = challenges.len();
    let size = 2_usize.pow(size as u32) as usize;
    let mut challenges = challenges
        .into_iter()
        .enumerate()
        .map(|(i, (challenge, inverse))| {
            let segment_size = size / (2_usize.pow(i as u32 + 1));
            let challenge_segment = repeat(challenge).take(segment_size);
            let inverse_segment = repeat(inverse).take(segment_size);
            let combined = inverse_segment.chain(challenge_segment);
            combined.cycle()
        })
        .collect::<Vec<_>>();
    let f = || {
        let elem = challenges
            .iter_mut()
            .filter_map(|iter| iter.next())
            .reduce(Mul::mul);
        elem
    };
    let s = std::iter::from_fn(f).take(size).collect::<Vec<_>>();
    debug_assert_eq!(size, s.len());
    s
}
