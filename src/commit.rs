use crate::{CoeffsOrEvals, Commitment, Fr, IpaScheme, UnsafeHidingCommitment};
use ark_ec::{AffineCurve, ProjectiveCurve, SWModelParameters};
use ark_ff::UniformRand;
use rand::Rng;
use std::ops::{Add, Mul, Sub};

pub trait CommitmentTrait<P, R>
where
    P: SWModelParameters,
    Self: Add<Self> + Sub<Self> + Mul<Fr<P>> + Sized,
    R: Rng,
{
    fn commit(scheme: &IpaScheme<P, R>, poly: impl Into<CoeffsOrEvals<P>>) -> Self;
}

impl<P, R> CommitmentTrait<P, R> for Commitment<P, false>
where
    P: SWModelParameters,
    R: Rng,
{
    fn commit(scheme: &IpaScheme<P, R>, poly: impl Into<CoeffsOrEvals<P>>) -> Self {
        let commitment = scheme.commit_simple(poly);
        Self(commitment.into_affine())
    }
}

impl<P, R> CommitmentTrait<P, R> for UnsafeHidingCommitment<P>
where
    P: SWModelParameters,
    R: Rng,
{
    fn commit(scheme: &IpaScheme<P, R>, poly: impl Into<CoeffsOrEvals<P>>) -> Self {
        let blinding_factor = {
            let rng = &mut *(scheme.rng.borrow_mut());
            Fr::<P>::rand(rng)
        };
        let commitment = scheme.commit_simple(poly);
        let commitment = commitment + scheme.blinding_basis.mul(blinding_factor);

        Self(commitment.into_affine(), blinding_factor)
    }
}
