use crate::{
    challenges::ChallengeGenerator,
    open::OpenTrait,
    utils::{compress, compress_basis, inner_product, scalar_inner_product, split},
    Fr, IpaScheme,
};
use ark_ec::{
    short_weierstrass_jacobian::GroupAffine, AffineCurve, ModelParameters, ProjectiveCurve,
    SWModelParameters,
};
use ark_ff::{Field, One, ToBytes, UniformRand};
use rand::Rng;
use std::fmt::Debug;

#[derive(Debug, Clone)]
pub struct Opening<P: SWModelParameters> {
    pub(crate) point: Fr<P>,
    pub(crate) eval: Fr<P>,
    pub(crate) rounds: Vec<(GroupAffine<P>, GroupAffine<P>)>,
    pub(crate) a: Fr<P>,
}
pub(crate) struct HidingRoundOutput<P: SWModelParameters> {
    pub(crate) lj: GroupAffine<P>,
    pub(crate) rj: GroupAffine<P>,
    a: Vec<Fr<P>>,
    b: Vec<Fr<P>>,
    basis: Vec<GroupAffine<P>>,
    blind: Fr<P>,
    challenges: Option<Vec<(Fr<P>, Fr<P>)>>,
}
pub(crate) struct RoundOutput<P: SWModelParameters> {
    pub(crate) lj: GroupAffine<P>,
    pub(crate) rj: GroupAffine<P>,
    a: Vec<Fr<P>>,
    b: Vec<Fr<P>>,
    basis: Vec<GroupAffine<P>>,
    challenges: Option<Vec<(Fr<P>, Fr<P>)>>,
}
#[derive(Clone, Copy, PartialEq, Eq)]
pub struct Commitment<T: SWModelParameters, const HIDING: bool>(pub(crate) GroupAffine<T>)
where
    GroupAffine<T>: Debug;

impl<T: SWModelParameters, const HIDING: bool> Debug for Commitment<T, HIDING>
where
    GroupAffine<T>: Debug,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_tuple("Commitment").field(&self.0).finish()
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct UnsafeHidingCommitment<T: SWModelParameters>(
    pub(crate) GroupAffine<T>,
    pub(crate) Fr<T>,
)
where
    GroupAffine<T>: Debug;

impl<T: SWModelParameters> UnsafeHidingCommitment<T>
where
    GroupAffine<T>: Debug,
{
    pub fn clean(self) -> Commitment<T, true> {
        Commitment(self.0)
    }
}
#[derive(Debug)]
pub struct HidingOpening<P: SWModelParameters> {
    pub(crate) point: Fr<P>,
    pub(crate) eval: Fr<P>,
    pub(crate) rounds: Vec<(GroupAffine<P>, GroupAffine<P>)>,
    pub(crate) a: Fr<P>,
    pub(crate) blinding_factor: Fr<P>,
}

impl<P: SWModelParameters> From<UnsafeHidingCommitment<P>> for Commitment<P, true> {
    fn from(unsafe_commitment: UnsafeHidingCommitment<P>) -> Self {
        Commitment(unsafe_commitment.0)
    }
}
impl<P, const HIDING: bool> From<Commitment<P, HIDING>> for Vec<u8>
where
    P: SWModelParameters,
{
    fn from(commit: Commitment<P, HIDING>) -> Self {
        let mut bytes = Vec::new();
        commit.0.write(&mut bytes).unwrap();
        bytes
    }
}

impl<P, R> IpaScheme<P, R>
where
    P: ModelParameters + SWModelParameters,
    Fr<P>: One,
    R: Rng,
{
    pub fn open<O>(&self, commitment: O::Commit, a: &[Fr<P>], point: Fr<P>, eval: Fr<P>) -> O
    where
        O: OpenTrait<P, R>,
    {
        O::open(self, commitment, a, point, eval)
    }

    pub(crate) fn open_recursive(
        prev: RoundOutput<P>,
        mut rounds: Vec<(GroupAffine<P>, GroupAffine<P>)>,
        point: Fr<P>,
        eval: Fr<P>,
        u: GroupAffine<P>,
    ) -> (Opening<P>, Option<Vec<(Fr<P>, Fr<P>)>>, GroupAffine<P>) {
        let RoundOutput {
            a,
            b,
            basis,
            challenges,
            ..
        } = prev;
        if a.len().is_one() {
            (
                Opening::<P> {
                    a: a[0],
                    rounds,
                    point,
                    eval,
                },
                challenges,
                basis[0],
            )
        } else {
            let prev = Self::round(&*basis, &*a, &*b, u, challenges);
            rounds.push((prev.lj, prev.rj));
            Self::open_recursive(prev, rounds, point, eval, u)
        }
    }
    pub(crate) fn round(
        basis: &[GroupAffine<P>],
        a: &[Fr<P>],
        b: &[Fr<P>],
        u: GroupAffine<P>,
        challenges: Option<Vec<(Fr<P>, Fr<P>)>>,
    ) -> RoundOutput<P> {
        assert_eq!(basis.len(), a.len());
        assert_eq!(basis.len(), b.len());
        assert!(basis.len() > 1);

        let ([lj, rj], [a, b], basis, blind, challenges) =
            Self::general_round(basis, a, b, u, None, None, None, challenges);
        debug_assert!(blind.is_none());
        RoundOutput {
            lj,
            rj,
            a,
            b,
            basis,
            challenges,
        }
    }

    fn general_round(
        basis: &[GroupAffine<P>],
        a: &[Fr<P>],
        b: &[Fr<P>],
        u: GroupAffine<P>,
        blinding_basis: Option<GroupAffine<P>>,
        blinding_factors: Option<[Fr<P>; 2]>,
        blind: Option<Fr<P>>,
        challenges: Option<Vec<(Fr<P>, Fr<P>)>>,
    ) -> (
        [GroupAffine<P>; 2],
        [Vec<Fr<P>>; 2],
        Vec<GroupAffine<P>>,
        Option<Fr<P>>,
        Option<Vec<(Fr<P>, Fr<P>)>>,
    ) {
        assert_eq!(basis.len(), a.len());
        assert_eq!(basis.len(), b.len());
        assert!(basis.len() > 1);
        let (a_l, a_r) = split(a);
        let (b_l, b_r) = split(b);
        let (g_l, g_r) = split(basis);

        let lj = inner_product(g_r, a_l) + u.mul(scalar_inner_product::<P>(a_l, b_r));
        let rj = inner_product(g_l, a_r) + u.mul(scalar_inner_product::<P>(a_r, b_l));
        let (lj, rj, factors) = match blinding_factors {
            Some(factors) => {
                let basis = blinding_basis.unwrap();
                (
                    lj + basis.mul(factors[0]),
                    rj + basis.mul(factors[1]),
                    Some(factors),
                )
            }
            None => (lj, rj, None),
        };

        let [lj, rj] = [lj, rj].map(|point| point.into_affine());

        let challenge = <ChallengeGenerator<P>>::round_challenge(&lj, &rj);
        let blind = factors.map(|[blind_l, blind_r]| {
            challenge.square() * blind_l
                + blind.unwrap()
                + challenge.inverse().unwrap().square() * blind_r
        });
        let inverse = challenge.inverse().unwrap();
        let a = compress::<P>(a_r, a_l, challenge, inverse);
        let b = compress::<P>(b_l, b_r, challenge, inverse);
        let basis = compress_basis(g_l, g_r, challenge);
        let challenges = challenges.map(|mut challenges| {
            challenges.push((challenge, inverse));
            challenges
        });
        ([lj, rj], [a, b], basis, blind, challenges)
    }

    pub(crate) fn open_recursive_hiding(
        prev: HidingRoundOutput<P>,
        mut rounds: Vec<(GroupAffine<P>, GroupAffine<P>)>,
        point: Fr<P>,
        eval: Fr<P>,
        u: GroupAffine<P>,
        blinding_basis: GroupAffine<P>,
        rng: &mut impl Rng,
    ) -> HidingOpening<P> {
        let HidingRoundOutput {
            a,
            b,
            basis,
            blind,
            challenges,
            ..
        } = prev;
        if a.len().is_one() {
            HidingOpening::<P> {
                a: a[0],
                rounds,
                point,
                eval,
                blinding_factor: blind,
            }
        } else {
            let prev =
                Self::hiding_round(&*basis, &*a, &*b, u, blinding_basis, rng, blind, challenges);
            rounds.push((prev.lj, prev.rj));
            Self::open_recursive_hiding(prev, rounds, point, eval, u, blinding_basis, rng)
        }
    }

    pub(crate) fn hiding_round(
        basis: &[GroupAffine<P>],
        a: &[Fr<P>],
        b: &[Fr<P>],
        u: GroupAffine<P>,
        blinding_basis: GroupAffine<P>,
        rng: &mut impl Rng,
        blind: Fr<P>,
        challenges: Option<Vec<(Fr<P>, Fr<P>)>>,
    ) -> HidingRoundOutput<P> {
        assert_eq!(basis.len(), a.len());
        assert_eq!(basis.len(), b.len());
        assert!(basis.len() > 1);

        let blinding_factors = [(); 2].map(|_| Fr::<P>::rand(rng));
        let ([lj, rj], [a, b], basis, blind, challenges) = Self::general_round(
            basis,
            a,
            b,
            u,
            Some(blinding_basis),
            Some(blinding_factors),
            Some(blind),
            challenges,
        );
        let blind = blind.unwrap();
        HidingRoundOutput {
            lj,
            rj,
            a,
            b,
            basis,
            blind,
            challenges,
        }
    }
}
