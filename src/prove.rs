use crate::{
    challenges::ChallengeGenerator,
    utils::{compress, compress_basis, inner_product, scalar_inner_product, split},
    Assert, Fr, IpaScheme, IsFalse, IsTrue,
};
use ark_ec::{
    short_weierstrass_jacobian::GroupAffine, AffineCurve, ModelParameters, ProjectiveCurve,
    SWModelParameters,
};
use ark_ff::{Field, One, UniformRand};
use rand::Rng;
use std::fmt::Debug;

#[derive(Debug)]
pub struct Opening<P: SWModelParameters> {
    pub(crate) point: Fr<P>,
    pub(crate) eval: Fr<P>,
    pub(crate) rounds: Vec<(GroupAffine<P>, GroupAffine<P>)>,
    pub(crate) a: Fr<P>,
}
struct HidingRoundOutput<P: SWModelParameters> {
    lj: GroupAffine<P>,
    rj: GroupAffine<P>,
    a: Vec<Fr<P>>,
    b: Vec<Fr<P>>,
    basis: Vec<GroupAffine<P>>,
    blind: Fr<P>,
}
struct RoundOutput<P: SWModelParameters> {
    lj: GroupAffine<P>,
    rj: GroupAffine<P>,
    a: Vec<Fr<P>>,
    b: Vec<Fr<P>>,
    basis: Vec<GroupAffine<P>>,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Commitment<T: SWModelParameters, const HIDING: bool>(pub(crate) GroupAffine<T>)
where
    GroupAffine<T>: Debug;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct UnsafeHidingCommitment<T: SWModelParameters>(
    pub(crate) GroupAffine<T>,
    pub(crate) Fr<T>,
)
where
    GroupAffine<T>: Debug;
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

impl<P, const HIDING: bool> IpaScheme<P, HIDING>
where
    P: ModelParameters + SWModelParameters,
    Fr<P>: One,
    Commitment<P, HIDING>: Clone + Copy,
{
    pub fn open(
        &self,
        commitment: Commitment<P, HIDING>,
        a: &[Fr<P>],
        point: Fr<P>,
        eval: Fr<P>,
    ) -> Opening<P>
    where
        Assert<HIDING>: IsFalse,
    {
        let u = ChallengeGenerator::inner_product_basis(&commitment, &point);
        let basis = &*self.basis;
        let b = self.b(point);
        //let mut rng = &self.rng;
        let first = Self::round(basis, a, &*b, u);
        let rounds = vec![(first.lj, first.rj)];
        let opening = Self::open_recursive(first, rounds, point, eval, u);
        opening
    }

    fn open_recursive(
        prev: RoundOutput<P>,
        mut rounds: Vec<(GroupAffine<P>, GroupAffine<P>)>,
        point: Fr<P>,
        eval: Fr<P>,
        u: GroupAffine<P>,
    ) -> Opening<P>
    where
        Assert<HIDING>: IsFalse,
    {
        let RoundOutput { a, b, basis, .. } = prev;
        if a.len().is_one() {
            Opening::<P> {
                a: a[0],
                rounds,
                point,
                eval,
            }
        } else {
            let prev = Self::round(&*basis, &*a, &*b, u);
            rounds.push((prev.lj, prev.rj));
            Self::open_recursive(prev, rounds, point, eval, u)
        }
    }
    fn round(
        basis: &[GroupAffine<P>],
        a: &[Fr<P>],
        b: &[Fr<P>],
        u: GroupAffine<P>,
    ) -> RoundOutput<P> {
        assert_eq!(basis.len(), a.len());
        assert_eq!(basis.len(), b.len());
        assert!(basis.len() > 1);

        let ([lj, rj], [a, b], basis, blind) =
            Self::general_round(basis, a, b, u, None, None, None);
        debug_assert!(blind.is_none());
        RoundOutput {
            lj,
            rj,
            a,
            b,
            basis,
        }
    }

    fn round_old(
        basis: &[GroupAffine<P>],
        a: &[Fr<P>],
        b: &[Fr<P>],
        u: GroupAffine<P>,
        blinding_basis: GroupAffine<P>,
        rng: &mut impl Rng,
        blind: Fr<P>,
    ) -> HidingRoundOutput<P> {
        assert_eq!(basis.len(), a.len());
        assert_eq!(basis.len(), b.len());
        assert!(basis.len() > 1);
        let (a_l, a_r) = split(a);
        let (b_l, b_r) = split(b);
        let (g_l, g_r) = split(basis);

        let [blind_r, blind_l] = [(); 2].map(|_| Fr::<P>::rand(rng));
        //.map(|field| blinding_basis.mul(field));
        let lj = inner_product(g_r, a_l)
            + blinding_basis.mul(blind_l)
            + u.mul(scalar_inner_product::<P>(a_l, b_r));
        let rj = inner_product(g_l, a_r)
            + blinding_basis.mul(blind_r)
            + u.mul(scalar_inner_product::<P>(a_r, b_l));

        let [lj, rj] = [lj, rj].map(|point| point.into_affine());

        let challenge = <ChallengeGenerator<P>>::round_challenge(&lj, &rj);
        let blind =
            challenge.square() * blind_l + blind + challenge.inverse().unwrap().square() * blind_r;
        let a = compress::<P>(a_r, a_l, challenge);
        let b = compress::<P>(b_l, b_r, challenge);
        let basis = compress_basis(g_l, g_r, challenge);
        HidingRoundOutput {
            lj,
            rj,
            a,
            b,
            basis,
            blind,
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
    ) -> (
        [GroupAffine<P>; 2],
        [Vec<Fr<P>>; 2],
        Vec<GroupAffine<P>>,
        Option<Fr<P>>,
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
        let a = compress::<P>(a_r, a_l, challenge);
        let b = compress::<P>(b_l, b_r, challenge);
        let basis = compress_basis(g_l, g_r, challenge);
        ([lj, rj], [a, b], basis, blind)
    }

    fn open_recursive_hiding(
        prev: HidingRoundOutput<P>,
        mut rounds: Vec<(GroupAffine<P>, GroupAffine<P>)>,
        point: Fr<P>,
        eval: Fr<P>,
        u: GroupAffine<P>,
        blinding_basis: GroupAffine<P>,
        rng: &mut impl Rng,
    ) -> HidingOpening<P>
    where
        Assert<HIDING>: IsTrue,
    {
        let HidingRoundOutput {
            a, b, basis, blind, ..
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
            let prev = Self::hiding_round(&*basis, &*a, &*b, u, blinding_basis, rng, blind);
            rounds.push((prev.lj, prev.rj));
            Self::open_recursive_hiding(prev, rounds, point, eval, u, blinding_basis, rng)
        }
    }

    pub fn open_hiding(
        &self,
        commitment: UnsafeHidingCommitment<P>,
        a: &[Fr<P>],
        point: Fr<P>,
        eval: Fr<P>,
        rng: &mut impl Rng,
    ) -> HidingOpening<P>
    where
        Assert<HIDING>: IsTrue,
    {
        let UnsafeHidingCommitment(commitment, blinding) = commitment;
        let commitment = Commitment::<_, true>(commitment);
        let u = ChallengeGenerator::inner_product_basis(&commitment, &point);
        let basis = &*self.basis;
        let blinding_basis = self.blinding_basis;
        let b = self.b(point);
        //let mut rng = &self.rng;
        let first = Self::hiding_round(basis, a, &*b, u, blinding_basis, rng, blinding);
        let rounds = vec![(first.lj, first.rj)];
        let opening =
            Self::open_recursive_hiding(first, rounds, point, eval, u, blinding_basis, rng);
        opening
    }

    fn hiding_round(
        basis: &[GroupAffine<P>],
        a: &[Fr<P>],
        b: &[Fr<P>],
        u: GroupAffine<P>,
        blinding_basis: GroupAffine<P>,
        rng: &mut impl Rng,
        blind: Fr<P>,
    ) -> HidingRoundOutput<P>
    where
        Assert<HIDING>: IsTrue,
    {
        assert_eq!(basis.len(), a.len());
        assert_eq!(basis.len(), b.len());
        assert!(basis.len() > 1);

        let blinding_factors = [(); 2].map(|_| Fr::<P>::rand(rng));
        let ([lj, rj], [a, b], basis, blind) = Self::general_round(
            basis,
            a,
            b,
            u,
            Some(blinding_basis),
            Some(blinding_factors),
            Some(blind),
        );
        let blind = blind.unwrap();
        HidingRoundOutput {
            lj,
            rj,
            a,
            b,
            basis,
            blind,
        }
    }
}
