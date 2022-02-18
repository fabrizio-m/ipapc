use crate::challenges::ChallengeGenerator;
use ark_ec::{
    short_weierstrass_jacobian::{GroupAffine, GroupProjective},
    AffineCurve, ModelParameters, ProjectiveCurve, SWModelParameters,
};
use ark_ff::{Field, One, UniformRand};
use rand::{prelude::StdRng, Rng, SeedableRng};
use std::{
    convert::identity,
    fmt::Debug,
    iter::{repeat, successors},
    ops::Add,
};
use utils::{compress, compress_basis, inner_product, scalar_inner_product, split};

mod challenges;
mod homomorphism;
#[cfg(test)]
mod tests;
mod utils;

//type Poly<Fr> = DensePolynomial<Fr>;
type Fr<P> = <GroupAffine<P> as AffineCurve>::ScalarField;
pub struct IpaScheme<P, const HIDING: bool>
where
    P: ModelParameters + SWModelParameters,
{
    basis: Vec<GroupAffine<P>>,
    blinding_basis: GroupAffine<P>,
    //rng: R,
    max_degree: usize,
}
#[derive(Debug)]
pub struct HidingOpening<P: SWModelParameters> {
    point: Fr<P>,
    eval: Fr<P>,
    rounds: Vec<(GroupAffine<P>, GroupAffine<P>)>,
    a: Fr<P>,
    blinding_factor: Fr<P>,
}
#[derive(Debug)]
pub struct Opening<P: SWModelParameters> {
    point: Fr<P>,
    eval: Fr<P>,
    rounds: Vec<(GroupAffine<P>, GroupAffine<P>)>,
    a: Fr<P>,
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
pub struct Commitment<T: SWModelParameters, const HIDING: bool>(GroupAffine<T>)
where
    GroupAffine<T>: Debug;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct UnsafeHidingCommitment<T: SWModelParameters>(GroupAffine<T>, Fr<T>)
where
    GroupAffine<T>: Debug;

impl<P: SWModelParameters> From<UnsafeHidingCommitment<P>> for Commitment<P, true> {
    fn from(unsafe_commitment: UnsafeHidingCommitment<P>) -> Self {
        Commitment(unsafe_commitment.0)
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Blinding<P: SWModelParameters>(Fr<P>);

pub enum Init<T: SWModelParameters> {
    Seed(u64),
    Elements(Vec<GroupAffine<T>>, GroupAffine<T>),
}

impl<P, const HIDING: bool> IpaScheme<P, HIDING>
where
    P: ModelParameters + SWModelParameters,
    Fr<P>: One,
    Commitment<P, HIDING>: Clone + Copy,
{
    pub fn init(init: Init<P>, max_size: u8) -> Self {
        let max_degree = 2_usize.pow(max_size as u32);
        let (basis, blinding_basis) = init.to_elements(max_degree);
        Self {
            basis,
            max_degree,
            blinding_basis,
            //rng,
        }
    }
    pub fn commit_simple<'a>(&self, coeffs: Vec<Fr<P>>) -> GroupProjective<P> {
        debug_assert_eq!(coeffs.len(), self.max_degree);
        self.basis
            .iter()
            .zip(coeffs.into_iter())
            .map(|(g, f)| g.mul(f))
            .reduce(Add::add)
            .unwrap()
    }
    pub fn commit_hiding<'a>(
        &self,
        coeffs: Vec<Fr<P>>,
        rng: &mut impl Rng,
    ) -> UnsafeHidingCommitment<P>
    where
        Assert<HIDING>: IsTrue,
    {
        assert_eq!(coeffs.len(), self.max_degree);
        let blinding_factor = Fr::<P>::rand(rng);
        let commitment = self.commit_simple(coeffs);
        let commitment = commitment + self.blinding_basis.mul(blinding_factor);

        UnsafeHidingCommitment(commitment.into_affine(), blinding_factor)
    }
    pub fn commit<'a>(&self, coeffs: Vec<Fr<P>>) -> Commitment<P, HIDING>
    where
        Assert<HIDING>: IsFalse,
    {
        assert_eq!(coeffs.len(), self.max_degree);
        let commitment = self.commit_simple(coeffs);
        Commitment(commitment.into_affine())
    }
    fn b(&self, z: Fr<P>) -> Vec<Fr<P>> {
        successors(Some(<Fr<P>>::one()), |previous| Some(*previous * z))
            .take(self.basis.len())
            .collect()
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
        let final_check = final_commit == basis[0].mul(a) + u.mul(a * b);
        if final_check {
            Some(eval)
        } else {
            None
        }
    }
}

impl<T: SWModelParameters> Init<T> {
    fn to_elements(self, size: usize) -> (Vec<GroupAffine<T>>, GroupAffine<T>) {
        match self {
            Init::Seed(seed) => {
                let mut rng = StdRng::seed_from_u64(seed);
                let mut elems = repeat(())
                    .map(|_| {
                        let bytes: [u8; 32] = rng.gen();
                        let x = <T::BaseField as Field>::from_random_bytes(&bytes)?;
                        GroupAffine::<T>::get_point_from_x(x, false)
                    })
                    .filter_map(identity);
                let blind = elems.next().unwrap();
                (elems.take(size).collect(), blind)
            }
            Init::Elements(mut elems, blinding) => {
                assert!(elems.len() >= size);
                elems.truncate(size);
                assert_eq!(elems.len(), size);
                (elems, blinding)
            }
        }
    }
}

pub struct Assert<const COND: bool> {}

pub trait IsTrue {}
pub trait IsFalse {}

impl IsTrue for Assert<true> {}
impl IsFalse for Assert<false> {}
