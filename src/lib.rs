use crate::challenges::ChallengeGenerator;
use ark_ec::{
    short_weierstrass_jacobian::GroupAffine, AffineCurve, ModelParameters, ProjectiveCurve,
    SWModelParameters,
};
use ark_ff::{Field, One};
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
pub struct IpaScheme<P>
where
    P: ModelParameters + SWModelParameters,
{
    basis: Vec<GroupAffine<P>>,
    max_degree: usize,
}
#[derive(Debug)]
pub struct Opening<P: SWModelParameters> {
    point: Fr<P>,
    eval: Fr<P>,
    rounds: Vec<(GroupAffine<P>, GroupAffine<P>)>,
    a: Fr<P>,
}
struct RoundOutput<P: SWModelParameters> {
    lj: GroupAffine<P>,
    rj: GroupAffine<P>,
    a: Vec<Fr<P>>,
    b: Vec<Fr<P>>,
    basis: Vec<GroupAffine<P>>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Commitment<T: SWModelParameters>(GroupAffine<T>)
where
    GroupAffine<T>: Debug;

pub enum Init<T: SWModelParameters> {
    Seed(u64),
    Elements(Vec<GroupAffine<T>>),
}

impl<P> IpaScheme<P>
where
    P: ModelParameters + SWModelParameters,
    Fr<P>: One,
    Commitment<P>: Clone + Copy,
{
    pub fn init(init: Init<P>, max_size: u8) -> Self {
        let max_degree = 2_usize.pow(max_size as u32);
        let basis = init.to_elements(max_degree);
        Self { basis, max_degree }
    }
    pub fn commit<'a>(&self, coeffs: Vec<Fr<P>>) -> Commitment<P> {
        assert_eq!(coeffs.len(), self.max_degree);
        let point = self
            .basis
            .iter()
            .zip(coeffs.into_iter())
            .map(|(g, f)| g.mul(f))
            .reduce(Add::add)
            .unwrap();
        Commitment(point.into_affine())
    }
    fn b(&self, z: Fr<P>) -> Vec<Fr<P>> {
        successors(Some(<Fr<P>>::one()), |previous| Some(*previous * z))
            .take(self.basis.len())
            .collect()
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
        let (a_l, a_r) = split(a);
        let (b_l, b_r) = split(b);
        let (g_l, g_r) = split(basis);

        let lj = inner_product(g_r, a_l) + u.mul(scalar_inner_product::<P>(a_l, b_r)).into_affine();
        let rj = inner_product(g_l, a_r) + u.mul(scalar_inner_product::<P>(a_r, b_l)).into_affine();

        let challenge = <ChallengeGenerator<P>>::round_challenge(&lj, &rj);
        let a = compress::<P>(a_r, a_l, challenge);
        let b = compress::<P>(b_l, b_r, challenge);
        let basis = compress_basis(g_l, g_r, challenge);
        RoundOutput {
            lj,
            rj,
            a,
            b,
            basis,
        }
    }
    fn open_recursive(
        prev: RoundOutput<P>,
        mut rounds: Vec<(GroupAffine<P>, GroupAffine<P>)>,
        point: Fr<P>,
        eval: Fr<P>,
        u: GroupAffine<P>,
    ) -> Opening<P> {
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
    pub fn open(
        &self,
        commitment: Commitment<P>,
        a: &[Fr<P>],
        point: Fr<P>,
        eval: Fr<P>,
    ) -> Opening<P> {
        let u = ChallengeGenerator::inner_product_basis(&commitment, &point);
        let basis = &*self.basis;
        let b = self.b(point);
        let first = Self::round(basis, a, &*b, u);
        let rounds = vec![(first.lj, first.rj)];
        let opening = Self::open_recursive(first, rounds, point, eval, u);
        opening
    }
    pub fn verify(&self, commitment: Commitment<P>, open: Opening<P>) -> Option<Fr<P>> {
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
    fn to_elements(self, size: usize) -> Vec<GroupAffine<T>> {
        match self {
            Init::Seed(seed) => {
                let mut rng = StdRng::seed_from_u64(seed);
                let elems = repeat(())
                    .map(|_| {
                        let bytes: [u8; 32] = rng.gen();
                        let x = <T::BaseField as Field>::from_random_bytes(&bytes)?;
                        GroupAffine::<T>::get_point_from_x(x, false)
                    })
                    .filter_map(identity);
                elems.take(size).collect()
            }
            Init::Elements(mut elems) => {
                assert!(elems.len() >= size);
                elems.truncate(size);
                assert_eq!(elems.len(), size);
                elems
            }
        }
    }
}
