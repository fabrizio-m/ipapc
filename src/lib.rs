use ark_ec::{
    short_weierstrass_jacobian::{GroupAffine, GroupProjective},
    AffineCurve, ModelParameters, ProjectiveCurve, SWModelParameters,
};
use ark_ff::{Field, One, PrimeField, UniformRand};
pub use prove::{Commitment, HidingOpening, Opening, UnsafeHidingCommitment};
use rand::{prelude::StdRng, Rng, SeedableRng};
use std::{
    convert::identity,
    fmt::Debug,
    iter::{repeat, successors},
};

pub mod amortization;
mod challenges;
pub mod fft;
mod homomorphism;
pub mod prove;
#[cfg(test)]
mod tests;
mod utils;
pub mod verify;

//type Poly<Fr> = DensePolynomial<Fr>;
type Fr<P> = <GroupAffine<P> as AffineCurve>::ScalarField;
pub struct IpaScheme<P>
where
    P: ModelParameters + SWModelParameters,
{
    basis: Vec<GroupAffine<P>>,
    blinding_basis: GroupAffine<P>,
    max_degree: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Blinding<P: SWModelParameters>(Fr<P>);

#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Init<T: SWModelParameters> {
    Seed(u64),
    Elements(Vec<GroupAffine<T>>, GroupAffine<T>),
}

impl<P> IpaScheme<P>
where
    P: ModelParameters + SWModelParameters,
    Fr<P>: One,
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
    fn commit_simple<'a>(&self, coeffs: Vec<Fr<P>>) -> GroupProjective<P> {
        debug_assert_eq!(coeffs.len(), self.max_degree);
        let bases = &*self.basis;
        let coeffs = coeffs
            .into_iter()
            .map(|e| e.into_repr())
            .collect::<Vec<_>>();
        let scalars = &*coeffs;
        let result = ark_ec::msm::VariableBaseMSM::multi_scalar_mul(bases, scalars);
        result
    }
    pub fn commit_hiding<'a>(
        &self,
        coeffs: Vec<Fr<P>>,
        rng: &mut impl Rng,
    ) -> UnsafeHidingCommitment<P> {
        assert_eq!(coeffs.len(), self.max_degree);
        let blinding_factor = Fr::<P>::rand(rng);
        let commitment = self.commit_simple(coeffs);
        let commitment = commitment + self.blinding_basis.mul(blinding_factor);

        UnsafeHidingCommitment(commitment.into_affine(), blinding_factor)
    }
    pub fn commit<'a>(&self, coeffs: Vec<Fr<P>>) -> Commitment<P, false> {
        assert_eq!(coeffs.len(), self.max_degree);
        let commitment = self.commit_simple(coeffs);
        Commitment(commitment.into_affine())
    }
    fn b(&self, z: Fr<P>) -> Vec<Fr<P>> {
        successors(Some(<Fr<P>>::one()), |previous| Some(*previous * z))
            .take(self.basis.len())
            .collect()
    }
    ///returns the reference string as commitments
    pub fn string(&self) -> Vec<Commitment<P, false>> {
        self.basis
            .iter()
            .cloned()
            .map(|elem| Commitment(elem))
            .collect()
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
