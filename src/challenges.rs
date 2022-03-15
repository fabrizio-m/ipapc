use crate::{prove::Commitment, Fr};
use ark_ec::{
    short_weierstrass_jacobian::GroupAffine, AffineCurve, ProjectiveCurve, SWModelParameters,
};
use ark_ff::UniformRand;
use ark_serialize::CanonicalSerialize;
use blake2::{Blake2s256, Digest};
use rand::{prelude::StdRng, SeedableRng};
use std::marker::PhantomData;

#[derive(Clone)]
pub struct ChallengeGenerator<P: SWModelParameters> {
    data: Vec<u8>,
    _model: PhantomData<P>,
}

impl<P: SWModelParameters> ChallengeGenerator<P> {
    pub fn new() -> Self {
        Self {
            data: vec![],
            _model: PhantomData::default(),
        }
    }
    fn digest(&mut self, element: &GroupAffine<P>) {
        element.serialize_unchecked(&mut self.data).unwrap()
    }
    fn digest_scalar(&mut self, element: &Fr<P>) {
        element.serialize_unchecked(&mut self.data).unwrap()
    }
    ///digests commitments to generate the elements for amortization
    pub fn digest_for_amortization(&mut self, commitment: GroupAffine<P>) {
        commitment.serialize_unchecked(&mut self.data).unwrap()
    }
    ///generates the element for the lineal combination and the evaluation point
    pub fn amortization_elements(self) -> (Fr<P>, Fr<P>) {
        let mut rng = self.generate_rng();
        (<Fr<P>>::rand(&mut rng), <Fr<P>>::rand(&mut rng))
    }

    pub fn inner_product_basis<const HIDING: bool>(
        commitment: &Commitment<P, HIDING>,
        point: &Fr<P>,
    ) -> GroupAffine<P> {
        let mut challenge_generator = Self::new();
        challenge_generator.digest(&commitment.0);
        challenge_generator.digest_scalar(point);
        let mut rng = challenge_generator.generate_rng();

        <GroupAffine<P>>::prime_subgroup_generator()
            .mul(<Fr<P>>::rand(&mut rng))
            .into_affine()
    }
    pub fn round_challenge(lj: &GroupAffine<P>, rj: &GroupAffine<P>) -> Fr<P> {
        let mut challenge_generator = Self::new();
        challenge_generator.digest(lj);
        challenge_generator.digest(rj);
        let mut rng = challenge_generator.generate_rng();

        <Fr<P>>::rand(&mut rng)
    }
    fn generate_rng(self) -> StdRng {
        let Self { data, .. } = self;
        let seed: [u8; 32] = Blake2s256::digest(data).try_into().unwrap();
        StdRng::from_seed(seed)
    }
}
