use crate::Commitment;
use crate::Fr;
use ark_ec::{AffineCurve, ProjectiveCurve, SWModelParameters};
use std::ops::{Add, Mul, Neg, Sub};

impl<P: SWModelParameters> Add for Commitment<P> {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        Self(self.0 + rhs.0)
    }
}
impl<P: SWModelParameters> Mul<Fr<P>> for Commitment<P> {
    type Output = Self;

    fn mul(self, rhs: Fr<P>) -> Self::Output {
        Self(self.0.mul(rhs).into_affine())
    }
}
impl<P: SWModelParameters> Neg for Commitment<P> {
    type Output = Self;

    fn neg(self) -> Self::Output {
        Self(-self.0)
    }
}
impl<P: SWModelParameters> Sub for Commitment<P> {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self::Output {
        self + (-rhs)
    }
}

#[test]
fn homomorphisms() {
    use crate::Init;
    use crate::IpaScheme;
    type P = ark_pallas::PallasParameters;

    let scheme = IpaScheme::init(Init::<P>::Seed(1), 2);
    let p1 = [0, 1, 2, 3].map(Fr::<P>::from).to_vec();
    let p2 = [1, 9, 2, 3].map(Fr::<P>::from).to_vec();
    let p3: Vec<_> = p1.iter().zip(p2.iter()).map(|(a, b)| *a + *b).collect();
    let scalar = Fr::<P>::from(9);
    let p4 = p1.iter().map(|e| *e * scalar).collect();
    println!("{:#?}", p3.len());
    let c1 = scheme.commit(p1);
    let c2 = scheme.commit(p2);
    let c3 = scheme.commit(p3);
    let c4 = scheme.commit(p4);

    assert!(c1 + c2 == c3);
    assert!(c4 == c1 * scalar);
}
