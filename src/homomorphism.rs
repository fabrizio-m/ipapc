use crate::{Commitment, Fr};
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
    use ark_poly::UVPolynomial;
    use ark_poly::{univariate::DensePolynomial, Polynomial};
    use rand::prelude::StdRng;
    use rand::SeedableRng;
    type P = ark_pallas::PallasParameters;

    let make_scheme = || {
        let rng = StdRng::seed_from_u64(1);
        IpaScheme::init(Init::<P>::Seed(1), 2, rng)
    };

    let p1 = [0, 1, 2, 3].map(Fr::<P>::from).to_vec();
    let p2 = [1, 9, 2, 3].map(Fr::<P>::from).to_vec();
    let p3: Vec<_> = p1.iter().zip(p2.iter()).map(|(a, b)| *a + *b).collect();
    let scalar = Fr::<P>::from(9);
    let p4 = p1.iter().map(|e| *e * scalar).collect::<Vec<_>>();

    let check = |poly: Vec<_>| {
        let mut scheme = make_scheme();
        let (c, b) = scheme.commit(poly.clone());
        let point = Fr::<P>::from(43);
        let eval = {
            let poly = DensePolynomial::<Fr<P>>::from_coefficients_slice(&*poly);
            poly.evaluate(&point)
        };
        let open = scheme.open(c, b, &*poly, point, eval);
        scheme.verify(c, open).unwrap()
    };
    let c1 = check(p1);
    let c2 = check(p2);
    let c3 = check(p3);
    let c4 = check(p4);

    assert!(c1 + c2 == c3);
    assert!(c4 == c1 * scalar);
}
