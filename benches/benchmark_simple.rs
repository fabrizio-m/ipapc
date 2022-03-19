use ark_ec::{short_weierstrass_jacobian::GroupAffine, AffineCurve};
use ark_pallas::PallasParameters;
use ark_poly::{Polynomial, UVPolynomial};
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use ipapc::{Init, IpaScheme};
use rand::{prelude::ThreadRng, thread_rng, Rng};
use std::iter::repeat;

type Scheme = IpaScheme<PallasParameters>;
type Fr = <GroupAffine<PallasParameters> as AffineCurve>::ScalarField;

const SIZE: u8 = 10;

fn sample(size: u8) -> (Scheme, Vec<Fr>, ThreadRng) {
    let scheme = IpaScheme::<PallasParameters>::init(Init::Seed(1), size);
    let mut rng = thread_rng();
    let poly: Vec<Fr> = repeat(())
        .map(|_| rng.gen())
        .take(2_usize.pow(size as u32))
        .collect();
    (scheme, poly, rng)
}

pub fn commit(c: &mut Criterion) {
    c.bench_function("commit_simple", |b| {
        let (scheme, poly, _rng) = sample(SIZE);

        let coeffs = black_box(poly.to_vec());
        b.iter(|| {
            let _commit = scheme.commit(black_box(coeffs.clone()));
        })
    });
}
pub fn open(c: &mut Criterion) {
    c.bench_function("open_simple", |b| {
        let (scheme, poly, mut rng) = sample(SIZE);
        b.iter_batched(
            || {
                let coeffs = black_box(poly.to_vec());
                let commit = scheme.commit(coeffs.clone());
                let point: Fr = rng.gen();
                let (eval, poly) = {
                    let poly = ark_poly::univariate::DensePolynomial::<Fr>::from_coefficients_slice(
                        &*coeffs,
                    );
                    let eval = poly.evaluate(&point);
                    (eval, poly.coeffs)
                };
                (point, commit, poly, eval)
            },
            |(point, commitment, a, eval)| {
                let open = scheme.open(commitment, &*a, point, eval);
                (open, a)
            },
            criterion::BatchSize::SmallInput,
        );
    });
}
pub fn verify(c: &mut Criterion) {
    c.bench_function("verify_simple", |b| {
        let (scheme, poly, mut rng) = sample(SIZE);
        b.iter_batched(
            || {
                let coeffs = black_box(poly.to_vec());
                let commit = scheme.commit(coeffs.clone());
                let point: Fr = rng.gen();
                let (eval, poly) = {
                    let poly = ark_poly::univariate::DensePolynomial::<Fr>::from_coefficients_slice(
                        &*coeffs,
                    );
                    let eval = poly.evaluate(&point);
                    (eval, poly.coeffs)
                };
                //(point, commit, poly, eval)
                let open = scheme.open(commit, &*poly, point, eval);
                (commit, open)
            },
            |(commit, open)| {
                let open = scheme.verify(commit, open);
                open
            },
            criterion::BatchSize::SmallInput,
        );
    });
}
criterion_group!(benches, commit, open, verify);
criterion_main!(benches);
