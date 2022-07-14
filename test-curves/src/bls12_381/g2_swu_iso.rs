use crate::bls12_381::*;
use ark_ec::models::{
    short_weierstrass::{Affine, SWCurveConfig},
    CurveConfig,
};
use ark_ff::MontFp;

use ark_ec::hashing::curve_maps::swu::SWUParams;

type G2Affine = Affine<SwuIsoParameters>;

#[derive(Clone, Default, PartialEq, Eq)]
pub struct SwuIsoParameters;

impl CurveConfig for SwuIsoParameters {
    type BaseField = Fq2;
    type ScalarField = Fr;

    //sage: g2_iso.codomain().order() == g2_iso.domain().order()
    //True
    /// COFACTOR = (x^8 - 4 x^7 + 5 x^6) - (4 x^4 + 6 x^3 - 4 x^2 - 4 x + 13) //
    /// 9
    /// = 305502333931268344200999753193121504214466019254188142667664032982267604182971884026507427359259977847832272839041616661285803823378372096355777062779109
    #[rustfmt::skip]
    const COFACTOR: &'static [u64] = &[
        0xbc69f08f2ee75b35,
        0x84c6a0ea91b35288,
        0x8e2a8e9145ad7689,
        0x986ff031508ffe13,
        0x29c2f178731db956,
        0xd82bf015d1212b02,
        0xec0ec69d7477c1ae,
        0x954cbc06689f6a35,
        0x9894c0adebbf6b4e,
        0x8020005aaa95551
    ];
    // const COFACTOR: &'static [u64] = &[
    //     0xcf1c38e31c7238e5,
    //     0x1616ec6e786f0c70,
    //     0x21537e293a6691ae,
    //     0xa628f1cb4d9e82ef,
    //     0xa68a205b2e5a7ddf,
    //     0xcd91de4547085aba,
    //     0x91d50792876a202,
    //     0x5d543a95414e7f1,
    // ];

    /// COFACTOR_INV = COFACTOR^{-1} mod r
    /// 26652489039290660355457965112010883481355318854675681319708643586776743290055
    #[rustfmt::skip]
    const COFACTOR_INV: Fr = MontFp!("26652489039290660355457965112010883481355318854675681319708643586776743290055");
}

// https://datatracker.ietf.org/doc/draft-irtf-cfrg-hash-to-curve/
// Hashing to Elliptic Curves
// 8.8.2.  BLS12-381 G2
//   *  E': y'^2 = x'^3 + A' * x' + B', where
//
//      -  A' = 240 * I
//
//      -  B' = 1012 * (1 + I)
//
//   * Z: -(2 + I)

impl SWCurveConfig for SwuIsoParameters {
    /// COEFF_A = [240, 0]
    const COEFF_A: Fq2 = Fq2::new(MontFp!("0"), MontFp!("240"));

    /// COEFF_B = [1012, 1012]
    const COEFF_B: Fq2 = Fq2::new(MontFp!("1012"), MontFp!("1012"));

    /// AFFINE_GENERATOR_COEFFS = (G2_GENERATOR_X, G2_GENERATOR_Y)
    const GENERATOR: G2Affine = G2Affine::new_unchecked(G2_GENERATOR_X, G2_GENERATOR_Y);
}

pub const G2_GENERATOR_X: Fq2 = Fq2::new(G2_GENERATOR_X_C0, G2_GENERATOR_X_C1);
pub const G2_GENERATOR_Y: Fq2 = Fq2::new(G2_GENERATOR_Y_C0, G2_GENERATOR_Y_C1);

// sage: gen_p = E2p.random_point()
// sage: gen_p *= 305502333931268344200999753193121504214466019254188142667664032982267604182971884026507427359259977847832272839041616661285803823378372096355777062779109
// sage: gen_p
// (3742867872338536021971219193407057626975566811775929228638115610845414823755734136602223916420662890679205212867815*X + 2008970337784971958599114573863246208442825790877271736434783886837659225024554344101831904846457859657000371883048 : 347020557560111175564150279686487416538085532798697633906066128669670325197422119639641125873165591545442034801173*X + 3002750085289950562156391971969303286031832683890522668239999505933394032768942051901168474820339778663917384775384 : 1)

#[rustfmt::skip]
pub const G2_GENERATOR_X_C0: Fq = MontFp!("2008970337784971958599114573863246208442825790877271736434783886837659225024554344101831904846457859657000371883048");

/// G2_GENERATOR_X_C1 =
#[rustfmt::skip]
pub const G2_GENERATOR_X_C1: Fq = MontFp!("3742867872338536021971219193407057626975566811775929228638115610845414823755734136602223916420662890679205212867815");

/// G2_GENERATOR_Y_C0 =
#[rustfmt::skip]
pub const G2_GENERATOR_Y_C0: Fq = MontFp!("3002750085289950562156391971969303286031832683890522668239999505933394032768942051901168474820339778663917384775384");

/// G2_GENERATOR_Y_C1 =
#[rustfmt::skip]
pub const G2_GENERATOR_Y_C1: Fq = MontFp!("347020557560111175564150279686487416538085532798697633906066128669670325197422119639641125873165591545442034801173");

impl SWUParams for SwuIsoParameters {
    // ZETA = -(2 + u) as per IETF draft.
    const ZETA: Fq2 = Fq2::new(MontFp!("-2"), MontFp!("-1"));
}
