use crate::{
    hashing::{
        curve_maps::{
            swu::{SWUMap, SWUParams},
            wb::WBParams,
        },
        map_to_curve_hasher::MapToCurve,
    },
    models::SWModelParameters,
    short_weierstrass_jacobian::GroupAffine,
    ModelParameters,
};
//use ark_ff::{
//    biginteger::BigInteger64,
//    field_new,
//    fields::{FftParameters, Fp64, Fp64Parameters, FpParameters},
//};

use ark_ff::SquareRootField;
use ark_std::vec::Vec;
use hashbrown::HashMap;

use ark_ff::{
    biginteger::{BigInteger256, BigInteger64},
    field_new,
    fields::{FftParameters, Fp384, Fp384Parameters, FpParameters, Fp256, Fp256Parameters, Fp64, Fp64Parameters},
};
use crate::{AffineCurve};

use ark_ff::{Zero, One, Field, PrimeField};



use ark_test_curves::bls12_381::fq::{Fq as Fq381};
use ark_test_curves::bls12_381::fr::{Fr as Fr381};
use ark_test_curves::bls12_381::g1::{Parameters as BLS381_Parameters};

#[cfg(all(feature = "default", feature = "std"))]
use crate::hashing::{
    curve_maps::wb::WBMap, field_hashers::DefaultFieldHasher,
    map_to_curve_hasher::MapToCurveBasedHasher, HashToCurve,
};

pub type F127 = Fp64<F127Parameters>;

pub struct F127Parameters;

impl Fp64Parameters for F127Parameters {}
impl FftParameters for F127Parameters {
    type BigInt = BigInteger64;

    // N = 126 => s = 1
    const TWO_ADICITY: u32 = 1;

    // sage: FF(3)^63
    // 126
    #[rustfmt::skip]
    const TWO_ADIC_ROOT_OF_UNITY: BigInteger64 = BigInteger64([126]);
}

impl FpParameters for F127Parameters {
    /// MODULUS = 127
    #[rustfmt::skip]
    const MODULUS: BigInteger64 = BigInteger64([127]);

    const MODULUS_BITS: u32 = 7;

    const CAPACITY: u32 = Self::MODULUS_BITS - 1;

    const REPR_SHAVE_BITS: u32 = 64 - Self::MODULUS_BITS;

    // Nearst greater power of 2^64 (2^64R) to 127 is 0 so R = 1
    // sage: FF(2^64)
    // 2
    #[rustfmt::skip]
    const R: BigInteger64 = BigInteger64([2]);

    /// R2 = R^2 % Self::MODULUS
    #[rustfmt::skip]
    const R2: BigInteger64 = BigInteger64([4]);

    /// INV = -MODULUS^{-1} mod 2^64
    // sage: R = Integers(2^64)
    // sage: R
    // Ring of integers modulo 18446744073709551616
    // sage: m = R(127)
    // sage: m^(-1)
    // 9150747060186627967
    // sage: -m^(-1)
    // 9295997013522923649
    const INV: u64 = 9295997013522923649;

    // sage: FF(3).multiplicative_order()
    // 126
    // Montgomery conversion 3 * 2 = 6 % 127
    /// GENERATOR = 3
    #[rustfmt::skip]
    const GENERATOR: BigInteger64 = BigInteger64([6]);

    /// (Self::MODULUS - 1) / 2
    #[rustfmt::skip]
    const MODULUS_MINUS_ONE_DIV_TWO: BigInteger64 = BigInteger64([63]);

    // T and T_MINUS_ONE_DIV_TWO, where MODULUS - 1 = 2^S * T
    // For T coprime to 2

    /// t for 2^s * t = MODULUS - 1
    // T = (MODULUS - 1) / 2^S =
    // 12208678567578594777604504606729831043093128246378069236549469339647
    // sage: factor(127-1)
    // 2 * 3^2 * 7
    // sage: (127-1)/2
    // 63
    #[rustfmt::skip]
    const T: BigInteger64 = BigInteger64([63]);

    // (T - 1) / 2 = (63 - 1)/2
    #[rustfmt::skip]
    const T_MINUS_ONE_DIV_TWO: BigInteger64 = BigInteger64([31]);
}

const F127_ZERO: F127 = field_new!(F127, "0");
const F127_ONE: F127 = field_new!(F127, "1");

struct TestSWUMapToCurveParams;

impl ModelParameters for TestSWUMapToCurveParams {
    type BaseField = F127;
    type ScalarField = F127;
}
/// just because not defining another field
///
/// from itertools import product
/// p = 127
/// FF = GF(p)
/// for a,b in product(range(0,p), range(0,p)):
///     try:
///         E = EllipticCurve([FF(a),FF(b)])
///         if E.order() == p:
///             print(E)
///     except:
///         pass
///
/// y^2 = x^3 + x + 63
impl SWModelParameters for TestSWUMapToCurveParams {
    /// COEFF_A = 1
    const COEFF_A: F127 = F127_ONE;

    /// COEFF_B = 1
    #[rustfmt::skip]
    const COEFF_B: F127 = field_new!(F127, "63");

    const COFACTOR: &'static [u64] = &[1];

    #[rustfmt::skip]
    const COFACTOR_INV: F127 = F127_ONE;

    /// AFFINE_GENERATOR_COEFFS = (G1_GENERATOR_X, G1_GENERATOR_Y)
    const AFFINE_GENERATOR_COEFFS: (Self::BaseField, Self::BaseField) =
        (field_new!(F127, "62"), field_new!(F127, "70"));
}

impl SWUParams for TestSWUMapToCurveParams {
    const XI: F127 = field_new!(F127, "-1");
    const ZETA: F127 = field_new!(F127, "3");
    const XI_ON_ZETA_SQRT: F127 = field_new!(F127, "13");
}

/// test that field_new make a none zero element out of 1
#[test]
fn test_field_element_construction() {
    let a1 = F127::from(1);
    let a2 = F127::from(2);
    let a3 = F127::from(125);

    assert!(F127::from(0) == a2 + a3);
    assert!(F127::from(0) == a2 * a1 + a3);
}

#[test]
fn test_field_division() {
    let num = F127::from(0x3d);
    let den = F127::from(0x7b);
    let num_on_den = F127::from(0x50);

    assert!(num / den == num_on_den);
}

/// Testing checking the hashing parameters are sane
/// check zeta is a non-square
#[test]
fn checking_the_hashing_parameters() {
    assert!(SquareRootField::legendre(&TestSWUMapToCurveParams::ZETA).is_qr() == false);
}

/// The point of the test is to get a  simpl SWU compatible curve
/// and make simple hash
#[cfg(all(feature = "default", feature = "std"))]
#[test]
fn hash_arbitary_string_to_curve_swu() {
    use blake2::VarBlake2b;

    let test_swu_to_curve_hasher = MapToCurveBasedHasher::<
        GroupAffine<TestSWUMapToCurveParams>,
        DefaultFieldHasher<VarBlake2b>,
        SWUMap<TestSWUMapToCurveParams>,
    >::new(&[1])
    .unwrap();

    let hash_result = test_swu_to_curve_hasher.hash(b"if you stick a Babel fish in your ear you can instantly understand anything said to you in any form of language.").expect("fail to hash the string to curve");

    #[cfg(feature = "std")]
    println!("{:?}, {:?}", hash_result, hash_result.x,);

    assert!(hash_result.x != F127_ZERO);
}

/// the test use a simple SWU compatible curve
/// and map the whole field to it. We observe the map behaviour. Specifically
/// The map is not constant and that everything can be mapped and nobody panics
#[test]
fn map_field_to_curve_swu() {
    let test_map_to_curve = SWUMap::<TestSWUMapToCurveParams>::new_map_to_curve(&[0]).unwrap();

    let mut map_range: Vec<GroupAffine<TestSWUMapToCurveParams>> = vec![];
    for current_field_element in 0..127 {
        map_range.push(
            test_map_to_curve
                .map_to_curve(F127::from(current_field_element as u64))
                .unwrap(),
        );
    }

    let mut counts = HashMap::new();

    let mode = map_range
        .iter()
        .copied()
        .max_by_key(|&n| {
            let count = counts.entry(n).or_insert(0);
            *count += 1;
            *count
        })
        .unwrap();

    #[cfg(feature = "std")]
    println!(
        "mode {} repeated {} times",
        mode,
        counts.get(&mode).unwrap()
    );

    assert!(*counts.get(&mode).unwrap() != 127);
}

// Testing wb19 on  small curvse
// E_isogenous : Elliptic Curve defined by y^2 = x^3 + 109*x + 124 over Finite
// Field of size 127 E : y^2 = x^3 + 3
// Isogeny map
// (10*x^18*y + 59*x^17*y + 41*x^16*y + 48*x^15*y - 7*x^14*y + 6*x^13*y +
// 5*x^12*y + 62*x^11*y + 12*x^10*y + 36*x^9*y - 49*x^8*y - 18*x^7*y - 63*x^6*y
// - 43*x^5*y - 60*x^4*y - 18*x^3*y + 30*x^2*y - 57*x*y - 34*y)/(x^18 + 44*x^17
// - 63*x^16 + 52*x^15 + 3*x^14 + 38*x^13 - 30*x^12 + 11*x^11 - 42*x^10 - 13*x^9
// - 46*x^8 - 61*x^7 - 16*x^6 - 55*x^5 + 18*x^4 + 23*x^3 - 24*x^2 - 18*x + 32)
struct TestSWU127MapToIsogenousCurveParams;

// first we define the isogenous curve
// sage: E_isogenous.order()
// 127
impl ModelParameters for TestSWU127MapToIsogenousCurveParams {
    type BaseField = F127;
    type ScalarField = F127;
}

// E_isogenous : Elliptic Curve defined by y^2 = x^3 + 109*x + 124 over Finite
// Field of size 127
impl SWModelParameters for TestSWU127MapToIsogenousCurveParams {
    /// COEFF_A = 1
    const COEFF_A: F127 = field_new!(F127, "109");

    /// COEFF_B = 1
    #[rustfmt::skip]
    const COEFF_B: F127 = field_new!(F127, "124");

    const COFACTOR: &'static [u64] = &[1];

    #[rustfmt::skip]
    const COFACTOR_INV: F127 = F127_ONE;

    /// AFFINE_GENERATOR_COEFFS = (G1_GENERATOR_X, G1_GENERATOR_Y)
    const AFFINE_GENERATOR_COEFFS: (Self::BaseField, Self::BaseField) =
        (field_new!(F127, "84"), field_new!(F127, "2"));
}

impl SWUParams for TestSWU127MapToIsogenousCurveParams {
    const XI: F127 = field_new!(F127, "-1");
    const ZETA: F127 = field_new!(F127, "3");
    const XI_ON_ZETA_SQRT: F127 = field_new!(F127, "13");
}

// The struct defining our parameters for the target curve
struct TestWBF127MapToCurveParams;

impl ModelParameters for TestWBF127MapToCurveParams {
    type BaseField = F127;
    type ScalarField = F127;
}

impl SWModelParameters for TestWBF127MapToCurveParams {
    /// COEFF_A = 0
    const COEFF_A: F127 = F127_ZERO;

    /// COEFF_B = 3
    #[rustfmt::skip]
    const COEFF_B: F127 = field_new!(F127, "3");

    const COFACTOR: &'static [u64] = &[1];

    #[rustfmt::skip]
    const COFACTOR_INV: F127 = F127_ONE;

    /// AFFINE_GENERATOR_COEFFS = (G1_GENERATOR_X, G1_GENERATOR_Y)
    const AFFINE_GENERATOR_COEFFS: (Self::BaseField, Self::BaseField) =
        (field_new!(F127, "62"), field_new!(F127, "70"));
}

// E_isogenous : Elliptic Curve defined by y^2 = x^3 + 109*x + 124 over Finite
// Field of size 127 psi_x: (-57*x^13 - 21*x^12 + 10*x^11 + 34*x^10 + 40*x^9 -
// 13*x^8 + 32*x^7 - 32*x^6 + 23*x^5 - 14*x^4 + 39*x^3 + 23*x^2 + 63*x +
// 4)/(x^12 - 13*x^11 + 11*x^10 - 33*x^9 - 30*x^8 + 30*x^7 + 34*x^6 - 44*x^5 +
// 63*x^4 - 20*x^3 - 10*x^2 + 31*x + 2)

// psi_y: (10*x^18*y + 59*x^17*y + 41*x^16*y + 48*x^15*y - 7*x^14*y + 6*x^13*y +
// 5*x^12*y + 62*x^11*y + 12*x^10*y + 36*x^9*y - 49*x^8*y - 18*x^7*y - 63*x^6*y
// - 43*x^5*y - 60*x^4*y - 18*x^3*y + 30*x^2*y - 57*x*y - 34*y)/(x^18 + 44*x^17
// - 63*x^16 + 52*x^15 + 3*x^14 + 38*x^13 - 30*x^12 + 11*x^11 - 42*x^10 - 13*x^9
// - 46*x^8 - 61*x^7 - 16*x^6 - 55*x^5 + 18*x^4 + 23*x^3 - 24*x^2 - 18*x + 32)
impl WBParams for TestWBF127MapToCurveParams {
    type IsogenousCurve = TestSWU127MapToIsogenousCurveParams;

    const PHI_X_NOM: &'static [<Self::IsogenousCurve as ModelParameters>::BaseField] = &[
        field_new!(F127, "4"),
        field_new!(F127, "63"),
        field_new!(F127, "23"),
        field_new!(F127, "39"),
        field_new!(F127, "-14"),
        field_new!(F127, "23"),
        field_new!(F127, "-32"),
        field_new!(F127, "32"),
        field_new!(F127, "-13"),
        field_new!(F127, "40"),
        field_new!(F127, "34"),
        field_new!(F127, "10"),
        field_new!(F127, "-21"),
        field_new!(F127, "-57"),
    ];

    const PHI_X_DEN: &'static [<Self::IsogenousCurve as ModelParameters>::BaseField] = &[
        field_new!(F127, "2"),
        field_new!(F127, "31"),
        field_new!(F127, "-10"),
        field_new!(F127, "-20"),
        field_new!(F127, "63"),
        field_new!(F127, "-44"),
        field_new!(F127, "34"),
        field_new!(F127, "30"),
        field_new!(F127, "-30"),
        field_new!(F127, "-33"),
        field_new!(F127, "11"),
        field_new!(F127, "-13"),
        field_new!(F127, "1"),
    ];

    const PHI_Y_NOM: &'static [<Self::IsogenousCurve as ModelParameters>::BaseField] = &[
        field_new!(F127, "-34"),
        field_new!(F127, "-57"),
        field_new!(F127, "30"),
        field_new!(F127, "-18"),
        field_new!(F127, "-60"),
        field_new!(F127, "-43"),
        field_new!(F127, "-63"),
        field_new!(F127, "-18"),
        field_new!(F127, "-49"),
        field_new!(F127, "36"),
        field_new!(F127, "12"),
        field_new!(F127, "62"),
        field_new!(F127, "5"),
        field_new!(F127, "6"),
        field_new!(F127, "-7"),
        field_new!(F127, "48"),
        field_new!(F127, "41"),
        field_new!(F127, "59"),
        field_new!(F127, "10"),
    ];

    const PHI_Y_DEN: &'static [<Self::IsogenousCurve as ModelParameters>::BaseField] = &[
        field_new!(F127, "32"),
        field_new!(F127, "-18"),
        field_new!(F127, "-24"),
        field_new!(F127, "23"),
        field_new!(F127, "18"),
        field_new!(F127, "-55"),
        field_new!(F127, "-16"),
        field_new!(F127, "-61"),
        field_new!(F127, "-46"),
        field_new!(F127, "-13"),
        field_new!(F127, "-42"),
        field_new!(F127, "11"),
        field_new!(F127, "-30"),
        field_new!(F127, "38"),
        field_new!(F127, "3"),
        field_new!(F127, "52"),
        field_new!(F127, "-63"),
        field_new!(F127, "44"),
        field_new!(F127, "1"),
    ];
}

/// The point of the test is to get a  simpl SWU compatible curve
/// and make simple hash
#[cfg(all(feature = "default", feature = "std"))]
#[test]
fn hash_arbitary_string_to_curve_wb() {
    use blake2::VarBlake2b;

    let test_wb_to_curve_hasher = MapToCurveBasedHasher::<
        GroupAffine<TestWBF127MapToCurveParams>,
        DefaultFieldHasher<VarBlake2b>,
        WBMap<TestWBF127MapToCurveParams>,
    >::new(&[1])
    .unwrap();

    let hash_result = test_wb_to_curve_hasher.hash(b"if you stick a Babel fish in your ear you can instantly understand anything said to you in any form of language.").expect("fail to hash the string to curve");

    #[cfg(feature = "std")]
    println!("the wb hash is: {:?}", hash_result);

    assert!(hash_result.x != F127_ZERO);
}
//////BLS12-381 Tests
struct TestSWUMapToCurveBLS12_381Params;

impl ModelParameters for TestSWUMapToCurveBLS12_381Params {
    type BaseField = Fq381;
    type ScalarField = Fr381;
}

impl SWModelParameters for TestSWUMapToCurveBLS12_381Params {
    //const COEFF_A: Fq = field_new!(Fq, "012190336318893619529228877361869031420615612348429846051986726275283378313155663745811710833465465981901188123677"); //sage doesn't approve of this
     const COEFF_A: Fq381 = field_new!(Fq381, "2858613208430792460670318198342879349494999260436483523154854961351063857243634726019465176474256126859776719994977");

    #[rustfmt::skip]
    const COEFF_B: Fq381 = field_new!(Fq381, "2906670324641927570491258158026293881577086121416628140204402091718288198173574630967936031029026176254968826637280");

    //sage: g1_iso.domain().order()/52435875175126190479447740508185965837690552500527637822603658699938581184513
    //76329603384216526031706109802092473003
    /// COFACTOR = (x - 1)^2 / 3  = 76329603384216526031706109802092473003
    const COFACTOR: &'static [u64] = &[0x8c00aaab0000aaab, 0x396c8c005555e156];

    /// COFACTOR_INV = COFACTOR^{-1} mod r
    /// = 52435875175126190458656871551744051925719901746859129887267498875565241663483
    #[rustfmt::skip]
    const COFACTOR_INV: Fr381 = field_new!(Fr381, "52435875175126190458656871551744051925719901746859129887267498875565241663483");

    /// AFFINE_GENERATOR_COEFFS = (G1_GENERATOR_X, G1_GENERATOR_Y)
    const AFFINE_GENERATOR_COEFFS: (Self::BaseField, Self::BaseField) =
        (field_new!(Fq381, "628127623378585612843095022119507708025289394540669560027004601611569871267541856210323712812047239723504248810248"), field_new!(Fq381, "344075650239127142968089520704786925624745533124141202280424406752399324209523628375922007963596482424831722220273"));

}

// impl SWModelParameters for TestSWUMapToCurveBLS12_381Params {
//     /// COEFF_A = 0
//     const COEFF_A: Fq381 = field_new!(Fq381, "0");

//     /// COEFF_B = 4
//     #[rustfmt::skip]
//     const COEFF_B: Fq381 = field_new!(Fq381, "4");

//     /// COFACTOR = (x - 1)^2 / 3  = 76329603384216526031706109802092473003
//     const COFACTOR: &'static [u64] = &[0x8c00aaab0000aaab, 0x396c8c005555e156];

//     /// COFACTOR_INV = COFACTOR^{-1} mod r
//     /// = 52435875175126190458656871551744051925719901746859129887267498875565241663483
//     #[rustfmt::skip]
//     const COFACTOR_INV: Fr381 = field_new!(Fr381, "52435875175126190458656871551744051925719901746859129887267498875565241663483");

//     /// AFFINE_GENERATOR_COEFFS = (G1_GENERATOR_X, G1_GENERATOR_Y)
//     const AFFINE_GENERATOR_COEFFS: (Self::BaseField, Self::BaseField) =
//         (G1_GENERATOR_X, G1_GENERATOR_Y);

//     #[inline(always)]
//     fn mul_by_a(_: &Self::BaseField) -> Self::BaseField {
//         Self::BaseField::zero()
//     }
// }

// /// G1_GENERATOR_X =
// /// 3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507
// #[rustfmt::skip]
// pub const G1_GENERATOR_X: Fq381 = field_new!(Fq381, "3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507");

// /// G1_GENERATOR_Y =
// /// 1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724655353366534992391756441569
// #[rustfmt::skip]
// pub const G1_GENERATOR_Y: Fq381 = field_new!(Fq381, "1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724655353366534992391756441569");

// Rust does not let us to re-use the BLS12 381 definition
// impl SWModelParameters for TestSWUMapToCurveBLS12_381Params {
//     /// COEFF_A = 1
//     const COEFF_A: Fq381 = <BLS381_Parameters as SWModelParameters>::COEFF_A;

//     /// COEFF_B = 1
//     #[rustfmt::skip]
//     const COEFF_B: Fq381 = BLS381_Parameters::COEFF_B; 

//     const COFACTOR: &'static [u64] = &[1];

//     #[rustfmt::skip]
//     const COFACTOR_INV: Fr381 = BLS381_Parameters::COFACTOR_INV;

//     // /// AFFINE_GENERATOR_COEFFS = (G1_GENERATOR_X, G1_GENERATOR_Y)
//     const AFFINE_GENERATOR_COEFFS: (Self::BaseField, Self::BaseField) = BLS381_Parameters::AFFINE_GENERATOR_COEFFS;


// }

impl SWUParams for TestSWUMapToCurveBLS12_381Params {    
    const XI : Fq381 = field_new!(Fq381, "11"); //a nonsquare in Fq ietf standard
    const ZETA: Fq381 = field_new!(Fq381, "2"); //arbitatry primitive root of unity (element)CHECK THE STANDARD
    const XI_ON_ZETA_SQRT: Fq381 = field_new!(Fq381, "4002409555221667393417789825735904156556882819939007885332058136124031650490837864442687629129015664037894272559787"); ////square root of THETA=Xi/Zet

}

 /// The point of the test is to get a  simpl SWU compatible curve
/// and make simple hash
#[cfg(all(feature = "default", feature = "std"))]
#[test]
fn hash_arbitary_string_to_curve_swu_BLS_ISOGENY() {
    use blake2::{VarBlake2b};

    let test_swu_to_curve_hasher = MapToCurveBasedHasher::<GroupAffine<TestSWUMapToCurveBLS12_381Params>, DefaultFieldHasher<VarBlake2b>, SWUMap<TestSWUMapToCurveBLS12_381Params>>::new(&[1]).unwrap();
    
    let hash_result = test_swu_to_curve_hasher.hash(b"if you stick a Babel fish in your ear you can instantly understand anything said to you in any form of language.").expect("fail to hash the string t1o curve");

    
    println!("{:?}, {:?}", hash_result, hash_result.x, );

    assert!(hash_result.x != field_new!(Fq381, "0"));

}
