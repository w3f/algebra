use ark_ff::{
    biginteger::{BigInteger256, BigInteger64},
    field_new,
    fields::{FftParameters, Fp384, Fp384Parameters, FpParameters, Fp256, Fp256Parameters, Fp64, Fp64Parameters},
};
use crate::{ModelParameters, models::SWModelParameters, AffineCurve};
use crate::short_weierstrass_jacobian::GroupAffine;
use crate::hashing::curve_maps::swu::{SWUParams, SWUMap};
use super::map_to_curve_hasher::{MapToCurveBasedHasher, MapToCurve};
use crate::hashing::{map_to_curve_hasher::HashToField, field_hashers::DefaultFieldHasher,HashToCurve};
use ark_ff::{Zero, One, Field, PrimeField, SquareRootField};
use ark_std::vec::Vec;
use std::collections::HashMap;

use ark_test_curves::bls12_381::fq::{Fq as Fq381};
use ark_test_curves::bls12_381::fr::{Fr as Fr381};
use ark_test_curves::bls12_381::g1::{Parameters as BLS381_Parameters};
    
pub type F127 = Fp64<F127Parameters>;

pub struct F127Parameters;

impl Fp64Parameters for F127Parameters {}
impl FftParameters for F127Parameters {
    type BigInt = BigInteger64;

    // N = 126 => s = 1
    const TWO_ADICITY: u32 = 1;

    //sage: FF(3)^63
    //126
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

    // Nearst power of 2^64 to 127 is 0 so R = 1 but maybe they mean larger
    // otherwise square root panics
    //sage: FF(2^64)
    //2
    #[rustfmt::skip]
    const R: BigInteger64 = BigInteger64([2]);

    #[rustfmt::skip]
    const R2: BigInteger64 = BigInteger64([4]);

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
    /// GENERATOR = 3
    #[rustfmt::skip]
    const GENERATOR: BigInteger64 = BigInteger64([3]);

    #[rustfmt::skip]
    const MODULUS_MINUS_ONE_DIV_TWO: BigInteger64 = BigInteger64([63]);

    // T and T_MINUS_ONE_DIV_TWO, where MODULUS - 1 = 2^S * T
    // For T coprime to 2

    // T = (MODULUS - 1) / 2^S =
    // 12208678567578594777604504606729831043093128246378069236549469339647
    //sage: factor(127-1)
    //2 * 3^2 * 7
    //sage: (127-1)/2
    //63
    #[rustfmt::skip]
    const T: BigInteger64 = BigInteger64([63]);

    // (T - 1) / 2 =
    #[rustfmt::skip]
    const T_MINUS_ONE_DIV_TWO: BigInteger64 = BigInteger64([31]);
}

const F127_ZERO: F127 = field_new!(F127, "0");
const F127_ONE: F127 = field_new!(F127, "1");

struct TestSWUMapToCurveF127Params;

impl ModelParameters for TestSWUMapToCurveF127Params {
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
///
impl SWModelParameters for TestSWUMapToCurveF127Params {
    /// COEFF_A = 1
    const COEFF_A: F127 = field_new!(F127, "1");

    /// COEFF_B = 1
    #[rustfmt::skip]
    const COEFF_B: F127 = field_new!(F127, "63");

    const COFACTOR: &'static [u64] = &[1];

    #[rustfmt::skip]
    const COFACTOR_INV: F127 = field_new!(F127, "1");

    /// AFFINE_GENERATOR_COEFFS = (G1_GENERATOR_X, G1_GENERATOR_Y)
    const AFFINE_GENERATOR_COEFFS: (Self::BaseField, Self::BaseField) =
        (field_new!(F127, "62"), field_new!(F127, "70"));


}
impl SWUParams for TestSWUMapToCurveF127Params {

    const XI : F127 = field_new!(F127, "-1");
    const ZETA: F127 = field_new!(F127, "3");
    const XI_ON_ZETA_SQRT: F127 = field_new!(F127, "14");

}

///test that field_new make a none zero element out of 1
#[test]
fn test_field_element_construction() {
    let a1 = F127::from(1);
    let a2 = F127::from(2);
    let a3 = F127::from(125);

    assert!(F127::from(0) == a2 + a3);
}
    

/// Testing checking the hashing parameters are sane
#[test]
fn chceking_the_hsahing_parameters() {
    ///check zeta is a non-square
    assert!(SquareRootField::legendre(&TestSWUMapToCurveF127Params::ZETA).is_qr() == false);
    
}

/// The point of the test is to get a  simpl SWU compatible curve
/// and make simple hash
#[test]
fn hash_arbitary_string_to_curve_swu() {
    use blake2::{VarBlake2b};

    let test_swu_to_curve_hasher = MapToCurveBasedHasher::<GroupAffine<TestSWUMapToCurveF127Params>, DefaultFieldHasher<VarBlake2b>, SWUMap<TestSWUMapToCurveF127Params>>::new(&[1]).unwrap();
    
    let hash_result = test_swu_to_curve_hasher.hash(b"if you stick a Babel fish in your ear you can instantly understand anything said to you in any form of language.").expect("fail to hash the string to curve");

    
    println!("{:?}, {:?}", hash_result, hash_result.x, );

    assert!(hash_result.x != field_new!(F127, "0"));

}

/// the test use a simple SWU compatible curve
/// and map the whole field to it. We observe the map behaviour. Specifically
/// The map is not constant and that everything can be mapped and nobody panics
#[test]
fn map_field_to_curve_swu() {
    let test_map_to_curve =  SWUMap::<TestSWUMapToCurveF127Params>::new_map_to_curve(&[0]).unwrap();

    let mut map_range : Vec<GroupAffine<TestSWUMapToCurveF127Params>> = vec![];
    for current_field_element in 0..127 {
        map_range.push(test_map_to_curve.map_to_curve(F127::from(current_field_element)).unwrap());
    }

    let mut counts = HashMap::new();

    let mode = map_range.iter().copied().max_by_key(|&n| {
        let count = counts.entry(n).or_insert(0);
        *count += 1;
        *count
    }).unwrap();

    println!("mode {} repeated {} times", mode, counts.get(&mode).unwrap());

    assert!(*counts.get(&mode).unwrap() != 127);
}

//////BLS12-381 Tests
struct TestSWUMapToCurveBLS12_381Params;

impl ModelParameters for TestSWUMapToCurveBLS12_381Params {
    type BaseField = Fq381;
    type ScalarField = Fr381;
}

impl SWModelParameters for TestSWUMapToCurveBLS12_381Params {
    /// COEFF_A = 0
    const COEFF_A: Fq381 = field_new!(Fq381, "0");

    /// COEFF_B = 4
    #[rustfmt::skip]
    const COEFF_B: Fq381 = field_new!(Fq381, "4");

    /// COFACTOR = (x - 1)^2 / 3  = 76329603384216526031706109802092473003
    const COFACTOR: &'static [u64] = &[0x8c00aaab0000aaab, 0x396c8c005555e156];

    /// COFACTOR_INV = COFACTOR^{-1} mod r
    /// = 52435875175126190458656871551744051925719901746859129887267498875565241663483
    #[rustfmt::skip]
    const COFACTOR_INV: Fr381 = field_new!(Fr381, "52435875175126190458656871551744051925719901746859129887267498875565241663483");

    /// AFFINE_GENERATOR_COEFFS = (G1_GENERATOR_X, G1_GENERATOR_Y)
    const AFFINE_GENERATOR_COEFFS: (Self::BaseField, Self::BaseField) =
        (G1_GENERATOR_X, G1_GENERATOR_Y);

    #[inline(always)]
    fn mul_by_a(_: &Self::BaseField) -> Self::BaseField {
        Self::BaseField::zero()
    }
}

/// G1_GENERATOR_X =
/// 3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507
#[rustfmt::skip]
pub const G1_GENERATOR_X: Fq381 = field_new!(Fq381, "3685416753713387016781088315183077757961620795782546409894578378688607592378376318836054947676345821548104185464507");

/// G1_GENERATOR_Y =
/// 1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724655353366534992391756441569
#[rustfmt::skip]
pub const G1_GENERATOR_Y: Fq381 = field_new!(Fq381, "1339506544944476473020471379941921221584933875938349620426543736416511423956333506472724655353366534992391756441569");

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

    const XI : Fq381 = field_new!(Fq381, "-1");
     const ZETA: Fq381 = field_new!(Fq381, "3");
     const XI_ON_ZETA_SQRT: Fq381 = field_new!(Fq381, "14");

 }

/// The point of the test is to get a  simpl SWU compatible curve
/// and make simple hash
#[test]
fn hash_arbitary_string_to_curve_swu() {
    use blake2::{VarBlake2b};

    let test_swu_to_curve_hasher = MapToCurveBasedHasher::<GroupAffine<TestSWUMapToCurveF127Params>, DefaultFieldHasher<VarBlake2b>, SWUMap<TestSWUMapToCurveF127Params>>::new(&[1]).unwrap();
    
    let hash_result = test_swu_to_curve_hasher.hash(b"if you stick a Babel fish in your ear you can instantly understand anything said to you in any form of language.").expect("fail to hash the string to curve");

    
    println!("{:?}, {:?}", hash_result, hash_result.x, );

    assert!(hash_result.x != field_new!(F127, "0"));

}
