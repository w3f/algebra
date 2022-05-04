use core::marker::PhantomData;

use crate::{models::SWModelParameters, ModelParameters};
use ark_ff::batch_inversion;
use ark_poly::{univariate::DensePolynomial, Polynomial, UVPolynomial};

use crate::{
    hashing::{map_to_curve_hasher::MapToCurve, HashToCurveError},
    models::short_weierstrass_jacobian::GroupAffine,
    AffineCurve,
};

use super::swu::{SWUMap, SWUParams};
type BaseField<MP> = <MP as ModelParameters>::BaseField;

/// Trait defining the necessary parameters for the WB hash-to-curve method
/// for the curves of Weierstrass form of:
/// of y^2 = x^3 + a*x + b where b != 0 but `a` can be zero like BLS-381 curve.
/// From [\[WB2019\]]
///
/// - [\[WB2019\]] <http://dx.doi.org/10.46586/tches.v2019.i4.154-179>
pub trait WBParams: SWModelParameters + Sized {
    // The isogenous curve should be defined over the same base field but it can have
    // different scalar field type IsogenousCurveScalarField :
    type IsogenousCurve: SWUParams<BaseField = BaseField<Self>>;

    const PHI_X_NOM: &'static [BaseField<Self>];
    const PHI_X_DEN: &'static [BaseField<Self>];

    const PHI_Y_NOM: &'static [BaseField<Self>];
    const PHI_Y_DEN: &'static [BaseField<Self>];

    fn isogeny_map(
        domain_point: GroupAffine<Self::IsogenousCurve>,
    ) -> Result<GroupAffine<Self>, HashToCurveError> {
        let x_num = DensePolynomial::from_coefficients_slice(&Self::PHI_X_NOM[..]);
        let x_den = DensePolynomial::from_coefficients_slice(&Self::PHI_X_DEN[..]);

        let y_num = DensePolynomial::from_coefficients_slice(&Self::PHI_Y_NOM[..]);
        let y_den = DensePolynomial::from_coefficients_slice(&Self::PHI_Y_DEN[..]);

        let mut v: [BaseField<Self>; 2] = [
            x_den.evaluate(&domain_point.x),
            y_den.evaluate(&domain_point.x),
        ];
        batch_inversion(&mut v);
        let img_x = x_num.evaluate(&domain_point.x) * v[0];
        let img_y = (y_num.evaluate(&domain_point.x) * domain_point.y) * v[1];

        Ok(GroupAffine::new(img_x, img_y, false))
    }
}

pub struct WBMap<P: WBParams> {
    swu_field_curve_hasher: SWUMap<P::IsogenousCurve>,
    curve_params: PhantomData<fn() -> P>,
}

impl<P: WBParams> MapToCurve<GroupAffine<P>> for WBMap<P> {
    /// Constructs a new map if `P` represents a valid map.
    fn new_map_to_curve() -> Result<Self, HashToCurveError> {
        // Verifying that the isogeny maps the generator of the SWU curve into us
        let isogenous_curve_generator = GroupAffine::<P::IsogenousCurve>::new(
            P::IsogenousCurve::AFFINE_GENERATOR_COEFFS.0,
            P::IsogenousCurve::AFFINE_GENERATOR_COEFFS.1,
            false,
        );

        match P::isogeny_map(isogenous_curve_generator) {
            Ok(point_on_curve) => {
                if !point_on_curve.is_on_curve() {
                    return Err(HashToCurveError::MapToCurveError(format!("the isogeny maps the generator of its domain: {} into {} which does not belong to its codomain.",isogenous_curve_generator, point_on_curve)));
                }
            },
            Err(e) => return Err(e),
        }

        Ok(WBMap {
            swu_field_curve_hasher: SWUMap::<P::IsogenousCurve>::new_map_to_curve().unwrap(),
            curve_params: PhantomData,
        })
    }

    /// Map random field point to a random curve point
    /// inspired from
    /// <https://github.com/zcash/pasta_curves/blob/main/src/hashtocurve.rs>
    fn map_to_curve(
        &self,
        element: <GroupAffine<P> as AffineCurve>::BaseField,
    ) -> Result<GroupAffine<P>, HashToCurveError> {
        // first we need to map the field point to the isogenous curve
        let point_on_isogenious_curve = self.swu_field_curve_hasher.map_to_curve(element).unwrap();
        P::isogeny_map(point_on_isogenious_curve)
    }
}

#[cfg(test)]
mod test {
    use crate::hashing::HashToCurve;
    use crate::{
        hashing::{
            curve_maps::{
                swu::SWUParams,
                wb::{WBMap, WBParams},
            },
            field_hashers::DefaultFieldHasher,
            map_to_curve_hasher::MapToCurveBasedHasher,
        },
        models::SWModelParameters,
        short_weierstrass_jacobian::GroupAffine,
        ModelParameters,
    };
    use ark_ff::{fields::Fp64, MontBackend, MontFp};

    #[derive(ark_ff::MontConfig)]
    #[modulus = "127"]
    #[generator = "6"]
    pub struct F127Config;
    pub type F127 = Fp64<MontBackend<F127Config, 1>>;

    const F127_ZERO: F127 = MontFp!(F127, "0");
    const F127_ONE: F127 = MontFp!(F127, "1");

    /// The struct defining our parameters for the target curve of hashing
    struct TestWBF127MapToCurveParams;

    impl ModelParameters for TestWBF127MapToCurveParams {
        const COFACTOR: &'static [u64] = &[1];

    #[rustfmt::skip]
        const COFACTOR_INV: F127 = F127_ONE;

        type BaseField = F127;
        type ScalarField = F127;
    }

    /// E: Elliptic Curve defined by y^2 = x^3 + 3 over Finite
    /// Field of size 127
    impl SWModelParameters for TestWBF127MapToCurveParams {
        /// COEFF_A = 0
        const COEFF_A: F127 = F127_ZERO;

        /// COEFF_B = 3
    #[rustfmt::skip]
        const COEFF_B: F127 = MontFp!(F127, "3");

        /// AFFINE_GENERATOR_COEFFS = (G1_GENERATOR_X, G1_GENERATOR_Y)
        const AFFINE_GENERATOR_COEFFS: (Self::BaseField, Self::BaseField) =
            (MontFp!(F127, "62"), MontFp!(F127, "70"));
    }

    /// Testing WB19 hashing on a small curve
    /// E_isogenous : Elliptic Curve defined by y^2 = x^3 + 109*x + 124 over Finite
    /// Field of size 127
    /// Isogenous to E : y^2 = x^3 + 3
    struct TestSWU127MapToIsogenousCurveParams;

    /// First we define the isogenous curve
    /// sage: E_isogenous.order()
    /// 127
    impl ModelParameters for TestSWU127MapToIsogenousCurveParams {
        const COFACTOR: &'static [u64] = &[1];

    #[rustfmt::skip]
        const COFACTOR_INV: F127 = F127_ONE;

        type BaseField = F127;
        type ScalarField = F127;
    }

    /// E_isogenous : Elliptic Curve defined by y^2 = x^3 + 109*x + 124 over Finite
    /// Field of size 127
    impl SWModelParameters for TestSWU127MapToIsogenousCurveParams {
        /// COEFF_A = 109
        const COEFF_A: F127 = MontFp!(F127, "109");

        /// COEFF_B = 124
    #[rustfmt::skip]
        const COEFF_B: F127 = MontFp!(F127, "124");

        /// AFFINE_GENERATOR_COEFFS = (G1_GENERATOR_X, G1_GENERATOR_Y)
        const AFFINE_GENERATOR_COEFFS: (Self::BaseField, Self::BaseField) =
            (MontFp!(F127, "84"), MontFp!(F127, "2"));
    }

    /// SWU parameters for E_isogenous
    impl SWUParams for TestSWU127MapToIsogenousCurveParams {
        /// NON-SQUARE = - 1
        const XI: F127 = MontFp!(F127, "-1");
        /// A Primitive Root of unity = 3
        const ZETA: F127 = MontFp!(F127, "3");
        /// sqrt(Xi/Zeta)
        const XI_ON_ZETA_SQRT: F127 = MontFp!(F127, "13");
    }

    /// E_isogenous : Elliptic Curve defined by y^2 = x^3 + 109*x + 124 over Finite
    /// Field of size 127
    /// With psi: E_isogenous -> E
    /// psi = (psi_x(x,y), psi_y(x,y))
    /// where
    /// psi_x: (-57*x^13 - 21*x^12 + 10*x^11 + 34*x^10 + 40*x^9 -
    /// 13*x^8 + 32*x^7 - 32*x^6 + 23*x^5 - 14*x^4 + 39*x^3 + 23*x^2 + 63*x +
    /// 4)/(x^12 - 13*x^11 + 11*x^10 - 33*x^9 - 30*x^8 + 30*x^7 + 34*x^6 - 44*x^5 +
    /// 63*x^4 - 20*x^3 - 10*x^2 + 31*x + 2)
    ///
    /// psi_y: (10*x^18*y + 59*x^17*y + 41*x^16*y + 48*x^15*y - 7*x^14*y + 6*x^13*y +
    /// 5*x^12*y + 62*x^11*y + 12*x^10*y + 36*x^9*y - 49*x^8*y - 18*x^7*y - 63*x^6*y
    /// - 43*x^5*y - 60*x^4*y - 18*x^3*y + 30*x^2*y - 57*x*y - 34*y)/(x^18 + 44*x^17
    /// - 63*x^16 + 52*x^15 + 3*x^14 + 38*x^13 - 30*x^12 + 11*x^11 - 42*x^10 - 13*x^9
    /// - 46*x^8 - 61*x^7 - 16*x^6 - 55*x^5 + 18*x^4 + 23*x^3 - 24*x^2 - 18*x + 32)
    impl WBParams for TestWBF127MapToCurveParams {
        type IsogenousCurve = TestSWU127MapToIsogenousCurveParams;

        const PHI_X_NOM: &'static [<Self::IsogenousCurve as ModelParameters>::BaseField] = &[
            MontFp!(F127, "4"),
            MontFp!(F127, "63"),
            MontFp!(F127, "23"),
            MontFp!(F127, "39"),
            MontFp!(F127, "-14"),
            MontFp!(F127, "23"),
            MontFp!(F127, "-32"),
            MontFp!(F127, "32"),
            MontFp!(F127, "-13"),
            MontFp!(F127, "40"),
            MontFp!(F127, "34"),
            MontFp!(F127, "10"),
            MontFp!(F127, "-21"),
            MontFp!(F127, "-57"),
        ];

        const PHI_X_DEN: &'static [<Self::IsogenousCurve as ModelParameters>::BaseField] = &[
            MontFp!(F127, "2"),
            MontFp!(F127, "31"),
            MontFp!(F127, "-10"),
            MontFp!(F127, "-20"),
            MontFp!(F127, "63"),
            MontFp!(F127, "-44"),
            MontFp!(F127, "34"),
            MontFp!(F127, "30"),
            MontFp!(F127, "-30"),
            MontFp!(F127, "-33"),
            MontFp!(F127, "11"),
            MontFp!(F127, "-13"),
            MontFp!(F127, "1"),
        ];

        const PHI_Y_NOM: &'static [<Self::IsogenousCurve as ModelParameters>::BaseField] = &[
            MontFp!(F127, "-34"),
            MontFp!(F127, "-57"),
            MontFp!(F127, "30"),
            MontFp!(F127, "-18"),
            MontFp!(F127, "-60"),
            MontFp!(F127, "-43"),
            MontFp!(F127, "-63"),
            MontFp!(F127, "-18"),
            MontFp!(F127, "-49"),
            MontFp!(F127, "36"),
            MontFp!(F127, "12"),
            MontFp!(F127, "62"),
            MontFp!(F127, "5"),
            MontFp!(F127, "6"),
            MontFp!(F127, "-7"),
            MontFp!(F127, "48"),
            MontFp!(F127, "41"),
            MontFp!(F127, "59"),
            MontFp!(F127, "10"),
        ];

        const PHI_Y_DEN: &'static [<Self::IsogenousCurve as ModelParameters>::BaseField] = &[
            MontFp!(F127, "32"),
            MontFp!(F127, "-18"),
            MontFp!(F127, "-24"),
            MontFp!(F127, "23"),
            MontFp!(F127, "18"),
            MontFp!(F127, "-55"),
            MontFp!(F127, "-16"),
            MontFp!(F127, "-61"),
            MontFp!(F127, "-46"),
            MontFp!(F127, "-13"),
            MontFp!(F127, "-42"),
            MontFp!(F127, "11"),
            MontFp!(F127, "-30"),
            MontFp!(F127, "38"),
            MontFp!(F127, "3"),
            MontFp!(F127, "52"),
            MontFp!(F127, "-63"),
            MontFp!(F127, "44"),
            MontFp!(F127, "1"),
        ];
    }

    /// The point of the test is to get a simple WB compatible curve
    /// and make simple hash
    #[test]
    fn hash_arbitrary_string_to_curve_wb() {
        use sha2::Sha256;
        let test_wb_to_curve_hasher = MapToCurveBasedHasher::<
            GroupAffine<TestWBF127MapToCurveParams>,
            DefaultFieldHasher<Sha256, 128>,
            WBMap<TestWBF127MapToCurveParams>,
        >::new(&[1])
        .unwrap();

        let hash_result = test_wb_to_curve_hasher.hash(b"if you stick a Babel fish in your ear you can instantly understand anything said to you in any form of language.").expect("fail to hash the string to curve");

        assert!(
            hash_result.x != F127_ZERO && hash_result.y != F127_ZERO,
            "we assume that not both a and b coefficienst are zero for the test curve"
        );

        assert!(
            hash_result.is_on_curve(),
            "hash results into a point off the curve"
        );
    }
}
