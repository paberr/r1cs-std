use crate::fields::{fp::FpVar, quadratic_extension::*};
use ark_ff::{
    fields::{Fp2Config, Fp2ConfigWrapper, QuadExtConfig},
    PrimeField,
};

use super::{nonnative::NonNativeFieldVar, FieldOpsBounds, FieldVar};

/// A quadratic extension field constructed over a prime field.
/// This is the R1CS equivalent of `ark_ff::Fp2<P>`.
pub type Fp2Var<P> = QuadExtVar<FpVar<<P as Fp2Config>::Fp>, Fp2ConfigWrapper<P>>;
pub type NonNativeFp2Var<P, ConstraintF> = GenericQuadExtVar<
    NonNativeFieldVar<<P as Fp2Config>::Fp, ConstraintF>,
    ConstraintF,
    Fp2ConfigWrapper<P>,
>;

pub type GenericFp2Var<P, ConstraintF, BF> =
    GenericQuadExtVar<BF, ConstraintF, Fp2ConfigWrapper<P>>;

impl<P: Fp2Config, ConstraintF: PrimeField, BF: FieldVar<P::Fp, ConstraintF>>
    QuadExtVarConfig<BF, ConstraintF> for Fp2ConfigWrapper<P>
where
    for<'a> &'a BF: FieldOpsBounds<'a, <P as Fp2Config>::Fp, BF>,
{
    fn mul_base_field_var_by_frob_coeff(fe: &mut BF, power: usize) {
        *fe *= Self::FROBENIUS_COEFF_C1[power % Self::DEGREE_OVER_BASE_PRIME_FIELD];
    }
}
