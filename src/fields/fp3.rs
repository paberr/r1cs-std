use crate::fields::{cubic_extension::*, fp::FpVar};
use ark_ff::{
    fields::{CubicExtConfig, Fp3ConfigWrapper},
    Fp3Config, PrimeField,
};

use super::{FieldOpsBounds, FieldVar};

/// A cubic extension field constructed over a prime field.
/// This is the R1CS equivalent of `ark_ff::Fp3<P>`.
pub type Fp3Var<P> = CubicExtVar<FpVar<<P as Fp3Config>::Fp>, Fp3ConfigWrapper<P>>;

pub type GenericFp3Var<P, ConstraintF, BF> = CubicExtVar<BF, Fp3ConfigWrapper<P>, ConstraintF>;

impl<P: Fp3Config, ConstraintF: PrimeField, BF: FieldVar<P::Fp, ConstraintF>>
    CubicExtVarConfig<BF, ConstraintF> for Fp3ConfigWrapper<P>
where
    for<'a> &'a BF: FieldOpsBounds<'a, <P as Fp3Config>::Fp, BF>,
{
    fn mul_base_field_vars_by_frob_coeff(c1: &mut BF, c2: &mut BF, power: usize) {
        *c1 *= Self::FROBENIUS_COEFF_C1[power % Self::DEGREE_OVER_BASE_PRIME_FIELD];
        *c2 *= Self::FROBENIUS_COEFF_C2[power % Self::DEGREE_OVER_BASE_PRIME_FIELD];
    }
}
