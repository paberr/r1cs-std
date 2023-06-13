use crate::fields::{fp3::Fp3Var, quadratic_extension::*};
use ark_ff::{fields::fp6_2over3::*, Fp3Config, PrimeField, QuadExtConfig};

use super::{fp3::GenericFp3Var, FieldOpsBounds, FieldVar};

/// A sextic extension field constructed as the tower of a
/// quadratic extension over a cubic extension field.
/// This is the R1CS equivalent of `ark_ff::fp6_2over3::Fp6<P>`.
pub type Fp6Var<P> = QuadExtVar<Fp3Var<<P as Fp6Config>::Fp3Config>, Fp6ConfigWrapper<P>>;
pub type GenericFp6Var<P, ConstraintF, BF> = GenericQuadExtVar<
    GenericFp3Var<<P as Fp6Config>::Fp3Config, ConstraintF, BF>,
    ConstraintF,
    Fp6ConfigWrapper<P>,
>;

impl<
        P: Fp6Config,
        ConstraintF: PrimeField,
        BF: FieldVar<<P::Fp3Config as Fp3Config>::Fp, ConstraintF>,
    > QuadExtVarConfig<GenericFp3Var<P::Fp3Config, ConstraintF, BF>, ConstraintF>
    for Fp6ConfigWrapper<P>
where
    for<'b> &'b BF: FieldOpsBounds<'b, <<P as Fp6Config>::Fp3Config as Fp3Config>::Fp, BF>,
{
    fn mul_base_field_var_by_frob_coeff(
        fe: &mut GenericFp3Var<P::Fp3Config, ConstraintF, BF>,
        power: usize,
    ) {
        fe.c0 *= Self::FROBENIUS_COEFF_C1[power % Self::DEGREE_OVER_BASE_PRIME_FIELD];
        fe.c1 *= Self::FROBENIUS_COEFF_C1[power % Self::DEGREE_OVER_BASE_PRIME_FIELD];
        fe.c2 *= Self::FROBENIUS_COEFF_C1[power % Self::DEGREE_OVER_BASE_PRIME_FIELD];
    }
}
