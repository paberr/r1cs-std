use crate::fields::{fp2::Fp2Var, quadratic_extension::*};
use ark_ff::{
    fields::{Fp4ConfigWrapper, QuadExtConfig},
    Fp2Config, Fp4Config, PrimeField,
};

use super::{fp2::GenericFp2Var, FieldOpsBounds, FieldVar};

/// A quartic extension field constructed as the tower of a
/// quadratic extension over a quadratic extension field.
/// This is the R1CS equivalent of `ark_ff::Fp4<P>`.
pub type Fp4Var<P> = QuadExtVar<Fp2Var<<P as Fp4Config>::Fp2Config>, Fp4ConfigWrapper<P>>;

pub type GenericFp4Var<P, ConstraintF, BF> = GenericQuadExtVar<
    GenericFp2Var<<P as Fp4Config>::Fp2Config, ConstraintF, BF>,
    ConstraintF,
    Fp4ConfigWrapper<P>,
>;

impl<
        P: Fp4Config,
        ConstraintF: PrimeField,
        BF: FieldVar<<P::Fp2Config as Fp2Config>::Fp, ConstraintF>,
    > QuadExtVarConfig<GenericFp2Var<<P as Fp4Config>::Fp2Config, ConstraintF, BF>, ConstraintF>
    for Fp4ConfigWrapper<P>
where
    for<'b> &'b BF: FieldOpsBounds<'b, <<P as Fp4Config>::Fp2Config as Fp2Config>::Fp, BF>,
{
    fn mul_base_field_var_by_frob_coeff(
        fe: &mut GenericFp2Var<<P as Fp4Config>::Fp2Config, ConstraintF, BF>,
        power: usize,
    ) {
        fe.c0 *= Self::FROBENIUS_COEFF_C1[power % Self::DEGREE_OVER_BASE_PRIME_FIELD];
        fe.c1 *= Self::FROBENIUS_COEFF_C1[power % Self::DEGREE_OVER_BASE_PRIME_FIELD];
    }
}
