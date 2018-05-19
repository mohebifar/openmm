{
const real invR2 = invR*invR;
const real invR3 = invR2*invR;
const real invR4 = invR3*invR;
const real invR5 = invR4*invR;
const real invR6 = invR5*invR;
const real invR7 = invR6*invR;
const real invR8 = invR7*invR;
const real invR9 = invR8*invR;
const real invR10 = invR9*invR;

// const real r6 = r2*r2*r2;

const real d = sigmaEpsilon1.y + sigmaEpsilon2.y;
const real d2 = d * d * 0.5f;
const real d3 = d2 * d * 0.3333333333f;
const real d4 = d3 * d * 0.25f;
const real d5 = d4 * d * 0.2f;
const real d6 = d5 * d * 0.1666666667f;
const real d7 = d6 * d * 0.1428571429f;
const real d8 = d7 * d * 0.125f;
const real d9 = d8 * d * 0.1111111111f;
const real d10 = d9 * d * 0.1f;

const real mdr = -d * r;
real expTerm = exp(mdr);

real c6Deriv = 6.0f * invR6 + expTerm * (
    invR6 * (mdr - 6.0f) +
    d * invR5  * (mdr - 5.0f) +
    d2 * invR4 * (mdr - 4.0f) +
    d3 * invR3 * (mdr - 3.0f) +
    d4 * invR2 * (mdr - 2.0f) +
    d5 * invR * (mdr - 1.0f) +
    d6 * mdr
);

real c8Deriv = 8.0f * invR8 + expTerm * (
    invR8 * (mdr - 8.0f) +
    d * invR7 * (mdr - 7.0f) +
    d2 * invR6 * (mdr - 6.0f) +
    d3 * invR5 * (mdr - 5.0f) +
    d4 * invR4 * (mdr - 4.0f) +
    d5 * invR3 * (mdr - 3.0f) +
    d6 * invR2 * (mdr - 2.0f) +
    d7 * invR * (mdr - 1.0f) +
    d8 * mdr
);

real c10Deriv = 10.0f * invR10 + expTerm * (
    invR10 * (mdr - 10.0f) +
    d * invR9 * (mdr - 9.0f) +
    d2 * invR8 * (mdr - 8.0f) +
    d3 * invR7 * (mdr - 7.0f) +
    d4 * invR6 * (mdr - 6.0f) +
    d5 * invR5 * (mdr - 5.0f) +
    d6 * invR4 * (mdr - 4.0f) +
    d7 * invR3 * (mdr - 3.0f) +
    d8 * invR2 * (mdr - 2.0f) +
    d9 * invR * (mdr - 1.0f) +
    d10 * mdr
);

real c6E = invR6 - expTerm * (
    invR6 +
    d * invR5 +
    d2 * invR4 +
    d3 * invR3 +
    d4 * invR2 +
    d5 * invR +
    d6
);

real c8E = invR8 - expTerm * (
    invR8 +
    d * invR7 +
    d2 * invR6 +
    d3 * invR5 +
    d4 * invR4 +
    d5 * invR3 +
    d6 * invR2 +
    d7 * invR +
    d8
);

real c10E = invR10 - expTerm * (
    invR10 +
    d * invR9 +
    d2 * invR8 +
    d3 * invR7 +
    d4 * invR6 +
    d5 * invR5 +
    d6 * invR4 +
    d7 * invR3 +
    d8 * invR2 +
    d9 * invR +
    d10
);

// real combinedA;
real combinedA = buckingham1.x * buckingham2.x;
real combinedB = buckingham1.y * buckingham2.y * 0.5f;

/*
if (combinedB == 0.0f) {
            combinedB = 0.0f;
            combinedA = 0.0f;
} else {
	combinedB = combinedB / (buckingham1.y + buckingham2.y);
	combinedA = pow(buckingham1.x * buckingham1.y, 1.0f/buckingham1.y) * pow(buckingham2.x * buckingham2.y, 1.0f/buckingham2.y);
	combinedA = pow(combinedA, combinedB) / (2.0f * combinedB);
}
*/

//const real rvdw = 0.376243f;
#if USE_EWALD
    bool needCorrection = hasExclusions && isExcluded && atom1 != atom2 && atom1 < NUM_ATOMS && atom2 < NUM_ATOMS;
    unsigned int includeInteraction = ((!isExcluded && r2 < CUTOFF_SQUARED) || needCorrection);
    const real alphaR = EWALD_ALPHA*r;
    const real expAlphaRSqr = EXP(-alphaR*alphaR);
    const real prefactor = 138.935456f*posq1.w*posq2.w*invR;

#ifdef USE_DOUBLE_PRECISION
    const real erfcAlphaR = erfc(alphaR);
#else
    // This approximation for erfc is from Abramowitz and Stegun (1964) p. 299.  They cite the following as
    // the original source: C. Hastings, Jr., Approximations for Digital Computers (1955).  It has a maximum
    // error of 1.5e-7.

    const real t = RECIP(1.0f+0.3275911f*alphaR);
    const real erfcAlphaR = (0.254829592f+(-0.284496736f+(1.421413741f+(-1.453152027f+1.061405429f*t)*t)*t)*t)*t*expAlphaRSqr;
#endif
    real tempForce = 0.0f;
#if HAS_LENNARD_JONES
    // The multiplicative term to correct for the multiplicative terms that are always
    // present in reciprocal space.  The real terms have an additive contribution
    // added in, but for excluded terms the multiplicative term is just subtracted.
    // These factors are needed in both clauses of the needCorrection statement, so
    // I declare them up here.
    #if DO_LJPME
        const real dispersionAlphaR = EWALD_DISPERSION_ALPHA*r;
        const real dar2 = dispersionAlphaR*dispersionAlphaR;
        const real dar4 = dar2*dar2;
        const real dar6 = dar4*dar2;
        const real expDar2 = EXP(-dar2);
        const real c6 = cACoefficients1.x * cACoefficients2.x;
        const real coef = invR6*c6;
        const real eprefac = 1.0f + dar2 + 0.5f*dar4;
        const real dprefac = eprefac + dar6/6.0f;
    #endif
#endif
    if (needCorrection) {
        // Subtract off the part of this interaction that was included in the reciprocal space contribution.

        if (1-erfcAlphaR > 1e-6) {
            real erfAlphaR = ERF(alphaR); // Our erfc approximation is not accurate enough when r is very small, which happens with Drude particles.
            tempForce = -prefactor*(erfAlphaR-alphaR*expAlphaRSqr*TWO_OVER_SQRT_PI);
            tempEnergy += -prefactor*erfAlphaR;
        }
        else {
            includeInteraction = false;
            tempEnergy -= TWO_OVER_SQRT_PI*EWALD_ALPHA*138.935456f*posq1.w*posq2.w;
        }
#if HAS_LENNARD_JONES
        #if DO_LJPME
            // The multiplicative grid term
            tempEnergy += coef*(1.0f - expDar2*eprefac);
            tempForce += 6.0f*coef*(1.0f - expDar2*dprefac);
        #endif
#endif
    }
    else {
#if HAS_LENNARD_JONES
        // real sig = sigmaEpsilon1.x + sigmaEpsilon2.x;
        // real sig2 = invR*sig;
        // sig2 *= sig2;
        // real sig6 = sig2*sig2*sig2;
        // real eps = sigmaEpsilon1.y*sigmaEpsilon2.y;
        // real epssig6 = sig6*eps;

        real c6 = cACoefficients1.x * cACoefficients2.x;
        real c8 = cACoefficients1.y * cACoefficients2.y;
        real c10 = cBCoefficients1.x * cBCoefficients2.x;
        real c12 = cBCoefficients1.y * cBCoefficients2.y * invR6 * invR6;

        real buckinghamExp = -2.0f * combinedB * r;
        real buckinghamRepulsion = combinedA * EXP(buckinghamExp);

    	// tempForce = -buckinghamExp * buckinghamRepulsion + 12.0f * c12 - 6.0f * invR6 * c6 - 8.0f * invR8 * c8 - 10.0f * invR10 * c10;
    	// real ljEnergy = buckinghamRepulsion + c12 - invR6 * c6 - invR8 * c8 - invR10 * c10;

    	tempForce = -buckinghamExp * buckinghamRepulsion + 12.0f * c12 - c6Deriv * c6 - c8Deriv * c8 - c10Deriv * c10;
    	real ljEnergy = buckinghamRepulsion + c12 - c6E * c6 - c8E * c8 -c10E * c10;

//	if (combinedB != 0.0f) {
//		printf("SEC, r: %f, Force: %f, Energy: %f, rep: %f \n", r, tempForce, ljEnergy, buckinghamRepulsion);
//	}

        #if USE_LJ_SWITCH
        if (r > LJ_SWITCH_CUTOFF) {
            real x = r-LJ_SWITCH_CUTOFF;
            real switchValue = 1+x*x*x*(LJ_SWITCH_C3+x*(LJ_SWITCH_C4+x*LJ_SWITCH_C5));
            real switchDeriv = x*x*(3*LJ_SWITCH_C3+x*(4*LJ_SWITCH_C4+x*5*LJ_SWITCH_C5));
            tempForce = tempForce*switchValue - ljEnergy*switchDeriv*r;
            ljEnergy *= switchValue;
        }
        #endif
#if DO_LJPME
        // The multiplicative grid term
        ljEnergy += coef*(1.0f - expDar2*eprefac);
        tempForce += 6.0f*coef*(1.0f - expDar2*dprefac);
        // The potential shift accounts for the step at the cutoff introduced by the
        // transition from additive to multiplicative combintion rules and is only
        // needed for the real (not excluded) terms.  By addin these terms to ljEnergy
        // instead of tempEnergy here, the includeInteraction mask is correctly applied.
        // sig2 = sig*sig;
        // sig6 = sig2*sig2*sig2*INVCUT6;
        // epssig6 = eps*sig6;

        // BUCKINGHAM - COMMENTED OUT THE POTENTIAL SHIFT

        // The additive part of the potential shift
        // ljEnergy += epssig6*(1.0f - sig6);
        // ljEnergy += buckinghamRepulsion + c12 - c6 - c8 - c10;
        // The multiplicative part of the potential shift
        // ljEnergy += MULTSHIFT6*c6;
#endif
        tempForce += prefactor*(erfcAlphaR+alphaR*expAlphaRSqr*TWO_OVER_SQRT_PI);
        tempEnergy += includeInteraction ? ljEnergy + prefactor*erfcAlphaR : 0;
#else
        tempForce = prefactor*(erfcAlphaR+alphaR*expAlphaRSqr*TWO_OVER_SQRT_PI);
        tempEnergy += includeInteraction ? prefactor*erfcAlphaR : 0;
#endif
    }
    dEdR += includeInteraction ? tempForce*invR*invR : 0;
#else
#ifdef USE_CUTOFF
    unsigned int includeInteraction = (!isExcluded && r2 < CUTOFF_SQUARED);
#else
    unsigned int includeInteraction = (!isExcluded);
#endif
    real tempForce = 0.0f;
  #if HAS_LENNARD_JONES
    // real sig = sigmaEpsilon1.x + sigmaEpsilon2.x;
    // real sig2 = invR*sig;
    // sig2 *= sig2;
    // real sig6 = sig2*sig2*sig2;
    // real epssig6 = sig6*(sigmaEpsilon1.y*sigmaEpsilon2.y);
    // tempForce = epssig6*(12.0f*sig6 - 6.0f);
    // real ljEnergy = includeInteraction ? epssig6*(sig6 - 1) : 0;

    real c6 = cACoefficients1.x * cACoefficients2.x;
    real c8 = cACoefficients1.y * cACoefficients2.y;
    real c10 = cBCoefficients1.x * cBCoefficients2.x;
    real c12 = cBCoefficients1.y * cBCoefficients2.y * invR6 * invR6;

    // printf("B: %f, A: %f\n", combinedB, combinedA);
    real buckinghamExp = -2.0f * combinedB * r;
    real buckinghamRepulsion = combinedA * EXP(buckinghamExp);

    // tempForce = -buckinghamExp * buckinghamRepulsion + 12.0f * c12 - 6.0f * invR6 * c6 - 8.0f * invR8 * c8 - 10.0f * invR10 * c10;
    // real ljEnergy = includeInteraction ? buckinghamRepulsion + c12 - invR6 * c6 - invR8 * c8 - invR10 * c10 : 0;

    tempForce = -buckinghamExp * buckinghamRepulsion + 12.0f * c12 - c6Deriv *c6 - c8Deriv * c8 - c10Deriv * c10;
    real ljEnergy = includeInteraction ? buckinghamRepulsion + c12 - c6E * c6 - c8E *c8 -c10E * c10 : 0;


    #if USE_LJ_SWITCH
    if (r > LJ_SWITCH_CUTOFF) {
        real x = r-LJ_SWITCH_CUTOFF;
        real switchValue = 1+x*x*x*(LJ_SWITCH_C3+x*(LJ_SWITCH_C4+x*LJ_SWITCH_C5));
        real switchDeriv = x*x*(3*LJ_SWITCH_C3+x*(4*LJ_SWITCH_C4+x*5*LJ_SWITCH_C5));
        tempForce = tempForce*switchValue - ljEnergy*switchDeriv*r;
        ljEnergy *= switchValue;
    }
    #endif
    tempEnergy += ljEnergy;
  #endif
#if HAS_COULOMB
  #ifdef USE_CUTOFF
    const real prefactor = 138.935456f*posq1.w*posq2.w;
    tempForce += prefactor*(invR - 2.0f*REACTION_FIELD_K*r2);
    tempEnergy += includeInteraction ? prefactor*(invR + REACTION_FIELD_K*r2 - REACTION_FIELD_C) : 0;
  #else
    const real prefactor = 138.935456f*posq1.w*posq2.w*invR;
    tempForce += prefactor;
    tempEnergy += includeInteraction ? prefactor : 0;
  #endif
#endif
    dEdR += includeInteraction ? tempForce*invR*invR : 0;
#endif
}
