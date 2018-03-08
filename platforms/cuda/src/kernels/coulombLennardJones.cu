{
const real invR2 = invR*invR;
const real invR6 = invR2*invR2*invR2;
const real r6 = r2*r2*r2;

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

	    real combinedA;
        real combinedB = buckingham1.y * buckingham2.y;

        if (combinedB == 0.0f) {
            combinedB = 0.0f;
            combinedA = 0.0f;
        } else {
            combinedB = combinedB / (buckingham1.y + buckingham2.y);
            combinedA = pow(buckingham1.x * buckingham1.y, 1.0f/buckingham1.y) * pow(buckingham2.x * buckingham2.y, 1.0f/buckingham2.y);
            combinedA = pow(combinedA, combinedB) / (2.0f * combinedB);
        }
        //	printf("SEC, B: %f, A: %f\n", combinedB, combinedA);
        real buckinghamExp = -2.0f * combinedB * r;
        real buckinghamRepulsion = combinedA * EXP(buckinghamExp);

    	real rvdw = sigmaEpsilon1.y + sigmaEpsilon2.y;
    	real rvdw2 = rvdw*rvdw;
    	real rvdw6 = rvdw2*rvdw2*rvdw2;

    	real c6D = 1.0/(r6+rvdw6);
    	real c8D = 1.0/(r6*r2 + rvdw6*rvdw2);
    	real c10D = 1.0/(r6*r2*r2 + rvdw6*rvdw2*rvdw2);

        // tempForce = buckinghamRepulsion * combinedB * r + 12.0f * c12 - 6.0f * c6 - 8.0f * c8 - 10.0f * c10;
        // real ljEnergy = buckinghamRepulsion + c12 - c6 - c8 - c10;

    	tempForce = buckinghamRepulsion * combinedB * r + 12.0f * c12 - (6.0f*c6*c6D*c6D*r6) - (8.0f*c8*c8D*c8D*r6*r2) - (10.0f*c10*c10D*c10D*r6*r2*r2);
    	real ljEnergy = buckinghamRepulsion + c12 - c6*c6D - c8*c8D - c10*c10D;

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

    real combinedB = buckingham1.y * buckingham2.y;
    combinedB = combinedB == 0.0f ? 0.0f : combinedB / (buckingham1.y + buckingham2.y);
    real combinedA = pow(buckingham1.x * buckingham1.y, 1.0f/buckingham1.y) * pow(buckingham2.x * buckingham2.y, 1.0f/buckingham2.y);
    combinedA = pow(combinedA, combinedB) / (2.0f * combinedB);
    // printf("B: %f, A: %f\n", combinedB, combinedA);
    real buckinghamExp = -2.0f * combinedB * r;
    real buckinghamRepulsion = combinedA * EXP(buckinghamExp);


    real rvdw = sigmaEpsilon1.y + sigmaEpsilon2.y;
    real rvdw2 = rvdw*rvdw;
    real rvdw6 = rvdw2*rvdw2*rvdw2;

    real c6D = 1.0/(r6+rvdw6);
    real c8D = 1.0/(r6*r2 + rvdw6*rvdw2);
    real c10D = 1.0/(r6*r2*r2 + rvdw6*rvdw2*rvdw2);

//    tempForce = buckinghamRepulsion * combinedB * r + 12.0f * c12 - 6.0f * c6 - 8.0f * c8 - 10.0f * c10;
//    real ljEnergy = includeInteraction ? buckinghamRepulsion + c12 - c6 - c8 - c10 : 0;

    tempForce = buckinghamRepulsion * combinedB * r + 12.0f * c12 - (6.0f*c6*c6D*c6D*r6) - (8.0f*c8*c8D*c8D*r6*r2) - (10.0f*c10*c10D*c10D*r6*r2*r2);
    real ljEnergy = includeInteraction ? buckinghamRepulsion + c12 - c6*c6D - c8*c8D - c10*c10D : 0;


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
