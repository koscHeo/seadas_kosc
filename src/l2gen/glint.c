#include "l12_proto.h"

/* ============================================================================

c       This corrects the sunglint contamination reflectances in the
c       SeaWiFS 8 bands using the Cox & Munk model for the ocean
c       sun glitter radiance distribution as function of the sea surface
c       wind.
c
c Inputs:
c       num_iter, I, --- iteration number in the atmospheric corrections.
c       nband, I, --- number of bands for the sensor (e.g., 8 for SeaWiFS).
c       glint_coef, R, --- glitter radiance (F0=1) from Cox & Munk.
c       air_mass, R, --- airmass value.
c       mu0, R, --- cosine of the solar zenith angle.
c       taur(nband), R, --- Rayleigh optical thicknesses.
c       taua(nband), R, --- retrieved aerosol optical thicknesses.
c       La(nband), R, --- aerosol reflectances at 8 SeaWiFS bands.
c Output:
c       TLg(nband), R, --- sunglint radiances at 8 SeaWiFS bands.
c Credits:
c       Written by
c        Menghua Wang
c        UMBC, Code 970.2, NASA/GSFC, Greenbelt, MD
c        10-7-1999
c
c       Modification History
c       - Modified to return glint radiance at TOA rather than to 
c         correct an input reflectance. Bryan Franz, 14 October 1999.
c       - Added low rhoa enhancements from M. Wang, 28 November 1999.
c	- Added a check and necessary correction to make sure there is 
c	  no over-correction, i.e., there is no pixel lost due to sun 
c	  glint corrections.  Menghua Wang, 2-18-00
c ============================================================================= */

static float taua_est(float x) { return(-0.8 - 0.4*log(x)); };

void glint_rad(int32_t iter, int32_t nband, int32_t nir_s, int32_t nir_l, 
               float glint_coef, float airmass, 
               float mu0, float F0[], float taur[], float  taua[], float La[], float TLg[])
{
    static float pi=PI;
    static float radeg=RADEG;
    static float glint_min=GLINT_MIN;
    static float taua_min=0.08;
    static float taua_ave=0.1;
    static float rhoa_min=0.01;
    static float rhoa_min2=0.008;
    static int32_t  iter_max=2;
    static float rfac=0.8;

    int32_t  ib; 
    float taua_c;
    float refl_test;
    float taua_ave2;
    float fac;


    /* If the number of iterations exceeds the maximum, we just return.  This assumes   */
    /* that the calling routine saved the TLg from the previous iteration. If the glint */
    /* coefficient is very low, return 0.                                               */

    if (iter > iter_max) 
        return;

    else if (glint_coef <= glint_min)
        for (ib=0; ib<nband; ib++)
	    TLg[ib] = 0.0;

    else {

        refl_test = pi/mu0*(La[nir_l]/F0[nir_l]-glint_coef*exp(-(taur[nir_l]+taua_ave)*airmass));

        if (refl_test <= rhoa_min)
	    taua_ave2 = taua_est(10.*MAX(refl_test,0.0001));
        else
  	    taua_ave2 = taua_ave;

	for (ib=nband-1; ib>=0; ib--) {

            if (iter <= 1)
	        taua_c = taua_ave2;
            else if (taua[nir_l] <= taua_min)
	        taua_c = taua_est(taua[nir_l]);
            else
	        taua_c = taua[ib];

            /* Check for over-correction */
            if (ib == nir_l) 
                refl_test = pi/mu0*(La[ib]/F0[ib] - glint_coef*exp(-(taur[ib]+taua_c)*airmass));

            if (refl_test <= rhoa_min2)
	        TLg[ib] = F0[ib]*glint_coef * exp(-(taur[ib]+1.5*taua_c)*airmass);
            else
	        TLg[ib] = F0[ib]*glint_coef * exp(-(taur[ib]+    taua_c)*airmass);
	}


        /* Make sure there is no over-correction */
        if (La[nir_l] > 0.0 && La[nir_s] > 0.0) {
	    fac = MAX(TLg[nir_l]/La[nir_l], TLg[nir_s]/La[nir_s]);
	    if (fac >= rfac) {
	        for (ib=0; ib<nband; ib++) {
		    TLg[ib] = rfac*TLg[ib]/fac;
		}
 
	    }

	}

    }
}





  
