/******************************************************************************
 * This routine computes daily Photosynthetically Active Radiation (PAR) just
 * above the surface.
 *
 * The atmosphere and surface are modeled as a 2-layer system. The first layer
 * contains molecules and aerosols and is located above the second layer, the
 * cloud/surface layer. The TOA radiances are transformed into reflectance
 * and corrected for gaseous absorption. The resulting reflectances are
 * further corrected for molecular and aerosol scattering and combined to
 * yield the instantaneous albedo of the cloud/surface layer in the PAR
 * interval, A. This albedo is then used to compute instantaneous PAR. Daily
 * PAR is finally obtained by integrating instantaneous PAR from sunrise to
 * sunset, taking into account statistically diurnal effects on A.
 *
 *      *** Important: This subroutine uses an external file containing
 *      aerosol optical properties.
 *
 *      Output:
 *      PAR: Photosynthetically Active Radiation (mW/cm^2/um)
 *
 *    Key changes in v 2.0:
 *      * wavelength integration weighting function added
 *
 *    Key changes in v 1.8:
 *      * AsatMax renamed to AcldMax to indicate "cloud" emphasis
 *      * Use (Asat-AsBar) rather than Asat in the test that determines if
 *        Asat is too large
 *      * Added 670 section to aerosols table since we need that to choose
 *        the right model in GetPhase subroutine, + modified code in that
 *        routine
 *      * Changed TauCons from 20 to 15
 *
 *    Key changes in v 1.7:
 *      * Corrected Tau term in spherical albedo expressions
 *      * Modified initial Asat definition with AsBar
 *      * Put 'cloudy' test after computation of new Asat
 *      * Replaced archaic intrinsic routines (COSD, SIND, ACOSD, JMOD, LONG)
 *
 *    Key changes in v 1.6:
 *      * Angstrom exponents are expected to be primarily positive,
 *        that is, in the range (-0.1, 2.0) ... initial test was removed.
 *      * Surface albedo tables were added, with interpolation subroutines
 *      * TauAer(500nm) calculation was added (for Surface albedo tables)
 *      * Replaced AsBar calculation with call to As interpolation subroutine
 *      * Added test to estimate if pixel is clear or cloudy
 *      * Replaced AsAvg calculation with call to As interpolation subroutine
 *      * Removed AsBar and AsAvg from calculation of Aavg, in integration
 *      * Added Abar vs. AsAvg test for "function" integrand term
 *
 *    Key changes in v 1.5:
 *      * addition of water vapor term in gaseous transmittance in integration
 *      * comparison test for Abar vs. AsAvg
 *
 *      Authors: Robert Frouin <RFrouin@ucsd.edu>, scientific algorithms
 *              John McPherson <JMcPherson@ucsd.edu>, program structure
 *
 *      Adapted from the original FORTRAN to C and better integrated with l2gen
 *      by Robert Lossing and Sean Bailey
 *
 * References:
 *
 * Frouin, R., and B. Chertock, 1992: A technique for global monitoring of net
 * solar irradiance at the ocean surface. Part I: Model. J. Appl. Meteo., 31,
 * 1056-1066.
 *
 * Frouin, R., D. W. Ligner, and C. Gautier, 1989: A Simple analytical formula
 * to compute clear sky total and photosynthetically available solar irradiance
 * at the ocean surface. J. Geophys. Res., 94, 9731-9742.
 *
 * Tanre, D., M. Herman, P.-Y. Deschamps, and A. De Leffe, 1979: Atmospheric
 * modeling for Space measurements of ground reflectances, including
 * bi-directional properties. Appl. Optics, 18, 21,3587-21,3597.
 *
 * Young, D. F., P. Minnis, D. R. Doelling, G. G. Gibson, and T. Wong:
 * Temporal interpolation methods for the Clouds and the Earth's Radiant
 * Energy System (CERES) Experiment. J. Appl. Meteo., 37, 572-590.
 *
 *****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "l2_struc.h"
#include "par_utils.h"
#include <timeutils.h>
#include "l12_proto.h"

float calc_par(l2str *l2rec, int ip, int nbands, float *Lt, float taua,
        float angstrom, float *wl, float *fo, float *ko3, float *taumolbar) {

    int i, j, band;
    float r[nbands], csolz, transg;
    float weight[nbands];
    float rp[nbands], transa, sa, ps0, del;
    float asymfac, beta, tau[nbands], taumol[nbands], tauaer[nbands];
    float gamma, cos2x, pmol;
    float csenz, sinsz, sinvz, cosra, transav;

    float par0, integtrans, delt, saavg;
    float tmstep, sz, cossz, tg, td;
    float sola;
    float asavg, func[11], ra[nbands], rc[nbands];

    float sumSa, sumrc, sumtdir, sumtdif, denom, tatemp, ravgsat;
    float tdir, tdif, tdiravg, tdifavg, astemp;
    float q, q4, taucons, q4tau, f, asat, asbar, abar, acldmax;
    float fcsolz, fcsenz, fcossz;
    float tg_oz, tg_wv;
    float trise, tset;
    float par;

    int widx490, widx510;
    float ang500, ta500, slope;

    int cloudy = FALSE;

    static float kavg_oz = 5.2e-05;
    static float kavg_wv = 0.002;

//	need some memory allocation
    float *aerphasefunc;
    float *omega;
    float *modelAngstrom;

    int16 month;
    int16 mday;

    aerphasefunc = (float*) allocateMemory(nbands * sizeof(float),
            "aerphasefunc");
    omega = (float*) allocateMemory(nbands * sizeof(float), "omega");
    modelAngstrom = (float*) allocateMemory(nbands * sizeof(float),
            "modelAngstrom");

    int asat_is_max = 0;

    /****************************************************************************
     * Perform data checks; assign defaults if necessary: 					    *
     * If data invalid and uncorrectable, PAR set to BAD_FLT missing value.	    *
     ***************************************************************************/

    if (l2rec->solz[ip] < 0.0f || l2rec->solz[ip] > 90.0f) {
        printf("Error: subroutine CalcPAR. computation failed.\n");
        printf("SolZen (outside valid range 0-90) = %f\n", l2rec->solz[ip]);
        par = BAD_FLT;
        return par;
    }

    if (l2rec->senz[ip] < 0.0f || l2rec->senz[ip] > 90.0f) {
        printf("Error: subroutine CalcPAR. computation failed.\n");
        printf("ViewZen (outside valid range 0-90) = %f\n", l2rec->senz[ip]);
        par = BAD_FLT;
        return par;
    }

    if (l2rec->delphi[ip] < -180.0f || l2rec->delphi[ip] > 180.0f) {
        printf("Error: subroutine CalcPAR. computation failed.\n");
        printf("delphi (outside valid range -180:180) = %f\n",
                l2rec->delphi[ip]);
        par = BAD_FLT;
        return par;
    }

    float dobson = l2rec->oz[ip];

    csolz = l2rec->csolz[ip];
    csenz = l2rec->csenz[ip];
    cos2x = powf(cosf((l2rec->scattang[ip] * M_PI / 180.0)), 2.0);

    if (l2rec->scattang[ip] < 0.0 || l2rec->scattang[ip] > 180.0) {
        printf("Error: subroutine CalcPAR. computation failed.\n");
        printf("scatang (outside valid range 0-180) = %f\n",
                l2rec->scattang[ip]);
        printf("scatang = -(cossz*zosvz)+(sinsz*sinvz*cosra)\n");
        par = BAD_FLT;
        return par;
    }

    if (taua <= 0.0f) {
        taua = 0.1f;
    }

    float surfpress = l2rec->pr[ip];
    if (surfpress <= 0.0f) {
        surfpress = 1013.15f;
    }

    /* Need to use conventional Angstrom exponent.
     if it's missing GetAerPhase should work and return a default value */

    GetAerPhase(l2rec, ip, nbands, angstrom, aerphasefunc, omega,
            modelAngstrom);

    /****************************************************************************
     * Angstrom exponents are conventionally > 0, but actually < 0 			    *
     ****************************************************************************/

    angstrom = -1 * angstrom;

    yd2md(*l2rec->year, *l2rec->day, &month, &mday);

    if (dobson < 0.0f) {
        dobson = EstimateDobson(*l2rec->year, month, mday, l2rec->lat[ip]);
    }

    float watvap = l2rec->wv[ip];
    if (watvap < 0.0f) {
        watvap = EstimateWatVap(*l2rec->year, month, mday, l2rec->lat[ip]);
    }

    /*****************************************************************************
     * Transform TOA radiances into reflectances:			                     *
     ****************************************************************************/
    // TODO: consider passing in rhot instead of Lt
    float Airmass = (1.f / csolz) + (1.f / csenz);

    for (band = 0; band < nbands; band++) {
        r[band] = (M_PI * (Lt[band])) / (fo[band] * l2rec->fsol * csolz);
        // Compute gaseous transmittance:
        transg = expf((-ko3[band] * dobson * Airmass));
        rp[band] = r[band] / transg;
    }

    /****************************************************************************
     * Compute reflectances of the cloud-surface subsystem:				        *
     ***************************************************************************/

    /* Setup for transmittance equation:*/

    ps0 = 1013.5f;
    del = 0.0095f;

    pmol = (2.0f * (1.0f - del) * 0.75f * (1.0f + cos2x) + 3.0f * del)
            / (2.0f + del);

    asymfac = 0.6666667f;
    beta = 0.5f * (1.0f + asymfac);
    gamma = 1.0f - asymfac;

    /******************************************************************************
     * Compute diffuse atmospheric transmittance:								  *
     ******************************************************************************/
    sumSa = 0.0f;
    sumrc = 0.0f;
    sumtdir = 0.0f;
    sumtdif = 0.0f;
    denom = 0.0f;

    weight[0] = (wl[0] - 400.0) + (wl[1] - wl[0]) / 2.0;
    for (band = 0; band < nbands - 1; band++) {
        weight[band] = (wl[band] - wl[band - 1]) / 2.0
                + (wl[band + 1] - wl[band]) / 2.0;
    }
    weight[nbands - 1] = (700.0 - wl[nbands - 1])
            + (wl[nbands - 1] - wl[nbands - 2]) / 2.0;

    for (band = 0; band < nbands; band++) {
        taumol[band] = taumolbar[band] * (surfpress / ps0);

        tauaer[band] = taua
                * powf(wl[band] / (float) l2rec->input->aer_wave_long,
                        modelAngstrom[band]);

        tau[band] = taumol[band] + gamma * tauaer[band];

        ra[band] = ((taumol[band] * pmol)
                + (omega[band] * tauaer[band] * aerphasefunc[band]))
                / (4.0f * csolz * csenz);

        /* Compute diffuse atmospheric transmittance:*/

        transa = expf(-(taumol[band] + tauaer[band]) / csolz)
                * expf((0.52f * taumol[band] + beta * tauaer[band]) / csolz);
        transav = expf(-(taumol[band] + tauaer[band]) / csenz)
                * expf((0.52f * taumol[band] + beta * tauaer[band]) / csenz);

        // Get direct and diffuse transmittances
        tdir = expf(-tau[band] / csolz);

        tdif = (expf(-tau[band] / csolz))
                * (expf((0.52f * taumol[band] + beta * tauaer[band]) / csolz)
                        - 1.f);
        sumtdir += (tdir * fo[band]);
        sumtdif += (tdif * fo[band]);

        /* Compute spherical albedo of the atmosphere:*/
        sa = expf(-(taumol[band] + gamma * tauaer[band]))
                * (0.92f * taumol[band] + gamma * tauaer[band]);

        sumSa += sa * fo[band] * weight[band];
        denom += fo[band] * weight[band];

        rc[band] = (rp[band] - ra[band])
                / ((transa * transav) + (sa * (rp[band] - ra[band])));
        sumrc += (rc[band] * fo[band] * weight[band]);
    }

    //Derive tauaer(500) from (469,555)
    // indexes into modelAngstrom need to use bindex
    switch (l2rec->sensorID) {
    case HMODIST:
    case HMODISA:
        widx490 = windex(469., wl, nbands);
        widx510 = windex(555., wl, nbands);
        break;
    case VIIRS:
        widx490 = windex(490., wl, nbands);
        widx510 = windex(555., wl, nbands);
        break;
    case OCTS:
        widx490 = windex(490., wl, nbands);
        widx510 = windex(520., wl, nbands);
        break;
    default:
        widx490 = windex(490., wl, nbands);
        widx510 = windex(510., wl, nbands);
        break;
    }

    slope = (modelAngstrom[widx510] - modelAngstrom[widx490])
            / (wl[widx510] - wl[widx490]);
    ang500 = modelAngstrom[widx490] + slope * (wl[widx510] - wl[widx490]);
    ta500 = taua * powf((500.0 / (float) l2rec->input->aer_wave_long), ang500);

    /******************************************************************************
     * Compute daily average PAR:												  *
     ******************************************************************************/

    /* Set up terms constant in the integration:*/
    // Sa_bar:
    saavg = sumSa / denom;
    // R_bar:
    ravgsat = sumrc / denom;

    /******************************************************************************
     *Get simple bi-directional reflectance factor:								*
     ******************************************************************************/

    q = 1.0f / (3.0f * (1.0f - 0.853f));
    q4 = 4.0f * q;

    taucons = 15.0f;
    q4tau = q4 / (taucons + q4);

    fcsenz = (3.f / 7.f) * (1.0f + 2.0f * csenz);
    fcsolz = (3.f / 7.f) * (1.0f + 2.0f * csolz);

    f = (1.0f - q4tau * fcsolz)
            / (0.49f * ((1.f + 4.f * csolz * csenz) / (csolz + csenz))
                    - q4tau * fcsolz * fcsenz);

    /******************************************************************************
     * Get As_bar:																 *
     ******************************************************************************/
    asbar = interp_as_taulow(csolz, ta500);

    /******************************************************************************
     * Get Albedo from R_bar:													  *
     ******************************************************************************/
    asat = (ravgsat - asbar) * f + asbar;

    /******************************************************************************
     * Test Asat to see if it's too large; adjust as necessary:					  *
     ******************************************************************************/

    acldmax = 1.0f - (q4 / (300.0f + q4)) * fcsolz;

    if ((asat - asbar) >= acldmax) {
        asat_is_max = 1;
        taucons = 300.0f;
        q4tau = q4 / (taucons + q4);
    } else {
        asat_is_max = 0;
        taucons = q4 * ((fcsolz / (1.0f - (asat - asbar))) - 1.0f);

        /***************************************************************************
         * Compute new Asat if needed:											   *
         ***************************************************************************/

        if (taucons > 15.0f) {
            q4tau = q4 / (taucons + q4);

            f = (1.0f - q4tau * fcsolz)
                    / (0.49f * ((1.f + 4.f * csolz * csenz) / (csolz + csenz))
                            - q4tau * fcsolz * fcsenz);
            asat = (ravgsat - asbar) * f + asbar;

            /* New Test:*/

            if ((asat - asbar) >= acldmax) {
                asat_is_max = 1;
                taucons = 300.0f;
                q4tau = q4 / (taucons + q4);
            }
        } else {
            taucons = 15.0f;
            q4tau = q4 / (taucons + q4);
        }
    }

    // Test to estimate if pixel is clear or cloudy:
    cloudy = FALSE;
    if (asat > asbar)
        cloudy = TRUE;

    /*****************************************************************************
     * Get parameters used in integration time-steps:							 *
     *****************************************************************************/

    /*  Get the rise and set time for this day (for integral):*/

    triseset(l2rec->day, l2rec->lon[ip], l2rec->lat[ip], &trise, &tset);

    /*  Set up for 10 intervals*/
    delt = (tset - trise) / 10.0f;

    /***************************************************************************
     *	Get parameters used in integration time-steps:					       *
     ***************************************************************************/

    for (i = 0; i < 11; i++) {
        tmstep = trise + i * delt;
        sunangs_(l2rec->year, l2rec->day, &tmstep, l2rec->lon + ip,
                l2rec->lat + ip, &sz, &sola);
        cossz = cosf(sz * (float) M_PI / 180.0f);
        fcossz = (3.f / 7.f) * (1.0f + 2.0f * cossz);

        /***************************************************************************
         * Get gaseous transmittance term:										   *
         ***************************************************************************/
        if (sz >= 89.5f) {
            tg_oz = 0.0f;
            tg_wv = 0.0f;
        } else {
            //kavg_oz (5.2E-5) is valid for dobson in true dobson units, not the ozone units in l2gen
            tg_oz = expf(-kavg_oz * powf((dobson * 1000. / cossz), 0.99f));
            tg_wv = expf(-kavg_wv * powf((watvap / cossz), 0.87f));
        }

        tg = tg_oz * tg_wv;
        /***************************************************************************
         * Get average diffuse atmospheric transmittance and As_bar terms:         *
         ***************************************************************************/
        td = 0.0f;

        if (sz < 89.5f) {
            for (band = 0; band < nbands; band++) {
                tatemp = expf(-(taumol[band] + tauaer[band]) / cossz)
                        * expf(
                                (0.52f * taumol[band] + beta * tauaer[band])
                                        / cossz);
                td += tatemp * fo[band] * weight[band];
            }
            td /= denom;
        }

        if (cossz < 0.05) {
            cossz = 0.05;
        }
        if (cloudy) {
            asavg = interp_as_tauhigh(cossz);
        } else {
            asavg = interp_as_taulow(cossz, ta500);
        }

        /*****************************************************************************
         *Get A_bar term:															 *
         *****************************************************************************/

        if (asat_is_max) {
            abar = 1.0f - q4tau * fcossz;
        } else {
            abar = asat * (1.0f - q4tau * fcossz) / (1.0f - q4tau * fcsolz);
        }

        /*new test:*/
        if (abar >= 1.0f) {
            taucons = 300.0f;
            q4tau = q4 / (taucons + q4);
            abar = 1.0f - q4tau * fcossz;
        }

        /****************************************************************************
         * Get "function" integrand term:											*
         ****************************************************************************/

        if (abar <= asavg) {
            func[i] = cossz * tg * td / (1.0 - abar * saavg);
        } else {
            func[i] = cossz * tg * td * (1.0f - abar) / (1.0f - asavg)
                    / (1.0f - abar * saavg);
        }
    }

    /*	Use trapezoidal rule for integration:*/

    integtrans = 0.0f;

    for (i = 0; i < 10; i++) {
        integtrans += delt * (func[i] + func[i + 1]) / 2.0f;
    }

    /* Finally, determine the PAR:*/
    // par0 determined as the mean value [400-700nm] for Thuillier 2003
    // from run/data/common/Thuillier_F0.dat
    par0 = 176.323f;
    par = par0 * l2rec->fsol * integtrans / 24.0f;

    free(aerphasefunc);
    free(omega);
    free(modelAngstrom);

    return par;
}
