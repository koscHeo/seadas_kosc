/* --------------------------------------------------------------- */
/* get_par.c - computes photosynthetically active radiation        */
/*                                                                 */
/* Inputs:                                                         */
/*     l2rec - level-2 structure containing one complete scan      */
/*             after atmospheric correction.                       */
/* Outputs:                                                        */
/*     par   - Photosynthetically Active Radiation just above the  */
/*             surface from SeaWiFS level 1b instantaneous         */
/*             radiances at 412, 443, 490, 510, 555, and 670 nm.   */
/*                                                                 */
/* Algorithm Provided By: Robert Frouin,                           */
/*                        Scripps Institution of Oceanography      */
/*                                                                 */
/* Written By: B. A. Franz, SAIC GSC, SIMBIOS, 30 July 1999        */
/* Modified:                                                       */
/*  Fall 2013, Robert Lossing, SAIC - Adapted underlying code to C */
/* --------------------------------------------------------------- */

#include <stdlib.h>
#include <math.h>
#include "l12_proto.h"
#include "smi_climatology.h"
#include "par_utils.h"

static int* alloc_bindx( int nbands,  int* bindx) {
    if ((bindx = (int *) calloc(nbands, sizeof(int))) == NULL) {
        printf("-E- : Error allocating memory to bindx in get_par\n");
        exit(FATAL_ERROR);
    }
    return bindx;
}

void get_par(l2str *l2rec, float par[]) {
    int32_t ip, ib, ipb, iw;

    float *Lt;

    float angst;
    float taua;
    float delphi;

    float hour = *l2rec->msec / 3.6E6;

    static int32_t lastDay = -999;
    static int32_t month = -999;
    static int32_t monthDay = -999;

    static int32_t mask = SEAICE | LAND | HIGLINT;

    static float *lambda;
    static float *Fobar;
    static float *Taur;
    static float *kO3;
    static int *bindx;
    static int nwave;
    static int firstCall = TRUE;

    if (firstCall) {

        firstCall = FALSE;


        /* Initialize climatologies */
        smi_climatology_init(l2rec->input->alphafile, *l2rec->day, ALPHA510);
        smi_climatology_init(l2rec->input->tauafile, *l2rec->day, TAUA865);

        switch (l2rec->sensorID) {

        case HMODISA:
        case HMODIST:
            nwave = 3;
            bindx = alloc_bindx(nwave, bindx);
            bindx[0] = 2;
            bindx[1] = 6;
            bindx[2] = 7;
            break;
        case SEAWIFS:
            nwave = 6;
            bindx = alloc_bindx(nwave, bindx);
            for (ib = 0; ib < nwave; ib++)
                bindx[ib] = ib;
            break;
        case VIIRS:
            nwave = 5;
            bindx = alloc_bindx(nwave, bindx);
            for (ib = 0; ib < nwave; ib++)
                bindx[ib] = ib;
            break;
        case MERIS:
            nwave = 7;
            bindx = alloc_bindx(nwave, bindx);
           for (ib = 0; ib < nwave; ib++)
                bindx[ib] = ib;
            break;
        case OCTS:
            nwave = 6;
            bindx = alloc_bindx(nwave, bindx);
            for (ib = 0; ib < nwave; ib++)
                bindx[ib] = ib;
            break;
        default:
            printf("PAR not supported for this sensor (%d).\n",
                    l2rec->sensorID);
            exit(1);
            break;
        }
        if ( (lambda = (float *)calloc(nwave,sizeof(float))) == NULL) {
            printf("-E- : Error allocating memory to lambda in get_par\n");
            exit(FATAL_ERROR);
        }
        if ( (Fobar = (float *)calloc(nwave,sizeof(float))) == NULL) {
            printf("-E- : Error allocating memory to Fobar in get_par\n");
            exit(FATAL_ERROR);
        }
        if ( (Taur = (float *)calloc(nwave,sizeof(float))) == NULL) {
            printf("-E- : Error allocating memory to Taur in get_par\n");
            exit(FATAL_ERROR);
        }
        if ( (kO3 = (float *)calloc(nwave,sizeof(float))) == NULL) {
            printf("-E- : Error allocating memory to kO3 in get_par\n");
            exit(FATAL_ERROR);
        }

        /* Get band-pass dependent quantities */
        // consider passing in an array of band indexes to index l2rec->(lambda,k03,Fobar,Taur) in calc_par
        for (iw = 0; iw < nwave; iw++) {
            ib = bindx[iw];
            lambda[iw] = l2rec->fwave[ib];
            kO3[iw] = l2rec->k_oz[ib];
            Fobar[iw] = l2rec->Fobar[ib];
            Taur[iw] = l2rec->Tau_r[ib];
        }
    }
    if ( (Lt = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("-E- : Error allocating memory to Lt in get_par\n");
        exit(FATAL_ERROR);
    }

    for (ip = 0; ip < l2rec->npix; ip++) {

        /* Grab radiances for this pixel */
        for (ib = 0; ib < nwave; ib++) {
            ipb = ip * l2rec->nbands + bindx[ib];
            Lt[ib] = l2rec->Lt[ipb];
        }

        /* Skip pixel if masked */
        if (Lt[0] <= 0.0 || (l2rec->flags[ip] & mask) != 0
                || l2rec->solz[ip] > 90.0) {
            par[ip] = BAD_FLT;
            l2rec->flags[ip] |= PRODFAIL;
            continue;
        }

        /* Get angstrom and AOT from climatology */
        angst = smi_climatology(l2rec->lon[ip], l2rec->lat[ip], ALPHA510);
        taua = smi_climatology(l2rec->lon[ip], l2rec->lat[ip], TAUA865);

        switch (l2rec->sensorID) {

        case SEAWIFS:
        case HMODIST:
        case HMODISA:
        case MERIS:
        case VIIRS:
        case OCTS:

            par[ip] = calc_par(l2rec, ip, nwave, Lt, taua, angst, lambda,
                    Fobar, kO3, Taur);
            break;

        default:
            printf("PAR not supported for this sensor (%d).\n",
                    l2rec->sensorID);
            exit(1);
            break;
        }
        /* Convert to E/D/m^2 */
        if (par[ip] != BAD_FLT)
            par[ip] *= 1.193;

    }
    free(Lt);
}

/*
 Subject:
 PAR routine
 Date:
 Fri, 23 Jul 1999 08:55:29 -0700
 From:
 Robert Frouin <rfrouin@ucsd.edu>
 To:
 chuck@seawifs, gfargion@simbios, gene@seawifs, wang@simbios, franz@seawifs
 CC:
 jmcpherson@ucsd.edu




 Greetings:

 A routine to compute daily PAR from SeaWiFS level 1b radiances is available
 at the following address: http://genius.ucsd.edu/~john/SeaWiFS_dir/ under
 the rubrique "PAR subroutine and test program".

 The routine requires as input year, month, day, time, latitude, longitude,
 SeaWiFS radiances in the first 6 spectral bands, solar zenith angle,
 viewing zenith angle, relative azimuth angle, aerosol optical thickness at
 865 nm, Angstrom coefficient, ozone amount, and surface pressure. Routine
 output is daily PAR.

 Thus a daily PAR value is computed for each instantaneous SeaWiFS
 observation, clear or cloudy. Diurnal variations are taken into account
 statistically. The algorithm is described succintly in the routine.

 During our discussion at GSFC in June, a first routine was supposed to be
 developed to provide a normalized cloud/surface albedo and then a second
 routine to compute daily PAR from the normalized albedo, and the second
 routine was to be applied when binning to the 9 km resolution. Now daily
 PAR is obtained using a single routine, which is more convenient.

 The binning to the 9 km resolution should be done as follows. First,
 weight-average the daily PAR estimates obtained from all SeaWiFS
 observations during the same day at each location (there might be several
 SeaWiFS observations of a surface target during the same day). The weight
 is the cosine of the sun zenith angle for the SeaWiFS observation. That is:

 PAR_avg = sum{cos[tetas(i)]*PAR(i)}/sum{cos[tetas(i)]}

 Second, simply average the values at all the locations within the 9 km bins.

 The routine requires aerosol data, ozone amount, surface pressure. If these
 parameters are missing (-999 or less), default values are used. If Eric
 Vermote's aerosol climatology (or Menghua's) is not available yet, please
 use default values for tests.

 At this time, the statistical diurnal function does not depend on latitude,
 longitude, and date, but will depend on these parameters in the second
 version of the code. Creating a date and location dependent function
 requires analysing several years of ISCCP data. We have the data, but a
 couple of weeks is needed to accomplish the task.

 I will present the algorithm at the SeaWiFS atmospheric correction meeting
 next week, and prepare a detailed documentation.

 Best, Robert.


 Robert Frouin
 Scripps Institution of Oceanography
 University of California San Diego
 9500 Gilman Drive
 La Jolla, CA 92093-0221
 Voice Tel.: 619/534-6243
 Fax Tel.: 619/534-7452


 */
