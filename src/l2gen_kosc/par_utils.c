/******************************************************************************
 *
 *			Photosynthetic Available Radiation Utilities
 *				Sensor-Independent Shared Subroutine
 *					by Robert J. Lossing
 *						July 30, 2013
 *				Ocean Color Processing Group
 *
 *******************************************************************************/

/******************************************************************************
 *
 *      Retrieve Phase function and Omega, given input scattering angle
 *        and Angstrom coefficient.
 *
 *      Authors: Robert Frouin <RFrouin@ucsd.edu>, scientific algorithms
 *              John McPherson <JMcPherson@ucsd.edu>, program structure
 *
 *      Note: Angstrom exponents in tables aren't in increasing order,
 *        need to determine proper order to use:
 *   			Model  omega(6)       alpha(6,8)     true_order
 *      			1   1.0000000E+00 -8.9877717E-02  2nd
 *     			2   9.9999797E-01 -9.1842167E-02  1st *min
 *     			3   9.8501903E-01  5.1216149E-01  8
 *      			4   9.8840600E-01  4.0654629E-01  7
 *     			5   9.9607497E-01  1.9749035E-01  4
 *      			6   9.9880499E-01  8.1851527E-02  3rd
 *      			7   9.7776997E-01  7.7916741E-01  9
 *     			8   9.9351001E-01  3.9742991E-01  6
 *    			9   9.9790698E-01  2.1781628E-01  5
 *     			10  9.5734900E-01  1.6471386E+00  12th *max
 *     			11  9.8160601E-01  1.4880587E+00  11th
 *     			12  9.9180502E-01  1.2930322E+00  10th
 *
 *******************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "par_utils.h"
#define NMODEL 12
#define NS 15
#define NT 9

static float *tablewavelengths;
static float *tablephaseangles;
static float *tablealphas;  //angstrom exponent
static float *tableomegas;  //single scattering albedo
static float *tableaerphasefunc;

void GetAerPhase(l2str *l2rec, int ip, int32_t nbands, float angstrom,
        float *aerphasefunc, float *omega, float *modelAngstrom) {

    float phaseangle = l2rec->scattang[ip];
    float distx, dista, p1, p2;
    int i, j, ix1, ix2, ia1, ia2, widx;

    static int firstCall = TRUE;

    // The data table is not in monotonically-increasing order, so here are the indices in proper order:

    int ii[12] = { 1, 0, 5, 4, 8, 7, 3, 2, 6, 11, 10, 9 };

    // Read the aerosol data file:

    if (firstCall) {
        firstCall = FALSE;
        tablewavelengths = (float*) allocateMemory(nbands * sizeof(float),
                "tablewavelengths");
        tablephaseangles = (float*) allocateMemory(NPHASE * sizeof(float),
                "tablephaseangles");
        tablealphas = (float*) allocateMemory(nbands * NMODEL * sizeof(float),
                "tablealphas");
        tableomegas = (float*) allocateMemory(nbands * NMODEL * sizeof(float),
                "tableomegas");
        tableaerphasefunc = (float*) allocateMemory(
                NPHASE * nbands * NMODEL * sizeof(float), "tableaerphasefunc");
        read_aerosol_par(l2rec, nbands, tablewavelengths, tablephaseangles,
                tablealphas, tableomegas, tableaerphasefunc);
    }

    // Find out which tablephaseangles-values in the table will be used:

    if (phaseangle < tablephaseangles[0]) {
//        printf("Input phase angle too small, setting to 0. Was %f", phaseangle);
        phaseangle = tablephaseangles[0];
        ix1 = 0;
    } else if (phaseangle >= tablephaseangles[NPHASE - 1]) {
//        printf("Input phase angle too large, setting to 180. Was %f",
//                phaseangle);
        phaseangle = tablephaseangles[NPHASE - 1];
        ix1 = NPHASE - 2;
    } else {
        for (i = NPHASE - 2; i > 0; i--) {
            if (phaseangle >= tablephaseangles[i]) {
                ix1 = i;
                break;
            }
        }
    }
    ix2 = ix1 + 1;

    distx = (tablephaseangles[ix2] - phaseangle)
            / (tablephaseangles[ix2] - tablephaseangles[ix1]);

    // Find out which models in the table to use:

    widx = windex(670.0, tablewavelengths, nbands);

    if (angstrom < tablealphas[widx * NMODEL + ii[0]]) {
        angstrom = tablealphas[widx * NMODEL + ii[0]];
        ia1 = 0;
    } else if (angstrom >= tablealphas[widx * NMODEL + ii[NMODEL - 1]]) {
        angstrom = tablealphas[widx * NMODEL + ii[NMODEL - 1]];
        ia1 = NMODEL - 2;
    } else {
        for (j = NMODEL - 2; j > 0; j--) {
            if (angstrom >= tablealphas[widx * NMODEL + ii[j]]) {
                ia1 = j;
                break;
            }
        }
    }
    ia2 = ia1 + 1;
    dista = (tablealphas[widx * NMODEL + ii[ia2]] - angstrom)
            / (tablealphas[widx * NMODEL + ii[ia2]]
                    - tablealphas[widx * NMODEL + ii[ia1]]);

    // Loop through the wavelengths and perform interpolations:

    for (i = 0; i < nbands; i++) {
        omega[i] = (tableomegas[i * NMODEL + ii[ia1]] * dista)
                + (tableomegas[i * NMODEL + ii[ia2]] * (1.0 - dista));

        modelAngstrom[i] = (tablealphas[i * NMODEL + ii[ia1]] * dista)
                + (tablealphas[i * NMODEL + ii[ia2]] * (1.0 - dista));

        p1 =
                (tableaerphasefunc[i * NMODEL * NPHASE + ii[ia1] * NPHASE + ix1]
                        * distx)
                        + (tableaerphasefunc[i * NMODEL * NPHASE
                                + ii[ia1] * NPHASE + ix2] * (1.0 - distx));

        p2 =
                (tableaerphasefunc[i * NMODEL * NPHASE + ii[ia2] * NPHASE + ix1]
                        * distx)
                        + (tableaerphasefunc[i * NMODEL * NPHASE
                                + ii[ia2] * NPHASE + ix2] * (1.0 - distx));

        aerphasefunc[i] = (p1 * dista) + (p2 * (1.0 - dista));
    }

    return;
}

/********************************************************************************
 *
 *	Description:
 *       This is to read the SeaWiFS aerosol data file for the PAR.
 *                               Menghua Wang 6/22/99
 *	Parameters:
 *       angl, R, ---- scattering angles.
 *       omega0, R, --- single scattering albedo.
 *       alpha, R, --- Angstrom coefficient alpha(lambda,865).
 *       s11, R, --- aerosol phase function.
 *
 *	SeaWiFS aerosol models:
 *       1-2: Oceanic RH=90 and 99%
 *       3-6: Maritime RH=50%, 70, 90, and 99%
 *       7-9: Coastal RH=50, 90, and 99%
 *       10-12: Tropospheric RH=50, 90, and 99%
 *
 *******************************************************************************/
void read_aerosol_par(l2str *l2rec, int32_t nbands, float *tablewavelengths,
        float *tablephaseangles, float *tablealphas, float *tableomegas,
        float *tableaerphasefunc) {

    int imodel[NMODEL];

    char aerosolfilename[FILENAME_MAX] = "";
    char *dataroot;
    FILE *aerosoldata;
    //Get OCDATAROOT path from user:
    if ((dataroot = getenv("OCDATAROOT")) == NULL ) {
        printf("OCDATAROOT environment variable is not defined.\n");
        return;
    }
    switch (l2rec->sensorID) {
    case HMODIST:
    case HMODISA:
        sprintf(aerosolfilename,
                "%s/modis/aerosol/modis_aerosol_par_3bands.dat", dataroot);
        break;
    default:
        sprintf(aerosolfilename,
                "%s/%s/aerosol/%s_aerosol_par.dat", dataroot, sensorDir[l2rec->sensorID],sensorDir[l2rec->sensorID]);
        break;
    }

    // Opens the aerosol data file for reading.
    if ((aerosoldata = fopen(aerosolfilename, "r")) == NULL) {
        printf("Error reading aerosol table for PAR from %s.\n",aerosolfilename);
        exit(EXIT_FAILURE);
    }
    printf("Loading aerosol properties for PAR from %s.\n",aerosolfilename);

    int j, jj, iph;
    int num_model;
    for (jj = 0; jj < NPHASE; jj++) {
        if ((fscanf(aerosoldata, "%f", &tablephaseangles[jj])) == 0) {
            printf("Problem reading phase angles from %s.\n", aerosolfilename);
            exit(EXIT_FAILURE);
        }
    }
    for (j = 0; j < nbands; j++) {
        for (iph = 0; iph < NMODEL; iph++) {
            if ((fscanf(aerosoldata, "%d %d %f", &imodel[iph], &num_model,
                    &tablewavelengths[j])) == 0) {
                printf("Problem reading model number and wavelength from %s.\n",
                        aerosolfilename);
                exit(EXIT_FAILURE);
            }
            if ((fscanf(aerosoldata, "%f %f", &tableomegas[j * NMODEL + iph],
                    &tablealphas[j * NMODEL + iph])) == 0) {
                printf("Problem reading SSA and Angstrom from %s.\n",
                        aerosolfilename);
                exit(EXIT_FAILURE);
            }
            for (jj = 0; jj < NPHASE; jj++) {
                if ((fscanf(aerosoldata, "%f",
                        &tableaerphasefunc[j * NMODEL * NPHASE + iph * NPHASE
                                + jj])) == 0) {
                    printf("Problem reading phase function values from %s.\n",
                            aerosolfilename);
                    exit(EXIT_FAILURE);
                }
            }
        }

    }
    fclose(aerosoldata);
}

/********************************************************************************
 *
 *	  Estimate Dobson units from climatology, given the month and
 *      latitude.  Table has 12 columns, one per month, and 35 rows,
 *      from 85 N to -85, in steps of 5 deg lat.
 *      Perform 2d interpolation (latitude,month_day) to get estimate.
 *
 *     Author: John McPherson <JMcPherson@ucsd.edu>
 *
 *      Dobson look-up table from Keating et al., Ozone reference models
 *      for the middle atmosphere, Handbook for MAP, Vol. 31, 1-36, 1989
 *
 */

float EstimateDobson(int32_t year, int32_t month, int32_t day, float lat) {

    int i1, i2, m1, m2;
    float daymid, d1, d2;
    float t, t1, t2, fac, fac2, difflat;
    float dobson;

    int tabdobson[12][35] = { { 395, 395, 395, 395, 395, 392, 390, 387, 376,
            354, 322, 292, 269, 254, 248, 246, 247, 251, 255, 260, 266, 271,
            277, 286, 295, 306, 319, 334, 344, 344, 338, 331, 324, 320, 316 },

    { 433, 433, 433, 436, 432, 428, 426, 418, 402, 374, 338, 303, 278, 261, 251,
            246, 248, 250, 254, 258, 262, 265, 270, 278, 286, 294, 303, 313,
            322, 325, 324, 317, 306, 299, 294 },

    { 467, 470, 460, 459, 451, 441, 433, 420, 401, 377, 347, 316, 291, 271, 260,
            254, 254, 255, 257, 259, 261, 264, 269, 277, 284, 289, 296, 305,
            312, 315, 317, 312, 305, 299, 295 },

    { 467, 465, 462, 455, 444, 431, 421, 410, 395, 373, 348, 325, 304, 287, 275,
            267, 261, 259, 258, 259, 260, 263, 271, 278, 284, 289, 297, 306,
            314, 318, 319, 313, 302, 302, 302 },

    { 411, 414, 416, 415, 410, 406, 402, 394, 382, 363, 342, 324, 307, 291, 279,
            271, 264, 260, 258, 257, 258, 264, 271, 281, 291, 303, 312, 318,
            322, 323, 322, 322, 322, 322, 322 },

    { 371, 371, 370, 368, 367, 372, 375, 372, 360, 341, 323, 311, 301, 290, 282,
            275, 268, 263, 259, 256, 258, 264, 273, 289, 306, 319, 327, 328,
            328, 337, 337, 337, 337, 337, 337 },

    { 333, 332, 332, 334, 338, 346, 350, 346, 335, 321, 310, 302, 296, 289, 284,
            280, 274, 268, 262, 259, 261, 268, 279, 295, 315, 331, 340, 342,
            338, 344, 340, 340, 340, 340, 340 },

    { 311, 308, 308, 313, 320, 327, 330, 326, 319, 310, 303, 298, 291, 286, 283,
            281, 277, 273, 268, 264, 266, 274, 288, 306, 327, 343, 353, 355,
            351, 339, 325, 307, 294, 294, 294 },

    { 283, 291, 302, 308, 312, 317, 318, 313, 307, 300, 295, 290, 284, 279, 279,
            279, 278, 276, 272, 270, 273, 282, 295, 313, 333, 348, 360, 367,
            368, 353, 324, 291, 267, 253, 230 },

    { 299, 299, 299, 309, 315, 317, 317, 312, 302, 291, 283, 280, 275, 270, 268,
            267, 263, 263, 265, 269, 277, 287, 301, 317, 336, 354, 371, 387,
            402, 402, 374, 333, 294, 274, 259 },

    { 314, 314, 314, 314, 332, 332, 327, 322, 311, 297, 284, 276, 270, 263, 261,
            260, 258, 259, 264, 270, 278, 286, 298, 311, 323, 335, 350, 366,
            381, 390, 388, 376, 357, 346, 341 },

    { 358, 358, 358, 358, 358, 358, 353, 349, 338, 320, 299, 281, 267, 256, 252,
            251, 251, 253, 257, 264, 272, 279, 287, 297, 307, 318, 332, 347,
            358, 365, 366, 364, 358, 356, 353 } };

    int days[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

    //   Set up for the time interpolation:

    if (year % 4 == 0)
        days[1] = days[1] + 1;

    daymid = days[month - 1] / 2.0;

    if (day >= daymid) {
        m1 = month - 1;
        m2 = m1 + 1;
        if (m2 > 12)
            m2 = 1;
        t = day;
    } else {
        m2 = month - 1;
        m1 = m2 - 1;
        if (m1 < 1)
            m1 = 12;
        t = days[m1] + day;
    }

    t1 = days[m1] / 2.0;
    t2 = days[m1] + (days[m2] / 2.0);

    // Perform the lat. interpolation, and further set up for time interp.:

    if (lat >= 85.0) {
        d1 = tabdobson[m1][1];
        d2 = tabdobson[m2][1];
    } else if (lat <= -85.0) {
        d1 = tabdobson[m1][35];
        d2 = tabdobson[m2][35];
    } else {
        i1 = 17 - (int) (lat / 5.0);
        i2 = i1 + 1;
        fac = (tabdobson[month][i2] - tabdobson[month][i1]) / (-5.0);
        difflat = lat - (90.0 - (i1 * 5.0));
        d1 = tabdobson[m1][i1] + (fac * difflat);
        d2 = tabdobson[m2][i1] + (fac * difflat);
    }

    // Complete the time interpolation to get final estimate:

    fac2 = (d2 - d1) / (t2 - t1);

    dobson = d1 + (t - t1) * fac2;
    return dobson;
}

/* ****************************************************************************
 *
 *       Interpolate water vapor table, as for dobson (ozone).
 *       Table has 12 columns, one per month, and 35 rows,
 *       from 85 N to 85 S, in steps of 5 deg lat.
 *       winter peak (Jan 15), summer peak (July 15), north hemis
 *       Perform 2d interpolation (latitude,month_day) to get estimate.
 *
 *      Tropical: if ( ABS( Lat ) .LE. 20.0 ) WatVap = 4.12
 *
 *      Mid-latitude case: if ( ABS( Lat ) .LE. 60.0 ) Then
 *        if ( (North.And.Summer) .Or. (South.And.Winter) ) WatVap = 2.93
 *        else WatVap = 0.85
 *
 *      Sub-arctic case:
 *         if ( (North.And.Summer) .Or. (South.And.Winter) ) WatVap = 2.10
 *         else WatVap = 0.42
 *
 *      after McClatchey, R.A., et al, _Optical Properties of the Atmosphere_,
 *      AFCRL 71-0279 354, 91 pp., 1971
 *
 */

float EstimateWatVap(int32_t year, int32_t month, int32_t day, float lat) {

    int i1, i2, m1, m2;
    month -= 1;
    float daymid, d1, d2;
    float t, t1, t2, fac, fac2, difflat;
    float watvap;

    float tabwv[12][35] = { { 0.42, 0.47, 0.52, 0.56, 0.61, 0.66, 0.71, 0.75,
            0.80, 0.85, 1.26, 1.67, 2.08, 2.48, 2.89, 3.30, 3.71, 4.12, 3.97,
            3.82, 3.67, 3.53, 3.38, 3.23, 3.08, 2.93, 2.84, 2.75, 2.65, 2.56,
            2.47, 2.38, 2.28, 2.19, 2.10 },

    { 0.70, 0.76, 0.81, 0.87, 0.92, 0.98, 1.03, 1.09, 1.14, 1.20, 1.56, 1.93,
            2.29, 2.66, 3.02, 3.39, 3.75, 4.12, 3.93, 3.74, 3.54, 3.35, 3.16,
            2.97, 2.78, 2.58, 2.50, 2.41, 2.33, 2.24, 2.16, 2.07, 1.99, 1.90,
            1.82 },

    { 0.98, 1.04, 1.11, 1.17, 1.23, 1.29, 1.36, 1.42, 1.48, 1.54, 1.87, 2.19,
            2.51, 2.83, 3.15, 3.48, 3.80, 4.12, 3.88, 3.65, 3.41, 3.18, 2.94,
            2.71, 2.47, 2.24, 2.16, 2.08, 2.00, 1.93, 1.85, 1.77, 1.69, 1.62,
            1.54 },

    { 1.26, 1.33, 1.40, 1.47, 1.54, 1.61, 1.68, 1.75, 1.82, 1.89, 2.17, 2.45,
            2.73, 3.01, 3.28, 3.56, 3.84, 4.12, 3.84, 3.56, 3.28, 3.01, 2.73,
            2.45, 2.17, 1.89, 1.82, 1.75, 1.68, 1.61, 1.54, 1.47, 1.40, 1.33,
            1.26 },

    { 1.54, 1.62, 1.69, 1.77, 1.85, 1.93, 2.00, 2.08, 2.16, 2.24, 2.47, 2.71,
            2.94, 3.18, 3.41, 3.65, 3.88, 4.12, 3.80, 3.48, 3.15, 2.83, 2.51,
            2.19, 1.87, 1.54, 1.48, 1.42, 1.36, 1.29, 1.23, 1.17, 1.11, 1.04,
            0.98 },

    { 1.82, 1.90, 1.99, 2.07, 2.16, 2.24, 2.33, 2.41, 2.50, 2.58, 2.78, 2.97,
            3.16, 3.35, 3.54, 3.74, 3.93, 4.12, 3.75, 3.39, 3.02, 2.66, 2.29,
            1.93, 1.56, 1.20, 1.14, 1.09, 1.03, 0.98, 0.92, 0.87, 0.81, 0.76,
            0.70 },

    { 2.10, 2.19, 2.28, 2.38, 2.47, 2.56, 2.65, 2.75, 2.84, 2.93, 3.08, 3.23,
            3.38, 3.53, 3.67, 3.82, 3.97, 4.12, 3.71, 3.30, 2.89, 2.49, 2.08,
            1.67, 1.26, 0.85, 0.80, 0.75, 0.71, 0.66, 0.61, 0.56, 0.52, 0.47,
            0.42 },

    { 1.82, 1.90, 1.99, 2.07, 2.16, 2.24, 2.33, 2.41, 2.50, 2.58, 2.78, 2.97,
            3.16, 3.35, 3.54, 3.74, 3.93, 4.12, 3.75, 3.39, 3.02, 2.66, 2.29,
            1.93, 1.56, 1.20, 1.14, 1.09, 1.03, 0.98, 0.92, 0.87, 0.81, 0.76,
            0.70 },

    { 1.54, 1.62, 1.69, 1.77, 1.85, 1.93, 2.00, 2.08, 2.16, 2.24, 2.47, 2.71,
            2.94, 3.18, 3.41, 3.65, 3.88, 4.12, 3.80, 3.48, 3.15, 2.83, 2.51,
            2.19, 1.87, 1.54, 1.48, 1.42, 1.36, 1.29, 1.23, 1.17, 1.11, 1.04,
            0.98 },

    { 1.26, 1.33, 1.40, 1.47, 1.54, 1.61, 1.68, 1.75, 1.82, 1.89, 2.17, 2.45,
            2.73, 3.01, 3.28, 3.56, 3.84, 4.12, 3.84, 3.56, 3.28, 3.01, 2.73,
            2.45, 2.17, 1.89, 1.82, 1.75, 1.68, 1.61, 1.54, 1.47, 1.40, 1.33,
            1.26 },

    { 0.98, 1.04, 1.11, 1.17, 1.23, 1.29, 1.36, 1.42, 1.48, 1.54, 1.87, 2.19,
            2.51, 2.83, 3.15, 3.48, 3.80, 4.12, 3.88, 3.65, 3.41, 3.18, 2.94,
            2.71, 2.47, 2.24, 2.16, 2.08, 2.00, 1.93, 1.85, 1.77, 1.69, 1.62,
            1.54 },

    { 0.70, 0.76, 0.81, 0.87, 0.92, 0.98, 1.03, 1.09, 1.14, 1.20, 1.56, 1.93,
            2.29, 2.66, 3.02, 3.39, 3.75, 4.12, 3.93, 3.74, 3.54, 3.35, 3.16,
            2.97, 2.78, 2.58, 2.50, 2.41, 2.33, 2.24, 2.16, 2.07, 1.99, 1.90,
            1.82 } };

    int days[12] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

    // Set up for the time interpolation:

    if (year % 4 == 0)
        days[1] = days[1] + 1;

    daymid = days[month] / 2.0;

    if (day >= daymid) {
        m1 = month;
        m2 = m1 + 1;
        if (m2 > 11)
            m2 = 0;
        t = day;
    } else {
        m2 = month;
        m1 = m2 - 1;
        if (m1 < 0)
            m1 = 11;
        t = days[m1] + day;
    }

    t1 = days[m1] / 2.0;
    t2 = days[m1] + (days[m2] / 2.0);

    // Perform the lat. interpolation, and further set up for time interp.:

    if (lat >= 85.0) {
        d1 = tabwv[m1][0];
        d2 = tabwv[m2][0];
    } else if (lat <= (-85.0)) {
        d1 = tabwv[m1][34];
        d2 = tabwv[m2][34];
    }

    else {
        i1 = 16 - (int) (lat / 5.0);
        i2 = i1 + 1;
        fac = (tabwv[month][i2] - tabwv[month][i1]) / (-5.0);
        difflat = lat - (90.0 - (i1 * 5.0));
        d1 = tabwv[m1][i1] + (fac * difflat);
        d2 = tabwv[m2][i1] + (fac * difflat);
    }

    // Complete the time interpolation to get final estimate:

    fac2 = (d2 - d1) / (t2 - t1);

    watvap = d1 + (t - t1) * fac2;

    return watvap;
}

/********************************************************************************
 *       		Get Earth-Sun distance correction	(VarSol)
 *
 *            Calculation of the variability of the solar constant during the year.
 *            jday is the number of the day in the month
 *            dsol is a multiplicative factor to apply to the mean value of
 *            solar constant.
 *
 ********************************************************************************/

float varsol(int32_t jday) {

    float om;
    float dsol;

    om = (0.9856f * (float) (jday - 4)) * M_PI / 180.f;
    dsol = 1.f / powf((1.f - .01673f * cosf(om)), 2.f);

    return dsol;
}

/*************************************************************************************
 *        Compute the rise and set time for the given Julian day of year	         *
 *************************************************************************************/

void triseset(int32_t *jday, float xlon, float xlat, float *trise, float *tset) {

    float fac, xlo, xla, xj, a1, a2, et, a3, delta, cosah, ah1, ah2, tsv1, tsv2,
            tsm1, tsm2;

    fac = (float) M_PI / 180.f;

    xlo = xlon * fac;
    xla = xlat * fac;
    xj = (float) (*jday - 1);

    a1 = (1.00554f * xj - 6.28306f) * fac;
    a2 = (1.93946f * xj + 23.35089f) * fac;
    et = -7.67825f * sinf(a1) - 10.09176f * sinf(a2);
    a3 = (0.9683f * xj - 78.00878f) * fac;

    delta = 23.4856f * sinf(a3) * fac;
    cosah = (-(sinf(xla) * sinf(delta))) / (cosf(xla) * cosf(delta));

    if ((cosah < -1.f) || (cosah > 1.f)) {
        *trise = 0.0f;
        *tset = 24.0f;
        return;
    }

    ah1 = acosf(cosah);
    ah2 = -ah1;

    tsv1 = ah1 / (15.f * fac);
    tsv2 = ah2 / (15.f * fac);
    tsm1 = tsv1 + 12.f - et / 60.f;
    tsm2 = tsv2 + 12.f - et / 60.f;

    *trise = tsm2 - (xlo / (15.f * fac));
    *tset = tsm1 - (xlo / (15.f * fac));
    return;

}

/********************************************************************************
 *     Get Sun position
 ********************************************************************************/

float get_solz(int32_t jday, float time, float lon, float lat) {
    float asol;
    double tsm, xla, xj, tet, a1, a2, a3, a4, a5, et, tsv, ah, b1, b2, b3, b4,
            b5, b6, b7, delta, amuzero, elev, az, caz, azim, pi2;

    //     solar position (zenithal angle asol,azimuthal angle phi0 in degrees)
    //    j is the day number in the year
    //
    //   mean solar time (heure decimale) (aka decimal time)
    double fac = M_PI / 180.;
    tsm = time + lon / 15.;
    xla = lat * fac;
    xj = (float) (jday - 1);
    tet = 2. * M_PI * xj / 365.;

    //    time equation (in mn.dec)
    a1 = 0.000075;
    a2 = 0.001868;
    a3 = 0.032077;
    a4 = 0.014615;
    a5 = 0.040849;
    et = a1 + a2 * cos(tet) - a3 * sin(tet) - a4 * cos(2. * tet)
            - a5 * sin(2. * tet);
    et = et * 12. * 60. / M_PI;

    //     true solar time

    tsv = tsm + et / 60.;
    tsv -= 12.;

    //     hour angle

    ah = tsv * 15. * fac;

    //     solar declination (in radian)

    b1 = 0.006918;
    b2 = 0.399912;
    b3 = 0.070257;
    b4 = 0.006758;
    b5 = 0.000907;
    b6 = 0.002697;
    b7 = 0 - .001480;
    delta = b1 - b2 * cos(tet) + b3 * sin(tet) - b4 * cos(2. * tet)
            + b5 * sin(2. * tet) - b6 * cos(3. * tet) + b7 * sin(3.f * tet);

    //     elevation,azimuth

    amuzero = sin(xla) * sin(delta) + cos(xla) * cos(delta) * cos(ah);
    elev = asin(amuzero);
    az = cos(delta) * sin(ah) / cos(elev);
    if ((fabs(az) - 1.000) > 0.00000)
        az = copysign(az, 1.);

    caz = (-cos(xla) * sin(delta) + sin(xla) * cos(delta) * cos(ah))
            / cos(elev);
    azim = asin(az);

    if (caz <= 0.)
        azim = M_PI - azim;

    if (caz > 0. && az <= 0.)
        azim = 2. * M_PI + azim;

    azim = azim + M_PI;
    pi2 = 2.f * M_PI;

    if (azim > pi2)
        azim -= pi2;

    elev /= fac;

    //    conversion in degrees

    asol = (float) (90. - elev);
//	phi0 = azim / fac;

    return asol;

}
/******************************************************************************
 *      Given inputs cos(SZ) and tauA500, interpolate the appropriate 		  *
 *       surface albedo from the table. This is for the case where tau < 1.0	  *
 *       Inputs are assumed to be in the valid ranges (the associated solar 	  *
 *       zenith angle should be in range 0.0 - 87.134 degrees, and tau		  *
 *      should be in the range 0.0 - 0.99)									  *
 ******************************************************************************/

float interp_as_taulow(float csz, float tau) {

    int s1, s2, t1, t2, i;

    float stot, sdist, ttot, tdist, slope, alb1[2], as;

    float cszt[NS] = { 0.05, 0.09, 0.15, 0.21, 0.27, 0.33, 0.39, 0.45, 0.52,
            0.60, 0.68, 0.76, 0.84, 0.92, 1.00 };

    float taut[NT] = { 0.00, 0.05, 0.10, 0.16, 0.24, 0.35, 0.50, 0.70, 0.99 };

    //       Store the table here rather than read a new file. Fortran stores
    //      2-d arrays internally by columns (col1, col2, col3 ...)

    float ast[NS][NT] = { { 0.1592253, 0.1258933, 0.1088253, 0.0968594,
            0.0880966, 0.0815492, 0.0766266, 0.0721790, 0.0679746 }, {
            0.1944972, 0.1637040, 0.1432702, 0.1240363, 0.1066057, 0.0922638,
            0.0826043, 0.0757601, 0.0700167 }, { 0.1978650, 0.1785707,
            0.1643212, 0.1482773, 0.1304572, 0.1119496, 0.0963537, 0.0839010,
            0.0740948 }, { 0.1785221, 0.1673987, 0.1589092, 0.1486409,
            0.1359193, 0.1209508, 0.1061981, 0.0920763, 0.0793196 }, {
            0.1531512, 0.1473233, 0.1426633, 0.1366985, 0.1289257, 0.1188203,
            0.1078114, 0.0957847, 0.0831087 }, { 0.1281988, 0.1255739,
            0.1234628, 0.1204020, 0.1161201, 0.1101539, 0.1030820, 0.0943791,
            0.0841410 }, { 0.1061354, 0.1052812, 0.1048179, 0.1036617,
            0.1016918, 0.0986726, 0.0948040, 0.0894514, 0.0822926 }, {
            0.0877530, 0.0881279, 0.0883850, 0.0884469, 0.0880587, 0.0870923,
            0.0854952, 0.0826607, 0.0782380 }, { 0.0708031, 0.0720173,
            0.0727886, 0.0736250, 0.0741367, 0.0746769, 0.0747173, 0.0741294,
            0.0722523 }, { 0.0567974, 0.0582123, 0.0593260, 0.0604172,
            0.0615151, 0.0627740, 0.0639047, 0.0647193, 0.0650727 }, {
            0.0472026, 0.0485713, 0.0495830, 0.0507434, 0.0519943, 0.0536504,
            0.0551176, 0.0566950, 0.0581553 }, { 0.0406631, 0.0419177,
            0.0429259, 0.0440614, 0.0451823, 0.0467889, 0.0483670, 0.0501810,
            0.0522433 }, { 0.0366926, 0.0377500, 0.0384530, 0.0393431,
            0.0404503, 0.0419569, 0.0434430, 0.0451987, 0.0473454 }, {
            0.0343793, 0.0352937, 0.0358116, 0.0364891, 0.0374246, 0.0385732,
            0.0398924, 0.0414707, 0.0435983 }, { 0.0331561, 0.0337733,
            0.0341567, 0.0346916, 0.0354239, 0.0364011, 0.0374280, 0.0387133,
            0.0405543 } };

    //     Set up interpolation:

    //      Locate subcells in array to use:

    s1 = NS - 2;
    s2 = NS - 1;
    for (i = 0; i < NS - 1; i++) {
        if (csz < cszt[i + 1]) {
            s1 = i;
            s2 = i + 1;
            break;
        }
    }
    t1 = NT - 2;
    t2 = NT - 1;
    for (i = 0; i < NT - 1; i++) {
        if (tau < taut[i + 1]) {
            t1 = i;
            t2 = i + 1;
            break;
        }
    }

    stot = cszt[s2] - cszt[s1];
    sdist = csz - cszt[s1];
    ttot = taut[t2] - taut[t1];
    tdist = tau - taut[t1];

    slope = (ast[s2][t1] - ast[s1][t1]) / stot;
    alb1[0] = ast[s1][t1] + (slope * sdist);
    slope = (ast[s2][t2] - ast[s1][t2]) / stot;
    alb1[1] = ast[s1][t2] + (slope * sdist);

    slope = (alb1[1] - alb1[0]) / ttot;
    as = alb1[0] + (slope * tdist);

    return as;
}

/********************************************************************************
 *      Given input cos(SZ), interpolate the appropriate surface albedo			*
 *       from the table. This is for the case where tau > 1.0					*
 *       Input is assumed to be in the valid range (the associated solar			*
 *       zenith angle should be in range 0.0 - 87.134 degrees)					*
 ********************************************************************************/

float interp_as_tauhigh(float csz) {

    int s1, s2, i;

    float stot, sdist, slope, as;

    float cszt[NS] = { 0.05, 0.09, 0.15, 0.21, 0.27, 0.33, 0.39, 0.45, 0.52,
            0.60, 0.68, 0.76, 0.84, 0.92, 1.00 };

    float ast[NS] = { 0.0623166, 0.0630070, 0.0640739, 0.0650397, 0.0657760,
            0.0661690, 0.0659415, 0.0651271, 0.0635101, 0.0610652, 0.0583470,
            0.0556146, 0.0530646, 0.0509498, 0.0490149 };

    //     Set up interpolation:

    //      Locate subcells in array to use:

    s1 = NS - 2;
    s2 = NS - 1;
    for (i = 0; i < NS; i++) {
        if (csz < cszt[i + 1]) {
            s1 = i;
            s2 = i + 1;
            break;
        }
    }

    stot = cszt[s2] - cszt[s1];
    sdist = csz - cszt[s1];

    slope = (ast[s2] - ast[s1]) / stot;
    as = ast[s1] + (slope * sdist);

    return as;
}

