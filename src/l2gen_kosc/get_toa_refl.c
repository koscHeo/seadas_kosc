/* --------------------------------------------------------------- */
/* get_toa_refl.c - compute top-of-atmosphere reflectance.         */
/*                                                                 */
/* Inputs:                                                         */
/*     l2rec - level-2 structure containing one complete scan      */
/*             after atmospheric correction.                       */
/*     band  = band number 0-7                                     */
/* Outputs:                                                        */
/*     rhot  - toa reflectance                                     */
/*                                                                 */
/* Written By: B. Franz, SAIC GSC, SIMBIOS Project, 11 April 2000  */
/*                                                                 */
/* --------------------------------------------------------------- */

#include <stdlib.h>
#include <math.h>
#include "l12_proto.h"

void get_toa_refl(l2str *l2rec, int band, float rhot[])
{
    static float pi    = 3.141592654;
    static float radeg = 57.29577951;

    float mu0;
    int32_t  ip, ipb;

    for (ip=0; ip<l2rec->npix; ip++) {
        ipb = ip*l2rec->nbands + band;
        mu0 = cos(l2rec->solz[ip]/radeg);
        rhot[ip] = pi * l2rec->Lt[ipb]/l2rec->Fo[band]/mu0;
    }

    return;
}
