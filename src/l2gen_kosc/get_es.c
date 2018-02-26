/* --------------------------------------------------------------- */
/* get_es.c - compute surface irradiance (Ed(0+))                  */
/*                                                                 */
/* Inputs:                                                         */
/*     l2rec - level-2 structure containing one complete scan      */
/*             after atmospheric correction.                       */
/*     band  - band number (0-7)                                   */
/*                                                                 */
/* Outputs:                                                        */
/*     Es for specified band                                       */
/*                                                                 */
/* Algorithm Provided By: M. Wang                                  */
/* Written By: B. A. Franz, SAIC GSC, SIMBIOS Project, 4 Aug 1999  */
/*                                                                 */
/* --------------------------------------------------------------- */

#include <stdlib.h>
#include <math.h>
#include "l12_proto.h"

/**
 * now returns in W/m2/nm
 */
void get_es(l2str *l2rec, int band, float Es[])
{
    static float radeg = 3.141592654/180.;
    int32_t  ip;
    int32_t  ipb;

    for (ip=0; ip<l2rec->npix; ip++) {
        ipb = ip*l2rec->nbands+band;
        if (l2rec->La[ipb] > 0.0) {
            Es[ip] = l2rec->Fo[band]
               * l2rec->tg_sol[ipb]
               * l2rec->t_sol   [ipb]
               * cos(l2rec->solz[ip] * radeg) 
               * 10.0;
	} else
            Es[ip] = 0.0;
    }

    return;
}

