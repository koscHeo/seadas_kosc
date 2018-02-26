/* --------------------------------------------------------------- */
/* get_ice_frac.c - dump the ice fraction ancillary data           */
/*                                                                 */
/* Inputs:                                                         */
/*     l2rec - level-2 structure containing one complete scan      */
/*             after atmospheric correction.                       */
/*                                                                 */
/* Outputs:                                                        */
/*     ice - array where ice data is written                       */
/*                                                                 */
/* Written By: Don Shea, OBPG,  23 Apr 2009                        */
/*                                                                 */
/* --------------------------------------------------------------- */

#include "l12_proto.h"

void get_ice_frac(l2str *l2rec, float ice[])
{
    int32_t  ip;

    for (ip=0; ip<l2rec->npix; ip++) {
        if(l2rec->flags[ip] & LAND)
            ice[ip] = 0.0;
        else
            ice[ip] = ice_fraction(l2rec->lon[ip], l2rec->lat[ip]);
    }

}

