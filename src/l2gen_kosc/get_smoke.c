/*---------------------------------------------------------------------*/
/* get_smoke.c - Eric Vermote's smoke index                            */
/*                                                                     */
/* Inputs:                                                             */
/*     l2rec - level-2 structure containing one complete scan after    */
/*             atmospheric correction.                                 */
/* Outputs:                                                            */
/*     smoke - smoke index, 1 value per pixel.                         */
/*                                                                     */
/* Written by: Bryan Franz, SAIC-GSC, February 2000                    */
/*                                                                     */
/*---------------------------------------------------------------------*/

#include <stdlib.h>
#include <math.h>
#include "l12_proto.h"

static float fill =  0.0;

void get_smoke(l2str *l2rec, float smoke[])
{
   int32_t   ip, ipb;

   float  rho1; 
   float  rho2;
   float  rho4;
   float  rho5;
   float  rho6;

   for (ip=0; ip<l2rec->npix; ip++) {

       ipb = l2rec->nbands*ip;

       rho1 = l2rec->rhos[ipb+0];
       rho2 = l2rec->rhos[ipb+1];
       rho4 = l2rec->rhos[ipb+3];
       rho5 = l2rec->rhos[ipb+4];
       rho6 = l2rec->rhos[ipb+5];

       if (rho1 < 0.0 || rho6 > 0.4)
           smoke[ip] = fill;

       else if (rho2 < rho4 || rho5 < rho6)
           smoke[ip] = 0.0;

       else
           smoke[ip] = rho1 * (rho2 - rho4) * (rho5 - rho6) * 10000.0;

   }

}

