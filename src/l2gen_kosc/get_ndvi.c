/*---------------------------------------------------------------------*/
/* get_ndvi.c - vegetation index classification for MSl12.             */
/*                                                                     */
/* Inputs:                                                             */
/*     l2rec - level-2 structure containing one complete scan after    */
/*             atmospheric correction.                                 */
/* Outputs:                                                            */
/*     ndvi  - vegetation index for land, 1 value per pixel.           */
/*                                                                     */
/* Written by: Bryan Franz, SAIC-GSC, February 2000                    */
/*                                                                     */
/*---------------------------------------------------------------------*/

#include <stdlib.h>
#include <math.h>
#include "l12_proto.h"

static float undef  = -2.0;

//static float minval = -2.0;
//static float maxval =  2.0;

// as per C. Tucker (11/2014), no cut-off at -2
static float minval = -1000.0;
static float maxval =  1000.0;

void get_ndvi(l2str *l2rec, float ndvi[])
{
   int32_t   ip, ipb;
   float  red; 
   float  nir;
   int    ib670 = windex(670,l2rec->fwave,l2rec->nbands);
   int    ib865 = windex(865,l2rec->fwave,l2rec->nbands);

   if (ib670 < 0 || ib865 < 0) {
       printf("NDVI requires bands near 670 and 865nm\n");
       exit(1);
   }

   for (ip=0; ip<l2rec->npix; ip++) {

       ndvi[ip] = undef;

       ipb = l2rec->nbands*ip;

       red = l2rec->rhos[ipb+ib670];
       nir = l2rec->rhos[ipb+ib865];

       if (red > 0.0 && nir > 0.0) {

           ndvi[ip] = (nir - red)
                    / (nir + red);

           ndvi[ip] = MIN(MAX(ndvi[ip],minval),maxval);
       }
   }
}


void get_evi(l2str *l2rec, float evi[])
{
   static float L=1.0, c1=6.0, c2=7.5;

   int32_t   ip, ipb;
   float  blu;
   float  red;
   float  nir;
   double val;
   int    ib412 = bindex_get(412);
   int    ib670 = bindex_get(670);
   int    ib865 = bindex_get(865);

   if (ib412 < 0 || ib670 < 0 || ib865 < 0) {
       printf("EVI requires bands near 412, 670, and 865nm\n");
       exit(1);
   }

   for (ip=0; ip<l2rec->npix; ip++) {

       evi[ip] = undef;

       ipb = l2rec->nbands*ip;

       blu = l2rec->rhos[ipb+ib412];
       red = l2rec->rhos[ipb+ib670];
       nir = l2rec->rhos[ipb+ib865];

       if ( red <= 0.0 || nir <= 0.0 )
           continue;

       else {
           if ( blu > 0.0  &&  ( blu <= red  ||  red <= nir ) ) {

               /* Most cases - EVI formula */

               if ( (val = L + nir + c1 * red - c2 * blu) == 0 )
                 continue;
               else
                 evi[ip] = 2.5 * (nir - red) / val;

           } else {   

               /* Backup - SAVI formula */

               if ( (val = 0.5 + nir + red) == 0 )
                 continue;
               else
                 evi[ip] = 1.5 * (nir - red) / val;
           }
           evi[ip] = MIN(MAX(evi[ip],minval),maxval);
       }
   }
}

