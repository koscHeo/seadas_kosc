#include <stdio.h>
#include <stdlib.h>
#include "l2_struc.h"
#include "l12_parms.h"
#include "l12_proto.h"

/* --------------------------------------------------------- */
/* init_l2() - initialize a Level-2 record                   */
/* --------------------------------------------------------- */
void init_l2( l2str *l2rec , int32_t nbands)
{
     char  *p = l2rec->data;

     int32_t  npix  = l2rec->npix;
     int32_t  nband = nbands;

     int32_t  ip, ib, ipb;

     for (ip=0; ip<npix; ip++) {

         for (ib=0; ib<nband; ib++) {
	     ipb = ip*nband+ib;
             l2rec->taua [ipb] = BAD_FLT;
             l2rec->TLg  [ipb] = BAD_FLT;
             l2rec->La   [ipb] = BAD_FLT;
             l2rec->Lw   [ipb] = BAD_FLT;
             l2rec->nLw  [ipb] = BAD_FLT;
             l2rec->nLw_unc[ipb] = 0.0;
             l2rec->brdf [ipb] = BAD_FLT;
             l2rec->Rrs  [ipb] = BAD_FLT;
             l2rec->Rrs_unc[ipb] = 0.0;
             l2rec->a    [ipb] = BAD_FLT;
             l2rec->bb   [ipb] = BAD_FLT;
	 }

         l2rec->chl[ip] = BAD_FLT;      
         l2rec->eps[ip] = BAD_FLT;         

         l2rec->aerratio  [ip] = BAD_FLT;    
         l2rec->aermodmin [ip] = BAD_INT;
         l2rec->aermodmax [ip] = BAD_INT;
         l2rec->aerratio2 [ip] = BAD_FLT;    
         l2rec->aermodmin2[ip] = BAD_INT;
         l2rec->aermodmax2[ip] = BAD_INT;
         l2rec->num_iter  [ip] = BAD_INT;
     }
}
