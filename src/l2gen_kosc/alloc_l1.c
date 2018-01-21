
#include <stdio.h>
#include <stdlib.h>
#include "l1_struc.h"
#include "l12_parms.h"
#include "l12_proto.h"

void free_l1( l1str *l1rec )
{
     free((void *) l1rec->data);
}


/* --------------------------------------------------------- */
/* alloc_l1() - allocates 1 level-1b record to hold data for */
/*              a single scan of "npix" pixels.              */
/* --------------------------------------------------------- */
int32_t alloc_l1( int32_t npix, int32_t nbands, int32_t nbandsir, int32_t n_inprods, l1str *l1rec )
{
     char  *p;
     int32_t len;
     int32_t i;

     l1rec->npix = npix;

     /*                                                      */
     /* allocate data block as contiguous bytes              */
     /*                                                      */
     len = sizeof(int32_t)*nbands
         + 8*sizeof(float)*nbands
         + 3*sizeof(int32_t)
         + 7*sizeof(float)*npix
         + 2*sizeof(float)*npix*nbands
         + 2*sizeof(float)*npix*nbandsir
         + 1*sizeof(float)*npix               /* prtemp */
         + 1*sizeof(float)*npix
         + 2*sizeof(int32_t)*npix
         +22*sizeof(float)*npix
         +20*sizeof(float)*npix*nbands
         + 1*sizeof(int32_t)*npix
         + 1*sizeof(short)*npix
         + 1*sizeof(short)*npix
	 + 1*sizeof(short)*npix               /* ssttype */
         + 1*sizeof(short)*npix
         +17*sizeof(char )*npix
         +n_inprods*sizeof(float)*npix
         + 1*sizeof(int32)*npix
         + 1*sizeof(int32_t)*npix          /* pixdet */
         + 1*sizeof(float)*npix*nbands; /* smile_delta */

     /* Force to 4-byte increments for good measure */
     len = (len/4+1)*4;

     if ((p = (char *) malloc(len)) == NULL) {
         fprintf(stderr,
             "-E- %s line %d: Memory allocation failure.\n",
             __FILE__,__LINE__);                                 
         return(0);
     }

     /* Note: positional allocation is in order of datatype size, to 
        ensure that all 4-byte words start on 4-byte boundaries. Some 
        machines seem to have trouble if this is not done. */

     l1rec->length = len;     
     l1rec->data   = p;

     /* This first block must not be changed, or flat binary output L1 won't work */

     l1rec->year   = (int32_t   *) p;  p += sizeof(int32_t);
     l1rec->day    = (int32_t   *) p;  p += sizeof(int32_t); 
     l1rec->msec   = (int32_t   *) p;  p += sizeof(int32_t);

     l1rec->lon    = (float  *) p;  p += sizeof(float)*npix;
     l1rec->lat    = (float  *) p;  p += sizeof(float)*npix;
     l1rec->solz   = (float  *) p;  p += sizeof(float)*npix;
     l1rec->sola   = (float  *) p;  p += sizeof(float)*npix;
     l1rec->senz   = (float  *) p;  p += sizeof(float)*npix;
     l1rec->sena   = (float  *) p;  p += sizeof(float)*npix;

     l1rec->Lt     = (float  *) p;  p += sizeof(float)*npix*nbands;
     l1rec->Lt_unc = (float  *) p;  p += sizeof(float)*npix*nbands;

     l1rec->Ltir   = (float  *) p;  p += sizeof(float)*npix*nbandsir;
     l1rec->Bt     = (float  *) p;  p += sizeof(float)*npix*nbandsir;

     l1rec->prtemp = (float  *) p;  p += sizeof(float)*npix;

     l1rec->delphi = (float  *) p;  p += sizeof(float)*npix;
     l1rec->csolz  = (float  *) p;  p += sizeof(float)*npix;
     l1rec->csenz  = (float  *) p;  p += sizeof(float)*npix;
     l1rec->scattang = (float  *) p;  p += sizeof(float)*npix;

     l1rec->pixnum = (int32_t   *) p;  p += sizeof(int32_t)*npix;
     l1rec->nobs   = (int32_t   *) p;  p += sizeof(int32_t)*npix;

     l1rec->alpha  = (float  *) p;  p += sizeof(float)*npix;
     l1rec->ws     = (float  *) p;  p += sizeof(float)*npix; 
     l1rec->wd     = (float  *) p;  p += sizeof(float)*npix; 
     l1rec->mw     = (float  *) p;  p += sizeof(float)*npix; 
     l1rec->zw     = (float  *) p;  p += sizeof(float)*npix; 
     l1rec->pr     = (float  *) p;  p += sizeof(float)*npix; 
     l1rec->oz     = (float  *) p;  p += sizeof(float)*npix; 
     l1rec->wv     = (float  *) p;  p += sizeof(float)*npix; 
     l1rec->rh     = (float  *) p;  p += sizeof(float)*npix; 
     l1rec->no2_tropo = (float  *) p;  p += sizeof(float)*npix; 
     l1rec->no2_strat = (float  *) p;  p += sizeof(float)*npix; 
     l1rec->no2_frac = (float  *) p;  p += sizeof(float)*npix; 
     l1rec->height = (float  *) p;  p += sizeof(float)*npix; 
     l1rec->glint_coef   = (float  *) p;  p += sizeof(float)*npix;
     l1rec->cloud_albedo = (float  *) p;  p += sizeof(float)*npix;
     l1rec->aerindex     = (float  *) p;  p += sizeof(float)*npix;
     l1rec->sstref       = (float  *) p;  p += sizeof(float)*npix;
     l1rec->sssref       = (float  *) p;  p += sizeof(float)*npix;
     l1rec->rho_cirrus   = (float  *) p;  p += sizeof(float)*npix;

     l1rec->t_h2o    = (float  *) p;  p += sizeof(float)*npix*nbands; 
     l1rec->t_o2     = (float  *) p;  p += sizeof(float)*npix*nbands; 
     l1rec->tg_sol   = (float  *) p;  p += sizeof(float)*npix*nbands; 
     l1rec->tg_sen   = (float  *) p;  p += sizeof(float)*npix*nbands; 
     l1rec->t_sol    = (float  *) p;  p += sizeof(float)*npix*nbands; 
     l1rec->t_sen    = (float  *) p;  p += sizeof(float)*npix*nbands; 
     l1rec->rhof     = (float  *) p;  p += sizeof(float)*npix*nbands; 
     l1rec->tLf      = (float  *) p;  p += sizeof(float)*npix*nbands; 
     l1rec->Lr       = (float  *) p;  p += sizeof(float)*npix*nbands; 
     l1rec->L_q      = (float  *) p;  p += sizeof(float)*npix*nbands; 
     l1rec->L_u      = (float  *) p;  p += sizeof(float)*npix*nbands; 
     l1rec->polcor   = (float  *) p;  p += sizeof(float)*npix*nbands; 
     l1rec->dpol     = (float  *) p;  p += sizeof(float)*npix*nbands; 
     l1rec->TLg      = (float  *) p;  p += sizeof(float)*npix*nbands;
     l1rec->sw_n     = (float  *) p;  p += sizeof(float)*npix*nbands;
     l1rec->sw_a     = (float  *) p;  p += sizeof(float)*npix*nbands;
     l1rec->sw_bb    = (float  *) p;  p += sizeof(float)*npix*nbands;
     l1rec->sw_a_avg  = (float  *) p;  p += sizeof(float)*npix*nbands;
     l1rec->sw_bb_avg = (float  *) p;  p += sizeof(float)*npix*nbands;
     l1rec->rhos     = (float  *) p;  p += sizeof(float)*npix*nbands;

     l1rec->flags    = (int32_t   *) p;  p += sizeof(int32_t )*npix;

     l1rec->elev     = (float  *) p;  p += sizeof(float)*npix;
     l1rec->ancqc    = (short  *) p;  p += sizeof(short)*npix; 

     l1rec->ssttype  = (short  *) p;  p += sizeof(short)*npix;

     l1rec->mask     = (char   *) p;  p += sizeof(char )*npix; 
     l1rec->hilt     = (char   *) p;  p += sizeof(char )*npix; 
     l1rec->cloud    = (char   *) p;  p += sizeof(char )*npix; 
     l1rec->glint    = (char   *) p;  p += sizeof(char )*npix; 
     l1rec->land     = (char   *) p;  p += sizeof(char )*npix; 
     l1rec->swater   = (char   *) p;  p += sizeof(char )*npix;
     l1rec->ice      = (char   *) p;  p += sizeof(char )*npix;
     l1rec->solzmax  = (char   *) p;  p += sizeof(char )*npix;
     l1rec->senzmax  = (char   *) p;  p += sizeof(char )*npix;
     l1rec->stlight  = (char   *) p;  p += sizeof(char )*npix;
     l1rec->absaer   = (char   *) p;  p += sizeof(char )*npix;
     l1rec->navfail  = (char   *) p;  p += sizeof(char )*npix;
     l1rec->navwarn  = (char   *) p;  p += sizeof(char )*npix;
     l1rec->darkpix  = (char   *) p;  p += sizeof(char )*npix;
     l1rec->filter   = (char   *) p;  p += sizeof(char )*npix;
     l1rec->cirrus   = (char   *) p;  p += sizeof(char )*npix;
     l1rec->slot     = (unsigned char*) p;  p += sizeof(char )*npix;

     l1rec->pixdet   = (int32_t   *) p;  p += sizeof(int32_t )*npix;
     l1rec->radcor   = (float  *) p;  p += sizeof(float)*npix*nbands;

     l1rec->n_inprods    = n_inprods; 
     l1rec->in_prods     = (float **) malloc(sizeof(float*)*n_inprods);

     for ( i = 0; i < n_inprods; i++ ) {
         l1rec->in_prods[i]  = (float  *) p;  p += sizeof(float)*npix;
     }

     l1rec->in_flags  = (int32_t  *) p;  p += sizeof(int32)*npix;

     l1rec->iwave  = (int32_t   *) p;  p += sizeof(int32_t)*nbands;
     l1rec->fwave  = (float    *) p;   p += sizeof(float)*nbands;
     l1rec->Fo     = (float    *) p;   p += sizeof(float)*nbands;
     l1rec->Fobar  = (float    *) p;   p += sizeof(float)*nbands;
     l1rec->Fonom  = (float    *) p;   p += sizeof(float)*nbands;
     l1rec->Tau_r  = (float    *) p;   p += sizeof(float)*nbands;
     l1rec->k_oz   = (float    *) p;   p += sizeof(float)*nbands;
     l1rec->aw     = (float    *) p;   p += sizeof(float)*nbands;
     l1rec->bbw    = (float    *) p;   p += sizeof(float)*nbands;

     if ((len - (int32)(p - l1rec->data)) < 0) {
         printf("%s Line %d: bad allocation on L1 record\n",__FILE__,__LINE__);
         exit(1);
     }

     if(want_verbose)
        printf("Allocated %d bytes in L1 record.\n",(int32)(p - l1rec->data));


     return(1);
}



