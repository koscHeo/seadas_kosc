
#include <stdio.h>
#include <stdlib.h>
#include "l1_struc.h"
#include "l12_parms.h"
#include "l12_proto.h"

void free_l1_n( l1str *l1rec_n )
{
     free((void *) l1rec_n[0].data);
}


/* --------------------------------------------------------- */
/* alloc_l1() - allocates 1 level-1b record to hold data for */
/*              a single scan of "npix" pixels.              */
/* --------------------------------------------------------- */
int32_t alloc_l1_n( int32_t nline, int32_t npix, int32_t nbands, int32_t nbandsir, int32_t n_inprods, l1str *l1rec_n )
{
     char  *p;
     int32_t len;
     int32_t i, j;


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
         + 1*sizeof(int32_t)*npix
         + 1*sizeof(int32_t)*npix          /* pixdet */
         + 1*sizeof(float)*npix*nbands; /* smile_delta */

     /* Force to 4-byte increments for good measure */
     len = (len/4+1)*4;

     int32_t total_len = len * nline;

     if ((p = (char *) malloc(total_len)) == NULL) {
         fprintf(stderr,
             "-E- %s line %d: Memory allocation failure.\n",
             __FILE__,__LINE__);                                 
         return(0);
     }

     /* Note: positional allocation is in order of datatype size, to 
        ensure that all 4-byte words start on 4-byte boundaries. Some 
        machines seem to have trouble if this is not done. */
        

     for (i=0; i<nline; i++) {

         l1rec_n[i].npix   = npix;
         l1rec_n[i].length = len;     
         l1rec_n[i].data   = p;

         /* This first block must not be changed, or flat binary output L1 won't work */

         l1rec_n[i].year   = (int32_t   *) p;  p += sizeof(int32_t);
         l1rec_n[i].day    = (int32_t   *) p;  p += sizeof(int32_t); 
         l1rec_n[i].msec   = (int32_t   *) p;  p += sizeof(int32_t);

         l1rec_n[i].lon    = (float  *) p;  p += sizeof(float)*npix;
         l1rec_n[i].lat    = (float  *) p;  p += sizeof(float)*npix;
         l1rec_n[i].solz   = (float  *) p;  p += sizeof(float)*npix;
         l1rec_n[i].sola   = (float  *) p;  p += sizeof(float)*npix;
         l1rec_n[i].senz   = (float  *) p;  p += sizeof(float)*npix;
         l1rec_n[i].sena   = (float  *) p;  p += sizeof(float)*npix;

         l1rec_n[i].Lt     = (float  *) p;  p += sizeof(float)*npix*nbands;
         l1rec_n[i].Lt_unc = (float  *) p;  p += sizeof(float)*npix*nbands;

         l1rec_n[i].Ltir   = (float  *) p;  p += sizeof(float)*npix*nbandsir;
         l1rec_n[i].Bt     = (float  *) p;  p += sizeof(float)*npix*nbandsir;

         l1rec_n[i].prtemp = (float  *) p;  p += sizeof(float)*npix;

         l1rec_n[i].delphi = (float  *) p;  p += sizeof(float)*npix;
         l1rec_n[i].csolz  = (float  *) p;  p += sizeof(float)*npix;
         l1rec_n[i].csenz  = (float  *) p;  p += sizeof(float)*npix;
         l1rec_n[i].scattang = (float  *) p;  p += sizeof(float)*npix;

         l1rec_n[i].pixnum = (int32_t   *) p;  p += sizeof(int32_t)*npix;
         l1rec_n[i].nobs   = (int32_t   *) p;  p += sizeof(int32_t)*npix;

         l1rec_n[i].alpha  = (float  *) p;  p += sizeof(float)*npix;
         l1rec_n[i].ws     = (float  *) p;  p += sizeof(float)*npix; 
         l1rec_n[i].wd     = (float  *) p;  p += sizeof(float)*npix; 
         l1rec_n[i].mw     = (float  *) p;  p += sizeof(float)*npix; 
         l1rec_n[i].zw     = (float  *) p;  p += sizeof(float)*npix; 
         l1rec_n[i].pr     = (float  *) p;  p += sizeof(float)*npix; 
         l1rec_n[i].oz     = (float  *) p;  p += sizeof(float)*npix; 
         l1rec_n[i].wv     = (float  *) p;  p += sizeof(float)*npix; 
         l1rec_n[i].rh     = (float  *) p;  p += sizeof(float)*npix; 
         l1rec_n[i].no2_tropo = (float  *) p;  p += sizeof(float)*npix; 
         l1rec_n[i].no2_strat = (float  *) p;  p += sizeof(float)*npix; 
         l1rec_n[i].no2_frac = (float  *) p;  p += sizeof(float)*npix; 
         l1rec_n[i].height = (float  *) p;  p += sizeof(float)*npix; 
         l1rec_n[i].glint_coef   = (float  *) p;  p += sizeof(float)*npix;
         l1rec_n[i].cloud_albedo = (float  *) p;  p += sizeof(float)*npix;
         l1rec_n[i].aerindex     = (float  *) p;  p += sizeof(float)*npix;
         l1rec_n[i].sstref       = (float  *) p;  p += sizeof(float)*npix;
         l1rec_n[i].sssref       = (float  *) p;  p += sizeof(float)*npix;
         l1rec_n[i].rho_cirrus   = (float  *) p;  p += sizeof(float)*npix;

         l1rec_n[i].t_h2o    = (float  *) p;  p += sizeof(float)*npix*nbands; 
         l1rec_n[i].t_o2     = (float  *) p;  p += sizeof(float)*npix*nbands; 
         l1rec_n[i].tg_sol   = (float  *) p;  p += sizeof(float)*npix*nbands; 
         l1rec_n[i].tg_sen   = (float  *) p;  p += sizeof(float)*npix*nbands; 
         l1rec_n[i].t_sol    = (float  *) p;  p += sizeof(float)*npix*nbands; 
         l1rec_n[i].t_sen    = (float  *) p;  p += sizeof(float)*npix*nbands; 
         l1rec_n[i].rhof     = (float  *) p;  p += sizeof(float)*npix*nbands; 
         l1rec_n[i].tLf      = (float  *) p;  p += sizeof(float)*npix*nbands; 
         l1rec_n[i].Lr       = (float  *) p;  p += sizeof(float)*npix*nbands; 
         l1rec_n[i].L_q      = (float  *) p;  p += sizeof(float)*npix*nbands; 
         l1rec_n[i].L_u      = (float  *) p;  p += sizeof(float)*npix*nbands; 
         l1rec_n[i].polcor   = (float  *) p;  p += sizeof(float)*npix*nbands; 
         l1rec_n[i].dpol     = (float  *) p;  p += sizeof(float)*npix*nbands; 
         l1rec_n[i].TLg      = (float  *) p;  p += sizeof(float)*npix*nbands;
         l1rec_n[i].sw_n     = (float  *) p;  p += sizeof(float)*npix*nbands;
         l1rec_n[i].sw_a     = (float  *) p;  p += sizeof(float)*npix*nbands;
         l1rec_n[i].sw_bb    = (float  *) p;  p += sizeof(float)*npix*nbands;
         l1rec_n[i].sw_a_avg  = (float  *) p;  p += sizeof(float)*npix*nbands;
         l1rec_n[i].sw_bb_avg = (float  *) p;  p += sizeof(float)*npix*nbands;
         l1rec_n[i].rhos     = (float  *) p;  p += sizeof(float)*npix*nbands;

         l1rec_n[i].flags    = (int32_t   *) p;  p += sizeof(int32_t )*npix;

         l1rec_n[i].elev     = (float  *) p;  p += sizeof(float)*npix;
         l1rec_n[i].ancqc    = (short  *) p;  p += sizeof(short)*npix; 

         l1rec_n[i].ssttype  = (short  *) p;  p += sizeof(short)*npix;

         l1rec_n[i].mask     = (char   *) p;  p += sizeof(char )*npix; 
         l1rec_n[i].hilt     = (char   *) p;  p += sizeof(char )*npix; 
         l1rec_n[i].cloud    = (char   *) p;  p += sizeof(char )*npix; 
         l1rec_n[i].glint    = (char   *) p;  p += sizeof(char )*npix; 
         l1rec_n[i].land     = (char   *) p;  p += sizeof(char )*npix; 
         l1rec_n[i].swater   = (char   *) p;  p += sizeof(char )*npix;
         l1rec_n[i].ice      = (char   *) p;  p += sizeof(char )*npix;
         l1rec_n[i].solzmax  = (char   *) p;  p += sizeof(char )*npix;
         l1rec_n[i].senzmax  = (char   *) p;  p += sizeof(char )*npix;
         l1rec_n[i].stlight  = (char   *) p;  p += sizeof(char )*npix;
         l1rec_n[i].absaer   = (char   *) p;  p += sizeof(char )*npix;
         l1rec_n[i].navfail  = (char   *) p;  p += sizeof(char )*npix;
         l1rec_n[i].navwarn  = (char   *) p;  p += sizeof(char )*npix;
         l1rec_n[i].darkpix  = (char   *) p;  p += sizeof(char )*npix;
         l1rec_n[i].filter   = (char   *) p;  p += sizeof(char )*npix;
         l1rec_n[i].cirrus   = (char   *) p;  p += sizeof(char )*npix;
         l1rec_n[i].slot     = (unsigned char*) p;  p += sizeof(char )*npix;

         l1rec_n[i].pixdet   = (int32_t   *) p;  p += sizeof(int32_t )*npix;
         l1rec_n[i].radcor   = (float  *) p;  p += sizeof(float)*npix*nbands;

         l1rec_n[i].n_inprods    = n_inprods; 
         l1rec_n[i].in_prods     = (float **) malloc(sizeof(float*)*n_inprods);
         for ( j = 0; j < n_inprods; j++ ) {
             l1rec_n[i].in_prods[j]  = (float  *) p;  p += sizeof(float)*npix;
         }

         l1rec_n[i].in_flags  = (int32_t  *) p;  p += sizeof(int32_t)*npix;

         l1rec_n[i].iwave  = (int32_t   *) p;  p += sizeof(int32_t)*nbands;
         l1rec_n[i].fwave  = (float    *) p;   p += sizeof(float)*nbands;
         l1rec_n[i].Fo     = (float    *) p;   p += sizeof(float)*nbands;
         l1rec_n[i].Fobar  = (float    *) p;   p += sizeof(float)*nbands;
         l1rec_n[i].Fonom  = (float    *) p;   p += sizeof(float)*nbands;
         l1rec_n[i].Tau_r  = (float    *) p;   p += sizeof(float)*nbands;
         l1rec_n[i].k_oz   = (float    *) p;   p += sizeof(float)*nbands;
         l1rec_n[i].aw     = (float    *) p;   p += sizeof(float)*nbands;
         l1rec_n[i].bbw    = (float    *) p;   p += sizeof(float)*nbands;

         if ((len - (int32_t)(p - l1rec_n[i].data)) < 0) {
             printf("%s Line %d: bad allocation on L1(n-lines) record\n",__FILE__,__LINE__);
             exit(1);
         }

     }

     if(want_verbose)
         printf("Allocated %d bytes in L1(n-lines) record.\n",total_len);


     return(1);
}



