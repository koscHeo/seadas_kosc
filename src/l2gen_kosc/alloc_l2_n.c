#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "l2_struc.h"
#include "l12_parms.h"

void free_l2_n( l2str *l2rec_n )
{
    free((void *) l2rec_n[0].data);
}


/* --------------------------------------------------------- */
/* alloc_l2() - allocates 1 level-2 record to hold data for  */
/*              a single scan of "npix" pixels.              */
/* --------------------------------------------------------- */
int alloc_l2_n( int32_t nline,  int32_t npix, int32_t nbands, l2str *l2rec_n )
{
     char *p;
     int32_t i;
     int32_t len =  3*sizeof(int32_t ) 
              +  4*sizeof(int32_t )*npix
              + 10*sizeof(float)*npix
              + 13*sizeof(float)*npix*nbands
              +  1*sizeof(int32_t)*npix;

     if (len % 4 != 0)
         len = len / 4 * 4 + 4;
    
    int32_t total_len = len * nline;

     if ((p = (char *) malloc(total_len)) == NULL) {
         fprintf(stderr,
             "-E- %s line %d: Memory allocation failure.\n",
             __FILE__,__LINE__);                                 
         return(0);
     }

     for (i=0; i<nline; i++) {

         l2rec_n[i].length = len;
         l2rec_n[i].npix   = npix;
         l2rec_n[i].data   = p;

         l2rec_n[i].year      = (int32_t   *) p;  p += sizeof(int32_t);
         l2rec_n[i].day       = (int32_t   *) p;  p += sizeof(int32_t); 
         l2rec_n[i].msec      = (int32_t   *) p;  p += sizeof(int32_t);
         l2rec_n[i].lon       = (float  *) p;  p += sizeof(float)*npix;
         l2rec_n[i].lat       = (float  *) p;  p += sizeof(float)*npix;
         l2rec_n[i].solz      = (float  *) p;  p += sizeof(float)*npix;
         l2rec_n[i].sola      = (float  *) p;  p += sizeof(float)*npix;
         l2rec_n[i].senz      = (float  *) p;  p += sizeof(float)*npix;
         l2rec_n[i].sena      = (float  *) p;  p += sizeof(float)*npix;
         l2rec_n[i].Lt        = (float  *) p;  p += sizeof(float)*npix*nbands;
         l2rec_n[i].Lt_unc    = (float  *) p;  p += sizeof(float)*npix*nbands;
         l2rec_n[i].aermodmin = (int32_t   *) p;  p += sizeof(int32_t)*npix;
         l2rec_n[i].aermodmax = (int32_t   *) p;  p += sizeof(int32_t)*npix;
         l2rec_n[i].aerratio  = (float  *) p;  p += sizeof(float)*npix;
         l2rec_n[i].aermodmin2= (int32_t   *) p;  p += sizeof(int32_t)*npix;
         l2rec_n[i].aermodmax2= (int32_t   *) p;  p += sizeof(int32_t)*npix;
         l2rec_n[i].aerratio2 = (float  *) p;  p += sizeof(float)*npix;
         l2rec_n[i].eps       = (float  *) p;  p += sizeof(float)*npix;
         l2rec_n[i].taua      = (float  *) p;  p += sizeof(float)*npix*nbands;
         l2rec_n[i].TLg       = (float  *) p;  p += sizeof(float)*npix*nbands;
         l2rec_n[i].La        = (float  *) p;  p += sizeof(float)*npix*nbands;
         l2rec_n[i].Lw        = (float  *) p;  p += sizeof(float)*npix*nbands;
         l2rec_n[i].nLw       = (float  *) p;  p += sizeof(float)*npix*nbands;
         l2rec_n[i].nLw_unc   = (float  *) p;  p += sizeof(float)*npix*nbands;
         l2rec_n[i].brdf      = (float  *) p;  p += sizeof(float)*npix*nbands;
         l2rec_n[i].Rrs       = (float  *) p;  p += sizeof(float)*npix*nbands;
         l2rec_n[i].Rrs_unc   = (float  *) p;  p += sizeof(float)*npix*nbands;
         l2rec_n[i].a         = (float  *) p;  p += sizeof(float)*npix*nbands;
         l2rec_n[i].bb        = (float  *) p;  p += sizeof(float)*npix*nbands;
         l2rec_n[i].chl       = (float  *) p;  p += sizeof(float)*npix;
         l2rec_n[i].num_iter  = (int32_t   *) p;  p += sizeof(int32_t )*npix;

         if ((len - (int32_t)(p - l2rec_n[i].data)) < 0) {
             printf("%s Line %d: bad allocation on L2(n-lines) record\n",__FILE__,__LINE__);
             exit(1);
         }

     }

     printf("Allocated %d bytes in L2(n-lines) record.\n", total_len);

     return(total_len);
}
