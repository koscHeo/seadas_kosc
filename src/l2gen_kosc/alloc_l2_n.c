#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "l2_struc_n.h"
#include "l12_parms.h"

void free_l2_n( l2str_n *l2rec )
{
     free((void *) l2rec->data);
}


/* --------------------------------------------------------- */
/* alloc_l2() - allocates 1 level-2 record to hold data for  */
/*              a single scan of "npix" pixels.              */
/* --------------------------------------------------------- */
int alloc_l2_n( int32_t nline,  int32_t npix, int32_t nbands, l2str_n *l2rec )
{
     char *p;
     int32_t len =  3*sizeof(int32_t ) 
              +  4*sizeof(int32_t )*nline*npix
              + 10*sizeof(float)*nline*npix
              + 13*sizeof(float)*nline*npix*nbands
              +  1*sizeof(int32_t)*nline*npix;

     if (len % 4 != 0)
         len = len / 4 * 4 + 4;

     if ((p = (char *) malloc(len)) == NULL) {
         fprintf(stderr,
             "-E- %s line %d: Memory allocation failure.\n",
             __FILE__,__LINE__);                                 
         return(0);
     }

     l2rec->length = len;
     l2rec->npix   = npix;
     l2rec->data   = p;

     l2rec->year      = (int32_t   *) p;  p += sizeof(int32_t);
     l2rec->day       = (int32_t   *) p;  p += sizeof(int32_t); 
     l2rec->msec      = (int32_t   *) p;  p += sizeof(int32_t);
     l2rec->lon       = (float  *) p;  p += sizeof(float)*nline*npix;
     l2rec->lat       = (float  *) p;  p += sizeof(float)*nline*npix;
     l2rec->solz      = (float  *) p;  p += sizeof(float)*nline*npix;
     l2rec->sola      = (float  *) p;  p += sizeof(float)*nline*npix;
     l2rec->senz      = (float  *) p;  p += sizeof(float)*nline*npix;
     l2rec->sena      = (float  *) p;  p += sizeof(float)*nline*npix;
     l2rec->Lt        = (float  *) p;  p += sizeof(float)*nline*npix*nbands;
     l2rec->Lt_unc    = (float  *) p;  p += sizeof(float)*nline*npix*nbands;
     l2rec->aermodmin = (int32_t   *) p;  p += sizeof(int32_t)*nline*npix;
     l2rec->aermodmax = (int32_t   *) p;  p += sizeof(int32_t)*nline*npix;
     l2rec->aerratio  = (float  *) p;  p += sizeof(float)*nline*npix;
     l2rec->aermodmin2= (int32_t   *) p;  p += sizeof(int32_t)*nline*npix;
     l2rec->aermodmax2= (int32_t   *) p;  p += sizeof(int32_t)*nline*npix;
     l2rec->aerratio2 = (float  *) p;  p += sizeof(float)*nline*npix;
     l2rec->eps       = (float  *) p;  p += sizeof(float)*nline*npix;
     l2rec->taua      = (float  *) p;  p += sizeof(float)*nline*npix*nbands;
     l2rec->TLg       = (float  *) p;  p += sizeof(float)*nline*npix*nbands;
     l2rec->La        = (float  *) p;  p += sizeof(float)*nline*npix*nbands;
     l2rec->Lw        = (float  *) p;  p += sizeof(float)*nline*npix*nbands;
     l2rec->nLw       = (float  *) p;  p += sizeof(float)*nline*npix*nbands;
     l2rec->nLw_unc   = (float  *) p;  p += sizeof(float)*nline*npix*nbands;
     l2rec->brdf      = (float  *) p;  p += sizeof(float)*nline*npix*nbands;
     l2rec->Rrs       = (float  *) p;  p += sizeof(float)*nline*npix*nbands;
     l2rec->Rrs_unc   = (float  *) p;  p += sizeof(float)*nline*npix*nbands;
     l2rec->a         = (float  *) p;  p += sizeof(float)*nline*npix*nbands;
     l2rec->bb        = (float  *) p;  p += sizeof(float)*nline*npix*nbands;
     l2rec->chl       = (float  *) p;  p += sizeof(float)*nline*npix;
     l2rec->num_iter  = (int32_t   *) p;  p += sizeof(int32_t )*nline*npix;

     if ((len - (int32)(p - l2rec->data)) < 0) {
         printf("%s Line %d: bad allocation on L2(n-lines) record\n",__FILE__,__LINE__);
         exit(1);
     }

     printf("Allocated %d bytes in L2(n-lines) record.\n",(int32)(p - l2rec->data));

     return(len);
}
