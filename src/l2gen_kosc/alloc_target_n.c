#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "target_struc_n.h"
#include "l12_parms.h"

/* --------------------------------------------------------- */
/* alloc_target() - allocates 1 target record to hold data   */
/*              for a single scan of "npix" pixels.          */
/* --------------------------------------------------------- */
int alloc_target_n( int32_t nline, int32_t npix, int32_t nbands, tgstr_n *rec )
{
     char *p;
     int32_t len = 3*sizeof(int32_t)               /* time   */
              + sizeof(float)*nline*npix           /* solz   */
              + 2*nbands*sizeof(float)*nline*npix; /* lw,nLw */

     if ((p = (char *) malloc(len)) == NULL) {
         fprintf(stderr,
             "-E- %s line %d: Memory allocation failure.\n",
             __FILE__,__LINE__);                                 
         return(0);
     }
     memset(p,'\0',len);

     rec->length = len;
     rec->npix   = npix;
     rec->data   = p;

     rec->year = (int32_t   *) p;  p += sizeof(int32_t);
     rec->day  = (int32_t   *) p;  p += sizeof(int32_t); 
     rec->msec = (int32_t   *) p;  p += sizeof(int32_t);
     rec->solz = (float  *) p;  p += sizeof(float)*nline*npix;
     rec->Lw   = (float  *) p;  p += sizeof(float)*nline*npix*nbands;
     rec->nLw  = (float  *) p;  p += sizeof(float)*nline*npix*nbands;

     return(len);
}
