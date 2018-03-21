#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "target_struc.h"
#include "l12_parms.h"

void free_target_n( tgstr *tgrec_n )
{
    free((void *) tgrec_n[0].data);
}


/* --------------------------------------------------------- */
/* alloc_target() - allocates 1 target record to hold data   */
/*              for a single scan of "npix" pixels.          */
/* --------------------------------------------------------- */
int alloc_target_n( int32_t nline, int32_t npix, int32_t nbands, tgstr *tgrec_n )
{
     char *p;
     int32_t i;
     int32_t len = 3*sizeof(int32_t)               /* time   */
              + sizeof(float)*npix                 /* solz   */
              + 2*nbands*sizeof(float)*npix;       /* lw,nLw */

     int32_t total_len = len * nline;

     if ((p = (char *) malloc(total_len)) == NULL) {
         fprintf(stderr,
                 "-E- %s line %d: Memory allocation failure.\n",
                 __FILE__,__LINE__);                                 
         return(0);
     }
     memset(p,'\0',total_len);

     for (i=0; i<nline; i++) {

         tgrec_n[i].length = len;
         tgrec_n[i].npix   = npix;
         tgrec_n[i].data   = p;

         tgrec_n[i].year = (int32_t   *) p;  p += sizeof(int32_t);
         tgrec_n[i].day  = (int32_t   *) p;  p += sizeof(int32_t); 
         tgrec_n[i].msec = (int32_t   *) p;  p += sizeof(int32_t);
         tgrec_n[i].solz = (float  *) p;  p += sizeof(float)*npix;
         tgrec_n[i].Lw   = (float  *) p;  p += sizeof(float)*npix*nbands;
         tgrec_n[i].nLw  = (float  *) p;  p += sizeof(float)*npix*nbands;

        if((len -  (int32_t)(p - tgrec_n[i].data)) < 0) {
            printf("%s Line %d: bad allocation on Target(n-lines) record\n", __FILE__, __LINE__);
            exit(1);
        }

     }

     printf("Allocated %d bytes in Target(n-lines) record. \n", total_len);

     return(total_len);
}
