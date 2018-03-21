#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "aer_struc.h"
#include "l12_parms.h"

void free_aer_n( aestr *aerec_n )
{
    free((void *) aerec_n[0].data);
}


/* --------------------------------------------------------- */
/* alloc_aer() - allocates 1 aerosol record to hold data     */
/*              for a single scan of "npix" pixels.          */
/* --------------------------------------------------------- */
int alloc_aer_n( int32_t nline, int32_t npix, int32_t nbands, aestr *aerec_n )
{
     char *p;
     int32_t i;
     int32_t len = 2*sizeof(int32_t)*npix       /* model #      */
              + 1*sizeof(float)*npix            /* model ratio  */
              + 1*nbands*sizeof(float)*npix;    /* taua         */
    
    int32_t total_len = len * nline;


     if ((p = (char *) malloc(total_len)) == NULL) {
         fprintf(stderr,
             "-E- %s line %d: Memory allocation failure.\n",
             __FILE__,__LINE__);                                 
         return(0);
     }
     memset(p,'\0',total_len);

     for (i=0; i<nline; i++) {

         aerec_n[i].length = len;
         aerec_n[i].npix   = npix;
         aerec_n[i].data   = p;

         aerec_n[i].mod_min = (int32_t   *) p;  p += sizeof(int32_t)*npix;
         aerec_n[i].mod_max = (int32_t   *) p;  p += sizeof(int32_t)*npix; 
         aerec_n[i].mod_rat = (float  *) p;  p += sizeof(float)*npix;
         aerec_n[i].taua    = (float  *) p;  p += sizeof(float)*npix*nbands;
        
        if((len -  (int32_t)(p - aerec_n[i].data)) < 0) {
            printf("%s Line %d: bad allocation on Aerosol(n-lines) record\n", __FILE__, __LINE__);
            exit(1);
        }

     }
     printf("Allocated %d bytes in Aerosol(n-lines) record. \n", total_len);

     return(total_len);
}
