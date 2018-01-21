/* -------------------------------------------------------------------- */
/* l1subpix() - sub-samples and/or crops a level-1 record in place      */
/*                                                                      */
/* Returns 0 on succss, 1 on error.                                     */
/*                                                                      */
/* Notes: Only record fields which are filled by the L1 read routines   */
/*        must be updated here.  All other fields will be filled later  */
/*        using the sub-sampled geometry and radiances.                 */
/*                                                                      */
/* Written By: B. A. Franz                                              */
/* -------------------------------------------------------------------- */

#include <string.h>
#include "filehandle.h"
#include "l1_struc.h"

int l1subpix(filehandle *l1file, l1str *l1rec)
{
    int32_t sp = l1file->spix;
    int32_t ep = l1file->epix;
    int32_t dp = l1file->dpix;

    if ( (sp == 0) && (ep == l1file->npix-1) && dp == 1)
        return (0);

    if (sp > ep || dp < 1 || dp > (ep - sp + 1)) {
        fprintf(stderr,
        "-E- %s Line %d: subpixel specification error (sp=%d,ep=%d,dp=%d).\n",
        __FILE__,__LINE__,sp,ep,dp);
        return (1);
    }

    if (dp == 1) {

        int32_t length;
        int32_t j;

        l1rec->npix = ep - sp + 1;
	
	length  = l1rec->npix*sizeof(float);

        memmove(l1rec->pixnum,&l1rec->pixnum[sp],l1rec->npix*sizeof(int32_t));
        memmove(l1rec->slot,  &l1rec->slot  [sp],l1rec->npix*sizeof(unsigned char));
        memmove(l1rec->nobs  ,&l1rec->nobs  [sp],l1rec->npix*sizeof(int32_t));

        memmove(l1rec->lon, &l1rec->lon [sp],length);
        memmove(l1rec->lat, &l1rec->lat [sp],length);
        memmove(l1rec->solz,&l1rec->solz[sp],length);
        memmove(l1rec->sola,&l1rec->sola[sp],length);
        memmove(l1rec->senz,&l1rec->senz[sp],length);
        memmove(l1rec->sena,&l1rec->sena[sp],length);

        memmove(l1rec->alpha,&l1rec->alpha[sp],length);

        memmove(l1rec->hilt   ,&l1rec->hilt   [sp],l1rec->npix*sizeof(char));
        memmove(l1rec->stlight,&l1rec->stlight[sp],l1rec->npix*sizeof(char));
        memmove(l1rec->navfail,&l1rec->navfail[sp],l1rec->npix*sizeof(char));
        memmove(l1rec->navwarn,&l1rec->navwarn[sp],l1rec->npix*sizeof(char));

        memmove(l1rec->Lt,     &l1rec->Lt    [sp*l1rec->nbands],length*l1rec->nbands);
        memmove(l1rec->Lt_unc, &l1rec->Lt_unc[sp*l1rec->nbands],length*l1rec->nbands);

        memmove(l1rec->sw_n,   &l1rec->sw_n   [sp*l1rec->nbands],length*l1rec->nbands);
        memmove(l1rec->sw_a,   &l1rec->sw_a   [sp*l1rec->nbands],length*l1rec->nbands);
        memmove(l1rec->sw_bb,  &l1rec->sw_bb  [sp*l1rec->nbands],length*l1rec->nbands);
        memmove(l1rec->sw_a_avg,   &l1rec->sw_a_avg   [sp*l1rec->nbands],length*l1rec->nbands);
        memmove(l1rec->sw_bb_avg,  &l1rec->sw_bb_avg  [sp*l1rec->nbands],length*l1rec->nbands);

        memmove(l1rec->Ltir, &l1rec->Ltir [sp*NBANDSIR],length*NBANDSIR);
        memmove(l1rec->Bt,   &l1rec->Bt   [sp*NBANDSIR],length*NBANDSIR);

        memmove(l1rec->rho_cirrus,&l1rec->rho_cirrus[sp],length);


	/* MERIS */
        memmove(l1rec->pixdet,&l1rec->pixdet[sp],l1rec->npix*sizeof(int32_t));
        memmove(l1rec->radcor, &l1rec->radcor[sp*l1rec->nbands],length*l1rec->nbands);
        memmove(l1rec->in_flags,&l1rec->in_flags[sp],l1rec->npix*sizeof(int32));
        for ( j = 0; j < l1rec->n_inprods; j++ )
            memmove(l1rec->in_prods[j], &l1rec->in_prods[j][sp],length);

	/* AVHRR */
        memmove(l1rec->prtemp,&l1rec->prtemp[sp],length);

    } else {

        int32_t i,j;

        l1rec->npix = (ep - sp)/dp + 1;

        for (i=0; i<l1rec->npix; i++) {

            l1rec->pixnum[i] = l1rec->pixnum[i*dp+sp];
            l1rec->slot[i] = l1rec->slot[i*dp+sp];
            l1rec->nobs[i] = l1rec->nobs[i*dp+sp];

            l1rec->lon [i] = l1rec->lon [i*dp+sp];
            l1rec->lat [i] = l1rec->lat [i*dp+sp];
            l1rec->solz[i] = l1rec->solz[i*dp+sp];
            l1rec->sola[i] = l1rec->sola[i*dp+sp];
            l1rec->senz[i] = l1rec->senz[i*dp+sp];
            l1rec->sena[i] = l1rec->sena[i*dp+sp];

            l1rec->alpha[i] = l1rec->alpha[i*dp+sp];

            l1rec->hilt   [i] = l1rec->hilt   [i*dp+sp];
            l1rec->stlight[i] = l1rec->stlight[i*dp+sp];
            l1rec->navfail[i] = l1rec->navfail[i*dp+sp];
            l1rec->navwarn[i] = l1rec->navwarn[i*dp+sp];

            l1rec->rho_cirrus[i] = l1rec->rho_cirrus[i*dp+sp];

            memmove(&l1rec->Lt   [i*l1rec->nbands],
                   &l1rec->Lt   [(i*dp+sp)*l1rec->nbands],l1rec->nbands*sizeof(float));

            memmove(&l1rec->sw_n [i*l1rec->nbands],
                   &l1rec->sw_n [(i*dp+sp)*l1rec->nbands],l1rec->nbands*sizeof(float));
            memmove(&l1rec->sw_a [i*l1rec->nbands],
                   &l1rec->sw_a [(i*dp+sp)*l1rec->nbands],l1rec->nbands*sizeof(float));
            memmove(&l1rec->sw_bb[i*l1rec->nbands],
                   &l1rec->sw_bb[(i*dp+sp)*l1rec->nbands],l1rec->nbands*sizeof(float));
            memmove(&l1rec->sw_a_avg [i*l1rec->nbands],
                   &l1rec->sw_a_avg [(i*dp+sp)*l1rec->nbands],l1rec->nbands*sizeof(float));
            memmove(&l1rec->sw_bb_avg[i*l1rec->nbands],
                   &l1rec->sw_bb_avg[(i*dp+sp)*l1rec->nbands],l1rec->nbands*sizeof(float));

            memmove(&l1rec->Ltir [i*NBANDSIR],
                   &l1rec->Ltir [(i*dp+sp)*NBANDSIR],NBANDSIR*sizeof(float));

            memmove(&l1rec->Bt   [i*NBANDSIR],
                   &l1rec->Bt   [(i*dp+sp)*NBANDSIR],NBANDSIR*sizeof(float));

	    /* MERIS */
            l1rec->pixdet[i] = l1rec->pixdet[i*dp+sp];
            memmove(&l1rec->radcor   [i*l1rec->nbands],
                   &l1rec->radcor   [(i*dp+sp)*l1rec->nbands],l1rec->nbands*sizeof(float));

            l1rec->in_flags[i] = l1rec->in_flags[i*dp+sp];
            for ( j = 0; j < l1rec->n_inprods; j++ )
                l1rec->in_prods[j][i] = l1rec->in_prods[j][i*dp+sp];

	    /* AVHRR */
            l1rec->prtemp[i] = l1rec->prtemp[i*dp+sp];
	}
    }

    l1rec->spix = 0;
    l1rec->epix = l1rec->npix-1;
    l1rec->dpix = 1;

    return (0);
}

