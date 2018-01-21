#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "l12_proto.h"

/* ---------------------------------------------------------- */
/* Transfers pointers and various pre-computed quantities     */
/* between the Level-1 abd Level-2 records.                   */ 
/* ---------------------------------------------------------- */

void cpl1l2(l1str *l1rec, l2str* l2rec) {

    int32_t npix = l1rec->npix;

    /*                                                      */
    /* Copy scan identifier info                            */
    /*                                                      */
    l2rec->iscan  = l1rec->iscan;
    l2rec->mside  = l1rec->mside;
    l2rec->detnum = l1rec->detnum;

    /*                                                      */
    /* Set pointers to ancillary data                       */
    /*                                                      */
    l2rec->ws = l1rec->ws;
    l2rec->wd = l1rec->wd;
    l2rec->mw = l1rec->mw;
    l2rec->zw = l1rec->zw;
    l2rec->pr = l1rec->pr;
    l2rec->oz = l1rec->oz;
    l2rec->wv = l1rec->wv;
    l2rec->rh = l1rec->rh;

    l2rec->no2_tropo = l1rec->no2_tropo;
    l2rec->no2_strat = l1rec->no2_strat;
    l2rec->no2_frac = l1rec->no2_frac;

    l2rec->height = l1rec->height;
    l2rec->elev   = l1rec->elev;
    l2rec->sstref = l1rec->sstref;      
    l2rec->ssttype = l1rec->ssttype;      
    l2rec->sssref = l1rec->sssref;      

    l2rec->rho_cirrus = l1rec->rho_cirrus;      

    l2rec->iwave = l1rec->iwave;
    l2rec->fwave = l1rec->fwave;
    l2rec->Fo    = l1rec->Fo;
    l2rec->Fobar = l1rec->Fobar;
    l2rec->Fonom = l1rec->Fonom;
    l2rec->Tau_r = l1rec->Tau_r;
    l2rec->k_oz  = l1rec->k_oz;
    l2rec->aw    = l1rec->aw;
    l2rec->bbw   = l1rec->bbw;

    l2rec->delphi  = l1rec->delphi;
    l2rec->alpha   = l1rec->alpha;
    l2rec->csolz   = l1rec->csolz;
    l2rec->csenz   = l1rec->csenz;
    l2rec->scattang = l1rec->scattang;

    l2rec->mask    = l1rec->mask;
    l2rec->pixnum  = l1rec->pixnum;
    l2rec->slot    = l1rec->slot;
    l2rec->nobs    = l1rec->nobs;

    /*                                                      */
    /* Copy scan time                                       */
    /*                                                      */
    *l2rec->year = *l1rec->year;
    *l2rec->day  = *l1rec->day;
    *l2rec->msec = *l1rec->msec;

    l2rec->tilt   = l1rec->tilt;
    l2rec->fsol   = l1rec->fsol;
    l2rec->flags  = l1rec->flags;
    l2rec->alt    = l1rec->alt;

    /*                                                      */
    /* Copy geolocation and view geometry and TOA radiance  */
    /*                                                      */
    memcpy(l2rec->lon ,l1rec->lon ,sizeof(float)*npix);
    memcpy(l2rec->lat ,l1rec->lat ,sizeof(float)*npix);
    memcpy(l2rec->solz,l1rec->solz,sizeof(float)*npix);
    memcpy(l2rec->sola,l1rec->sola,sizeof(float)*npix);
    memcpy(l2rec->senz,l1rec->senz,sizeof(float)*npix);
    memcpy(l2rec->sena,l1rec->sena,sizeof(float)*npix);
    memcpy(l2rec->Lt  ,l1rec->Lt  ,sizeof(float)*npix*l1rec->nbands);
    memcpy(l2rec->Lt_unc,l1rec->Lt_unc,sizeof(float)*npix*l1rec->nbands);

    /*                                                      */
    /* Copy pointers to precomputed atmospheric quantities  */
    /*                                                      */
    l2rec->t_h2o        = l1rec->t_h2o;
    l2rec->t_o2         = l1rec->t_o2;
    l2rec->tg_sol       = l1rec->tg_sol;
    l2rec->tg_sen       = l1rec->tg_sen;
    l2rec->t_sol        = l1rec->t_sol;
    l2rec->t_sen        = l1rec->t_sen;
    l2rec->rhof         = l1rec->rhof;
    l2rec->tLf          = l1rec->tLf;
    l2rec->Lr           = l1rec->Lr;
    l2rec->L_q          = l1rec->L_q;
    l2rec->L_u          = l1rec->L_u;
    l2rec->sw_n         = l1rec->sw_n;
    l2rec->sw_a         = l1rec->sw_a;
    l2rec->sw_bb        = l1rec->sw_bb;
    l2rec->sw_a_avg     = l1rec->sw_a_avg;
    l2rec->sw_bb_avg    = l1rec->sw_bb_avg;
    l2rec->polcor       = l1rec->polcor;
    l2rec->dpol         = l1rec->dpol;
    l2rec->glint_coef   = l1rec->glint_coef;
    l2rec->cloud_albedo = l1rec->cloud_albedo;
    l2rec->aerindex     = l1rec->aerindex;      
    l2rec->rhos         = l1rec->rhos;
    l2rec->Ltir         = l1rec->Ltir;
    l2rec->Bt           = l1rec->Bt;
    l2rec->pixdet       = l1rec->pixdet;
    l2rec->radcor       = l1rec->radcor;
    l2rec->in_prods     = l1rec->in_prods;

    return;
}


/* ---------------------------------------------------------- */
/* Converts a sensor-generic level-1b record to level-2       */ 
/*                                                            */
/* B. A. Franz, GSC, SIMBIOS Project, March 1998              */
/* W. Robinson, SAIC  15 Dec 2006  fix Western, Eastern most  */
/*     long for CZCS                                          */
/* ---------------------------------------------------------- */

int convl12( l1str *l1rec, l2str *l2rec, int32_t spix, int32_t epix, 
             instr *input, aestr *aerec)
{
    int32_t ip;                          /* Pixel index        */
    int32_t ib;                          /* Band index         */
    int32_t status;                      /* 0=OK, 1=bad        */
    int32_t npix = l1rec->npix;
    int32_t ipb;
    int day, year;

    /*                                                      */
    /* Clear the L2 record                                  */
    /*                                                      */
    init_l2(l2rec, l1rec->nbands);

    /*                                                      */
    /* Redefine npix of l2 record, as we may have subsamp'd */
    /*                                                      */
    l2rec->npix = npix;

    /*                                                         */
    /* Tranfer precomputed quantities and common-data pointers */ 
    /*                                                         */
    cpl1l2(l1rec,l2rec);

    /* Point L2 rec to computed SST, if requested (before atmcor) */
    if (input->proc_sst)
       l2rec->sst = get_sst(l2rec);
    else
        l2rec->sst = NULL;

    /* if glint_opt = 2 (use simple glint) set bring l1rec TLg to l2rec */
    if (input->glint_opt == 2)
      memcpy( l2rec->TLg, l1rec->TLg, sizeof(float)*npix*l1rec->nbands );

    /*                                                      */
    /* Loop through each pixel and do atmospheric correction*/ 
    /*                                                      */
    for (ip=spix; ip<=epix; ip++) {

        /* ------------------------------------------------ */
        /* Ocean processing                                 */
        /* ------------------------------------------------ */
        if ((input->proc_ocean != 0) && 
            !l1rec->mask[ip] && 
            l1rec->solz[ip] < SOLZNIGHT) {

	    
            if(input->atmocor) {

                /* set aerosol values from input rec, if supplied */
                if (input->aer_opt == AERSOA || input->aer_opt == AERSMA)
                    status = run_soa_sma(l2rec, ip);
                else
                    status = atmocor2(l2rec,aerec,ip);

                /*                                                           */
                /* If the atmospheric correction failed, flag and mask. Else,*/
                /* set flags which depend on complete atmospheric correction.*/
                /*                                                           */
                if (status != 0) {
                    l2rec->flags[ip] |= ATMFAIL;
                    l2rec->mask[ip] = 1;
                } else {
                    setflagbits(2,NULL,l2rec,ip);
                }
	    
            }

        } // if ocean

    } // for ip


    /* Load L2 rec with inherent optical properties */
    if (input->iop_opt > 0 && (input->proc_ocean != 0) && input->atmocor)
        get_iops(l2rec,input->iop_opt);

    return(0);
}



/* --------------------------------------------------------------- */
/* get_iops.c - load IOP (a & bb) fields in L2 rec.                */
/*                                                                 */
/* Inputs:                                                         */
/*     l2rec - level-2 structure containing one complete scan      */
/*             after atmospheric correction.                       */
/* Outputs:                                                        */
/*     iop_opt - algorithm selector                                */
/*                                                                 */
/* Written By: B. Franz, NASA/OBPG/SAIC, 25 Feb 2005               */
/*                                                                 */
/* --------------------------------------------------------------- */
void get_iops(l2str *l2rec, int32_t iop_opt)
{
    switch (iop_opt) {
      case IOPCARDER:
        iops_carder(l2rec);
        break;
      case IOPGSM:
        iops_gsm(l2rec);
        break;
      case IOPQAA:
        iops_qaa(l2rec);
        break;
      case IOPPML:
        iops_pml(l2rec);
        break;
      case IOPLAS:
        iops_las(l2rec);
        break;
      case IOPNIWA:
        iops_niwa(l2rec);
        break;
      case IOPGIOP:
        iops_giop(l2rec);
        break;
      case IOPSWIM:
    	iops_swim(l2rec);
    	break;
      default:
        break;
    }   

    return;
}
