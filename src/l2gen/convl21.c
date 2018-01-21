#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "l12_proto.h"
#include "filehandle.h"

/* ---------------------------------------------------------- */
/* Converts a sensor-generic level-2 record to level-1        */ 
/*                                                            */
/* B. A. Franz, GSC, SIMBIOS Project, August 1998             */
/* ---------------------------------------------------------- */

int convl21( l2str *l2rec, tgstr *tgrec, int32_t spix, int32_t epix, instr *input, 
        float *vLt, vcstr *vrec)
{
    static int    firstCall =  1;
    static int    vcal_opt  = -1;
    static float  *brdfsensor;
    static float  *brdfinsitu;
    static float  *F0        ;
    static float  *wave      ;
    static int32_t   nwvis;

    int32_t   ip;                          /* Pixel index       */
    int32_t   ib;                          /* Band index        */
    int32_t   ipb;                         /* Combined index    */
    float  tLw;
    float  solz_insitu;
    float  mu0_sensor;
    float  mu0_insitu;
    float  chl;
    float  tau;
    float  *Rrs;
    float  *nLw;
    float  *Lw;
    int    foundneg;

    /* Initialize static vars */
    if (firstCall) {
        firstCall = 0;
        if ((brdfsensor = (float *) calloc(l2rec->nbands,sizeof(float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for brdfsensor in convl21.\n",
                    __FILE__,__LINE__);
            exit(1);
        }
        if ((brdfinsitu = (float *) calloc(l2rec->nbands,sizeof(float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for brdfinsitu in convl21.\n",
                    __FILE__,__LINE__);
            exit(1);
        }
        if ((F0 = (float *) calloc(l2rec->nbands,sizeof(float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for F0 in convl21.\n",
                    __FILE__,__LINE__);
            exit(1);
        }
        if ((wave = (float *) calloc(l2rec->nbands,sizeof(float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for wave in convl21.\n",
                    __FILE__,__LINE__);
            exit(1);
        }

        for (ib=0; ib<l2rec->nbands; ib++) {
            brdfsensor[ib] = 1.0;
            brdfinsitu[ib] = 1.0;
            /* Disabling. Assume full-band target data.
            if (input->outband_opt >= 2) 
                F0[ib] = l2rec->Fonom[ib];
            else
                F0[ib] = l2rec->Fobar[ib];
             */
            F0[ib] = l2rec->Fobar[ib];
        }
        for (ib=0; ib<l2rec->nbands; ib++)
            wave[ib] = l2rec->fwave[ib];
        nwvis = rdsensorinfo(l2rec->sensorID,l2rec->input->evalmask,"NbandsVIS",  NULL);

        if (input->mode == INVERSE_LW || input->vcal_opt == INVERSE_LW)
            vcal_opt = INVERSE_LW;
        else if (input->mode == INVERSE_NLW || input->vcal_opt == INVERSE_NLW)
            vcal_opt = INVERSE_NLW;
        else if (input->mode == INVERSE_ZERO || input->vcal_opt == INVERSE_ZERO)
            vcal_opt = INVERSE_ZERO;
        else {
            printf("%s: Unknown calibration inversion option %d %d\n",
                    __FILE__, input->mode, input->vcal_opt);
            exit(1);
        }
    }
    if ((Rrs = (float *) calloc(l2rec->nbands,sizeof(float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for Rrs in convl21.\n",
                __FILE__,__LINE__);
        exit(1);
    }
    if ((nLw = (float *) calloc(l2rec->nbands,sizeof(float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for nLw in convl21.\n",
                __FILE__,__LINE__);
        exit(1);
    }
    if ((Lw = (float *) calloc(l2rec->nbands,sizeof(float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for Lw in convl21.\n",
                __FILE__,__LINE__);
        exit(1);
    }

    /* Initialize radiances to zero */
    memset(vLt,0,sizeof(float)*l2rec->npix*l2rec->nbands);

    /*                                                      */
    /* Loop through each pixel and reconstruct the TOA      */
    /* radiances at each band from the components of the    */
    /* atmospheric correction and the nLw                   */ 
    /*                                                      */
    for (ip=spix; ip<=epix; ip++) {

        if (vrec != NULL) {
            for (ib=0; ib<l2rec->nbands; ib++) {
                ipb = ip*l2rec->nbands + ib;
                vrec->vLt [ipb] = BAD_FLT;
                vrec->tLw [ipb] = BAD_FLT;
                vrec->Lw  [ipb] = BAD_FLT;
                vrec->nLw [ipb] = BAD_FLT;
                vrec->brdfsat[ipb] = BAD_FLT;
                vrec->brdftgt[ipb] = BAD_FLT;
            }
        }

        /* If atmospheric corr failed, go to next pixel     */
        if (l2rec->mask[ip] || (l2rec->flags[ip] & ATMFAIL) != 0)
            continue;

        /* If any target radiances are negative, go to next */
        foundneg=0;
        if (tgrec != NULL) {
            for (ib=0; ib<l2rec->nbands; ib++) {
                ipb = ip*l2rec->nbands + ib;
                if (vcal_opt == INVERSE_LW) {
                    if (tgrec->Lw[ipb] < 0.0)
                        foundneg=1;
                } else if (vcal_opt == INVERSE_NLW) {
                    if (tgrec->nLw[ipb] < 0.0)
                        foundneg=1;
                }
            }
        }
        if (foundneg) continue;

        /* Compute cos(solz) for sensor and in situ */
        mu0_sensor = cos(l2rec->solz[ip]/RADEG);
        if (tgrec != NULL)
            if (tgrec->solz[ip] >= 0.0)
                solz_insitu = tgrec->solz[ip];
            else
                solz_insitu = l2rec->solz[ip];
        else {
            if (input->vcal_solz >= 0.0)
                solz_insitu = input->vcal_solz;
            else
                solz_insitu = l2rec->solz[ip];
        }
        mu0_insitu = cos(solz_insitu/RADEG);


        /* Build target nLw from target Lw, using retrieved atmosphere */
        for (ib=0; ib<l2rec->nbands; ib++) {
            ipb  = ip*l2rec->nbands+ib;
            if (vcal_opt == INVERSE_LW) {
                tau  = -log(l2rec->tg_sol[ipb]*l2rec->t_sol[ipb])
                             * mu0_sensor;
                if (tgrec != NULL) 
                    Lw[ib] = tgrec->Lw[ipb];
                else
                    Lw[ib] = input->vcal_Lw[ib];
                nLw[ib] = Lw[ib]
                             / mu0_insitu
                             / exp(-tau/mu0_insitu)
                             / l2rec->fsol;
            } else if (vcal_opt == INVERSE_NLW) {
                if (tgrec != NULL) 
                    nLw[ib] = tgrec->nLw[ipb];
                else
                    nLw[ib] = input->vcal_nLw[ib];
            } else {
                nLw[ib] = 0.0;
            }
            Rrs[ib] = nLw[ib]/F0[ib];
        }

        /* Get chlorophyll from target nLw (if not specified) */
        if (vcal_opt != INVERSE_ZERO) {
            if (input->vcal_chl >= 0.0) 
                chl = input->vcal_chl;
            else {
                chl = get_default_chl(l2rec,Rrs);
                if (chl < 0.0) continue;
            }
        } else {
            chl = 0.0;
        }

        /* Compute f/Q corrections, if requested */ 
        if (input->brdf_opt != NOBRDF) {

            /* correction from sensor geometry to nadir view, zero sun angle */
            ocbrdf(l2rec,ip,l2rec->input->brdf_opt,wave,nwvis,
                    l2rec->solz[ip],l2rec->senz[ip],l2rec->delphi[ip],l2rec->ws[ip],
                    chl,nLw,F0,brdfsensor);

            /* correction from in situ geometry to nadir view, zero sun angle */
            if (vcal_opt == INVERSE_LW) {
                ocbrdf(l2rec,ip,l2rec->input->brdf_opt,wave,nwvis,
                        solz_insitu,0.0,0.0,l2rec->ws[ip],
                        chl,nLw,F0,brdfinsitu);
            }
        }

        /*                                                  */
        /* OK, loop through each band                       */
        /*                                                  */
        for (ib=0; ib<l2rec->nbands; ib++) {

            ipb  = ip*l2rec->nbands+ib;

            /*                                                     */
            /* The target type controls how tLw is reconstructed   */
            /*                                                     */
            switch (vcal_opt) {
            case INVERSE_ZERO:
                tLw = 0.0;
                break;
            case INVERSE_NLW:
                tLw = nLw[ib]
                          / brdfsensor[ib]
                                       * l2rec->polcor[ipb]
                                                       * l2rec->t_sol[ipb]
                                                                      * l2rec->t_sen[ipb]
                                                                                     * l2rec->tg_sol[ipb]
                                                                                                     * l2rec->tg_sen[ipb]
                                                                                                                     * mu0_sensor
                                                                                                                     * l2rec->fsol;
                break;
            case INVERSE_LW:
                tLw = nLw[ib] * brdfinsitu[ib]
                                           / brdfsensor[ib]
                                                        * l2rec->polcor[ipb]
                                                                        * l2rec->t_sol[ipb]
                                                                                       * l2rec->t_sen[ipb]
                                                                                                      * l2rec->tg_sol[ipb]
                                                                                                                      * l2rec->tg_sen[ipb]
                                                                                                                                      * mu0_sensor
                                                                                                                                      * l2rec->fsol;
                break;
            }


            vLt[ipb] = tLw
                    + ((( l2rec->TLg[ipb]
                                     + l2rec->La[ipb] )
                            * l2rec->t_o2[ipb]
                                          + l2rec->tLf[ipb]
                                                       + l2rec->Lr [ipb] )
                            * l2rec->polcor[ipb]
                                            - l2rec->radcor[ipb])
                                            * l2rec->tg_sol[ipb]
                                                            * l2rec->tg_sen[ipb];

            if (vrec != NULL) {
                vrec->tLw [ipb] = tLw;
                vrec->Lw  [ipb] = Lw[ib];
                vrec->nLw [ipb] = nLw[ib];
                vrec->brdfsat[ipb] = brdfsensor[ib];
                vrec->brdftgt[ipb] = brdfinsitu[ib];
            }
        }
    }

    free(Rrs);
    free(nLw);
    free(Lw);

    return(0);
}
