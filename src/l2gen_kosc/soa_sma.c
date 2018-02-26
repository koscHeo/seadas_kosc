#include "l12_proto.h"

void get_soa_prod_(int32_t *prodID, int32_t *pix, float *prod);
void get_sma_prod_(int32_t *prodID, int32_t *pix, float *prod);
void clear_soa_prod_(void);
void clear_sma_prod_(void);


/* ------------------------------------------------------------------ */
/* Run sprectral optimization algorithm on one pixel                  */
/* ------------------------------------------------------------------ */
int run_soa_sma(l2str *l2rec, int32_t ip)
{
    static float pi        = PI;
    static float radeg     = RADEG;
    static float badval    = BAD_FLT;
    static int   firstCall = 1; 

    static float *aphstar;

    float adg_s   = l2rec->input->gsm_adg_s;
    float bbp_s   = l2rec->input->gsm_bbp_s;

    int32_t  nwave   = l2rec->nbands;
    float *wave   = l2rec->fwave;
    float *Fo     = l2rec->Fo;
    float fsol    = l2rec->fsol;
    float *aw     = l2rec->aw;
    float *bbw    = l2rec->bbw;
    float *Tau_r  = l2rec->Tau_r;

    float solz    = l2rec->solz   [ip];
    float senz    = l2rec->senz   [ip];
    float delphi  = l2rec->delphi [ip];
    float lon     = l2rec->lon    [ip];
    float lat     = l2rec->lat    [ip];
    float wv      = l2rec->wv     [ip];

    float *Lw     = &l2rec->Lw    [ip*nwave];
    float *nLw    = &l2rec->nLw   [ip*nwave];
    float *Rrs    = &l2rec->Rrs   [ip*nwave];
    float *La     = &l2rec->La    [ip*nwave];
    float *taua   = &l2rec->taua  [ip*nwave];
    float *t_sol  = &l2rec->t_sol [ip*nwave];
    float *t_sen  = &l2rec->t_sen [ip*nwave];
    float *TLg    = &l2rec->TLg   [ip*nwave];
    float *tg_sol = &l2rec->tg_sol[ip*nwave];
    float *tg_sen = &l2rec->tg_sen[ip*nwave];

    int32_t  ib, ipb;
    int32_t  status;
    float mu0;
    int32_t  nneg;
    float *Ltemp  ;

    float *Rhor_in, *Ra_out,*Rw_out,*tstar,*tstar0,*optTaua,*optW0,\
           *optChl,*optAcdm,*optBbp,*optcdmtotot,*optMr,*optMi,*optV;


    /* Set static GSM parameters */
    if (firstCall) {

        int32_t   iw;
        int32_t   taphn = -1;
        float *taphw = NULL;
        float *taphs = NULL;

        firstCall = 0;
        if ( (aphstar = (float *)calloc(nwave,sizeof(float))) == NULL) {
            printf("-E- : Error allocating memory to dray_for_i\n");
            exit(FATAL_ERROR);
        }

        /* interpolate aphstar */
        taphw = l2rec->input->gsm_aphw;  /* aphstar table wavelengths  */
        taphs = l2rec->input->gsm_aphs;  /* aphstar table coefficients */
        taphn = 0;                       /* aphstar table size         */
        for (taphn=0; taphn<nwave; taphn++)
            if (taphw[taphn] < 0) break;
        for (iw=0; iw<nwave; iw++) {
            if (wave[iw] < (taphw[taphn-1]+10))
                aphstar[iw] = linterp(taphw,taphs,taphn,wave[iw]);
            else
                aphstar[iw] = 0.0;
        }
    }

    /* To be safe: disable out-of-band correction, since we get modeled rhow */
    /* from the SOA procedure, presumably at nominal center wavelengths (BAF)*/
    l2rec->input->outband_opt=0;

    /* To be safe: disable brdf correction.  It should be handled within SOA */
    l2rec->input->brdf_opt=0;


    /* Initialize output values. If any input radiances are negative, just return. */

    if (ip == 0) {
        clear_soa_prod_();
        clear_sma_prod_();
    }

    status = 0;
    nneg   = 0;
/* Allocate memory
 *
 */

    if ( (Rhor_in = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("-E- : Error allocating memory to Rhor_in\n");
        exit(FATAL_ERROR);
    }
    if ( (Ra_out = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("-E- : Error allocating memory to RA_out\n");
        exit(FATAL_ERROR);
    }
    if ( (Rw_out = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("-E- : Error allocating memory to Rw_out\n");
        exit(FATAL_ERROR);
    }
    if ( (tstar = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("-E- : Error allocating memory to tstar\n");
        exit(FATAL_ERROR);
    }
    if ( (tstar0 = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("-E- : Error allocating memory to tstar0\n");
        exit(FATAL_ERROR);
    }
    if ( (optTaua = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("-E- : Error allocating memory to optTaua\n");
        exit(FATAL_ERROR);
    }
    if ( (optW0 = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("-E- : Error allocating memory to optW0\n");
        exit(FATAL_ERROR);
    }
    if ( (optChl = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("-E- : Error allocating memory to dray_for_i\n");
        exit(FATAL_ERROR);
    }
    if ( (optAcdm = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("-E- : Error allocating memory to optAcdm\n");
        exit(FATAL_ERROR);
    }
    if ( (optBbp = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("-E- : Error allocating memory to optBdp\n");
        exit(FATAL_ERROR);
    }
    if ( (optcdmtotot = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("-E- : Error allocating memory to optcdmtotot\n");
        exit(FATAL_ERROR);
    }
    if ( (optMr = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("-E- : Error allocating memory to optMr\n");
        exit(FATAL_ERROR);
    }
    if ( (optMi = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("-E- : Error allocating memory to optMi\n");
        exit(FATAL_ERROR);
    }
    if ( (optV = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("-E- : Error allocating memory to optV\n");
        exit(FATAL_ERROR);
    }

    for (ib=0; ib<nwave; ib++) {

        ipb = ip*nwave + ib;

        t_sol [ib]  = 1.0;
        t_sen [ib]  = 1.0;
        TLg   [ib]  = 0.0;
        La    [ib]  = badval;
        Lw    [ib]  = badval;
        nLw   [ib]  = badval;
        Rrs   [ib]  = badval;
        taua  [ib]  = badval;

        if (l2rec->Lt[ipb] <= 0.0)
	    nneg++;
    }
    l2rec->eps[ip] = badval;


    /* If any expected channels are negative */
    if (nneg > 0) {
        status = 1;
        free(Rhor_in);
        free(Ra_out);
        free(Rw_out);
        free(tstar);
        free(tstar0);
        free(optTaua);
        free(optW0);
        free(optChl);
        free(optAcdm);
        free(optBbp);
        free(optcdmtotot);
        free(optMr);
        free(optMi);
        free(optV);
 
        return(status);
    }

    mu0 = cos( solz/radeg );
     
    if ((Ltemp = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("-E- : Error allocating memory to run_soa_sma\n");
        exit(FATAL_ERROR);
    }

    /* Remove pre-computed atmospheric effects */
    /* --------------------------------------- */
    for (ib=0; ib<nwave; ib++) {

        ipb = ip*nwave + ib;

	/* Copy TOA radiance to temp var, eliminate missing bands */
 	Ltemp[ib]  = l2rec->Lt[ipb];

        /* Remove whitecap radiance */
        Ltemp[ib] -= l2rec->tLf[ipb];

        /* Correct for ozone absorption.  We correct for inbound and outbound here, */
        /* then we put the inbound back when computing Lw.                          */
        Ltemp[ib] = Ltemp[ib]/l2rec->tg_sol[ipb]/l2rec->tg_sen[ipb];

        /*  Apply polarization correction */
        Ltemp[ib] /= l2rec->polcor[ipb];

        /*  Remove glint (currently set to zero) */
        Ltemp[ib] -= l2rec->TLg[ipb];

        /*  Solar irradiance and Rayleigh reflectance */
        Rhor_in[ib] = pi*l2rec->Lr[ipb]/Fo[ib]/mu0;
    }
      

    if (l2rec->input->aer_opt == AERSOA) {

      atmcor_soa_(  
	&(l2rec->sensorID),     
	sensorDir[l2rec->sensorID],
	&nwave,
	wave,
        &(l2rec->iscan),
        &ip,
        &solz,
        &senz,
        &delphi,
        &lat,
        &lon,
        Ltemp,
        Rhor_in,
        Fo,
        Tau_r,
        aw,
        bbw,
        aphstar,
        &adg_s,
        &bbp_s,
        Ra_out, 
        Rw_out,
        &wv,
        tstar,
        tstar0,
        optTaua,
        optW0,
        optChl, 
        optAcdm,
        optcdmtotot,
        optBbp, 
        optMr, 
        optMi,
        optV,
        &status);

    } else if (l2rec->input->aer_opt == AERSMA) {

      /* baf
      atmcor_sma_(
	Ltemp1,       
        &(l2rec->iscan),
        &ip,
        &solz,
        &senz,
        &delphi,
        &lat,
        &lon,
        Ltemp,
        Rhor_in,
        Fo,
        Ra_out, 
        Rw_out,
        &wv,
        tstar,
        tstar0,
        optTaua,
        optW0,
        optChl, 
        optAcdm, 
        optBbp, 
        optdom,
        &status);
      */   
      status = 1;
    }

    if (status == 0) {

        for (ib=0; ib<nwave; ib++) {

            ipb = ip*nwave + ib;

            t_sen[ib] = tstar[ib];
	        t_sol[ib] = tstar0[ib];
            taua [ib] = optTaua[ib];
            La   [ib] = Ra_out[ib]*Fo[ib]/pi;
            nLw  [ib] = Rw_out[ib]*Fo[ib]/pi;
            Rrs  [ib] = Rw_out[ib]/pi;
	    Lw   [ib] = nLw[ib]*t_sol[ib]*tg_sol[ib]*mu0*fsol;
        }

        l2rec->eps[ip] = 1.0;

    } else {

        l2rec->flags[ip] |= ATMFAIL;
    }

    /* Compute final chl from final nLw (needed for flagging) */
    l2rec->chl[ip] = get_default_chl(l2rec,Rrs);

    free(Rhor_in);
    free(Ra_out);
    free(Rw_out);
    free(tstar);
    free(tstar0);
    free(optTaua);
    free(optW0);
    free(optChl);
    free(optAcdm);
    free(optBbp);
    free(optcdmtotot);
    free(optMr);
    free(optMi);
    free(optV);
    free(Ltemp);
    return(status);
}

   


/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
void get_soa(l2str *l2rec, int32_t prodID, float prod[])
{
    int32_t   ip;

    for (ip=1; ip<=l2rec->npix; ip++)
        if (l2rec->mask[ip-1] == 0)
            get_soa_prod_(&prodID,&ip,&(prod[ip-1]));
	else
	    prod[ip-1] = -1.0;
      
    return;
}

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
void get_sma(l2str *l2rec, int32_t prodID, float prod[])
{
    int32_t   ip;

    for (ip=1; ip<=l2rec->npix; ip++)
        if (l2rec->mask[ip-1] == 0)
            get_sma_prod_(&prodID,&ip,&(prod[ip-1]));
	else
	    prod[ip-1] = -1.0;
      
    return;
}

/*    printf("%10.4f",optChl) */
