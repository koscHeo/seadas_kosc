#include "l12_proto.h"

static int ib412;
static int ib443;
static int ib490;
static int ib510;
static int ib555;
static int ib670;
static int ib865;
static int ibcloud;
static float *wave;

char    isCoccolith ( l2str *l2rec, int32_t ip );
char    isTurbid    ( l2str *l2rec, int32_t ip );
float   aerindex    ( l2str *l2rec, int32_t ip );


void l1_mask_set(l1str *l1rec, int32_t ip)
{
      l1rec->mask[ip] = (l1rec->landMaskOn    && l1rec->land   [ip]) ||
                        (l1rec->bathMaskOn    && l1rec->swater [ip]) || 
                        (l1rec->cloudMaskOn   && l1rec->cloud  [ip]) || 
                        (l1rec->glintMaskOn   && l1rec->glint  [ip]) ||
                        (l1rec->stlightMaskOn && l1rec->stlight[ip]) ||
                        (l1rec->senzMaskOn    && l1rec->senzmax[ip]) ||
                        (l1rec->solzMaskOn    && l1rec->solzmax[ip]) ||
                        (l1rec->hiltMaskOn    && l1rec->hilt   [ip]) ||
                        (l1rec->filter[ip])   ||
                        (l1rec->navfail[ip]);
}

int setflags( instr *input, l1str *l1rec )
{
    static float pi    = PI;
    static float radeg = RADEG;
    static int checkCirrus = 0;
    static int firstCall   = 1;
    int32 evalmask = input->evalmask;
    float mu0;
    int32_t  ip, ipb, iw;
    float albedo;
    int32_t  nwave = l1rec->nbands;

    if (firstCall) {
        firstCall = 0;
        if ((wave = (float *)calloc(nwave,sizeof(float))) == NULL) {
            printf("-E- : Error allocating memory to setflags\n");
            exit(FATAL_ERROR);
        }

        for (iw=0; iw<nwave; iw++) {
            wave[iw] = l1rec->fwave[iw];
        }
        ib412 = windex(412.,wave,nwave);
        ib443 = windex(443.,wave,nwave);
        ib490 = windex(490.,wave,nwave);
        ib510 = windex(510.,wave,nwave);
        ib555 = windex(550.,wave,nwave);
        ib670 = windex(670.,wave,nwave);
        ib865 = windex(865.,wave,nwave);
        ibcloud = windex(input->cloud_wave,wave,nwave);

        printf("\nUsing %6.1f nm channel for cloud flagging over water.\n",wave[ibcloud]);
        printf("Using %6.1f nm channel for cloud flagging over land.\n\n",wave[ib412]);

        if (input->cirrus_thresh[0] >= 0.0 && input->cirrus_thresh[1] > input->cirrus_thresh[0]) {
            printf("\nApplying cirrus test for cloud flagging.\n\n");
            checkCirrus = 1;
        }
    }

    /*                                                     */
    /* Now set flags for each pixel                        */
    /*                                                     */
    for (ip=0; ip<l1rec->npix; ip++) {

      mu0     = cos( l1rec->solz[ip] / radeg );

      /* Check view angle limits */
      if (l1rec->senz[ip] > input->satzen) {
          l1rec->senzmax[ip] = ON;
      } 
      if (l1rec->solz[ip] > input->sunzen) {
          l1rec->solzmax[ip] = ON;
      } 
      

      /* Check for glint */
      if (!l1rec->land[ip] && l1rec->solz[ip] < SOLZNIGHT) {
          if (l1rec->glint_coef[ip] > input->glint) {
              l1rec->glint[ip] = ON;
          }
      }


      /* Check for clouds (daytime only) */
      if (l1rec->solz[ip] < SOLZNIGHT) {

          if (!l1rec->land[ip]) { 
              l1rec->cloud_albedo[ip] = l1rec->rhos[ip*nwave+ibcloud]
                  - l1rec->TLg[ip*nwave+ibcloud]*pi/mu0/l1rec->Fo[ibcloud];
              albedo = input->albedo;
          } else {
              l1rec->cloud_albedo[ip] = l1rec->rhos[ip*nwave+ib412];
              albedo = input->albedo * 4.0;
          }  
          if ((evalmask & MODCLOUD) != 0) {
              l1rec->cloud[ip] = get_cldmask(l1rec,ip);
          }else{
              if (l1rec->cloud_albedo[ip] > albedo)
                 l1rec->cloud[ip] = ON;
          }
          /* Eval switch to enable spectral cloud test, recover bright water */
          if ((input->cloud_eps > 0.0) != 0) {
              float min_cld_rhos   = 1.0;
              float max_cld_rhos   = 1.0;
              float cld_rhos_ratio = 0.0;
              if (l1rec->cloud[ip] == ON && !l1rec->land[ip] && l1rec->cloud_albedo[ip] < albedo*10) {
                  min_cld_rhos = MIN(MIN(MIN(
                      l1rec->rhos[ip*nwave+ib865],
                      l1rec->rhos[ip*nwave+ib670]),
                      l1rec->rhos[ip*nwave+ib555]),
                      l1rec->rhos[ip*nwave+ib412]);
                  max_cld_rhos = MAX(MAX(MAX(
                      l1rec->rhos[ip*nwave+ib865],
                      l1rec->rhos[ip*nwave+ib670]),
                      l1rec->rhos[ip*nwave+ib555]),
                      l1rec->rhos[ip*nwave+ib412]);
                  if (min_cld_rhos > 0.0) {
                      cld_rhos_ratio = max_cld_rhos/min_cld_rhos;
	              if (cld_rhos_ratio > input->cloud_eps  && (evalmask & MODCLOUD) == 0)
                          l1rec->cloud[ip] = OFF;
	          }
              }
          }

      }

      /* use MODIS cloud mask instead */
      /*
      if ((input->evalmask & MODCLDMSK) != 0) {
	 l1rec->cloud[ip] = (char) modis_cloud_mask(l1rec,ip);
      }
      */

      /* expand cloud flag with MODIS cirrus mask */
      if ((input->evalmask & MODCIRRUS) != 0 && (evalmask & MODCLOUD) == 0) {
	 l1rec->cloud[ip] |= (char) modis_cirrus_mask(l1rec,ip);
      }

      /* expand cloud flag with internal cirrus reflectance test */
      if (checkCirrus) {
          if (l1rec->rho_cirrus[ip] >= input->cirrus_thresh[0] && l1rec->rho_cirrus[ip] < input->cirrus_thresh[1]  && (evalmask & MODCLOUD) == 0)
             l1rec->cloud[ip] = 1;
      }

      /* Check for dark pixel */
      if (!l1rec->land[ip] && l1rec->solz[ip] < SOLZNIGHT) {
        for (iw=0; iw<nwave; iw++) {
          ipb = ip*nwave+iw;
          if ( l1rec->Lt[ipb] / 
               l1rec->tg_sol[ipb] /
               l1rec->tg_sen[ipb] -
               l1rec->Lr[ipb] < 0.0 ) {
              l1rec->darkpix[ip] = 1;
              break;
          }
	}
      }


      /* Set masking */
      l1_mask_set(l1rec,ip);

    }

    return(1);
}


void setflagbits(int level, l1str *l1rec, l2str *l2rec, int32_t ipix)
{
    int32_t  npix, spix, epix;
    int32_t  ip, ib, iw;
    int32_t  nwave;

    if (l1rec != NULL) {
        npix  = l1rec->npix;
        nwave = l1rec->nbands;
    } else if (l2rec != NULL) {
        npix  = l2rec->npix;
        nwave = l2rec->nbands;
    } else {
        printf("-E- %s line %d: attempt to set flags from NULL record.\n",
               __FILE__,__LINE__);
        exit(FATAL_ERROR);
    }

    if (ipix < 0) {
        spix = 0;
        epix = npix-1;
    } else {
        spix = ipix;
        epix = ipix;
    }

    /* Level-0, Level-1A flag bits (require readl1) */
    if (level == 0) {

      for (ip=spix; ip<=epix; ip++) {

        if (l1rec == NULL) {
            printf("-E- %s line %d: attempt to set flags from NULL record.\n",
               __FILE__,__LINE__);
            exit(FATAL_ERROR);
	}

        if (l1rec->hilt   [ip]) l1rec->flags[ip] |= HILT;
        if (l1rec->stlight[ip]) l1rec->flags[ip] |= STRAYLIGHT;
        if (l1rec->navfail[ip]) l1rec->flags[ip] |= NAVFAIL;
        if (l1rec->navwarn[ip]) l1rec->flags[ip] |= NAVWARN;

      }
    }

    /* Level-1B flag bits (require loadl1) */
    else if (level == 1) {

      for (ip=spix; ip<=epix; ip++) {

        if (l1rec == NULL) {
            printf("-E- %s line %d: attempt to set flags from NULL record.\n",
               __FILE__,__LINE__);
            exit(FATAL_ERROR);
	}

        if (l1rec->land   [ip]) l1rec->flags[ip] |= LAND;
        if (l1rec->swater [ip]) l1rec->flags[ip] |= COASTZ;
        if (l1rec->cloud  [ip]) l1rec->flags[ip] |= CLOUD;
        if (l1rec->ice    [ip]) l1rec->flags[ip] |= SEAICE;
        if (l1rec->glint  [ip]) l1rec->flags[ip] |= HIGLINT;
        if (l1rec->solzmax[ip]) l1rec->flags[ip] |= HISOLZEN;
        if (l1rec->senzmax[ip]) l1rec->flags[ip] |= HISATZEN;
        if (l1rec->filter [ip]) l1rec->flags[ip] |= FILTER;
        if (l1rec->glint_coef[ip] > GLINT_MIN) 
            l1rec->flags[ip] |= MODGLINT;
        for (iw=0; iw<nwave; iw++) {
  	    ib = iw;
	    if (l1rec->dpol[ip*nwave+ib] > l1rec->input->hipol) {
                l1rec->flags[ip] |= HIPOL;
                break;
	    }
	}
      }

    }

    /* Level-2 flag bits (require atmocor2) */
    else if (level == 2) {
      
      for (ip=spix; ip<=epix; ip++) {

        if (l2rec == NULL) {
            printf("-E- %s line %d: attempt to set flags from NULL record.\n",
               __FILE__,__LINE__);
            exit(FATAL_ERROR);
        }

        if (l2rec->eps[ip] < l2rec->input->epsmin || 
            l2rec->eps[ip] > l2rec->input->epsmax) {
            l2rec->flags[ip] |= ATMWARN;
	}

        if ( (l2rec->Lw[ip*nwave+ib490] < 0.0) ||
             (l2rec->Lw[ip*nwave+ib510] < 0.0) ||
             (l2rec->Lw[ip*nwave+ib555] < 0.0) ) {
            l2rec->flags[ip] |= ATMWARN;
        }

        if (l2rec->nLw[ip*l2rec->nbands+ib555] < l2rec->input->nlwmin)
            l2rec->flags[ip] |= LOWLW;

        if (isCoccolith(l2rec,ip))
            l2rec->flags[ip] |= COCCOLITH;

        if (l2rec->input->aer_opt != AERWANGSWIR && l2rec->input->aer_opt != AERRHSWIR) {
            if (isTurbid(l2rec,ip))
                l2rec->flags[ip] |= TURBIDW;
	}

        if (l2rec->chl[ip] == BAD_FLT)
            l2rec->flags[ip] |= CHLFAIL;

        else if (l2rec->chl[ip] > CHL_MAX || l2rec->chl[ip] < CHL_MIN)
            l2rec->flags[ip] |= CHLWARN;

        if (l2rec->num_iter[ip] > l2rec->input->aer_iter_max) {
            l2rec->flags[ip] |= MAXAERITER;
            l2rec->flags[ip] |= ATMWARN;
	}

        if (l2rec->input->absaer_opt >= 0) {
            l2rec->aerindex[ip] = aerindex(l2rec,ip);
            if (l2rec->aerindex[ip] < l2rec->input->absaer)
                l2rec->flags[ip] |= ABSAER;
	}
      }
    }

    else {
        printf("-E- %s line %d: attempt to set flags at bogus level.\n",
               __FILE__,__LINE__);
        exit(FATAL_ERROR);
    }
}


#define Between(a,x,b)(x >= a && x <= b)
char isCoccolith( l2str *l2rec, int32_t ip )
{
    static float firstCall = 1;
    static float c1,c2,c3,c4,c5,c6,c7,c8;

    float  nLw443 = l2rec->nLw[ip*l2rec->nbands+ib443];
    float  nLw510 = l2rec->nLw[ip*l2rec->nbands+ib510];
    float  nLw555 = l2rec->nLw[ip*l2rec->nbands+ib555];
    float   La670 = l2rec->La [ip*l2rec->nbands+ib670];

    if (firstCall) {
        firstCall = 0;
        c1 = l2rec->input->coccolith[0];
        c2 = l2rec->input->coccolith[1];
        c3 = l2rec->input->coccolith[2];
        c4 = l2rec->input->coccolith[3];
        c5 = l2rec->input->coccolith[4];
        c6 = l2rec->input->coccolith[5];
        c7 = l2rec->input->coccolith[6];
        c8 = l2rec->input->coccolith[7];
    }

    if (nLw443 >= c1                 &&
        nLw510 >= 0.0                &&
        nLw555 >= c2                 &&
         La670 <= 1.1                &&
        Between(c3,nLw443/nLw555,c4) &&
        Between(c5,nLw510/nLw555,c6) &&
        Between(c7,nLw443/nLw510,c8) )

        return(1);
    else
        return(0);
}


char isTurbid( l2str *l2rec, int32_t ip )
{
    
    static float Fo;     
    static int   firstCall = 1;

    if (firstCall) {
        firstCall = 0;
        if (l2rec->input->outband_opt >= 2)
            Fo = l2rec->Fonom[ib670];
	else
            Fo = l2rec->Fobar[ib670];
    }

    if (l2rec->nLw[ip*l2rec->nbands+ib670]/Fo > 0.0012)
        return(1);
    else
        return(0);
}


/* ---------------------------------------------------------------------------------------- */
/* aerindx() - compute aerosol index for absorbing aerosol test                             */
/* ---------------------------------------------------------------------------------------- */
float aerindex(l2str *l2rec, int32_t ip)
{
    static int    firstCall = 1;
    static int32_t   mask = ATMFAIL | LAND | CHLFAIL;
    static int    i510;
    static double poly_coeff[7];
    static double poly_coeff_412[7] = {0.73172938,-0.54565918,0.20022312,-0.12036241,0.11687968,-0.018732825,-0.0095574674}; 

    double poly_coeff_510[7] = {0.58543905,-0.013125745,-0.059568208,-0.016141849,0.0035106655,0.0012957265,1.4235445e-05};
    double poly_coeff_531[7] = {0.51483158,0.15893415,-0.051696975,-0.055474007,-0.0029635773,0.0053882411,0.0010106600};
    double poly_coeff_551[7] = {0.47507366,0.25216739,-0.0096908094,-0.070882408,-0.012501495,0.0061436085,0.0015798359};
    double poly_coeff_555[7] = {0.43681192,0.26663018,0.016592559,-0.068132662,-0.015470602,0.0051694309,0.0015132129};
    double poly_coeff_412_S[7] = {0.72884183,-0.54380414,0.20225533,-0.12180914,0.11508257,-0.017784535,-0.0095387145}; 
    double poly_coeff_412_M[7] = {0.73461692,-0.54751422,0.19819091,-0.11891567,0.11867679,-0.019681115,-0.0095762203}; 
    float  tLw_pred;
    float  Lt_pred;
    float  Lt_meas;
    float  mu0;
    float  index = 100.0, index510, index412;
    int32_t   ipb, i;
    double logchl, lchl, nLw_510, nLw_412;
    
    
    /* load parameters */
    if (firstCall) {
        
        /* determine index of the bands for the absorbing aerosol analysis */
	ib412  = windex(412.,l2rec->fwave,l2rec->nbands);
	ib510  = windex(510.,l2rec->fwave,l2rec->nbands);
	i510   = l2rec->iwave[ib510];
	
	switch (i510) {
	    case 510: 
		for (i=0; i<7; i++) poly_coeff[i] = poly_coeff_510[i];
		break;
	    case 531: 
		for (i=0; i<7; i++) poly_coeff[i] = poly_coeff_531[i];
		break;
	    case 551:
		for (i=0; i<7; i++) poly_coeff[i] = poly_coeff_551[i];
		break;
	    case 555: 
		for (i=0; i<7; i++) poly_coeff[i] = poly_coeff_555[i];
		break;
	    default:
		for (i=0; i<7; i++) poly_coeff[i] = poly_coeff_510[i];
		break;
        }
/*
	switch (l2rec->sensorID) {
	    case SEAWIFS: 
		for (i=0; i<7; i++) poly_coeff_412[i] = poly_coeff_412_S[i];
		break;
	    case MODISA: 
	    case HMODISA: 
		for (i=0; i<7; i++) poly_coeff_412[i] = poly_coeff_412_M[i];
		break;
	    default:
		break;
        }
*/	
	firstCall = 0;
       
    }
    

    if ((l2rec->flags[ip] & mask) != 0) 
        return(index);
	
    if (l2rec->chl[ip] < 0.1) 
        return(index);


    logchl = lchl = log10(l2rec->chl[ip]);
    nLw_510 = poly_coeff[0]+logchl*poly_coeff[1];
    nLw_412 = poly_coeff_412[0]+logchl*poly_coeff_412[1];
    for (i=0; i<5; i++) {
	 logchl *= lchl;
	 nLw_510 += poly_coeff[i+2]*logchl;
	 nLw_412 += poly_coeff_412[i+2]*logchl;
    }
    
    /* Convert water-leaving radiances from band-centered to band-averaged. */
    /*
    if (l2rec->input->outband_opt >= 2) {
          nLw_510 /= l2rec->f_outband[ip*l2rec->nbands + ib510];
          nLw_412 /= l2rec->f_outband[ip*l2rec->nbands + ib412];
    }
    */


    /* cos of solz */
    mu0 = cos(l2rec->solz[ip]/RADEG);        


    /* bring water-leaving radiance at 510nm to the TOA */
    ipb = ip*l2rec->nbands + ib510;
    tLw_pred = (float)nLw_510 / l2rec->brdf[ipb]
                      * l2rec->t_sol[ipb]
                      * l2rec->t_sen[ipb]
                      * l2rec->tg_sol[ipb]
                      * l2rec->tg_sen[ipb]
              	      * l2rec->polcor[ipb] 
              	      * l2rec->t_o2[ipb] 
 	              * mu0
                      * l2rec->fsol;
		      
    /* calculate TOA predicted radiance  */
    Lt_pred = tLw_pred
              + ( (l2rec->TLg[ipb] 
              +  l2rec->La[ipb] ) 
              *  l2rec->t_o2[ipb] 
              +  l2rec->Lr[ipb] 
              +  l2rec->tLf[ipb])
              *  l2rec->polcor[ipb] 
              *  l2rec->tg_sol[ipb] 
              *  l2rec->tg_sen[ipb];
    
    /* get measured TOA radiance */
    Lt_meas = l2rec->Lt[ipb];
		      
    /* obtain percent difference between measured and predicted TOA radiance */
    index510 = 100.0*(Lt_meas - Lt_pred)/Lt_meas;
    
    if (index510 >= 0.0) 
        return(index);
    

    /* bring water-leaving radiance to the TOA */
    ipb = ip*l2rec->nbands + ib412;
    tLw_pred = (float)nLw_412 / l2rec->brdf[ipb]
                      * l2rec->t_sol[ipb]
                      * l2rec->t_sen[ipb]
                      * l2rec->tg_sol[ipb]
                      * l2rec->tg_sen[ipb]
              	      * l2rec->polcor[ipb] 
              	      * l2rec->t_o2[ipb] 
 	              * mu0
                      * l2rec->fsol;
		      
    /* calculate predicted TOA radiance  */
    Lt_pred = tLw_pred
              + ( (l2rec->TLg[ipb] 
              +  l2rec->La[ipb] ) 
              *  l2rec->t_o2[ipb] 
              +  l2rec->Lr[ipb] 
              +  l2rec->tLf[ipb])
              *  l2rec->polcor[ipb] 
              *  l2rec->tg_sol[ipb] 
              *  l2rec->tg_sen[ipb];
    
    /* get measured TOA radiance */
    Lt_meas = l2rec->Lt[ipb];
		      
    /* obtain percent difference between measured and predicted TOA radiance */
    index412 = 100.0*(Lt_meas - Lt_pred)/Lt_meas + 2.0;

    return (index412);
    
}


