/* =========================================================== */
/* Module loadl1.c                                             */
/*                                                             */
/* Functions to fill a level-1b file with precomputed          */
/* atmospheric and masking data.                               */
/*                                                             */
/* Note: due to filtering capabilities of MSl12, it can not be */
/* assumed that the same l1rec pointer will be passed in later */
/* calls to this function, so all fields must be reloaded. In  */
/* addition, it is possible for MSl12 to process a L1B file    */
/* containing data from different time periods, so earth-sun   */
/* distance can not be assumed to be constant.                 */
/*                                                             */
/* Written By:                                                 */
/*                                                             */
/*     B. A. Franz                                             */
/*     SAIC General Sciences Corp.                             */
/*     NASA/SIMBIOS Project                                    */
/*     February 1999                                           */
/*                                                             */
/* Perturbed By:                                               */
/* E.Karakoylu (SAIC)                                          */
/* Summer 2015                                                 */         
/* =========================================================== */

#include <stdio.h>
#include "l12_proto.h"
#include "smi_climatology.h"

// noise option ------------//
#include <gsl/gsl_rng.h>    
#include <gsl/gsl_randist.h>
#include <sys/time.h>
//--------------------------//

unsigned long int random_seed(){
    /* Seed generator for gsl. */
    unsigned int seed;
    struct timeval tv;
    gettimeofday(&tv,0);
    return (tv.tv_sec + tv.tv_usec);
}   


float make_noise(float sigma){
    unsigned long randSeed = random_seed();
    float noise;
    gsl_rng *rng;
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,randSeed);
    noise = gsl_ran_gaussian(rng,sigma);
    gsl_rng_free(rng);
    return (noise);
}

float noise_model_swf(float lt,int32_t iw,float snr_mult){

    /* 
        The present model is of the form SNR = f(Lt)
        following the formula SNR = Lt/sigma.
        sigma is then obtained as  % Lt by 1/SNR.
        !! NEED TO CAP THE FIT.
     */
    
    static int firstpassperband=0;
    float sigma;
    int i,fitPolyOrder = 4;
    float snr=0;
    float ltSnrFitCoef[8][5] = {/*
        {-8.28726301e-07,3.85425664e-04,-9.10776926e-02,1.65881862e+01,4.54351582e-01},
        {-1.21871258e-06,5.21579320e-04,-1.14574109e-01,1.96509056e+01,4.18921861e-01},
        {-2.99068165e-06,1.05225457e-03,-1.90591166e-01,2.66343986e+01,6.67187489e-01},
        {-5.68939986e-06,1.67950509e-03,-2.56915149e-01,3.05832773e+01,9.34468454e-01},
        {-1.31635902e-05,3.09617393e-03,-3.73473556e-01,3.52394751e+01,3.54105899e-01},
        {-8.65458303e-04,1.18857306e-02,-8.37771886e-01,4.64496430e+01,4.14633422e-02},
        {-4.96827099e-04,4.50239057e-02,-2.10425126e-00,7.75862055e+01,5.18893137e-02},
        {-1.30487418e-03,9.35407901e-02,-3.40988182e-00,9.43414239e+01,7.84956550e-01}*/
        {-8.28726301e-03,3.85425664e-01,-9.10776926e+00,1.65881862e+02,4.54351582e-01},
        {-1.21871258e-02,5.21579320e-01,-1.14574109e+01,1.96509056e+02,4.18921861e-01},
        {-2.99068165e-02,1.05225457e-00,-1.90591166e+01,2.66343986e+02,6.67187489e-01},
        {-5.68939986e-02,1.67950509e-00,-2.56915149e+01,3.05832773e+02,9.34468454e-01},
        {-1.31635902e-01,3.09617393e-00,-3.73473556e+01,3.52394751e+02,3.54105899e-01},
        {-8.65458303e-01,1.18857306e+01,-8.37771886e+01,4.64496430e+02,4.14633422e-02},
        {-4.96827099e+00,4.50239057e+01,-2.10425126e+02,7.75862055e+02,5.18893137e-02},
        {-1.30487418e+01,9.35407901e+01,-3.40988182e+02,9.43414239e+02,7.84956550e-01}
    };
    float sigmaLim[8][2] = {   {0.00088,0.059},{0.00078,0.05012},
                                {0.00073,0.03689},{0.00074,0.03199},{0.0008, 0.05593},
                                {0.00108,0.04337},{0.00093,0.0261},{0.00109,0.02123}
    };


        
    for (i=0;i<fitPolyOrder;i++)
        snr += ltSnrFitCoef[iw][i] * pow(lt,(fitPolyOrder - i));
    snr += ltSnrFitCoef[iw][fitPolyOrder];
    snr *= snr_mult;
    sigma  = 1/snr;
    /*if (sigma < sigmaLim[iw][0])
        sigma = sigmaLim[iw][0];
    
    else if (sigma > sigmaLim[iw][1])
        sigma = sigmaLim[iw][1];
    */
    if (firstpassperband==iw){
        printf("Band: %d - Lt=%f - Sigma=%f\n",iw,lt,sigma);
        firstpassperband++;
    }
    return sigma;
}

float get_Lt_noise(float lt, int32_t iw,int32_t sensorID,float snr_fac){

    float sigma;
    float noise;

    switch (sensorID){
        
        case SEAWIFS:
            sigma = noise_model_swf(lt,iw,snr_fac);
            break; 
        
        case HMODISA:
        default:
            printf("-E-: %s line %d: Sensor not  supported.\n",__FILE__,__LINE__);
            exit(1);
            break;
    }    
    
    noise = make_noise(sigma);
    return noise;
}
int loadl1( filehandle *l1file, instr *input, l1str *l1rec, loadl1str *loadl1rec)
{
    double  radeg       = loadl1rec->radeg;
    int32_t sensorID    = loadl1rec->sensorID; 
    int32_t *Lambda     = loadl1rec->Lambda;
    float   *Fobar      = loadl1rec->Fobar;
    float   *Tau_r      = loadl1rec->Tau_r;
    float   *k_oz       = loadl1rec->k_oz;
    float   *aw         = loadl1rec->aw;
    float   *bbw        = loadl1rec->bbw;
    printf("radeg: %lf \n", radeg);
    printf("sensorID: %d \n", sensorID);
    printf("Lambda : %p \n", Lambda);
    printf("Fobar : %p \n", Fobar);
    printf("Tau_r : %p \n", Tau_r);
    printf("k_oz : %p \n", k_oz);
    printf("aw : %p \n", aw);
    printf("bbw : %p \n", bbw);
    
    int     navfail_cnt = 0;

    int32_t ip, ipb, ib, iw, ix;
    double  esdist;
    int32_t nbands      = l1rec->nbands;

    double  *rvs;
    double  temp;
    float   sigma_noise = l1rec->input->add_noise_sigma;                 //
    float   snr_factor  = l1rec->input->noise_scale;
    /* ------------------SINGLE BAND PERTURBATION--------------*/
    int     wiggle_band = l1rec->input->wiggle_band;
    float   wiggle_by   = l1rec->input->wiggle_by;
    /*-----------------------------------------------------*/


    if (sensorID != l1rec->sensorID) {

        loadl1rec->sensorID = l1rec->sensorID;

        if (rdsensorinfo(l1rec->sensorID,l1rec->input->evalmask,NULL,NULL) != l1rec->nbands) {
            printf("-E- %s : Unable to read sensor file\n",__FILE__);
            exit(1);
        }
        rdsensorinfo(l1rec->sensorID,l1rec->input->evalmask,"Lambda",(void **) &(loadl1rec->Lambda));
        rdsensorinfo(l1rec->sensorID,l1rec->input->evalmask,"Fobar", (void **) &(loadl1rec->Fobar));
        rdsensorinfo(l1rec->sensorID,l1rec->input->evalmask,"Tau_r", (void **) &(loadl1rec->Tau_r));
        rdsensorinfo(l1rec->sensorID,l1rec->input->evalmask,"k_oz",  (void **) &(loadl1rec->k_oz));
        rdsensorinfo(l1rec->sensorID,l1rec->input->evalmask,"aw",    (void **) &(loadl1rec->aw));
        rdsensorinfo(l1rec->sensorID,l1rec->input->evalmask,"bbw",   (void **) &(loadl1rec->bbw));

        printf("Loading land mask file from %s\n", input->land);
        if (land_mask_init(input->land) != 0) {
            printf("-E- %s : Unable to initialize land mask\n",__FILE__);
            exit(1);
        }

        printf("Loading bathymetry mask file from %s\n", input->water);
        if (bath_mask_init(input->water) != 0) {
            printf("-E- %s : Unable to initialize bath mask\n",__FILE__);
            exit(1);
        }

        printf("Loading ice mask file from %s\n", input->icefile);
        if (ice_mask_init(input->icefile,(int)(*l1rec->year),
                          (int)(*l1rec->day), input->ice_threshold) != 0) {
            printf("-E- %s : Unable to initialize ice mask\n",__FILE__);
            exit(1);
        }

        if(input->elevfile[0])
            printf("Loading elevation file from %s\n", input->elevfile);
        if(input->elev_auxfile[0])
            printf("Loading auxiliary elevation file from %s\n", input->elev_auxfile);
        elev_init(input->elevfile, input->elev_auxfile);

    }

    /* Get correction for Earth-Sun distance and apply to Fo  */
    esdist = esdist_(l1rec->year,l1rec->day,l1rec->msec);
    l1rec->fsol = pow(1.0/esdist,2);

    for (iw=0; iw<nbands; iw++) {

        l1rec->iwave [iw] = Lambda[iw];
        l1rec->fwave [iw] = (float) Lambda[iw];
        l1rec->Tau_r [iw] = Tau_r [iw];
        l1rec->k_oz  [iw] = k_oz  [iw];
        l1rec->aw    [iw] = aw    [iw]; // to be deleted
        l1rec->bbw   [iw] = bbw   [iw]; // to be deleted
        l1rec->Fobar [iw] = Fobar [iw];
        l1rec->Fo    [iw] = l1rec->Fobar[iw] * l1rec->fsol;

        get_f0_thuillier_ext(l1rec->iwave[iw],BANDW,&l1rec->Fonom[iw]);
    }

    /* Apply vicarious cross-calibration gains */

    for (ix=0; ix<input->xcal_nwave; ix++) {
        if ((input->xcal_opt[ix] & XCALRVS) != 0) {
            if ((ib = bindex_get(input->xcal_wave[ix])) < 0) {
                printf("-E- %sline %d: xcal wavelength %f does not match sensor\n",
                       __FILE__,__LINE__,input->xcal_wave[ix]);
                exit(1);
            };
            rvs = get_xcal(l1rec, XRVS, l1rec->iwave[ib]);
            for (ip=0;ip<l1rec->npix; ip++) {
                ipb  = ip*nbands+ib;
                if (rvs == 0x0) {
                  l1rec->Lt[ipb] = -999.0;
                  continue;
                }
                if (l1rec->Lt[ipb] > 0.0 && l1rec->Lt[ipb] < 1000.0)
                    l1rec->Lt[ipb] /= rvs[ip];
            }
        }
    }
    for (ip=0; ip<l1rec->npix; ip++) {
        /* Apply vicarious calibration */

        for (iw=0; iw<nbands; iw++) {
            ipb = ip*nbands+iw;

            if (l1rec->Lt[ipb] > 0.0 && l1rec->Lt[ipb] < 1000.0) {
                l1rec->Lt[ipb] *= input->gain[iw];
                l1rec->Lt[ipb] += input->offset [iw];

                if (sigma_noise){ // add noise if sigma >0
                    float noise;
                    //noise = make_noise(sigma_noise);
                    snr_factor = input->noise_scale;
                    if (snr_factor < 0)
                        noise = 0;
                    else
                    noise = get_Lt_noise(l1rec->Lt[ipb],iw,sensorID,snr_factor);
                    l1rec->Lt[ipb] *= (1 + noise);
                }

                else if (iw == wiggle_band)
                    l1rec->Lt[ipb] *= (1 + wiggle_by);
            }
        }
        
        /*** Geolocation-based lookups ***/
        if (!l1rec->navfail[ip]) {

            /* Enforce longitude convention */
            if (l1rec->lon[ip] < -180.)
                l1rec->lon[ip] += 360.0;

            else if (l1rec->lon[ip] > 180.0)
                l1rec->lon[ip] -= 360.0;

            /* Get terrain height */
            if (input->proc_land) {
                if (get_height(input->demfile,l1rec,ip,
                               l1file->terrain_corrected) != 0) {
                    printf("-E- %s line %d: Error getting terrain height.\n",
                           __FILE__,__LINE__);
                    return(1);
                }
            } else
                l1rec->height[ip] = 0.0;

            /* Set land, bathymetry and ice flags */
            if (input->proc_ocean != 2 &&
                input->format != FMT_L3BIN &&
                land_mask(l1rec->lat[ip],l1rec->lon[ip]) != 0) {
                l1rec->land[ip] = ON;
            }
            if (!l1rec->land[ip] &&
                bath_mask(l1rec->lat[ip],l1rec->lon[ip]) != 0) {
                l1rec->swater[ip] = ON;
            }
            if (!l1rec->land[ip] &&
                ice_mask(l1rec->lon[ip],l1rec->lat[ip]) != 0) {
                l1rec->ice[ip] = ON;
            }
        } else{
            navfail_cnt++;
        }
        /*** end Geolocation-based lookups ***/

        /* Set sea surface temperature and salinity, and seawater optical properties */
        for (iw=0; iw<nbands; iw++) {
            ipb = ip*nbands+iw;
            l1rec->sw_n [ipb] = 1.334;
            // center band
            l1rec->sw_a [ipb] = aw_spectra (l1rec->fwave[iw],BANDW);
            l1rec->sw_bb[ipb] = bbw_spectra(l1rec->fwave[iw],BANDW);
            // band-averaged
            l1rec->sw_a_avg [ipb] = aw [iw];
            l1rec->sw_bb_avg[ipb] = bbw[iw];
        }

        l1rec->sstref[ip] = BAD_FLT;
        l1rec->sssref[ip] = BAD_FLT;

        if (!l1rec->land[ip]) {
            float bbw_fac;
            l1rec->sstref[ip] = get_sstref(input->sstreftype,input->sstfile,l1rec,ip);
            l1rec->sssref[ip] = get_sssref(input->sssfile,l1rec->lon[ip],l1rec->lat[ip],(int)(*l1rec->day));
            if (l1rec->sstref[ip] > BAD_FLT && l1rec->sssref[ip] > BAD_FLT && input->seawater_opt > 0) {
                for (iw=0; iw<nbands; iw++) {
                    ipb = ip*nbands+iw;
                    l1rec->sw_n [ipb] = seawater_nsw(l1rec->fwave[iw],l1rec->sstref[ip],l1rec->sssref[ip],NULL);
                    // scale bbw based on ratio of center-band model results for actual sea state versus
                    // conditions of Morel measurements used to derive center and band-averaged bbw 
                    bbw_fac = seawater_bb(l1rec->fwave[iw],l1rec->sstref[ip],l1rec->sssref[ip])/seawater_bb (l1rec->fwave[iw],20.0,38.4);
		    l1rec->sw_bb[ipb] *= bbw_fac;
                    l1rec->sw_bb_avg[ipb] *= bbw_fac;
                }
            }
        }
        seawater_set(l1rec);

        /* Compute relative azimuth */
	/* CLASS AVHRR files contain relative azimuth so don't overwrite it */
	if (l1rec->sensorID != AVHRR) {
        l1rec->delphi[ip] = l1rec->sena[ip] - 180.0 - l1rec->sola[ip];
	}
        if (l1rec->delphi[ip] < -180.)
            l1rec->delphi[ip] += 360.0;
        else if (l1rec->delphi[ip] > 180.0)
            l1rec->delphi[ip] -= 360.0;

        /* Precompute frequently used trig relations */
        l1rec->csolz[ip] = cos(l1rec->solz[ip]/radeg);
        l1rec->csenz[ip] = cos(l1rec->senz[ip]/radeg);

        /* Scattering angle */
        temp   = sqrt((1.0-l1rec->csenz[ip]*l1rec->csenz[ip])*(1.0-l1rec->csolz[ip]*l1rec->csolz[ip]))
            * cos(l1rec->delphi[ip]/radeg);
        l1rec->scattang[ip] = acos(MAX(-l1rec->csenz[ip]*l1rec->csolz[ip] + temp,-1.0))*radeg;


    }

    /* add ancillary data */
    // the navfail_cnt test is a kludge to prevent processing failures for scans that are entirely invalid.
    if (navfail_cnt != l1rec->npix){
        if ( setanc(l1rec, input) != 0 )
            return(1);
    }
    if (input->windspeed > -999)
        for (ip=0; ip<l1rec->npix; ip++)
            l1rec->ws[ip] = input->windspeed;
    if (input->windangle > -999)
        for (ip=0; ip<l1rec->npix; ip++)
            l1rec->wd[ip] = input->windangle;
    if (input->pressure > -999)
        for (ip=0; ip<l1rec->npix; ip++)
            l1rec->pr[ip] = input->pressure;
    if (input->ozone > -999)
        for (ip=0; ip<l1rec->npix; ip++)
            l1rec->oz[ip] = input->ozone;
    if (input->watervapor > -999)
        for (ip=0; ip<l1rec->npix; ip++)
            l1rec->wv[ip] = input->watervapor;
    if (input->relhumid > -999)
        for (ip=0; ip<l1rec->npix; ip++)
            l1rec->rh[ip] = input->relhumid;

    /* SWIM bathymetry */
    for (ip=0; ip<l1rec->npix; ip++)
        if (!l1rec->navfail[ip])
            l1rec->elev[ip] = get_elev(l1rec->lat[ip], l1rec->lon[ip]);

    /* add atmospheric cnomponents that do not depend on Lt */
    for (ip=0; ip<l1rec->npix; ip++) {

        /* ------------------------------------------------ */
        /* Ocean processing                                 */
        /* ------------------------------------------------ */
        if ((input->proc_ocean != 0) && !l1rec->land[ip] && !l1rec->navfail[ip]) {

            if (input->ocrvc_opt == 0)
                atmocor1(l1rec,ip);

            /* set the smile_delta field */
            radcor(l1rec,ip,0);

            /* set polarization correction */
            polcor(l1rec,ip);

            /* add surface reflectance */
            get_rhos(l1rec,ip);

            /* assign uncertainty on Lt if not already set by sensor-specific i/o */
            for (iw=0; iw<nbands; iw++) {
                ipb = ip*nbands+iw;
                if (l1rec->Lt_unc[ipb] < BAD_FLT+1) {

                    l1rec->Lt_unc[ipb] = (l1rec->Lt[ipb]-l1rec->Lr[ipb])*input->gain_unc[iw];
                }
            }
        }

        /* ------------------------------------------------ */
        /* Land Processing                                  */
        /* ------------------------------------------------ */
        else if (input->proc_land && l1rec->land[ip] && !l1rec->navfail[ip]) {
            atmocor1_land(input,l1rec,ip);
            radcor(l1rec,ip,1);
            get_rhos(l1rec,ip);
        }

        /* ------------------------------------------------ */
        /* General Processing                               */
        /* ------------------------------------------------ */
        else {
            for (ib=0; ib<nbands; ib++) {
                ipb = ip*nbands+ib;
                l1rec->Lr      [ipb] = 0.0;
                l1rec->t_sol   [ipb] = 1.0;
                l1rec->t_sen   [ipb] = 1.0;
                l1rec->tg_sol[ipb] = 1.0;
                l1rec->tg_sen[ipb] = 1.0;
                l1rec->t_o2    [ipb] = 1.0;
                l1rec->t_h2o   [ipb] = 1.0;
                l1rec->polcor  [ipb] = 1.0;
                l1rec->radcor  [ipb] = 0.0;
            }
        }

    }
    /* set masks and flags */
    if ( setflags(input, l1rec) == 0 )
        return(1);

    setflagbits(1,l1rec,NULL,-1);

    return (0);
}
