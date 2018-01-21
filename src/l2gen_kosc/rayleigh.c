#include "l12_proto.h"

/* -------------------------------------------------------------------------------
 * rayleigh_iqu() - compute Rayleigh stokes vectors
 *
 * Description:
 *
 *   This computes the Rayleigh scattering radiance (including 
 *   polarization) for a given solar, viewing, and relative azimuthal
 *   angles, as well as the wind speed for the ocean surface roughness
 *   effects.  This is a modification of M. Wangs rayleigh_i subroutine, 
 *   incorporating the code to read the MODIS formatted hdf files.
 *
 * Inputs:
 *
 *   wspeed    ocean surface wind speed in (m/s).
 *   theta0    solar zenith angle (deg).
 *   theta     viewing zenith angle (deg).
 *   phi       relative azimuthal angle (deg).
 *             note: sign of phi matters for U component, sena - 180 - sola assumed
 *
 * Outputs:
 *
 *   ray_i(n)  Rayleigh radiances (F0=1) for n sensor bands.  	
 *   ray_q(n)  Rayleigh radiances (F0=1) for n sensor bands.  	
 *   ray_u(n)  Rayleigh radiances (F0=1) for n sensor bands.  	
 *
 * Hacked by:  B. Franz, December 2002.
 * C Version:  B. Franz, November 2004
 * W. Robinson, SAIC 20 May 2015  prevent solar zonith table index from 
 *               exceeding the table size
 *
 * ------------------------------------------------------------------------------- */

void rayleigh_iqu(int32_t sensorID, int32_t evalmask, int pol_opt,
        float ws, float solz, float senz, float phi,
        float ray_i[], float ray_q[], float ray_u[], int32_t nwave)
{

#define NSOL   45
#define NSEN   41
#define NORDER  3
#define NWIND   8

    static int  firstCall = 1;

    static float *ray_tau;
    static float *ray_dep;
    static float ray_sen  [NSEN];
    static float ray_sol  [NSOL];
    static float ray_sigma[NWIND];
    typedef float ray_array[NWIND][NSOL][NORDER][NSEN];
    static ray_array *ray_for_i;
    static ray_array *ray_for_q;
    static ray_array *ray_for_u;

    static int nsol    = NSOL;
    static int nsen    = NSEN;
    static int norder  = NORDER; 
    static int nwind   = NWIND;

    /* indices */
    int  i,j,k,l,m;
    int iwind;
    int isigma;
    int isigma1;
    int isigma2; 
    int isol1;
    int isol2;
    int isen1;
    int isen2;

    char *tmp_str;
    char  name   [H4_MAX_NC_NAME]  = "";
    char  sdsname[H4_MAX_NC_NAME]  = "";
    char  file   [FILENAME_MAX] = "";
    char  path   [FILENAME_MAX] = "";
    char  wavestr[10]           = "";

    int32 sd_id;
    int32 sds_id; 
    int32 status;
    int32 sds_index;
    int32 rank; 
    int32 nt; 
    int32 dims[H4_MAX_VAR_DIMS]; 
    int32 nattrs;
    int32 start[4] = {0,0,0,0}; 
    int32 edges[4] = {0,0,0,0};

    /* working variables */
    float sigma_m;
    float p, q, h;
    float f00, f10, f01, f11;
    float cosd_phi[NORDER];
    float sind_phi[NORDER];
    float ray_i_sig[2];
    float ray_q_sig[2];
    float ray_u_sig[2];

    int32_t iw;

    setvbuf(stdout, NULL, _IOLBF, 0);
    setvbuf(stderr, NULL, _IOLBF, 0);

    if (firstCall) { 

        int32_t *wave;

        nwave = rdsensorinfo(sensorID, evalmask, "Lambda", (void **) &wave);

        if ((tmp_str = getenv("OCDATAROOT")) == NULL) {
            printf("OCDATAROOT environment variable is not defined.\n");
            exit(1);
        }
        if ((evalmask & NEWRAYTAB) != 0) {
            strcpy(path,tmp_str); strcat(path,"/eval/"); strcat(path,sensorDir[sensorID]); strcat(path,"/");
        } else {
            strcpy(path,tmp_str); strcat(path,"/"); strcat(path,sensorDir[sensorID]); strcat(path,"/");
        }
        printf("\n");
        if ( (ray_for_i = (ray_array *)calloc(nwave,sizeof(ray_array))) == NULL) {
            printf("-E- : Error allocating memory to ray_for_i\n");
            exit(FATAL_ERROR);
        }

        if ( (ray_tau = (float *)calloc(nwave,sizeof(float))) == NULL) {
            printf("-E- : Error allocating memory to ray_for_q\n");
            exit(FATAL_ERROR);
        }
        if ( (ray_dep = (float *)calloc(nwave,sizeof(float))) == NULL) {
            printf("-E- : Error allocating memory to ray_for_q\n");
            exit(FATAL_ERROR);
        }

       if ( (ray_for_q = (ray_array *)calloc(nwave,sizeof(ray_array))) == NULL) {
            printf("-E- : Error allocating memory to ray_for_q\n");
            exit(FATAL_ERROR);
        }
        if ( (ray_for_u = (ray_array *)calloc(nwave,sizeof(ray_array))) == NULL) {
            printf("-E- : Error allocating memory to ray_for_u\n");
            exit(FATAL_ERROR);
        }

        for (iw=0; iw<nwave; iw++) {

            sprintf(wavestr,"%d",(int) wave[iw]);
            strcpy(file,path); strcat(file,"rayleigh/rayleigh_"); 
            strcat(file,sensorDir[sensorID]); 
            strcat(file,"_"); 
            strcat(file,wavestr); 
            strcat(file,"_iqu.hdf");


            /* Open the file and initiate the SD interface */
            sd_id = SDstart(file, DFACC_RDONLY);
            if (sd_id == -1) {
                printf("-E- %s:  Error opening file %s.\n",__FILE__,file);
                exit(1);
            } else
                printf("Loading Rayleigh LUT %s\n",file);

            strcpy(sdsname,"taur");
            sds_id = SDselect(sd_id, SDnametoindex(sd_id,sdsname));
            status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
            status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) ray_tau);
            if (status != 0) {
                printf("-E- %s:  Error reading SDS %s from %s.\n",__FILE__,sdsname,file);
                SDend(sd_id);
                exit(1);
            } else {
                status = SDendaccess(sds_id);
            }

            strcpy(sdsname,"depol");
            sds_id = SDselect(sd_id, SDnametoindex(sd_id,sdsname));
            status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
            status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) ray_dep);
            if (status != 0) {
                printf("-E- %s:  Error reading SDS %s from %s.\n",__FILE__,sdsname,file);
                SDend(sd_id);
                exit(1);
            } else {
                status = SDendaccess(sds_id);
            }

            strcpy(sdsname,"senz");
            sds_id = SDselect(sd_id, SDnametoindex(sd_id,sdsname));
            status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
            status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) ray_sen);
            if (status != 0) {
                printf("-E- %s:  Error reading SDS %s from %s.\n",__FILE__,sdsname,file);
                SDend(sd_id);
                exit(1);
            } else {
                status = SDendaccess(sds_id);
            }

            strcpy(sdsname,"solz");
            sds_id = SDselect(sd_id, SDnametoindex(sd_id,sdsname));
            status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
            status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) ray_sol);
            if (status != 0) {
                printf("-E- %s:  Error reading SDS %s from %s.\n",__FILE__,sdsname,file);
                SDend(sd_id);
                exit(1);
            } else {
                status = SDendaccess(sds_id);
            }

            strcpy(sdsname,"sigma");
            sds_id = SDselect(sd_id, SDnametoindex(sd_id,sdsname));
            status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
            status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) ray_sigma);
            if (status != 0) {
                printf("-E- %s:  Error reading SDS %s from %s.\n",__FILE__,sdsname,file);
                SDend(sd_id);
                exit(1);
            } else {
                status = SDendaccess(sds_id);
            }

            strcpy(sdsname,"i_ray");
            sds_id = SDselect(sd_id, SDnametoindex(sd_id,sdsname));
            status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
            status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) ray_for_i[iw]);
            if (status != 0) {
                printf("-E- %s:  Error reading SDS %s from %s.\n",__FILE__,sdsname,file);
                SDend(sd_id);
                exit(1);
            } else {
                status = SDendaccess(sds_id);
            }

            if (pol_opt > 0) {

                strcpy(sdsname,"q_ray");
                sds_id = SDselect(sd_id, SDnametoindex(sd_id,sdsname));
                status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
                status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) ray_for_q[iw]);
                if (status != 0) {
                    printf("-E- %s:  Error reading SDS %s from %s.\n",__FILE__,sdsname,file);
                    SDend(sd_id);
                    exit(1);
                } else {
                    status = SDendaccess(sds_id);
                }

                strcpy(sdsname,"u_ray");
                sds_id = SDselect(sd_id, SDnametoindex(sd_id,sdsname));
                status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
                status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) ray_for_u[iw]);
                if (status != 0) {
                    printf("-E- %s:  Error reading SDS %s from %s.\n",__FILE__,sdsname,file);
                    SDend(sd_id);
                    exit(1);
                } else {
                    status = SDendaccess(sds_id);
                }
            }
            SDend(sd_id);
        }

        firstCall = 0;
    }


    /* windspeed index into tables */

    sigma_m = 0.0731*sqrt(ws);
    isigma2 = 1;
    while (sigma_m > ray_sigma[isigma2] && isigma2 < nwind)
        isigma2++;
    isigma1 = isigma2 - 1;



    /* interpolate Rayleigh Stokes vectors for each band */

    for (iw=0; iw<nwave; iw++) {

        ray_i[iw] = 0.0;
        ray_q[iw] = 0.0;
        ray_u[iw] = 0.0;

        /* geometry indices into tables */

        if (iw == 0) {

            for (m=0; m<norder; m++) {
                cosd_phi[m] = cos(phi*m/RADEG);
                sind_phi[m] = sin(phi*m/RADEG);
            }

            /* solar zenith indices */

            isol1 = ((int)solz)/2;
            isol2 = isol1 + 1;

            /* sensor zenith indices */

            for (isen2=0; isen2<nsen; isen2++) {
                if ( senz < ray_sen[isen2] )
                    break;
            }
            isen1 = isen2 - 1;

            /* interpolation coefficients */
            /*  for solz > 88 or isol1 >= nsol - 1, use coeffs at end index */
            if( isol1 >= nsol - 1 )  {
              isol1 = nsol - 1;
              isol2 = isol1;
              p = 1.;
            } else 
              p = (solz-ray_sol[isol1])/(ray_sol[isol2]-ray_sol[isol1]);
            q = (senz-ray_sen[isen1])/(ray_sen[isen2]-ray_sen[isen1]);
        }


        /* interpolate the Rayleigh coefficients for each windspeed */

        for (isigma=isigma1; isigma<=isigma2; isigma++) {

            iwind = isigma-isigma1;  /* 0 or 1 */

            ray_i_sig [iwind] = 0.;
            ray_q_sig [iwind] = 0.;
            ray_u_sig [iwind] = 0.;

            /* I component */

            for (m=0; m<norder; m++) {

                f00 = ray_for_i[iw][isigma][isol1][m][isen1];
                f10 = ray_for_i[iw][isigma][isol2][m][isen1];
                f01 = ray_for_i[iw][isigma][isol1][m][isen2];
                f11 = ray_for_i[iw][isigma][isol2][m][isen2];

                ray_i_sig[iwind] = ray_i_sig[iwind] +
                        ((1.-p)*(1.-q)*f00 + p*q*f11 + p*(1.-q)*f10 + q*(1.-p)*f01) * cosd_phi[m];
            }

            if (pol_opt <= 0)
                continue;

            /* Q component */

            for (m=0; m<norder; m++) {

                f00 = ray_for_q[iw][isigma][isol1][m][isen1];
                f10 = ray_for_q[iw][isigma][isol2][m][isen1];
                f01 = ray_for_q[iw][isigma][isol1][m][isen2];
                f11 = ray_for_q[iw][isigma][isol2][m][isen2];

                ray_q_sig[iwind] = ray_q_sig[iwind] +
                        ((1.-p)*(1.-q)*f00 + p*q*f11 + p*(1.-q)*f10 + q*(1.-p)*f01) * cosd_phi[m];
            }


            /* U component */

            for (m=0; m<norder; m++) {

                f00 = ray_for_u[iw][isigma][isol1][m][isen1];
                f10 = ray_for_u[iw][isigma][isol2][m][isen1];
                f01 = ray_for_u[iw][isigma][isol1][m][isen2];
                f11 = ray_for_u[iw][isigma][isol2][m][isen2];

                ray_u_sig[iwind] = ray_u_sig[iwind] +
                        ((1.-p)*(1.-q)*f00 + p*q*f11 + p*(1.-q)*f10 + q*(1.-p)*f01) * sind_phi[m];
            }

        }


        /* do linear interpolation between wind speeds */

        if (isigma1 == isigma2) {

            ray_i[iw] = ray_i_sig[0];

            if (pol_opt > 0) {
                ray_q[iw] = ray_q_sig[0];
                ray_u[iw] = ray_u_sig[0];
            }

        } else {

            h = (sigma_m-ray_sigma[isigma1])/(ray_sigma[isigma2]-ray_sigma[isigma1]);

            ray_i[iw] = ray_i_sig[0] + (ray_i_sig[1] - ray_i_sig[0]) * h;

            if (pol_opt > 0) {
                ray_q[iw] = ray_q_sig[0] + (ray_q_sig[1] - ray_q_sig[0]) * h;
                ray_u[iw] = ray_u_sig[0] + (ray_u_sig[1] - ray_u_sig[0]) * h;
            }
        }

    }
}


/* ----------------------------------------------------------------- */
/* ray_press_wang() - Wang (2004) pressure correction                */
/* ----------------------------------------------------------------- */
float ray_press_wang(float taur, float airmass, float pr)
{
    static float p0=STDPR;
    float x = (-(0.6543-1.608*taur) + (0.8192-1.2541*taur)*log(airmass))*taur*airmass;
    return( (1.0-exp(-x*pr/p0))/(1.0-exp(-x)) );
}


/* ----------------------------------------------------------------- */
/* rayleigh() - rayleigh radiances with polarization, per band.      */
/* ----------------------------------------------------------------- */
void rayleigh(int32_t sensorID, int32_t evalmask, int32_t nwave, int pol_opt,
        float solz, float senz, float raz, float taur[],
        float Fo[], float pr, float ws,
        float Lr_i[], float Lr_q[], float Lr_u[])
{

    float mu0 = cos(solz/RADEG);
    float mu  = cos(senz/RADEG);
    float airmass;
    float fac;
    float *ray_i;
    float *ray_q;
    float *ray_u;
    int32_t  ib;

    if ((ray_i = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("-E- : Error allocating memory to rayleigh\n");
        exit(FATAL_ERROR);
    }
    if ((ray_q = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("-E- : Error allocating memory to rayleigh\n");
        exit(FATAL_ERROR);
    }
    if ((ray_u = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("-E- : Error allocating memory to rayleigh\n");
        exit(FATAL_ERROR);
    }

    /* Get Rayleigh multiple scattering coefficients */
    rayleigh_iqu(sensorID,evalmask,pol_opt,ws,solz,senz,raz,ray_i,ray_q,ray_u,nwave);

    /* Compute Rayleigh radiance for each band.  Include correction   */
    /* for variation in optical depth with atmospheric pressure.      */

    airmass = 1/mu0 + 1/mu;

    for (ib=0; ib<nwave; ib++) {

        fac = ray_press_wang(taur[ib],airmass,pr);

        Lr_i[ib] = Fo[ib]*ray_i[ib]*fac;
        Lr_q[ib] = Fo[ib]*ray_q[ib]*fac;
        Lr_u[ib] = Fo[ib]*ray_u[ib]*fac;
    }
    
    free(ray_i);
    free(ray_q);
    free(ray_u);
}



