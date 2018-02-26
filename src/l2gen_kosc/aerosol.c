/* ======================================================================================== */
/* module aerosol.c  - functions to facilitate aerosol model selection and application      */
/*                                                                                          */
/* Description:                                                                             */
/*                                                                                          */
/* This code replaces the original set of fortran subroutines developed by M.Wang, H.Gordon,*/
/* and others (e.g., rho_a_sub_quad, linear_abc, funct_eps, load_aer, load_ss11) as well as */
/* the original get_aerosol() developed for MSl12.                                          */
/*                                                                                          */
/* The functions herein read and interpolate the aerosol model tables, which are now stored */
/* as individual HDF files per model.  Where sensor wavelengths differ from tabulated model */
/* wavelengths, interpolation is performed. Efficiencies are gained by remembering and      */
/* reusing computational results when applicable.                                           */
/*                                                                                          */
/* Primary Function:                                                                        */
/* ----------------                                                                         */
/* aerosol() - compute aerosol reflectance using specified algorithm (main interface)       */
/*                                                                                          */
/* Secondary Functions:                                                                     */
/* -------------------                                                                      */
/* wangaer() - compute aerosol reflectance using Gordon & Wang 1994 algorithm               */
/* fixedaot() - compute aerosol reflectance for fixed aot(lambda)                           */
/* fixedaer() - compute aerosol reflectance for fixed aerosol model                         */
/* get_angstrom() - compute angstrom coefficient (interface to l2_hdf_generic)              */
/* diff_tran() - compute Rayleigh-aerosol diffuse trans for selected model pair, both paths */
/*                                                                                          */
/* Supporting Functions:                                                                    */
/* --------------------                                                                     */
/* load_aermod() - loads the entire aerosol model table for the specified model list        */
/* ss_to_ms_coef() - for one model, return coefficients of function relating single         */
/*                   scattering to multiple scattering at the input geometry.               */
/* rhoas_to_rhoa() - SS aerosol reflectance to MS aerosol reflectance                       */
/* rhoa_to_rhoas() - MS aerosol reflectance to SS aerosol reflectance                       */
/* model_select_wang() - M. Wang aerosol model selection process                 .          */
/* model_select_angst() - Select model pair based on input Angstrom coefficient             */
/* model_phase() - return phase function at model wavelengths at the input geometry.        */
/* model_epsilon() - return model epsilon at input wavelengths for input geometry.          */
/* model_transmittance() - compute path Rayleigh-aerosol diffuse trans for specified model  */
/* model_taua() - compute AOT at input wavelengths using specified model                    */
/* aeroob - out-of-band water-vapor correction                                              */
/*                                                                                          */
/*                                                                                          */
/* Written By:  B. Franz, SAIC, NASA/GSFC Ocean Biology Processing Group, Summer 2004.      */
/*                                                                                          */
/* ======================================================================================== */

#include <float.h>

#include "l12_proto.h"

#define MAXMODEL MAXAERMOD 
#define MAXSOLZ  33
#define MAXSENZ  35
#define MAXPHI   19
#define MAXSCATT 75
#define DTNTHETA 33

static float pi    = PI;
static float radeg = RADEG;
static float p0    = STDPR;

static int have_ms = 0;
static int have_rh = 0;
static int have_sd = 0;
static int use_rh  = 0;

static int32_t Nbands ;
static int32_t Maxband;                    /* must be >= NBANDS */

typedef struct aermod_struct {

    char  name[32];
    float rh;
    int   sd;

    /* angstrom exponent (nbands+1)*/
    float  *angstrom;

    /* single-scattering albedo(nbands+1), extinction coefficient(nbands+1), phase function */
    float  *albedo;
    float  *extc;
    float  **phase;

    /* quadratic coefficients for SS to MS relationship */
    float  *acost;
    float  *bcost;
    float  *ccost;
  
    /* cubic coefficients for ms_epsilon atmospheric correction ..ZA */
    float  *ams_all;
    float  *bms_all;
    float  *cms_all;
    float  *dms_all;
    float  *ems_all;

  
    /* Rayleigh-aerosol diffuse transmittance coeffs */  
    float  **dtran_a;
    float  **dtran_b;

    /* derived quantities */
    float  **lnphase;
    float  **d2phase;

} aermodstr;

aermodstr* alloc_aermodstr( int nbands, int nscatt, int nphi, int nsolz, int nsenz, int ntheta )
{
    aermodstr  *model;

    model = (aermodstr *) malloc(sizeof(aermodstr));
    model->angstrom = (float*) malloc(nbands*sizeof(float));
    model->albedo   = (float*) malloc(nbands*sizeof(float));
    model->extc     = (float*) malloc(nbands*sizeof(float));
    model->phase    = alloc2d_float( nscatt, nbands );
    model->lnphase  = alloc2d_float( nscatt, nbands );
    model->d2phase  = alloc2d_float( nscatt, nbands );
    model->dtran_a  = alloc2d_float( ntheta, nbands );
    model->dtran_b  = alloc2d_float( ntheta, nbands );
    model->acost    = (float*) malloc( nbands*nsolz*nphi*nsenz*sizeof(float));
    model->bcost    = (float*) malloc( nbands*nsolz*nphi*nsenz*sizeof(float));
    model->ccost    = (float*) malloc( nbands*nsolz*nphi*nsenz*sizeof(float));
    model->ams_all  = (float*) malloc( nbands*nsolz*nphi*nsenz*sizeof(float));
    model->bms_all  = (float*) malloc( nbands*nsolz*nphi*nsenz*sizeof(float));
    model->cms_all  = (float*) malloc( nbands*nsolz*nphi*nsenz*sizeof(float));
    model->dms_all  = (float*) malloc( nbands*nsolz*nphi*nsenz*sizeof(float)); 
    model->ems_all  = (float*) malloc( nbands*nsolz*nphi*nsenz*sizeof(float));   

    return model;
}

static void free_aermodstr( aermodstr *model )
{
    free(model->extc);
    free(model->albedo);
    free(model->angstrom);
    free(model->acost);
    free(model->bcost);
    free(model->ccost);
    free2d_float(model->phase);
    free2d_float(model->lnphase);
    free2d_float(model->d2phase);
    free2d_float(model->dtran_a);
    free2d_float(model->dtran_b);
    free(model->ams_all);
    free(model->bms_all);
    free(model->cms_all);
    free(model->dms_all);
    free(model->ems_all);

    free(model);
}

typedef struct aermodtab_struct {

    int32_t   sensorID;

    /* table dimensions */
    int32_t   nwave;
    int32_t   nmodel;
    int32_t   nsolz;
    int32_t   nsenz;
    int32_t   nphi;
    int32_t   nscatt;
    int32_t   dtran_nwave;
    int32_t   dtran_ntheta;

    /* table spectral bands and angles */
    float  *wave;
    float  *solz;
    float  *senz;
    float  *phi;
    float  *scatt;

    /* diffuse transmittance spectral bands and angles */
    float  *dtran_wave;
    float  *dtran_theta;
    float  *dtran_airmass;

    aermodstr **model;

} aermodtabstr;

static void free_aermodtabstr( aermodtabstr *aertab )
{
    int im;

    for ( im=0; im<aertab->nmodel+1; im++ )
        if ( aertab->model[im] )
            free_aermodstr(aertab->model[im]);

    free(aertab->wave);
    free(aertab->solz);
    free(aertab->senz);
    free(aertab->phi);
    free(aertab->scatt);
    free(aertab->dtran_wave);
    free(aertab->dtran_theta);
    free(aertab->dtran_airmass);
}

typedef struct alphaT_struct { 
    int32_t   modnum; 
    float  angstrom;
} alphaTstr;

typedef struct epsilonT_struct {
    int32_t   modnum;
    float  eps_obs;
} epsilonTstr;

/* aerosol table */
static aermodtabstr *aertab = NULL;

void aerosol_free( void )
{
    if ( aertab )
        free_aermodtabstr(aertab);
}


/* global variable declarations */
static int   loaded = 0; 
static int   interp = 0;
static int32_t  *iwatab;
static int32_t  *iwdtab;

static int32_t  iwnir_s=-1;
static int32_t  iwnir_l=-1;

static int   imm50 = -1;
static int   imc50 = -1;
static int   imc70 = -1;
static int   imt90 = -1;
static int   imt99 = -1;
static int   wang_modx = 0;

static float mu0;
static float mu;
static float airmass;

static int32_t  evalmask = 0;
static int32_t  aer_opt = 0;
static float airmass_plp;
static float airmass_sph;



/* ---------------------------------------------------------------------------------------- */
/* first_deriv() - returns first derivative (dy/dx) of 1st or last array indices using a    */
/*                 4-pt Lagrangian interpolation.  NOTE: It is assumed that 4 points exist. */
/* ---------------------------------------------------------------------------------------- */
float first_deriv(float x[],float y[], int n)
{
    float a1,a2,a3,a4,a5,a6,d1;
    
    if (n == 0) {

        a1 = x[0]-x[1];
        a2 = x[0]-x[2];
        a3 = x[0]-x[3];
        a4 = x[1]-x[2];
        a5 = x[1]-x[3];
        a6 = x[2]-x[3];

        d1 = y[0]*(1.0/a1+1.0/a2+1.0/a3) 
           - a2*a3*y[1]/(a1*a4*a5)
           + a1*a3*y[2]/(a2*a4*a6)
           - a1*a2*y[3]/(a3*a5*a6);

    } else {

        a1 = x[n-1]-x[n-4];
	a2 = x[n-1]-x[n-3];
	a3 = x[n-1]-x[n-2];
	a4 = x[n-2]-x[n-4];
	a5 = x[n-2]-x[n-3];
	a6 = x[n-3]-x[n-4];

        d1 = -a2*a3*y[n-4]/(a6*a4*a1)
           +  a1*a3*y[n-3]/(a6*a5*a2)
           -  a1*a2*y[n-2]/(a4*a5*a3)
 	   +  y[n-1]*(1.0/a1+1.0/a2+1.0/a3);
    }

    return(d1);   
}


/* ---------------------------------------------------------------------------------------- */
/* load_aermod() - loads the entire aerosol model table for the specified model list        */
/* ---------------------------------------------------------------------------------------- */
int load_aermod(int32_t sensorID,float wave[],int32_t nwave, char *aermodfile, char models[MAXAERMOD][32], int32_t nmodels)
{
    char *tmp_str;
    int32 sd_id;
    int32 sds_id; 
    int32 status;
    int32 rank; 
    int32 nt; 
    int32 dims[H4_MAX_VAR_DIMS]; 
    int32 nattrs;
    int32 start[4] = {0,0,0,0}; 

    int32_t   mwave, msolz, msenz, mphi, mscatt, dtran_mwave, dtran_mtheta;

    float d1phase1;
    float d1phaseN;
    float rh;
    int16 sd;

    char  name   [H4_MAX_NC_NAME]  = "";
    char  sdsname[H4_MAX_NC_NAME]  = "";
    char  file   [FILENAME_MAX] = "";
    char  path   [FILENAME_MAX] = "";

    int iw, im, is, iwbase, i;
    static int firstCall=1;

    if (firstCall==1) {
        if ((iwatab = (int32_t *)calloc(nwave,sizeof(int32_t))) == NULL) {
            printf("Unable to allocate space for iwatab.\n");
            exit(1);
        }
        if ((iwdtab = (int32_t *)calloc(nwave,sizeof(int32_t))) == NULL) {
            printf("Unable to allocate space for iwdtab.\n");
            exit(1);
        }
        firstCall = 0;
    }

    printf("Loading aerosol models from %s\n",aermodfile);

    for (im=0; im<nmodels+1; im++) {

        if (im < nmodels) {

            strcpy(file,path); strcat(file,aermodfile); 
            strcat(file,"_"); strcat(file,models[im]); strcat(file,".hdf");

            /* Open the file and initiate the SD interface */
            sd_id = SDstart(file, DFACC_RDONLY);
            if (sd_id == -1) {
                printf("-E- %s:  Error opening file %s.\n",
                   __FILE__,file);
                exit(1);
            }

	} else {

            strcpy(file,path); strcat(file,aermodfile); 
            strcat(file,"_default.hdf");

            /* Open the file and initiate the SD interface */
            sd_id = SDstart(file, DFACC_RDONLY);
            if (sd_id == -1) {
                printf("-E- %s:  Error opening file %s.\n",
                   __FILE__,file);
                exit(1);
            }
	}

       /* read attributes which should be constant between models */

        status = SDreadattr(sd_id,SDfindattr(sd_id,"Number of Wavelengths"),
           &mwave);
        status = SDreadattr(sd_id,SDfindattr(sd_id,
          "Number of Solar Zenith Angles"), &msolz);
        status = SDreadattr(sd_id,SDfindattr(sd_id,
          "Number of View Zenith Angles"), &msenz);
        status = SDreadattr(sd_id,SDfindattr(sd_id,
          "Number of Relative Azimuth Angles"), &mphi);
        status = SDreadattr(sd_id,SDfindattr(sd_id,
          "Number of Scattering Angles"), &mscatt);
        status = SDreadattr(sd_id,SDfindattr(sd_id,
          "Number of Diffuse Transmittance Wavelengths"  ),&dtran_mwave);
        status = SDreadattr(sd_id,SDfindattr(sd_id,
          "Number of Diffuse Transmittance Zenith Angles"),&dtran_mtheta);

        if (im == 0) {

            int32_t   nwave = mwave;
            int32_t   nsolz = msolz;
            int32_t   nsenz = msenz;
            int32_t   nphi = mphi;
            int32_t   nscatt = mscatt;
            int32_t   dtran_nwave = dtran_mwave;
            int32_t   dtran_ntheta = dtran_mtheta;

            printf("Number of Wavelengths                          %d\n", nwave);
            printf("Number of Solar Zenith Angles                  %d\n", nsolz);
            printf("Number of View Zenith Angles                   %d\n", nsenz);
            printf("Number of Relative Azimuth Angles              %d\n", nphi);
            printf("Number of Scattering Angles                    %d\n", nscatt);

            printf("Number of Diffuse Transmittance Wavelengths    %d\n", dtran_nwave);
            printf("Number of Diffuse Transmittance Zenith Angles  %d\n", dtran_ntheta);

            // allocate the aerosol table

            if ((aertab = (aermodtabstr *)calloc(1,sizeof(aermodtabstr))) == NULL) {
                printf("Unable to allocate space for aerosol table.\n");
                exit(1);
            }

            aertab->nmodel = nmodels;
            aertab->nwave  = nwave;
            aertab->nsolz  = nsolz;
            aertab->nsenz  = nsenz;
            aertab->nphi   = nphi;
            aertab->nscatt = nscatt;

            aertab->wave  = (float *) malloc(nwave*sizeof(float));
            aertab->solz  = (float *) malloc(nsolz*sizeof(float));
            aertab->senz  = (float *) malloc(nsenz*sizeof(float));
            aertab->phi   = (float *) malloc(nphi*sizeof(float));
            aertab->scatt = (float *) malloc(nscatt*sizeof(float));

            aertab->dtran_nwave  = dtran_nwave;
            aertab->dtran_ntheta = dtran_ntheta;

            aertab->dtran_wave    = (float *) malloc(dtran_nwave*sizeof(float));
            aertab->dtran_theta   = (float *) malloc(dtran_ntheta*sizeof(float));
            aertab->dtran_airmass = (float *) malloc(dtran_ntheta*sizeof(float));

            // allocate the model tables

            if ((aertab->model = (aermodstr **)calloc(1,(nmodels+1)*sizeof(aermodstr*))) == NULL) {
                printf("Unable to allocate space for %d aerosol models.\n",nmodels+1);
                exit(1);
            }
            for (i=0; i<nmodels+1; i++) {
                if ((aertab->model[i] = alloc_aermodstr( nwave+1, nscatt, MAXPHI, MAXSOLZ, MAXSENZ, dtran_ntheta)) == NULL) {
                    printf("Unable to allocate space for aerosol model %d.\n",im);
                    exit(1);
                }
            }

            /* read SDSes which are constant between models */
      
            strcpy(sdsname,"wave");
            sds_id = SDselect(sd_id, SDnametoindex(sd_id,sdsname));
            status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
            status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->wave);
            if (status != 0) {
                printf("-E- %s:  Error reading SDS %s from %s.\n",__FILE__,sdsname,file);
                exit(1);
            } else {
                status = SDendaccess(sds_id);
	    }

            strcpy(sdsname,"solz");
            sds_id = SDselect(sd_id, SDnametoindex(sd_id,sdsname));
            status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
            status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->solz);
            if (status != 0) {
                printf("-E- %s:  Error reading SDS %s from %s.\n",__FILE__,sdsname,file);
                exit(1);
            } else {
                status = SDendaccess(sds_id);
	    }

            strcpy(sdsname,"senz");
            sds_id = SDselect(sd_id, SDnametoindex(sd_id,sdsname));
            status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
            status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->senz);
            if (status != 0) {
                printf("-E- %s:  Error reading SDS %s from %s.\n",__FILE__,sdsname,file);
                exit(1);
            } else {
                status = SDendaccess(sds_id);
	    }

            strcpy(sdsname,"phi");
            sds_id = SDselect(sd_id, SDnametoindex(sd_id,sdsname));
            status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
            status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->phi);
            if (status != 0) {
                printf("-E- %s:  Error reading SDS %s from %s.\n",__FILE__,sdsname,file);
                exit(1);
            } else {
                status = SDendaccess(sds_id);
	    }

            strcpy(sdsname,"scatt");
            sds_id = SDselect(sd_id, SDnametoindex(sd_id,sdsname));
            status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
            status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->scatt);
            if (status != 0) {
                printf("-E- %s:  Error reading SDS %s from %s.\n",__FILE__,sdsname,file);
                exit(1);
            } else {
                status = SDendaccess(sds_id);
	    }

            strcpy(sdsname,"dtran_wave");
            sds_id = SDselect(sd_id, SDnametoindex(sd_id,sdsname));
            status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
            status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->dtran_wave);
            if (status != 0) {
                printf("-E- %s:  Error reading SDS %s from %s.\n",__FILE__,sdsname,file);
                exit(1);
            } else {
                status = SDendaccess(sds_id);
            }

            strcpy(sdsname,"dtran_theta");
            sds_id = SDselect(sd_id, SDnametoindex(sd_id,sdsname));
            status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
            status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->dtran_theta);
            if (status != 0) {
                printf("-E- %s:  Error reading SDS %s from %s.\n",__FILE__,sdsname,file);
                exit(1);
            } else {
                status = SDendaccess(sds_id);
            }

	} else {
         /*  check that all the aerosol files at least have the same 
             main dimensions  */
          if( ( aertab->nwave != mwave ) || ( aertab->nsolz != msolz ) ||
              ( aertab->nsenz != msenz ) || ( aertab->nphi != mphi ) ||
              ( aertab->nscatt != mscatt ) || 
              ( aertab->dtran_nwave != dtran_mwave ) ||
              ( aertab->dtran_ntheta != dtran_mtheta ) )
            {
            printf("-E- %s, %d:  Error, Aerosol table %s\n", 
              __FILE__, __LINE__, file );
            printf( "    has different dimensions from previous tables\n"  );
            exit(1);
            }
        }

        if (im < nmodels)
	    strncpy(aertab->model[im]->name,models[im],32);
        else
	    strncpy(aertab->model[im]->name,"default",32);

        status = SDreadattr(sd_id,SDfindattr(sd_id,"Relative Humidity"),&rh);
        if (status == 0) {
	    aertab->model[im]->rh = rh;
            have_rh = 1;
	} else
	    aertab->model[im]->rh = -1.0;

        status = SDreadattr(sd_id,SDfindattr(sd_id,"Size Distribution"),&sd);
        if (status == 0) {
	    aertab->model[im]->sd = sd;
            have_sd = 1;
	} else
	    aertab->model[im]->sd = -1;

        strcpy(sdsname,"albedo");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id,sdsname));
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->model[im]->albedo);
        if (status != 0) {
            printf("-E- %s:  Error reading SDS %s from %s.\n",__FILE__,sdsname,file);
            exit(1);
        } else {
            status = SDendaccess(sds_id);
        }

        strcpy(sdsname,"extc");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id,sdsname));
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->model[im]->extc);
        if (status != 0) {
            printf("-E- %s:  Error reading SDS %s from %s.\n",__FILE__,sdsname,file);
            exit(1);
        } else {
            status = SDendaccess(sds_id);
        }

	/* now computed from model parameters
        strcpy(sdsname,"angstrom");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id,sdsname));
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->model[im]->angstrom);
        if (status != 0) {
            printf("-E- %s:  Error reading SDS %s from %s.\n",__FILE__,sdsname,file);
            exit(1);
        } else {
            status = SDendaccess(sds_id);
        }
	*/

        strcpy(sdsname,"phase");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id,sdsname));
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->model[im]->phase[0]);
        if (status != 0) {
            printf("-E- %s:  Error reading SDS %s from %s.\n",__FILE__,sdsname,file);
            exit(1);
        } else {
            status = SDendaccess(sds_id);
        }

        strcpy(sdsname,"acost");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id,sdsname));
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->model[im]->acost);
        if (status != 0) {
            printf("-E- %s:  Error reading SDS %s from %s.\n",__FILE__,sdsname,file);
            exit(1);
        } else {
            status = SDendaccess(sds_id);
        }

        strcpy(sdsname,"bcost");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id,sdsname));
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->model[im]->bcost);
        if (status != 0) {
            printf("-E- %s:  Error reading SDS %s from %s.\n",__FILE__,sdsname,file);
            exit(1);
        } else {
            status = SDendaccess(sds_id);
        }

        strcpy(sdsname,"ccost");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id,sdsname));
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->model[im]->ccost);
        if (status != 0) {
            printf("-E- %s:  Error reading SDS %s from %s.\n",__FILE__,sdsname,file);
            exit(1);
        } else {
            status = SDendaccess(sds_id);
        }

        // check for multi-scattering epsilon tables
        if (SDnametoindex(sd_id,"ams_all") != FAIL) {

            have_ms = 1;

            strcpy(sdsname,"ams_all");
            sds_id = SDselect(sd_id, SDnametoindex(sd_id,sdsname));
            status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
            status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->model[im]->ams_all);
            if (status != 0) {
                printf("-E- %s:  Error reading SDS %s from %s.\n",__FILE__,sdsname,file);
                exit(1);
            } else {
                status = SDendaccess(sds_id);
            }

            strcpy(sdsname,"bms_all");
            sds_id = SDselect(sd_id, SDnametoindex(sd_id,sdsname));
            status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
            status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->model[im]->bms_all);
            if (status != 0) {
                printf("-E- %s:  Error reading SDS %s from %s.\n",__FILE__,sdsname,file);
                exit(1);
            } else {
                status = SDendaccess(sds_id);
            }

            strcpy(sdsname,"cms_all");
            sds_id = SDselect(sd_id, SDnametoindex(sd_id,sdsname));
            status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
            status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->model[im]->cms_all);
            if (status != 0) {
                printf("-E- %s:  Error reading SDS %s from %s.\n",__FILE__,sdsname,file);
                exit(1);
            } else {
                status = SDendaccess(sds_id);
            }

            strcpy(sdsname,"dms_all");
            sds_id = SDselect(sd_id, SDnametoindex(sd_id,sdsname));
            status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
            status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->model[im]->dms_all);
            if (status != 0) {
                printf("-E- %s:  Error reading SDS %s from %s.\n",__FILE__,sdsname,file);
                exit(1);
            } else {
                status = SDendaccess(sds_id);
            }

            strcpy(sdsname,"ems_all");
            sds_id = SDselect(sd_id, SDnametoindex(sd_id,sdsname));
            status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
            status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->model[im]->ems_all);
            if (status != 0) {
                printf("-E- %s:  Error reading SDS %s from %s.\n",__FILE__,sdsname,file);
                exit(1);
            } else {
                status = SDendaccess(sds_id);
            }

	}

        strcpy(sdsname,"dtran_a");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id,sdsname));
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->model[im]->dtran_a[0]);
        if (status != 0) {
            printf("-E- %s:  Error reading SDS %s from %s.\n",__FILE__,sdsname,file);
            exit(1);
        } else {
            status = SDendaccess(sds_id);
        }

        strcpy(sdsname,"dtran_b");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id,sdsname));
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) aertab->model[im]->dtran_b[0]);
        if (status != 0) {
            printf("-E- %s:  Error reading SDS %s from %s.\n",__FILE__,sdsname,file);
            exit(1);
        } else {
            status = SDendaccess(sds_id);
        }

        /* terminate access to the SD interface and close the file */
        status = SDend(sd_id);


        /* compute angstrom exponent for each model wavelength relative to max wavelength */
        iwbase = windex(865,aertab->wave,aertab->nwave);
        for (iw=0; iw<aertab->nwave; iw++) {
	  if (iw != iwbase) 
	    aertab->model[im]->angstrom[iw] = -log(aertab->model[im]->extc[iw]/aertab->model[im]->extc[iwbase])/
                                             log(aertab->wave[iw]/aertab->wave[iwbase]);
	}
        aertab->model[im]->angstrom[iwbase] = aertab->model[im]->angstrom[iwbase-1];

        /* precompute log of phase function and 2nd derivative (for cubic interp) */
        for (iw=0; iw<aertab->nwave; iw++) {
	    for (is=0; is<aertab->nscatt; is++) {
	        aertab->model[im]->lnphase[iw][is] = log(aertab->model[im]->phase[iw][is]);
	    }
            d1phase1 = first_deriv(aertab->scatt,&aertab->model[im]->lnphase[iw][0],0);
            d1phaseN = first_deriv(aertab->scatt,&aertab->model[im]->lnphase[iw][0],aertab->nscatt);
            spline(aertab->scatt,
                   &aertab->model[im]->lnphase[iw][0],
                   aertab->nscatt,
		   d1phase1,
                   d1phaseN,
		   &aertab->model[im]->d2phase[iw][0]);
	}

        /* precompute airmass for diffuse transmittance */
        for (iw=0; iw<aertab->dtran_nwave; iw++) {
	    for (is=0; is<aertab->dtran_ntheta; is++) {
	        aertab->dtran_airmass[is] = 1.0/cos(aertab->dtran_theta[is]/radeg);
	    }
	}
    }
    
    aertab->nmodel   = nmodels;
    aertab->sensorID = sensorID;

    /* map sensor wavelengths to table wavelengths */
    for (iw=0; iw<nwave; iw++) {
        iwatab[iw] = windex(wave[iw],aertab->wave,aertab->nwave);
        iwdtab[iw] = windex(wave[iw],aertab->dtran_wave,aertab->dtran_nwave);
        if ( abs(wave[iw]-aertab->wave[iwatab[iw]]) > 0.5 ) {
            printf("Aerosol model coefficients will be interpolated for %5.1f nm band.\n",wave[iw]);
            interp = 1;
        }
    }

    /* in case of more aertab wavelengths then sensor wavelengths */
    if (aertab->nwave != nwave)
        interp = 1;

    loaded = 1;

    /* set-up for Wang model cross-over correction */
    for (im=0; im<nmodels; im++) {
        if (strcmp(models[im],"m50") == 0) imm50 = im;
        else if (strcmp(models[im],"c50") == 0) imc50 = im;
        else if (strcmp(models[im],"c70") == 0) imc70 = im;
        else if (strcmp(models[im],"t90") == 0) imt90 = im;
        else if (strcmp(models[im],"t99") == 0) imt99 = im;
    }
    if (imm50 >= 0 && imc50 >= 0 && imc70 >= 0 && imt90 >= 0 && imt99 >= 0) {
        wang_modx = 1;
        printf("\nM. Wang model cross-over correction enabled.\n");
    }

    return(0);
}



#define INDEX(iw,isol,iphi,isen) (iw*aertab->nsolz*aertab->nphi*aertab->nsenz + isol*aertab->nphi*aertab->nsenz + iphi*aertab->nsenz + isen)

/* ---------------------------------------------------------------------------------------- */
/* ss_to_ms_coef() - for one model, return coefficients of function relating single         */
/*                   scattering to multiple scattering at the input geometry.               */
/*                                                                                          */
/* This is effectively a C version of M. Wangs linear_a_b_c.f.  The program optimizes for   */
/* multiple calls at the same geometry by computing for all models on the first call with a */
/* new geometry.  It returns pointers to the internal static arrays of coefficients.        */
/*                                                                                          */
/* B. Franz, 1 June 2004.                                                                   */
/* ---------------------------------------------------------------------------------------- */
void ss_to_ms_coef(int modnum, float solz, float senz, float phi, float **a, float **b, float **c)
{
    static float lastsolz = -999.;
    static float lastsenz = -999.;
    static float lastphi  = -999.;

    static int   computed[MAXMODEL];

    static float *a_coef[MAXMODEL];
    static float *b_coef[MAXMODEL];
    static float *c_coef[MAXMODEL];

    static float p, q, r;
    static float as000, as100, as010, as110, as001, as011, as101, as111;
    static float ai000, ai100, ai010, ai110, ai001, ai011, ai101, ai111;
    static float ac000, ac100, ac010, ac110, ac001, ac011, ac101, ac111;

    static int isolz1, isolz2;
    static int isenz1, isenz2;
    static int iphi1 , iphi2;

    float aphi;
    int   im, iw, i;
    static int firstCall=1;

    if (firstCall == 1) {
        firstCall = 0;
        for (i = 0; i < MAXMODEL; i ++) {
            if ((a_coef[i] = (float *)calloc(aertab->nwave,sizeof(float))) == NULL) {
                printf("Unable to allocate space for rhoa.\n");
                exit(1);
            }
            if ((b_coef[i] = (float *)calloc(aertab->nwave,sizeof(float))) == NULL) {
                printf("Unable to allocate space for rhoa.\n");
                exit(1);
            }
            if ((c_coef[i] = (float *)calloc(aertab->nwave,sizeof(float))) == NULL) {
                printf("Unable to allocate space for rhoa.\n");
                exit(1);
            }
        }
    }

    if (solz != lastsolz || senz != lastsenz || phi != lastphi) {

        for (im=0; im<aertab->nmodel; im++)
            computed[im] = 0;

        /* find bracketing solar indices */
        for (i=0; i<aertab->nsolz; i++) {
            if (solz < aertab->solz[i])
	      break;
        }
        isolz1 = MAX(i-1,0);
        isolz2 = MIN(i,aertab->nsolz-1);
        if (isolz2 != isolz1)
            r = (solz-aertab->solz[isolz1])/(aertab->solz[isolz2]-aertab->solz[isolz1]);
	else 
            r = 0.0;

        /* find bracketing view indices */
        for (i=0; i<aertab->nsenz; i++) {
            if (senz < aertab->senz[i])
	      break;
        }
        isenz1 = MAX(i-1,0);
        isenz2 = MIN(i,aertab->nsenz-1);
        if (isenz2 != isenz1)
            p = (senz-aertab->senz[isenz1])/(aertab->senz[isenz2]-aertab->senz[isenz1]);
        else
            p = 0.0;

        /* find bracketing azimuth indices */
        aphi = fabs(phi);
        for (i=0; i<aertab->nphi; i++) {
            if (aphi < aertab->phi[i])
	      break;
        }
        iphi1 = MAX(i-1,0);
        iphi2 = MIN(i,aertab->nphi-1);
        if (iphi2 != iphi1)
            q = (aphi-aertab->phi[iphi1])/(aertab->phi[iphi2]-aertab->phi[iphi1]);
        else
            q = 0.0;

        /* remember last geometry */
        lastsolz = solz;
        lastsenz = senz;
        lastphi  = phi;
    }

    if (!computed[modnum]) {

        im = modnum; computed[modnum] = 1;

        if (isolz2 == 0) {
 
            for (iw=0; iw<aertab->nwave; iw++) {

	      as000 = aertab->model[im]->acost[INDEX(iw,isolz1,0,isenz1)];
	      as100 = aertab->model[im]->acost[INDEX(iw,isolz1,0,isenz2)];
	      as001 = aertab->model[im]->acost[INDEX(iw,isolz2,0,isenz1)];
	      as101 = aertab->model[im]->acost[INDEX(iw,isolz2,0,isenz2)];

	      ai000 = aertab->model[im]->bcost[INDEX(iw,isolz1,0,isenz1)];
	      ai100 = aertab->model[im]->bcost[INDEX(iw,isolz1,0,isenz2)];
	      ai001 = aertab->model[im]->bcost[INDEX(iw,isolz2,0,isenz1)];
	      ai101 = aertab->model[im]->bcost[INDEX(iw,isolz2,0,isenz2)];

	      ac000 = aertab->model[im]->ccost[INDEX(iw,isolz1,0,isenz1)];
	      ac100 = aertab->model[im]->ccost[INDEX(iw,isolz1,0,isenz2)];
	      ac001 = aertab->model[im]->ccost[INDEX(iw,isolz2,0,isenz1)];
	      ac101 = aertab->model[im]->ccost[INDEX(iw,isolz2,0,isenz2)];

	      a_coef[im][iw] = (1.-p)*(1.-r)*as000 + p*r*as101
                             + (1.-p)*r*as001 + p*(1.-r)*as100;

	      b_coef[im][iw] = (1.-p)*(1.-r)*ai000 + p*r*ai101
	                     + (1.-p)*r*ai001 + p*(1.-q)*(1.-r)*ai100;

	      c_coef[im][iw] = (1.-p)*(1.-r)*ac000 + p*r*ac101
	                     + (1.-p)*r*ac001 + p*(1.-q)*(1.-r)*ac100;
	    }

        } else {

            for (iw=0; iw<aertab->nwave; iw++) {

	      as000 = aertab->model[im]->acost[INDEX(iw,isolz1,iphi1,isenz1)];
	      as100 = aertab->model[im]->acost[INDEX(iw,isolz1,iphi1,isenz2)];
	      as010 = aertab->model[im]->acost[INDEX(iw,isolz1,iphi2,isenz1)];
	      as110 = aertab->model[im]->acost[INDEX(iw,isolz1,iphi2,isenz2)];
	      as001 = aertab->model[im]->acost[INDEX(iw,isolz2,iphi1,isenz1)];
	      as011 = aertab->model[im]->acost[INDEX(iw,isolz2,iphi2,isenz1)];
	      as101 = aertab->model[im]->acost[INDEX(iw,isolz2,iphi1,isenz2)];
	      as111 = aertab->model[im]->acost[INDEX(iw,isolz2,iphi2,isenz2)];

	      ai000 = aertab->model[im]->bcost[INDEX(iw,isolz1,iphi1,isenz1)];
	      ai100 = aertab->model[im]->bcost[INDEX(iw,isolz1,iphi1,isenz2)];
	      ai010 = aertab->model[im]->bcost[INDEX(iw,isolz1,iphi2,isenz1)];
	      ai110 = aertab->model[im]->bcost[INDEX(iw,isolz1,iphi2,isenz2)];
	      ai001 = aertab->model[im]->bcost[INDEX(iw,isolz2,iphi1,isenz1)];
	      ai011 = aertab->model[im]->bcost[INDEX(iw,isolz2,iphi2,isenz1)];
	      ai101 = aertab->model[im]->bcost[INDEX(iw,isolz2,iphi1,isenz2)];
	      ai111 = aertab->model[im]->bcost[INDEX(iw,isolz2,iphi2,isenz2)];

	      ac000 = aertab->model[im]->ccost[INDEX(iw,isolz1,iphi1,isenz1)];
	      ac100 = aertab->model[im]->ccost[INDEX(iw,isolz1,iphi1,isenz2)];
	      ac010 = aertab->model[im]->ccost[INDEX(iw,isolz1,iphi2,isenz1)];
	      ac110 = aertab->model[im]->ccost[INDEX(iw,isolz1,iphi2,isenz2)];
	      ac001 = aertab->model[im]->ccost[INDEX(iw,isolz2,iphi1,isenz1)];
	      ac011 = aertab->model[im]->ccost[INDEX(iw,isolz2,iphi2,isenz1)];
	      ac101 = aertab->model[im]->ccost[INDEX(iw,isolz2,iphi1,isenz2)];
	      ac111 = aertab->model[im]->ccost[INDEX(iw,isolz2,iphi2,isenz2)];

	      a_coef[im][iw] = (1.-p)*(1.-q)*(1.-r)*as000 + p*q*r*as111
                             + p*(1.-q)*r*as101 + (1.-p)*q*(1.-r)*as010
                             + p*q*(1.-r)*as110 + (1.-p)*(1.-q)*r*as001
 	                     + (1.-p)*q*r*as011 + p*(1.-q)*(1.-r)*as100;

	      b_coef[im][iw] = (1.-p)*(1.-q)*(1.-r)*ai000 + p*q*r*ai111
                             + p*(1.-q)*r*ai101 + (1.-p)*q*(1.-r)*ai010
                             + p*q*(1.-r)*ai110 + (1.-p)*(1.-q)*r*ai001
	                     + (1.-p)*q*r*ai011 + p*(1.-q)*(1.-r)*ai100;

  	      c_coef[im][iw] = (1.-p)*(1.-q)*(1.-r)*ac000 + p*q*r*ac111
                             + p*(1.-q)*r*ac101 + (1.-p)*q*(1.-r)*ac010
                             + p*q*(1.-r)*ac110 + (1.-p)*(1.-q)*r*ac001
   	                     + (1.-p)*q*r*ac011 + p*(1.-q)*(1.-r)*ac100;
	    }
        }
    }

    /* return pointers to coeffs for this geometry */
    *a = &a_coef[modnum][0];
    *b = &b_coef[modnum][0];
    *c = &c_coef[modnum][0];

    return;
}

/* ---------------------------------------------------------------------------------------- */
/* fresnel_coef() - computes Fresnel reflectance coefficient for specified index of refr.   */
/* ---------------------------------------------------------------------------------------- */
float fresnel_coef(float mu, float index)
{
    float sq, r2, q1;

    sq = sqrt(pow(index,2.0)-1.0+pow(mu,2.0));
    r2 = pow((mu-sq)/(mu+sq),2.0);
    q1 = (1.0-pow(mu,2.0)-mu*sq)/(1.0-pow(mu,2.0)+mu*sq);

    return(r2*(q1*q1+1.0)/2.0);
}



/* ---------------------------------------------------------------------------------------- */
/* ms_eps_coef() -   for a given model, returns  ams, bms, cms, dms and ems coefficients to */
/*                   compute ms_reflectance at the input geometry.                          */
/*                   Also, the program optimizes for multiple calls at the same geometry    */
/*                   by computing for all models on the first call with a new geometry.     */
/*                   It returns pointers to the internal static arrays of coefficients.     */
/*                                                                                          */
/* Z. Ahmad July 08, 2014                                                                   */
/* ---------------------------------------------------------------------------------------- */
void ms_eps_coef(int modnum, int32_t iwnir_l, float wave[], float solz, float senz, float phi, 
                 float **a, float **b, float **c, float **d, float **e)


{
    static float lastsolz = -999.;
    static float lastsenz = -999.;
    static float lastphi  = -999.;

    static int   computed[MAXMODEL];

    static float *a_coef[MAXMODEL];
    static float *b_coef[MAXMODEL];
    static float *c_coef[MAXMODEL];
    static float *d_coef[MAXMODEL];
    static float *e_coef[MAXMODEL];

    static float p, q, r;
    static float as000, as100, as010, as110, as001, as011, as101, as111;
    static float ai000, ai100, ai010, ai110, ai001, ai011, ai101, ai111;
    static float ac000, ac100, ac010, ac110, ac001, ac011, ac101, ac111;
    static float ad000, ad100, ad010, ad110, ad001, ad011, ad101, ad111;
    static float ae000, ae100, ae010, ae110, ae001, ae011, ae101, ae111;


    static int isolz1, isolz2;
    static int isenz1, isenz2;
    static int iphi1 , iphi2;

    float aphi;
    int   im, iw, i;
    int   nwavez;
    static int firstCall = 1;

    if (firstCall == 1) {
        firstCall = 0;
        for (i = 0; i < MAXMODEL; i ++) {
            if ((a_coef[i] = (float *)calloc(aertab->nwave,sizeof(float))) == NULL) {
                printf("Unable to allocate space for a_coef.\n");
                exit(1);
            }
            if ((b_coef[i] = (float *)calloc(aertab->nwave,sizeof(float))) == NULL) {
                printf("Unable to allocate space for b_coef.\n");
                exit(1);
            }
            if ((c_coef[i] = (float *)calloc(aertab->nwave,sizeof(float))) == NULL) {
                printf("Unable to allocate space for c_coef.\n");
                exit(1);
            }
            if ((d_coef[i] = (float *)calloc(aertab->nwave,sizeof(float))) == NULL) {
                printf("Unable to allocate space for d_coef.\n");
                exit(1);
            }

            if ((e_coef[i] = (float *)calloc(aertab->nwave,sizeof(float))) == NULL) {
                printf("Unable to allocate space for e_coef.\n");
                exit(1);
            }

       }
    }

    if (solz != lastsolz || senz != lastsenz || phi != lastphi) {

        for (im=0; im<aertab->nmodel; im++)
            computed[im] = 0;

        /* find bracketing solar indices */
        for (i=0; i<aertab->nsolz; i++) {
            if (solz < aertab->solz[i])
              break;
        }
        isolz1 = MAX(i-1,0);
        isolz2 = MIN(i,aertab->nsolz-1);
        if (isolz2 != isolz1)
            r = (solz-aertab->solz[isolz1])/(aertab->solz[isolz2]-aertab->solz[isolz1]);
        else
            r = 0.0;

        /* find bracketing view indices */
        for (i=0; i<aertab->nsenz; i++) {
            if (senz < aertab->senz[i])
              break;
        }
        isenz1 = MAX(i-1,0);
        isenz2 = MIN(i,aertab->nsenz-1);
        if (isenz2 != isenz1)
            p = (senz-aertab->senz[isenz1])/(aertab->senz[isenz2]-aertab->senz[isenz1]);
        else
            p = 0.0;

        /* find bracketing azimuth indices */
        aphi = fabs(phi);
        for (i=0; i<aertab->nphi; i++) {
            if (aphi < aertab->phi[i])
              break;
        }
        iphi1 = MAX(i-1,0);
        iphi2 = MIN(i,aertab->nphi-1);
        if (iphi2 != iphi1)
            q = (aphi-aertab->phi[iphi1])/(aertab->phi[iphi2]-aertab->phi[iphi1]);
        else
            q = 0.0;

        /* remember last geometry */
        lastsolz = solz;
        lastsenz = senz;
        lastphi  = phi;
    
    }


    if (!computed[modnum]) {

        im = modnum; computed[modnum] = 1;

        if (isolz2 == 0) {

            for (iw=0; iw<aertab->nwave; iw++) {

              as000 = aertab->model[im]->ams_all[INDEX(iw,isolz1,0,isenz1)];
              as100 = aertab->model[im]->ams_all[INDEX(iw,isolz1,0,isenz2)];
              as001 = aertab->model[im]->ams_all[INDEX(iw,isolz2,0,isenz1)];
              as101 = aertab->model[im]->ams_all[INDEX(iw,isolz2,0,isenz2)];

              ai000 = aertab->model[im]->bms_all[INDEX(iw,isolz1,0,isenz1)];
              ai100 = aertab->model[im]->bms_all[INDEX(iw,isolz1,0,isenz2)];
              ai001 = aertab->model[im]->bms_all[INDEX(iw,isolz2,0,isenz1)];
              ai101 = aertab->model[im]->bms_all[INDEX(iw,isolz2,0,isenz2)];

              ac000 = aertab->model[im]->cms_all[INDEX(iw,isolz1,0,isenz1)];
              ac100 = aertab->model[im]->cms_all[INDEX(iw,isolz1,0,isenz2)];
              ac001 = aertab->model[im]->cms_all[INDEX(iw,isolz2,0,isenz1)];
              ac101 = aertab->model[im]->cms_all[INDEX(iw,isolz2,0,isenz2)];

              ad000 = aertab->model[im]->dms_all[INDEX(iw,isolz1,0,isenz1)];
              ad100 = aertab->model[im]->dms_all[INDEX(iw,isolz1,0,isenz2)];
              ad001 = aertab->model[im]->dms_all[INDEX(iw,isolz2,0,isenz1)];
              ad101 = aertab->model[im]->dms_all[INDEX(iw,isolz2,0,isenz2)];

              ae000 = aertab->model[im]->ems_all[INDEX(iw,isolz1,0,isenz1)];
              ae100 = aertab->model[im]->ems_all[INDEX(iw,isolz1,0,isenz2)];
              ae001 = aertab->model[im]->ems_all[INDEX(iw,isolz2,0,isenz1)];
              ae101 = aertab->model[im]->ems_all[INDEX(iw,isolz2,0,isenz2)];


              a_coef[im][iw] = (1.-p)*(1.-r)*as000 + p*r*as101
                             + (1.-p)*r*as001 + p*(1.-r)*as100;

              b_coef[im][iw] = (1.-p)*(1.-r)*ai000 + p*r*ai101
                             + (1.-p)*r*ai001 + p*(1.-q)*(1.-r)*ai100;

              c_coef[im][iw] = (1.-p)*(1.-r)*ac000 + p*r*ac101
                             + (1.-p)*r*ac001 + p*(1.-q)*(1.-r)*ac100;

              d_coef[im][iw] = (1.-p)*(1.-r)*ad000 + p*r*ad101
                             + (1.-p)*r*ad001 + p*(1.-q)*(1.-r)*ad100;

              e_coef[im][iw] = (1.-p)*(1.-r)*ae000 + p*r*ae101
                             + (1.-p)*r*ae001 + p*(1.-q)*(1.-r)*ae100;

            }

        } else {

   /*       printf("coeffs: as000,ai000,ac000,ad000,ae000\n");   */

            for (iw=0; iw<aertab->nwave; iw++) {

              as000 = aertab->model[im]->ams_all[INDEX(iw,isolz1,iphi1,isenz1)];
              as100 = aertab->model[im]->ams_all[INDEX(iw,isolz1,iphi1,isenz2)];
              as010 = aertab->model[im]->ams_all[INDEX(iw,isolz1,iphi2,isenz1)];
              as110 = aertab->model[im]->ams_all[INDEX(iw,isolz1,iphi2,isenz2)];
              as001 = aertab->model[im]->ams_all[INDEX(iw,isolz2,iphi1,isenz1)];
              as011 = aertab->model[im]->ams_all[INDEX(iw,isolz2,iphi2,isenz1)];
              as101 = aertab->model[im]->ams_all[INDEX(iw,isolz2,iphi1,isenz2)];
              as111 = aertab->model[im]->ams_all[INDEX(iw,isolz2,iphi2,isenz2)];

              ai000 = aertab->model[im]->bms_all[INDEX(iw,isolz1,iphi1,isenz1)];
              ai100 = aertab->model[im]->bms_all[INDEX(iw,isolz1,iphi1,isenz2)];
              ai010 = aertab->model[im]->bms_all[INDEX(iw,isolz1,iphi2,isenz1)];
              ai110 = aertab->model[im]->bms_all[INDEX(iw,isolz1,iphi2,isenz2)];
              ai001 = aertab->model[im]->bms_all[INDEX(iw,isolz2,iphi1,isenz1)];
              ai011 = aertab->model[im]->bms_all[INDEX(iw,isolz2,iphi2,isenz1)];
              ai101 = aertab->model[im]->bms_all[INDEX(iw,isolz2,iphi1,isenz2)];
              ai111 = aertab->model[im]->bms_all[INDEX(iw,isolz2,iphi2,isenz2)];

              ac000 = aertab->model[im]->cms_all[INDEX(iw,isolz1,iphi1,isenz1)];
              ac100 = aertab->model[im]->cms_all[INDEX(iw,isolz1,iphi1,isenz2)];
              ac010 = aertab->model[im]->cms_all[INDEX(iw,isolz1,iphi2,isenz1)];
              ac110 = aertab->model[im]->cms_all[INDEX(iw,isolz1,iphi2,isenz2)];
              ac001 = aertab->model[im]->cms_all[INDEX(iw,isolz2,iphi1,isenz1)];
              ac011 = aertab->model[im]->cms_all[INDEX(iw,isolz2,iphi2,isenz1)];
              ac101 = aertab->model[im]->cms_all[INDEX(iw,isolz2,iphi1,isenz2)];
              ac111 = aertab->model[im]->cms_all[INDEX(iw,isolz2,iphi2,isenz2)];

              ad000 = aertab->model[im]->dms_all[INDEX(iw,isolz1,iphi1,isenz1)];
              ad100 = aertab->model[im]->dms_all[INDEX(iw,isolz1,iphi1,isenz2)];
              ad010 = aertab->model[im]->dms_all[INDEX(iw,isolz1,iphi2,isenz1)];
              ad110 = aertab->model[im]->dms_all[INDEX(iw,isolz1,iphi2,isenz2)];
              ad001 = aertab->model[im]->dms_all[INDEX(iw,isolz2,iphi1,isenz1)];
              ad011 = aertab->model[im]->dms_all[INDEX(iw,isolz2,iphi2,isenz1)];
              ad101 = aertab->model[im]->dms_all[INDEX(iw,isolz2,iphi1,isenz2)];
              ad111 = aertab->model[im]->dms_all[INDEX(iw,isolz2,iphi2,isenz2)];

              ae000 = aertab->model[im]->ems_all[INDEX(iw,isolz1,iphi1,isenz1)];
              ae100 = aertab->model[im]->ems_all[INDEX(iw,isolz1,iphi1,isenz2)];
              ae010 = aertab->model[im]->ems_all[INDEX(iw,isolz1,iphi2,isenz1)];
              ae110 = aertab->model[im]->ems_all[INDEX(iw,isolz1,iphi2,isenz2)];
              ae001 = aertab->model[im]->ems_all[INDEX(iw,isolz2,iphi1,isenz1)];
              ae011 = aertab->model[im]->ems_all[INDEX(iw,isolz2,iphi2,isenz1)];
              ae101 = aertab->model[im]->ems_all[INDEX(iw,isolz2,iphi1,isenz2)];
              ae111 = aertab->model[im]->ems_all[INDEX(iw,isolz2,iphi2,isenz2)];


              a_coef[im][iw] = (1.-p)*(1.-q)*(1.-r)*as000 + p*q*r*as111
                             + p*(1.-q)*r*as101 + (1.-p)*q*(1.-r)*as010
                             + p*q*(1.-r)*as110 + (1.-p)*(1.-q)*r*as001
                             + (1.-p)*q*r*as011 + p*(1.-q)*(1.-r)*as100;

              b_coef[im][iw] = (1.-p)*(1.-q)*(1.-r)*ai000 + p*q*r*ai111
                             + p*(1.-q)*r*ai101 + (1.-p)*q*(1.-r)*ai010
                             + p*q*(1.-r)*ai110 + (1.-p)*(1.-q)*r*ai001
                             + (1.-p)*q*r*ai011 + p*(1.-q)*(1.-r)*ai100;

              c_coef[im][iw] = (1.-p)*(1.-q)*(1.-r)*ac000 + p*q*r*ac111
                             + p*(1.-q)*r*ac101 + (1.-p)*q*(1.-r)*ac010
                             + p*q*(1.-r)*ac110 + (1.-p)*(1.-q)*r*ac001
                             + (1.-p)*q*r*ac011 + p*(1.-q)*(1.-r)*ac100;

              d_coef[im][iw] = (1.-p)*(1.-q)*(1.-r)*ad000 + p*q*r*ad111
                             + p*(1.-q)*r*ad101 + (1.-p)*q*(1.-r)*ad010
                             + p*q*(1.-r)*ad110 + (1.-p)*(1.-q)*r*ad001
                             + (1.-p)*q*r*ad011 + p*(1.-q)*(1.-r)*ad100;

              e_coef[im][iw] = (1.-p)*(1.-q)*(1.-r)*ae000 + p*q*r*ae111
                             + p*(1.-q)*r*ae101 + (1.-p)*q*(1.-r)*ae010
                             + p*q*(1.-r)*ae110 + (1.-p)*(1.-q)*r*ae001
                             + (1.-p)*q*r*ae011 + p*(1.-q)*(1.-r)*ae100;


            }
        }
    }

    // root finding is quadratic, but table includes cubic term, make sure it's zero
    if (fabs(d_coef[modnum][iwnir_l]) > 1e-9) {
        printf("non zero cubic term found in longest NIR wavelength of aerosol table. Zia!!\n");
        exit(1);
    }

    /* return pointers to coeffs for this geometry */
    *a = &a_coef[modnum][0];    
    *b = &b_coef[modnum][0];    
    *c = &c_coef[modnum][0];    
    *d = &d_coef[modnum][0];    
    *e = &e_coef[modnum][0];

    return;
}


/* ---------------------------------------------------------------------------------------- */
/* model_select_ahmad() - select two aerosol models whose epsilon values bracket the        */
/*                        the observed ms epsilon, eps_obs                                  */
/*                                                                                          */
/* Z Ahmad July 2014.                                                                       */
/* ---------------------------------------------------------------------------------------- */
int comp_epsilonT(epsilonTstr *x, epsilonTstr *y) {return(x->eps_obs < y->eps_obs ? -1 : 1);}

void  model_select_ahmad(int32_t nmodels, int32_t *mindx, float eps_pred[], float eps_obs, int32_t *modmin,
                         int32_t *modmax, float *modrat)
{
    static epsilonTstr  epsilonT[MAXAERMOD];

    int im, im1, im2;
    static float weight;

    // set-up table to keep model epsilon and model number pairs together

    for (im=0; im < nmodels; im++) {
         epsilonT[im].eps_obs = eps_pred[im];
         epsilonT[im].modnum = im;
    /*     printf("%d %7.4f %7.4f\n",im,eps_pred[im],eps_obs);              */
    }

    // sort into ascending model epsilon order

    qsort(epsilonT,nmodels,sizeof(epsilonTstr),
         (int (*)(const void *,const void *)) comp_epsilonT);

    // find bounding epsilon indices in table

      for (im=0; im < nmodels; im++) {                                     
           if (eps_obs < epsilonT[im].eps_obs)                            
              break;                                                     
      }                                                                 

      im1 = MAX(MIN(im-1,nmodels-1),0);                                
      im2 = MAX(MIN(im  ,nmodels-1),0);                               

    // convert table indices to model indices of the input order             
    // compute model weighting

      *modmin = epsilonT[im1].modnum;                                
      *modmax = epsilonT[im2].modnum;                               
      *modrat =(eps_obs-epsilonT[im1].eps_obs)/(epsilonT[im2].eps_obs-epsilonT[im1].eps_obs);   
  /*    printf("%d %d %7.4f %7.4f %7.4f\n",im1,im2,eps_obs,epsilonT[im1].eps_obs,epsilonT[im2].eps_obs); */


    return;
}


/*------------------------------------------------------------------------------------------ */
/* comp_rhoa_ms_eps() -  for a given model, computes the AOT and aerosol reflectance         */
/*                                                                                           */
/* Z ahmad July2014                                                                          */
/*------------------------------------------------------------------------------------------ */
int comp_rhoa_ms_eps( int32_t nwave, float wave[], float solz, float senz, float phi,
                      float tau_iwnir_l, int32_t modl, float tau_pred[], float rho_pred[])

{
    float *ac,*bc,*cc,*dc,*ec;
    float ext_modl[nwave];
    float lg_tau_pred[nwave];
    float lg_rho_pred[nwave];
    int iw;

    /* get the coefficients for lg_rho vs lg_aot  */

     ms_eps_coef(modl,iwnir_l,wave,solz,senz,phi,&ac,&bc,&cc,&dc,&ec);

    /* get the extinction coefficients and compute AOT at all wavelengths */

    for (iw=0; iw<= nwave; iw++) {
        ext_modl[iw] = aertab->model[modl]->extc[iw];
    }

 /*   printf("tau_pred[iw],tau_iwnir_l\n");    */
    for (iw=0; iw<= nwave; iw++) {
        tau_pred[iw]=(ext_modl[iw]/ext_modl[iwnir_l])*tau_iwnir_l;
        lg_tau_pred[iw]=log(tau_pred[iw]);
    }

    /* compute rho_pred */

    for (iw=0; iw<=nwave; iw++) {
        lg_rho_pred[iw]=ac[iw] + 
                        bc[iw]*lg_tau_pred[iw] + 
                        cc[iw]*pow(lg_tau_pred[iw],2) +
                        dc[iw]*pow(lg_tau_pred[iw],3) +
                        ec[iw]*pow(lg_tau_pred[iw],4);
        rho_pred[iw]=exp(lg_rho_pred[iw]);
    }


    return(0);
}

/*------------------------------------------------------------------------------------------*/
/* ahmad_atm_corr() - compute aerosol reflectance at all wavelengths                        */
/*                                                                                          */
/* Z Ahmad. July 2014                                                                       */
/*----------------------------------------------------------------------------------------  */
int ahmad_atm_corr(int32_t sensorID, float wave[], int32_t nwave, int32_t iwnir_s, int32_t iwnir_l,
                   int32_t nmodels, int32_t mindx[],
                   float solz, float senz, float phi, float wv, float rhoa[],
                   int32_t *modmin, int32_t *modmax, float *modrat, float *epsnir,
                   float tau_pred_max[], float tau_pred_min[], float rho_pred_max[], float rho_pred_min[],
                   float tau_aer[], float rho_aer[] )
{


    float *ac,*bc,*cc,*dc,*ec;
    float ext_iwnir_l, ext_iwnir_s;
    double ax,bx,cx,dx,fx;
    float tau_iwnir_l[nmodels];
    float lg_tau_iwnir_s[nmodels],tau_iwnir_s[nmodels];
    float lg_rho_iwnir_s_pred[nmodels],rho_iwnir_s_pred[nmodels];
    float eps_pred[nmodels];
    int   iwtab,iwtab2,im,modl;
    int   im1, im2;
    float mwt;
    float lg_tau_iwnir_l;
    

    static float eps_obs;

    int iw;
    int km;
    
    int status = 0.0;


    /* compute the observed epsilon */

    eps_obs=rhoa[iwnir_s]/rhoa[iwnir_l];

/*             printf("rho_869,rho_748,eps_obs\n");                    */
/*    printf("%10.5f %10.5f %10.5f\n",rhoa[iwnir_l],rhoa[iwnir_s],eps_obs);   */


    // compute MS epsilon for all nmodels.  note that nmodels may be only a subset of
    // the full model suite.  mindx maps from subset index to full index in suite


    for (im=0; im < nmodels; im++) {

        modl = mindx[im];  // index in full model suite, needed for coefficient look-up

        /* compute AOT at longest aerosol wavelength (iwnir_l) */

        ms_eps_coef(modl,iwnir_l,wave,solz,senz,phi,&ac,&bc,&cc,&dc,&ec);
  
        ax = (double) ac[iwnir_l]-log((double)rhoa[iwnir_l]) ;
        bx = (double) bc[iwnir_l];
        cx = (double) cc[iwnir_l];

        fx = bx*bx-4.0*ax*cx;
        if (fx > 0.0 && cx != 0.0) {
            lg_tau_iwnir_l = 0.5*(-bx+sqrt(fx))/cx;
            tau_iwnir_l[im] = exp(lg_tau_iwnir_l);
        }  else  {
            status = 1;
            break;
        }
  
        /* compute AOT at shortest aerosol wavelength (iwnir_s) */

        ext_iwnir_l=aertab->model[modl]->extc[iwnir_l];
        ext_iwnir_s=aertab->model[modl]->extc[iwnir_s];

        tau_iwnir_s[im]=(ext_iwnir_s/ext_iwnir_l)*tau_iwnir_l[im];
        lg_tau_iwnir_s[im]=log(tau_iwnir_s[im]);


        /* compute reflectance at (iwnir_s) */

        lg_rho_iwnir_s_pred[im]=ac[iwnir_s] + bc[iwnir_s]*lg_tau_iwnir_s[im]  +
                                cc[iwnir_s]*pow(lg_tau_iwnir_s[im],2)  +
                                dc[iwnir_s]*pow(lg_tau_iwnir_s[im],3) +
                                ec[iwnir_s]*pow(lg_tau_iwnir_s[im],4);

        rho_iwnir_s_pred[im]=exp(lg_rho_iwnir_s_pred[im]);

        /* compute model epsilon */

        eps_pred[im]=rho_iwnir_s_pred[im]/rhoa[iwnir_l];

    }


    *epsnir = eps_obs;


    /* now, find the bounding models, but skip this if we only have one */
    /* model (as when vicariously calibrating)                          */

    if (nmodels > 1) {

        /* locate two model_epsilons that bracket the observed epsilon */
        /* this will return the model numbers for the subset of models */

        model_select_ahmad(nmodels,mindx,eps_pred,eps_obs,&im1,&im2,&mwt);

    } else {

        /* single model case (allows fixed model by limiting model suite) */

        im1 = 0;
        im2 = 0;
        mwt = 0.0;
    }

    /* map selection to full suite for subsequent processing and return */

    *modmin = mindx[im1];
    *modmax = mindx[im2];
    *modrat = mwt;

    /* compute tau_pred and rho_predicted */

    comp_rhoa_ms_eps(nwave,wave,solz,senz,phi,tau_iwnir_l[im1],*modmin,tau_pred_min,rho_pred_min);
    comp_rhoa_ms_eps(nwave,wave,solz,senz,phi,tau_iwnir_l[im2],*modmax,tau_pred_max,rho_pred_max);

    /* compute weighted tau_aer and rho_aer */

    for (iw=0; iw < nwave; iw++) {
       tau_aer[iw]=(1.0-mwt)*tau_pred_min[iw] + mwt*tau_pred_max[iw];
       rho_aer[iw]=(1.0-mwt)*rho_pred_min[iw] + mwt*rho_pred_max[iw];
    }


    return (status);
}


/* ---------------------------------------------------------------------------------------- */
/* model_phase() - return phase function at model wavelengths at the input geometry.        */
/*                                                                                          */
/* This is effectively a C version of load_ss.f by M. Wang.  The program optimizes for      */
/* multiple calls at the same geometry by computing for all models on the first call with a */
/* new geometry.  It returns a pointer to the internal static array for the requested model.*/
/*                                                                                          */
/* B. Franz, 1 June 2004.                                                                   */
/* ---------------------------------------------------------------------------------------- */
float *model_phase(int modnum, float solz, float senz, float phi)
{
    static float nw = 1.334;

    static int   computed[MAXMODEL];
    static float lastsolz = -999.;
    static float lastsenz = -999.;
    static float lastphi  = -999.;
    static float *phase[MAXMODEL];
    static float fres1,fres2;
    static float scatt1,scatt2;
    static int   firstCall = 1;

    float phase1, phase2;
    int   im, iw;

    if (firstCall == 1) {
        firstCall = 0;

        for (im = 0; im<MAXMODEL; im++) {
            if ((phase[im] = (float *)calloc(aertab->nwave,sizeof(float))) == NULL) {
                printf("Unable to allocate space for phase.\n");
                exit(1);
            }
        }

    }

    /* recalculate only if geometry changes */

    if (solz != lastsolz || senz != lastsenz || phi != lastphi) {

        float csolz = cos(solz/radeg);    
        float csenz = cos(senz/radeg);    
        float cphi  = cos(phi /radeg);    
        float temp;

        for (im=0; im<aertab->nmodel; im++)
            computed[im] = 0;

        /* determine scattering angles (direct and surface reflected) */
        temp   = sqrt((1.0-csenz*csenz)*(1.0-csolz*csolz)) * cphi;
        scatt1 = acos(MAX(-csenz*csolz + temp,-1.0))*radeg;
        scatt2 = acos(MIN( csenz*csolz + temp, 1.0))*radeg;

        /* compute Fresnel coefficients */
        fres1 = fresnel_coef(csenz,nw);
        fres2 = fresnel_coef(csolz,nw);

        lastsolz = solz;
        lastsenz = senz;
        lastphi  = phi;
    }

    if (!computed[modnum]) {

        im = modnum; computed[modnum] = 1;

        /* compute phase function for this geometry, all models */
        for (iw=0; iw<aertab->nwave; iw++) {
            splint(aertab->scatt,
                   &aertab->model[im]->lnphase[iw][0],
                   &aertab->model[im]->d2phase[iw][0],
                   aertab->nscatt,scatt1,&phase1);
            splint(aertab->scatt,
                   &aertab->model[im]->lnphase[iw][0],
                   &aertab->model[im]->d2phase[iw][0],
                   aertab->nscatt,scatt2,&phase2);
	    //          incident diffuse   reflected   diff  dir
            phase[im][iw] = exp(phase1) + exp(phase2)*(fres1+fres2);
	}
    }

    return(&phase[modnum][0]);
}


/* ---------------------------------------------------------------------------------------- */
/* aeroob_cf() - out-of-band water-vapor scale factor                                       */
/* ---------------------------------------------------------------------------------------- */
float aeroob_cf(int modnum,float solz, float senz, float phi)
{
    static int firstCall = 1;
    static int iw1;
    static int iw2;

    float *phase;  
    float  rhoas1, rhoas2;
    float  eps;
    float  cf;

    if (firstCall) {
        iw1 = windex(765,aertab->wave,aertab->nwave);
        iw2 = windex(865,aertab->wave,aertab->nwave);
        if (iw1 == iw2) iw1--;
        firstCall = 0;
    } 

    phase  = model_phase(modnum,solz,senz,phi);
    rhoas1 = aertab->model[modnum]->albedo[iw1] * phase[iw1] * aertab->model[modnum]->extc[iw1];
    rhoas2 = aertab->model[modnum]->albedo[iw2] * phase[iw2] * aertab->model[modnum]->extc[iw2];
    eps    = rhoas1/rhoas2;
    cf     = log(eps)/(aertab->wave[iw2]-aertab->wave[iw1]);

    return(cf);
}


/* ---------------------------------------------------------------------------------------- */
/* aeroob - out-of-band water-vapor correction                                              */
/* ---------------------------------------------------------------------------------------- */
float aeroob(int32_t sensorID, int32_t iw, float cf, float wv)
{
    static float *a01;
    static float *a02;
    static float *a03;
    static float *a04;
    static float *a05;
    static float *a06;
    static float *a07;
    static float *a08;
    static float *a09;
    static float *a10;
    static float *a11;
    static float *a12;

    static int   firstCall = 1;

    float  f;

    if (firstCall) {
        firstCall = 0;
        rdsensorinfo(sensorID,evalmask,"oobwv01",(void **) &a01);   /* coeff #1 per sensor wave */
        rdsensorinfo(sensorID,evalmask,"oobwv02",(void **) &a02);
        rdsensorinfo(sensorID,evalmask,"oobwv03",(void **) &a03);
        rdsensorinfo(sensorID,evalmask,"oobwv04",(void **) &a04);
        rdsensorinfo(sensorID,evalmask,"oobwv05",(void **) &a05);
        rdsensorinfo(sensorID,evalmask,"oobwv06",(void **) &a06);
        rdsensorinfo(sensorID,evalmask,"oobwv07",(void **) &a07);
        rdsensorinfo(sensorID,evalmask,"oobwv08",(void **) &a08);
        rdsensorinfo(sensorID,evalmask,"oobwv09",(void **) &a09);
        rdsensorinfo(sensorID,evalmask,"oobwv10",(void **) &a10);
        rdsensorinfo(sensorID,evalmask,"oobwv11",(void **) &a11);
        rdsensorinfo(sensorID,evalmask,"oobwv12",(void **) &a12);
        printf("\nLoading water-vapor correction coefficients.\n");
    }

    f = (a01[iw] + a02[iw]*airmass + cf*(a03[iw] + a04[iw]*airmass))
      + (a05[iw] + a06[iw]*airmass + cf*(a07[iw] + a08[iw]*airmass))*wv
      + (a09[iw] + a10[iw]*airmass + cf*(a11[iw] + a12[iw]*airmass))*wv*wv;

    return(f);
}


/* ---------------------------------------------------------------------------------------- */
/* rhoa_to_rhoas() - MS aerosol reflectance to SS aerosol reflectance                       */
/*                                                                                          */
/* B. Franz, 1 June 2004.                                                                   */
/* ---------------------------------------------------------------------------------------- */
int rhoa_to_rhoas(int32_t sensorID, int modnum, float solz, float senz, float phi, float wv,
                  float rhoa[], float wave[], int32_t nwave, int iw1, int iw2, float rhoas[])
{
    float *ac,*bc,*cc;
    double a,b,c;
    double f;
    int   iw;
    int   iwtab;
    float cf;
    int status = 0;

    ss_to_ms_coef(modnum,solz,senz,phi,&ac,&bc,&cc);

    cf = aeroob_cf(modnum,solz,senz,phi);

    for (iw=iw1; iw<=iw2; iw++) {
      if (rhoa[iw] < 1.e-20)
        rhoas[iw] = rhoa[iw];
      else {
        iwtab = iwatab[iw];
        a = (double) ac[iwtab];
        b = (double) bc[iwtab];
        c = (double) cc[iwtab];
        f = b*b - 4*c*( a - log((double)rhoa[iw]) );
	if (f > 0.00001) {     // this was 0.0, but small values caused segfault (BAF, 9/2014)
          if (fabs(c) > 1.e-20) {
            rhoas[iw] = exp(0.5*(-b+sqrt(f))/c);
          } else if (fabs(a) > 1.e-20 && fabs(b) > 1.e-20) {
 	    rhoas[iw] = pow(rhoa[iw]/a, 1./b);
          } else {
            status = 1;
            break;
	  }
          rhoas[iw] = rhoas[iw]/aeroob(sensorID,iw,cf,wv);
          if (!finite(rhoas[iw]) || rhoas[iw] < 1.e-20) {
            status=1;
            break;
	  }
	} else {
          status=1;
          break;
	}
      }
    }

    // return input values and failure status if any wavelengths failed

    if (status != 0) {
      for (iw=iw1; iw<=iw2; iw++)
        rhoas[iw] = rhoa[iw];
    }

    return(status);
}



/* ---------------------------------------------------------------------------------------- */
/* rhoas_to_rhoa() - SS aerosol reflectance to MS aerosol reflectance                       */
/*                                                                                          */
/* B. Franz, 1 June 2004.                                                                   */
/* ---------------------------------------------------------------------------------------- */
void rhoas_to_rhoa(int32_t sensorID, int modnum, float solz, float senz, float phi, float wv, 
                  float rhoas[], float wave[], int32_t nwave, int iw1, int iw2, float rhoa[])
{
    float *ac,*bc,*cc;
    float a,b,c;
    float lnrhoas;
    int   iw;
    int   iwtab;
    float cf;

    ss_to_ms_coef(modnum,solz,senz,phi,&ac,&bc,&cc);

    cf = aeroob_cf(modnum,solz,senz,phi);

    for (iw=iw1; iw<=iw2; iw++) {


	/* these changes ensure that rhoa1 -> rhoas -> rhoa2 == rhoa1 */
        /* but, that changes everything, slightly (tau)               */

        if (rhoas[iw] < 1.e-20)
             rhoa[iw] = rhoas[iw];
        else {
             iwtab = iwatab[iw];
             a = ac[iwtab];
             b = bc[iwtab];
             c = cc[iwtab];
             lnrhoas = log(rhoas[iw]*aeroob(sensorID,iw,cf,wv));
             rhoa[iw] = exp(a + b*lnrhoas + c*lnrhoas*lnrhoas);
        }

	/*
          iwtab = iwatab[iw];
          a = ac[iwtab];
          b = bc[iwtab];
          c = cc[iwtab];
          lnrhoas = log(rhoas[iw]);
 	  rhoa[iw] = exp(a + b*lnrhoas + c*lnrhoas*lnrhoas) 
                   * aeroob(sensorID,iw,cf,wv);
	*/
    }

    return;
}


/* ---------------------------------------------------------------------------------------- */
/* model_epsilon() - return model epsilon at input wavelengths for input geometry.          */
/*                                                                                          */
/* If the input wavelengths are not equal to the model wavelengths, the model epsilon will  */
/* be interpolated to the input wavelengths.  It is assumed that the longest (last) input   */
/* wavelength is equivalent to the longest wavelength of the model table.  Hence,           */
/* the function should always be called with the full array of input sensor wavelengths.    */
/*                                                                                          */
/* The program optimizes for multiple calls at the same geometry by computing for all       */
/* models on the first call with a new geometry.  It returns a pointer to the internal      */
/* static arrays of epsilon for the requested model.                               .        */
/*                                                                                          */
/* B. Franz, 1 June 2004.                                                                   */
/* ---------------------------------------------------------------------------------------- */
float *model_epsilon(int modnum, int32_t iwnir_l, float wave[], int32_t nwave, float solz, float senz, float phi)
{
    static float lastsolz = -999.;
    static float lastsenz = -999.;
    static float lastphi  = -999.;
    static int32_t  lastiwl  = -999;
    static int   lastmod  = -999;
    static float *epsilon[MAXMODEL];
    static int firstCall = 1;
    int i;
    float maxwave;
    /* recalculate only if geometry changes */

    if (firstCall == 1) {
        firstCall = 0;
       maxwave=MAX(aertab->nwave,nwave);
        for (i=0; i<MAXMODEL; i++) {
            if ((epsilon[i] = (float *)calloc(maxwave,sizeof(float))) == NULL) {
                printf("Unable to allocate space for epsilon.\n");
                exit(1);
            }

        }
    }

    /* recalculate only if geometry changes */

    if (modnum != lastmod || solz != lastsolz || senz != lastsenz || phi != lastphi || iwnir_l != lastiwl) {

        int    iwnir  = iwatab[iwnir_l];
        float *phase;  
        float  *lneps;
        float  rhoas1, rhoas2;
        int    im, iw, iwtab;

        if ((lneps = (float *)calloc(aertab->nwave,sizeof(float))) == NULL) {
            printf("Unable to allocate space for lneps.\n");
            exit(1);
        }

        im = modnum;
        phase = model_phase(im,solz,senz,phi);
        for (iw=0; iw<aertab->nwave; iw++) {
            rhoas1 = aertab->model[im]->albedo[iw]    * phase[iw]    * aertab->model[im]->extc[iw];
            rhoas2 = aertab->model[im]->albedo[iwnir] * phase[iwnir] * aertab->model[im]->extc[iwnir];
            epsilon[im][iw] = rhoas1/rhoas2;
            if (interp)
                lneps[iw] = log(epsilon[im][iw]);
        }
        if (interp) {
            for (iw=0; iw<nwave; iw++) {
                iwtab = iwatab[iw];
                if (aertab->wave[iwtab] != wave[iw] && wave[iw] > 0)
                    epsilon[im][iw] = exp(linterp(aertab->wave,lneps,aertab->nwave,wave[iw]));
                else
                    epsilon[im][iw] = exp(lneps[iwtab]);
            }
        }

        lastsolz = solz;
        lastsenz = senz;
        lastphi  = phi;
        lastiwl  = iwnir_l;
        lastmod  = modnum;
        free(lneps);

    }

     return(&epsilon[modnum][0]);
}


/* ---------------------------------------------------------------------------------------- */
/* model_select_wang() - M. Wang aerosol model selection process                 .          */
/*                                                                                          */
/* B. Franz, 1 June 2004.                                                                   */
/* ---------------------------------------------------------------------------------------- */
int model_select_wang(int32_t sensorID, float wave[], int32_t nwave, int32_t nmodel, int32_t mindx[],
                      float solz, float senz, float phi,
                      float wv, float rhoa[], int32_t iwnir_s, int32_t iwnir_l, 
                      int32_t *modmin, int32_t *modmax, float *modrat, float *epsnir)
{
    float *rhoas;
    float eps_ret  [MAXMODEL];
    float eps_mod  [MAXMODEL];
    float eps_err  [MAXMODEL];
    int   imod     [MAXMODEL];
    int   nmod = nmodel;
    float eps_ave;
    float *eps;
    float err_m;
    float err_p;
    int   jm, im, iim;
    int   eps_flg = 0;
    float wt;
    float tot_err;
    int   itmp;

    *modmin =  -1;
    *modmax =  -1;
    *modrat = 0.0;
    *epsnir = 0.0;

    if ((rhoas = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for rhoas.\n");
        exit(1);
    }

    /* get theoretical and retrieved epsilon for each model, and */
    /* compute average retrieved epsilon                         */

    eps_ave = 0.0;
    for (jm=0; jm<nmod; jm++) {

        im = mindx[jm];

        /* retrieve epsilon in NIR assuming each model */
        rhoa_to_rhoas(sensorID,im,solz,senz,phi,wv,rhoa,wave,nwave,iwnir_s,iwnir_l,rhoas);

        if (rhoas[iwnir_l] > 0.0000001)
        	eps_ret[jm] = rhoas[iwnir_s]/rhoas[iwnir_l];
        else
            eps_ret[jm] = 0;

        /* get model epsilon for each model at this geometry */
        eps = model_epsilon(im,iwnir_l,wave,nwave,solz,senz,phi);
        eps_mod[jm] = eps[iwnir_s];

        eps_ave += eps_ret[jm];
    }
    if(isfinite(eps_ave))
        eps_ave /= nmod;
    else
        eps_ave = 1.0;


    /* determine average retrieved epsilon for the four retrievals which most */
    /* closely match the theoretical model epsilons. the model set is reduced */
    /* from the full suite to the final four by iterative outlier rejection.  */

    while (nmod > 4) {

        /* compute differences between retrieved and model epsilon */
        for (im=0; im<nmodel; im++) {
            imod[im]    = im;
            eps_err[im] = eps_ave - eps_mod[im];
        }

        /* sort model indices by smallest to largest absolute differences */
        for (im=0; im<nmodel-1; im++) {
          for (iim=im+1; iim<nmodel; iim++)
            if ( fabs(eps_err[imod[im]]) > fabs(eps_err[imod[iim]]) ) {
              itmp      = imod[im ];
              imod[im ] = imod[iim];
              imod[iim] = itmp;
            }
	}

        /* recompute average retrieved epsilon over the n-2 closest models  */
        /* averaging is done as a weighted mean with wt=1-|eps_err|/tot_err */
        /* note that the sum of the weights is equal to n-1                 */

        nmod = nmod - 2;

        tot_err = 0.0;
        for (iim=0; iim<nmod; iim++) {
            im = imod[iim];
            tot_err += fabs(eps_err[im]);
        }

        eps_ave = 0.0;
        for (iim=0; iim<nmod; iim++) {
            im = imod[iim];
            wt = 1.0-fabs(eps_err[im])/tot_err;
            eps_ave += eps_ret[im]*wt;
        }
        eps_ave /= (nmod-1);
    }

    /* now select the two models which bracket eps_ave  */
    err_m   = 0-FLT_MAX;
    err_p   = FLT_MAX;
    for (im=0; im<nmodel; im++) {
        eps_err[im] = eps_ave - eps_mod[im];
        if (eps_err[im] >= 0.0) {
            if (eps_err[im] < err_p) {
                err_p   = eps_err[im];
                *modmin = im;
            }
        } else {
            if (eps_err[im] > err_m) {
                err_m   = eps_err[im];
                *modmax = im;
            }
        }
    }

    /* M. Wang's model cross-over correction */
    if (wang_modx && eps_mod[imc50] > eps_mod[imt99]) {
        if (*modmax == imt90)
            *modmin = imt99;
        else if (*modmax == imc50 && *modmin == imt99)
            *modmin = imc70;
        else if (*modmin == imm50)
            *modmax = imc50;
    }

    /* compute interpolation ratio */
    if (*modmin < 0) {
        /* no lower-bounding model */
        *modmin = *modmax;
        *modrat = 0.0;
        eps_flg = -1;
    } else if (*modmax < 0) {
        /* no upper-bounding model */
        *modmax = *modmin;
        *modrat = 0.0;
        eps_flg =  1;
    } else
        *modrat = (eps_ave-eps_mod[*modmin])/(eps_mod[*modmax]-eps_mod[*modmin]);

    *modmin = mindx[*modmin];
    *modmax = mindx[*modmax];

    /* return retrieved epsilon */
    *epsnir = eps_ave;

    free(rhoas);

    return(eps_flg);
}


/* ---------------------------------------------------------------------------------------- */
/* model_select_angst() - Select model pair based on input Angstrom coefficient             */
/*                                                                                          */
/* B. Franz, 1 August 2004.                                                                 */
/* ---------------------------------------------------------------------------------------- */
int compalphaT(alphaTstr *x, alphaTstr *y) {return(x->angstrom < y->angstrom ? -1 : 1);}

void model_select_angstrom(float angstrom, int32_t *modmin, int32_t *modmax, float *modrat)
{
    static alphaTstr alphaT[MAXAERMOD];
    static int firstCall = 1;

    int im, im1, im2;

    if (firstCall) {
 
        int ib = windex(520,aertab->wave,aertab->nwave);

        /* grab angstrom coefficients and sort in ascending order */ 
        for (im=0; im<aertab->nmodel; im++) {
            alphaT[im].angstrom = aertab->model[im]->angstrom[ib];
            alphaT[im].modnum = im;
        }
        qsort(alphaT,aertab->nmodel,sizeof(alphaTstr),
            (int (*)(const void *,const void *)) compalphaT);

        firstCall = 0;
    }
  
    for (im=0; im<aertab->nmodel; im++) {
         if (angstrom < alphaT[im].angstrom)
	    break;
    }
    im1 = MAX(MIN(im-1,aertab->nmodel-1),0);
    im2 = MAX(MIN(im  ,aertab->nmodel-1),0);

    *modmin = alphaT[im1].modnum;
    *modmax = alphaT[im2].modnum;

    if (im1 == im2) 
        *modrat = 1.0;
    else
        *modrat = (angstrom-alphaT[im1].angstrom) /
                  (alphaT[im2].angstrom-alphaT[im1].angstrom);

    return;
    
}

/* ---------------------------------------------------------------------------------------- */
/* model_taua() - compute AOT at input wavelengths using specified model                    */
/* ---------------------------------------------------------------------------------------- */
void model_taua(int32_t sensorID, int modnum, float wave[], int32_t nwave, int32_t iwnir_l, float rhoa[], 
                float solz, float senz, float phi, float wv, float taua[])
{
    float  *aot  ;
    float  *lnaot;
    float  *rhoas;

    int    iwnir  = iwatab[iwnir_l];
    float *phase  = model_phase(modnum,solz,senz,phi);
    float  csolz  = cos(solz/radeg);    
    float  csenz  = cos(senz/radeg);    
    int    iw, iwtab;
    float maxwave;

    maxwave=MAX(aertab->nwave,nwave);

    if ((rhoas = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for rhoas.\n");
        exit(1);
    }
    if ((aot = (float *)calloc(maxwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for aot.\n");
        exit(1);
    }
    if ((lnaot = (float *)calloc(maxwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for lnaot.\n");
        exit(1);
    }

    /* get SS aerosol reflectance at longest sensor wavelength */
    rhoa_to_rhoas(sensorID,modnum,solz,senz,phi,wv,rhoa,wave,nwave,iwnir_l,iwnir_l,rhoas);

    /* get aerosol optical thickness at longest sensor wavelength */
    aot[iwnir] = rhoas[iwnir_l]*(4.0*csolz*csenz)/(phase[iwnir]*aertab->model[modnum]->albedo[iwnir]); 

    /* get aerosol optical thickness at all other table wavelengths */
    for (iw=0; iw<aertab->nwave; iw++) {
        /* note to self: i actually have aot(865) for czcs at this point */
        aot[iw] = aot[iwnir] * aertab->model[modnum]->extc[iw] / aertab->model[modnum]->extc[iwnir];
        if (interp)
            lnaot[iw] = log(aot[iw]);
    }

    /* interpolate aot from table to sensor wavelengths */ 
    if (interp) {
        for (iw=0; iw<nwave; iw++) {
	    iwtab = iwatab[iw];
	    if (aertab->wave[iwtab] != wave[iw] && wave[iw] > 0)
	        taua[iw] = exp(linterp(aertab->wave,lnaot,aertab->nwave,wave[iw]));
            else
	        taua[iw] = aot[iwtab];
	}
    } else
        for (iw=0; iw<nwave; iw++)
            taua[iw] = aot[iw];

    free(rhoas);
    free(aot);
    free(lnaot);
    return;
}


/* ---------------------------------------------------------------------------------------- */
/* model_transmittance() - compute path Rayleigh-aerosol diffuse trans for specified model  */
/* ---------------------------------------------------------------------------------------- */
void model_transmittance(int modnum, float wave[], int32_t nwave, 
                         float theta, float taua[], float dtran[])
{
    static int firstCall = 1;
    static float *intexp;
    static int   *inttst;

    int   i1, i2, i;
    int   iw, iwtab;
    float a1, b1;
    float a2, b2;
    float x1,x2;
    float y1,y2;
    float xbar;
    float wt;

    if (firstCall) {
        float taur1, um1;
        float taur2, um2;
        firstCall = 0;

        intexp = (float *) malloc(nwave*sizeof(float));
        inttst =   (int *) malloc(nwave*sizeof(int));

        for (iw=0; iw<nwave; iw++) {
	    intexp[iw] = 1.0;
	    inttst[iw] = 0;
            iwtab = iwdtab[iw];
            if (fabs(wave[iw] - aertab->dtran_wave[iwtab]) > 0.51) {
	        um1 = aertab->dtran_wave[iwtab]/1000.0;
	        um2 = wave[iw]/1000.0;
                taur1 = 0.008569*pow(um1,-4)*(1.0+(0.0113*pow(um1,-2))+(0.00013*pow(um1,-4)));
                taur2 = 0.008569*pow(um2,-4)*(1.0+(0.0113*pow(um2,-2))+(0.00013*pow(um2,-4)));
		intexp[iw] = taur2/taur1;
                inttst[iw] = 1;
                printf("Interpolating diffuse transmittance for %d from %f by %f\n",
		       (int) wave[iw],aertab->dtran_wave[iwtab],intexp[iw]);
	    } 
	}
    }

    /* find bracketing zenith angle indices */
    for (i=0; i<aertab->dtran_ntheta; i++) {
        if (theta < aertab->dtran_theta[i])
	    break;
    }
    if (i == aertab->dtran_ntheta) {
        i1 = i-1;
        i2 = i-1;
        wt = 0.0;
    } else {
        i1 = MIN(MAX(i-1,0),aertab->dtran_ntheta-2);
        i2 = i1+1;
        x1   = aertab->dtran_airmass[i1];
        x2   = aertab->dtran_airmass[i2];
        xbar = 1.0/cos(theta/radeg);
        wt   = (xbar-x1)/(x2-x1);
    }

    /* use coefficients of nearest wavelength */
    for (iw=0; iw<nwave; iw++) {

        iwtab = iwdtab[iw];

        a1  = aertab->model[modnum]->dtran_a[iwtab][i1];
        b1  = aertab->model[modnum]->dtran_b[iwtab][i1];

        a2  = aertab->model[modnum]->dtran_a[iwtab][i2];
        b2  = aertab->model[modnum]->dtran_b[iwtab][i2];

        if (inttst[iw]) {
	    a1 = pow(a1,intexp[iw]);
	    a2 = pow(a2,intexp[iw]);
	}

        y1 = a1 * exp(-b1*taua[iw]);
        y2 = a2 * exp(-b2*taua[iw]);

        dtran[iw] = MAX(MIN( (1.0-wt)*y1 + wt*y2, 1.0), 1e-5);
    }

    return;
}


/* ---------------------------------------------------------------------------------------- */
/* diff_tran() - compute Rayleigh-aerosol diffuse trans for selected model pair, both paths */
/* ---------------------------------------------------------------------------------------- */
void diff_tran(int32_t sensorID,float wave[], int32_t nwave, int32_t iwnir_l, float solz, float senz, float phi, 
               float wv, float pr, float taur[], int32_t modmin, int32_t modmax, float modrat, float rhoa[], 
               float taua[], float tsol[], float tsen[], float tauamin[],float tauamax[], int tauaopt)
{
    int iw;
    float *tsolmin;
    float *tsolmax;
    float *tsenmin;
    float *tsenmax;

    if ((tsolmin = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for tsolmin.\n");
        exit(1);
    }
    if ((tsolmax = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for tsolmax.\n");
        exit(1);
    }
    if ((tsenmin = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for tsenmin.\n");
        exit(1);
    }
    if ((tsenmax = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for tsenmax.\n");
        exit(1);
    }

    /* get AOT per band for each model */
    if (tauaopt != 1) {
        model_taua(sensorID,modmin,wave,nwave,iwnir_l,rhoa,solz,senz,phi,wv,tauamin);
        model_taua(sensorID,modmax,wave,nwave,iwnir_l,rhoa,solz,senz,phi,wv,tauamax);
    }

    /* get diff trans sun to ground, per band for each model */
    model_transmittance(modmin,wave,nwave,solz,tauamin,tsolmin);
    model_transmittance(modmax,wave,nwave,solz,tauamax,tsolmax);

    /* get diff trans ground to sensor, per band for each model */
    model_transmittance(modmin,wave,nwave,senz,tauamin,tsenmin);
    model_transmittance(modmax,wave,nwave,senz,tauamax,tsenmax);

    /* interpolate and pressure correct */
    for (iw=0; iw<nwave; iw++) {
        taua[iw] = tauamin[iw]*(1.0-modrat) + tauamax[iw]*modrat;
        tsol[iw] = tsolmin[iw]*(1.0-modrat) + tsolmax[iw]*modrat;
        tsen[iw] = tsenmin[iw]*(1.0-modrat) + tsenmax[iw]*modrat;

        /* correct for pressure difference from standard pressure */
        tsol[iw] = tsol[iw]*exp(-0.5*taur[iw]/mu0*(pr/p0-1));
        tsen[iw] = tsen[iw]*exp(-0.5*taur[iw]/mu *(pr/p0-1));

        if ((evalmask & TRANSSPHER) != 0) {
	    /* correct for airmass difference, plane-parallel to spherical atmosphere */
            tsol[iw] = pow(tsol[iw],airmass_sph/airmass_plp);
            tsen[iw] = pow(tsen[iw],airmass_sph/airmass_plp);
	}
    }

    free(tsolmin);
    free(tsolmax);
    free(tsenmin);
    free(tsenmax);

    return;
}


/* ---------------------------------------------------------------------------------------- */
/* ahmadaer() - compute aerosol reflectance using MSEPS approach of Ahmad                   */
/*                                                                                          */
/* Z. Ahmad, August 2014.                                                                   */
/* ---------------------------------------------------------------------------------------- */
int ahmadaer(int32_t sensorID, float wave[], int32_t nwave, int32_t iwnir_s, int32_t iwnir_l,
            int32_t nmodels, int32_t mindx[],             
            float solz, float senz, float phi, float wv, float rhoa[], 
            int32_t *modmin, int32_t *modmax, float *modrat, float *epsnir,
            float tau_pred_min[], float tau_pred_max[])
{
    int   iw;

    float rho_pred_min[nwave],rho_pred_max[nwave];
    float rho_aer[nwave],tau_aer[nwave];

    if (!have_ms) {
        printf("\nThe multi-scattering epsilon atmospheric correction method requires\n");
        printf("ams_all, bms_all, cms_all, dms_all, ems in the aerosol model tables.\n");
        exit(1);
    }

    /* use the ms_epsilon method to get rhoa */
    if (ahmad_atm_corr(sensorID,wave,nwave,iwnir_s,iwnir_l,nmodels,mindx,solz,senz,phi,wv,
                       rhoa,modmin,modmax,modrat,epsnir,
                       tau_pred_max,tau_pred_min,rho_pred_max,rho_pred_min,tau_aer,rho_aer) !=0)
        return(1);

   for (iw=0; iw < nwave; iw++) {
       rhoa[iw] = rho_aer[iw]; 
     }

    return(0);
}


/* ---------------------------------------------------------------------------------------- */
/* wangaer() - compute aerosol reflectance using Gordon & Wang 1994 algorithm               */
/*                                                                                          */
/* B. Franz, 1 June 2004.                                                                   */
/* ---------------------------------------------------------------------------------------- */
int wangaer(int32_t sensorID, float wave[], int32_t nwave, int32_t iwnir_s, int32_t iwnir_l,
            int32_t nmodels, int32_t mindx[],             
            float solz, float senz, float phi, float wv, float rhoa[],
            int32_t *modmin, int32_t *modmax, float *modrat, float *epsnir,
            float tauamin[], float tauamax[])
{
    int    modflg;
    float *epsmin;
    float *epsmax;

    float epsmin1;
    float epsmax1;

    float *rhoasmin;
    float *rhoasmax;
    float *rhoamin;
    float *rhoamax;

    float cc = 0.0;
    int   iw;

    if ((rhoasmin = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for rhoasmin.\n");
        exit(1);
    }
    if ((rhoasmax = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for rhoasmax.\n");
        exit(1);
    }
    if ((rhoamin = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for rhoamin.\n");
        exit(1);
    }
    if ((rhoamax = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for rhoasmax.\n");
        exit(1);
    }

    /* find upper and lower-bounding models */
    modflg = model_select_wang(sensorID,wave,nwave,nmodels,mindx,solz,senz,phi,wv,rhoa,iwnir_s,iwnir_l,
                               modmin,modmax,modrat,epsnir);

    /* if no lower-bounding aerosol model, set-up for extrapolation */
    if (modflg < 0) 
        cc = log(*epsnir)/(wave[iwnir_l]-wave[iwnir_s]);

    /* get model epsilon for each bounding model, all wavelengths */
    epsmin = model_epsilon(*modmin,iwnir_l,wave,nwave,solz,senz,phi);
    epsmax = model_epsilon(*modmax,iwnir_l,wave,nwave,solz,senz,phi);

    /* get SS aerosol reflectance at longest wavelength for the two models */
    if (rhoa_to_rhoas(sensorID,*modmin,solz,senz,phi,wv,rhoa,wave,nwave,iwnir_l,iwnir_l,rhoasmin) != 0) {
        free(rhoamin);
        free(rhoasmin);
        free(rhoamax);
        free(rhoasmax);
        return(1);
    }
    if (rhoa_to_rhoas(sensorID,*modmax,solz,senz,phi,wv,rhoa,wave,nwave,iwnir_l,iwnir_l,rhoasmax) != 0) {
        free(rhoamin);
        free(rhoasmin);
        free(rhoamax);
        free(rhoasmax);
      return(1);
    }
    /* compute SS aerosol reflectance in all bands */
    for (iw=0; iw<nwave; iw++) {

        epsmin1 = epsmin[iw];
        epsmax1 = epsmax[iw];

        if (modflg < 0) {
	  epsmin1 = exp(cc*(wave[iwnir_l]-wave[iw]));
	  epsmax1 = epsmin1;
	}

        rhoasmin[iw] = rhoasmin[iwnir_l]*epsmin1;
        rhoasmax[iw] = rhoasmax[iwnir_l]*epsmax1;
    }

    /* compute MS aerosol reflectance in visible bands */
    rhoas_to_rhoa(sensorID,*modmin,solz,senz,phi,wv,rhoasmin,wave,nwave,0,nwave-1,rhoamin);
    rhoas_to_rhoa(sensorID,*modmax,solz,senz,phi,wv,rhoasmax,wave,nwave,0,nwave-1,rhoamax);

    /* interpolate between upper and lower-bounding models */
    for (iw=0; iw<nwave; iw++) {
        rhoa[iw] = rhoamin[iw]*(1.0-(*modrat)) + rhoamax[iw]*(*modrat);
    }

    model_taua(sensorID,*modmin,wave,nwave,iwnir_l,rhoa,solz,senz,phi,wv,tauamin);
    model_taua(sensorID,*modmax,wave,nwave,iwnir_l,rhoa,solz,senz,phi,wv,tauamax);

    free(rhoamin);
    free(rhoasmin);
    free(rhoamax);
    free(rhoasmax);

    return(0);
}

/* ---------------------------------------------------------------------------------------- */
/* model_select_franz() - Franz aerosol model selection process.                            */
/*                                                                                          */
/* B. Franz, 1 February 2009.                                                               */
/* ---------------------------------------------------------------------------------------- */

typedef struct rhoaT_struct {
    int32_t   modnum;
    float  rhoa;
    float  eps;
} rhoaTstr;

int comp_rhoaT(rhoaTstr *x, rhoaTstr *y) {return(x->rhoa < y->rhoa ? -1 : 1);}

int model_select_franz(int32_t sensorID, float wave[], int32_t nwave, int32_t nmodel, int32_t mindx[],
                      float solz, float senz, float phi,
                      float wv, float rhoa[], int32_t iwnir_s, int32_t iwnir_l,
                      int32_t *modmin, int32_t *modmax, float *modrat, float *epsnir)
{

    float *rhoas;
    float *rhoa_tmp;
    rhoaTstr rhoa_tab[MAXMODEL];

    float *eps;
    int   jm, im;
    int   jm1, jm2;
    float wt;

    *modmin =  -1;
    *modmax =  -1;
    *modrat = 0.0;
    *epsnir = 0.0;

    if ((rhoas = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for rhoas.\n");
        exit(1);
    }
    if ((rhoa_tmp = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for rhoas_tmp.\n");
        exit(1);
    }

    // predict MS aerosol reflectance assuming each model

    for (jm=0; jm<nmodel; jm++) {

        im = mindx[jm];

        // get model epsilon at this geometry 
        eps = model_epsilon(im,iwnir_l,wave,nwave,solz,senz,phi);

        // get SS aerosol reflectance at iwnir_l 
        if (rhoa_to_rhoas(sensorID,im,solz,senz,phi,wv,rhoa,wave,nwave,iwnir_l,iwnir_l,rhoas) != 0) {
            free(rhoas);
            free(rhoa_tmp);
            return(1);
        }
        // get SS aerosol reflectance at iwnir_s
        rhoas[iwnir_s] = rhoas[iwnir_l] * eps[iwnir_s];

        // get MS aerosol reflectance at iwnir_s and save
        rhoas_to_rhoa(sensorID,im,solz,senz,phi,wv,rhoas,wave,nwave,iwnir_s,iwnir_s,rhoa_tmp);

        rhoa_tab[jm].modnum = im;
        rhoa_tab[jm].rhoa = rhoa_tmp[iwnir_s];
        rhoa_tab[jm].eps  = eps[iwnir_s];
    }

    // put results in ascending order of predicted rhoa[iwnir_s]
    qsort(rhoa_tab,nmodel,sizeof(rhoaTstr),(int (*)(const void *,const void *)) comp_rhoaT);

    // compare observed rhoa with model predictions at iwnir_s to select models
    for (jm=0; jm<nmodel; jm++) {
        if (rhoa_tab[jm].rhoa > rhoa[iwnir_s])
            break;
    }
    if (jm == 0) {
        jm1 = 0;
        jm2 = 1;
    } else if (jm == nmodel) {
        jm1 = nmodel-2;
        jm2 = nmodel-1;
    } else {
        jm1 = jm  - 1;
        jm2 = jm1 + 1;
    }
    wt  = (rhoa[iwnir_s] - rhoa_tab[jm1].rhoa)/(rhoa_tab[jm2].rhoa - rhoa_tab[jm1].rhoa);

    *modmin = rhoa_tab[jm1].modnum;
    *modmax = rhoa_tab[jm2].modnum;
    *modrat = wt;
    *epsnir = rhoa_tab[jm1].eps*(1.0-wt) + rhoa_tab[jm2].eps*wt;

    free(rhoas);
    free(rhoa_tmp);

    return(0);
}

/* ---------------------------------------------------------------------------------------- */
/* franzaer() - compute aerosol reflectance using modified G&W algorithm                    */
/*                                                                                          */
/* B. Franz, February 2009.                                                                 */
/* ---------------------------------------------------------------------------------------- */
int franzaer(int32_t sensorID, float wave[], int32_t nwave, int32_t iwnir_s, int32_t iwnir_l,
            int32_t nmodels, int32_t mindx[],
            float solz, float senz, float phi, float wv, float rhoa[],
	    int32_t *modmin, int32_t *modmax, float *modrat, float *epsnir,
            float tauamin[], float tauamax[])
{
    float *epsmin;
    float *epsmax;

    float epsmin1;
    float epsmax1;

    float *rhoasmin;
    float *rhoasmax;
    float *rhoamin;
    float *rhoamax;

    int   iw;

    if ((rhoasmin = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for rhoasmin.\n");
        exit(1);
    }
    if ((rhoasmax = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for rhoasmax.\n");
        exit(1);
    }
    if ((rhoamin = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for rhoamin.\n");
        exit(1);
    }
    if ((rhoamax = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for rhoamax.\n");
        exit(1);
    }


    /* find upper and lower-bounding models */
    if (model_select_franz(sensorID,wave,nwave,nmodels,mindx,solz,senz,phi,wv,rhoa,iwnir_s,iwnir_l,
                           modmin,modmax,modrat,epsnir) != 0)
        return(1);

    /* get model epsilon for each bounding model, all wavelengths */
    epsmin = model_epsilon(*modmin,iwnir_l,wave,nwave,solz,senz,phi);
    epsmax = model_epsilon(*modmax,iwnir_l,wave,nwave,solz,senz,phi);

    /* get SS aerosol reflectance at longest wavelength for the two models */
    if (rhoa_to_rhoas(sensorID,*modmin,solz,senz,phi,wv,rhoa,wave,nwave,iwnir_l,iwnir_l,rhoasmin) != 0) {
        free(rhoasmin);
        free(rhoasmax);
        free(rhoamin);
        free(rhoamax);
       return(1);
    }
    if (rhoa_to_rhoas(sensorID,*modmax,solz,senz,phi,wv,rhoa,wave,nwave,iwnir_l,iwnir_l,rhoasmax) != 0) {
        free(rhoasmin);
        free(rhoasmax);
        free(rhoamin);
        free(rhoamax);
        return(1);
    }
    /* compute SS aerosol reflectance in all bands */
    for (iw=0; iw<nwave; iw++) {

        epsmin1 = epsmin[iw];
        epsmax1 = epsmax[iw];

        rhoasmin[iw] = rhoasmin[iwnir_l]*epsmin1;
        rhoasmax[iw] = rhoasmax[iwnir_l]*epsmax1;
    }

    /* compute MS aerosol reflectance in visible bands */
    rhoas_to_rhoa(sensorID,*modmin,solz,senz,phi,wv,rhoasmin,wave,nwave,0,nwave-1,rhoamin);
    rhoas_to_rhoa(sensorID,*modmax,solz,senz,phi,wv,rhoasmax,wave,nwave,0,nwave-1,rhoamax);

    /* interpolate between upper and lower-bounding models */
    for (iw=0; iw<nwave; iw++) {
        rhoa[iw] = rhoamin[iw]*(1.0-(*modrat)) + rhoamax[iw]*(*modrat);
    }

    model_taua(sensorID,*modmin,wave,nwave,iwnir_l,rhoa,solz,senz,phi,wv,tauamin);
    model_taua(sensorID,*modmax,wave,nwave,iwnir_l,rhoa,solz,senz,phi,wv,tauamax);

    free(rhoasmin);
    free(rhoasmax);
    free(rhoamin);
    free(rhoamax);

    return(0);
}

/* ---------------------------------------------------------------------------------------- */
/* ---------------------------------------------------------------------------------------- */
static int order_models(const void *p1, const void *p2)
{
    aermodstr *x = *(aermodstr **)p1;
    aermodstr *y = *(aermodstr **)p2;

    if (x->rh == y->rh) {
      if (x->sd > y->sd) 
        return( 1);
      else
        return(-1);
    } else {
      if (x->rh > y->rh) 
        return( 1);
      else
        return(-1);
    }
}


/* ---------------------------------------------------------------------------------------- */
/* rhaer() - compute aerosol reflectance using RH descrimination + desired selection scheme */
/* ---------------------------------------------------------------------------------------- */
int rhaer(int32_t sensorID, float wave[], int32_t nwave, int32_t iwnir_s, int32_t iwnir_l,
	  float solz, float senz, float phi, float wv, float rh, float pr, float taur[], float rhoa[],
          int32_t *modmin1, int32_t *modmax1, float *modrat1, int32_t *modmin2, int32_t *modmax2, float *modrat2, 
          float *eps, float taua[], float tsol[], float tsen[])
{
    static int   firstCall = 1;
    static int   nrh;
    static float rhtab[MAXAERMOD];
    static int   nsd;
    static int   sdtab[MAXAERMOD];

    float *rhoa1;
    float *rhoa2;
    float *taua1;
    float *taua2;
    float *tsol1;
    float *tsol2;
    float *tsen1;
    float *tsen2;
    float eps1;
    float eps2;
    //float modrat1;
    //float modrat2;

    float *tau_pred_min1;
    float *tau_pred_max1;
    float *tau_pred_min2;
    float *tau_pred_max2;

    //int32_t  modmin1;
    //int32_t  modmin2;
    //int32_t  modmax1;
    //int32_t  modmax2;
    int32_t  mindx1[MAXAERMOD];
    int32_t  mindx2[MAXAERMOD];
    int   irh1, irh2, irh;
    int   isd;
    float wt;

    int   iw, im;

    if (firstCall) {
        firstCall = 0;
        float lastrh = -1.0;
        int   lastsd = -1;
        if (!have_rh || !have_sd) {
	    printf("-E- %s line %d: This aerosol selection method requires models with a Relative Humidity attribute.\n",
		   __FILE__,__LINE__);
            exit(1);
	}

        // need in order of rh and sd within rh
        qsort(aertab->model,aertab->nmodel,sizeof(aermodstr*),(int (*)(const void *,const void *)) order_models);

        // count the number of model humidities and the number of model size distributions
        // note that use of a single model suite will yield nrh=1, which inherently avoids RH weighting that case

        nsd = 0;
        nrh = 0;

        for (im=0; im<aertab->nmodel; im++) {
            if (aertab->model[im]->rh != lastrh) {
	        rhtab[nrh] = aertab->model[im]->rh;
                lastrh     = rhtab[nrh];
                nrh++;
	    }
            if (nrh == 1 && aertab->model[im]->sd != lastsd) {
	        sdtab[nsd] = aertab->model[im]->sd;
                lastsd     = sdtab[nsd]; 
	        nsd++;
	    }
        }
        if (nrh*nsd != aertab->nmodel) {
	    printf("-E- %s line %d: number of humidities (%d) x number of size distributions (%d) must equal number of models (%d).\n",
		   __FILE__,__LINE__,nrh,nsd,aertab->nmodel);
            exit(1);
	} else {
  	    printf("%d aerosol models: %d humidities x %d size fractions\n",aertab->nmodel,nrh,nsd);
            for (irh=0; irh<nrh; irh++) {
              for (isd=0; isd<nsd; isd++) {
                im = irh*nsd+isd;
		printf("model %d, rh=%f, sd=%d, alpha=%f, name=%s\n", 
		     im,aertab->model[im]->rh, aertab->model[im]->sd, aertab->model[im]->angstrom[1],aertab->model[im]->name);
	      }
	    }
	}
    }

    // initialize
    if ((taua1 = (float *)malloc(nwave*sizeof(float))) == NULL) {
        printf("Unable to allocate space for taua1.\n");
        exit(1);
    }
    if ((taua2 = (float *)malloc(nwave*sizeof(float))) == NULL) {
        printf("Unable to allocate space for taua2.\n");
        exit(1);
    }
    if ((tsol1 = (float *)malloc(nwave*sizeof(float))) == NULL) {
        printf("Unable to allocate space for tsol1.\n");
        exit(1);
    }
    if ((tsol2 = (float *)malloc(nwave*sizeof(float))) == NULL) {
        printf("Unable to allocate space for tsol2.\n");
        exit(1);
    }
    if ((rhoa1 = (float *)malloc(nwave*sizeof(float))) == NULL) {
        printf("Unable to allocate space for rhoa1.\n");
        exit(1);
    }
    if ((rhoa2 = (float *)malloc(nwave*sizeof(float))) == NULL) {
        printf("Unable to allocate space for rhoa2.\n");
        exit(1);
    }
    if ((tsen1 = (float *)malloc(nwave*sizeof(float))) == NULL) {
        printf("Unable to allocate space for tsen1.\n");
        exit(1);
    }
    if ((tsen2 = (float *)malloc(nwave*sizeof(float))) == NULL) {
        printf("Unable to allocate space for tsen2.\n");
        exit(1);
    }
    if ((tau_pred_min1 = (float *)malloc(nwave*sizeof(float))) == NULL) {
        printf("Unable to allocate space for tau_pred_min1.\n");
        exit(1);
    }
    if ((tau_pred_min2 = (float *)malloc(nwave*sizeof(float))) == NULL) {
        printf("Unable to allocate space for tau_pred_min2.\n");
        exit(1);
    }
    if ((tau_pred_max1 = (float *)malloc(nwave*sizeof(float))) == NULL) {
        printf("Unable to allocate space for tau_pred_max1.\n");
        exit(1);
    }
    if ((tau_pred_max2 = (float *)malloc(nwave*sizeof(float))) == NULL) {
        printf("Unable to allocate space for tau_pred_max2.\n");
        exit(1);
    }

    //*modmin  = -1;
    //*modmax  = -1;
    //*modrat  =  0;
    //*modmin2 = -1;
    //*modmax2 = -1;
    //*modrat2 =  0;
 
    for (iw=0; iw<nwave; iw++) {
        taua[iw] = -1.0;
        tsol[iw] = -1.0;
        tsen[iw] = -1.0;
        rhoa1[iw] = rhoa[iw];
        rhoa2[iw] = rhoa[iw];
        rhoa [iw] = BAD_FLT;
    }

    // find RH index and wts

    if (nrh == 1 || rhtab[0] > rh) { // actual RH < smallest model RH or only one model RH
        irh1 = 0;
        irh2 = 0;
        wt = 0.0;
    } else if (rhtab[nrh-1] < rh) { // actual RH > smallestlargest model RH
        irh1 = nrh-1;
        irh2 = nrh-1;
        wt = 0.0;
    } else {
        for (irh=0; irh<nrh; irh++) {
            if (rhtab[irh] > rh)
                break;
        }
        irh1 = MIN(MAX(0,irh-1),nrh-2);
        irh2 = irh1 + 1;
        wt = (rh - rhtab[irh1])/(rhtab[irh2]-rhtab[irh1]);
    }

    // set indices of active model sets

    for (im=0; im<nsd; im++) {
        mindx1[im] = irh1*nsd + im;
        mindx2[im] = irh2*nsd + im;
    }

    // compute aerosol reflectances, aot, diffuse trans, eps from first model set
 
    if (aer_opt == AERRHMSEPS) {
        if ( ahmadaer(sensorID,wave,nwave,iwnir_s,iwnir_l,nsd,mindx1,
                     solz,senz,phi,wv,rhoa1,modmin1,modmax1,modrat1,&eps1,tau_pred_min1,tau_pred_max1) != 0) {
            free(taua1);
            free(taua2);
            free(tsol1);
            free(tsol2);
            free(tsen1);
            free(tsen2);
            free(rhoa1);
            free(rhoa2);
            free(tau_pred_min1);
            free(tau_pred_max1);
            free(tau_pred_min2);
            free(tau_pred_max2);
            return(1);
        }
    } else if (aer_opt == AERRHFRNIR) {
        if ( franzaer(sensorID,wave,nwave,iwnir_s,iwnir_l,nsd,mindx1,
		      solz,senz,phi,wv,rhoa1,modmin1,modmax1,modrat1,&eps1,tau_pred_min1,tau_pred_max1) != 0) {
            free(taua1);
            free(taua2);
            free(tsol1);
            free(tsol2);
            free(tsen1);
            free(tsen2);
            free(rhoa1);
            free(rhoa2);
            free(tau_pred_min1);
            free(tau_pred_max1);
            free(tau_pred_min2);
            free(tau_pred_max2);
            return(1);
        }
    } else {
        if ( wangaer(sensorID,wave,nwave,iwnir_s,iwnir_l,nsd,mindx1,
                     solz,senz,phi,wv,rhoa1,modmin1,modmax1,modrat1,&eps1,tau_pred_min1,tau_pred_max1) != 0) {
            free(taua1);
            free(taua2);
            free(tsol1);
            free(tsol2);
            free(tsen1);
            free(tsen2);
            free(rhoa1);
            free(rhoa2);
            free(tau_pred_min1);
            free(tau_pred_max1);
            free(tau_pred_min2);
            free(tau_pred_max2);
           return(1);
        }
    }
 
    diff_tran(sensorID,wave,nwave,iwnir_l,solz,senz,phi,wv,pr,taur,
              *modmin1,*modmax1,*modrat1,rhoa1,taua1,tsol1,tsen1,tau_pred_min1,tau_pred_max1,1);

    // compute aerosol reflectances, aot, diffuse trans, eps from second model set (if needed)

    if (irh2 != irh1) {

        if (aer_opt == AERRHMSEPS) {
            if ( ahmadaer(sensorID,wave,nwave,iwnir_s,iwnir_l,nsd,mindx2,
                         solz,senz,phi,wv,rhoa2,modmin2,modmax2,modrat2,&eps2,tau_pred_min2,tau_pred_max2) != 0) {
                free(taua1);
                free(taua2);
                free(tsol1);
                free(tsol2);
                free(tsen1);
                free(tsen2);
                free(rhoa1);
                free(rhoa2);
                free(tau_pred_min1);
                free(tau_pred_max1);
                free(tau_pred_min2);
                free(tau_pred_max2);
               return(1);
            }
      } else if (aer_opt == AERRHFRNIR) {
            if ( franzaer(sensorID,wave,nwave,iwnir_s,iwnir_l,nsd,mindx2,
			  solz,senz,phi,wv,rhoa2,modmin2,modmax2,modrat2,&eps2,tau_pred_min2,tau_pred_max2) != 0) {
                free(taua1);
                free(taua2);
                free(tsol1);
                free(tsol2);
                free(tsen1);
                free(tsen2);
                free(rhoa1);
                free(rhoa2);
                free(tau_pred_min1);
                free(tau_pred_max1);
                free(tau_pred_min2);
                free(tau_pred_max2);
              return(1);
            }
        } else {
            if ( wangaer(sensorID,wave,nwave,iwnir_s,iwnir_l,nsd,mindx2,
                         solz,senz,phi,wv,rhoa2,modmin2,modmax2,modrat2,&eps2,tau_pred_min2,tau_pred_max2) != 0) {
                free(taua1);
                free(taua2);
                free(tsol1);
                free(tsol2);
                free(tsen1);
                free(tsen2);
                free(rhoa1);
                free(rhoa2);
                free(tau_pred_min1);
                free(tau_pred_max1);
                free(tau_pred_min2);
                free(tau_pred_max2);
               return(1);
            }
        }

        diff_tran(sensorID,wave,nwave,iwnir_l,solz,senz,phi,wv,pr,taur,
                  *modmin2,*modmax2,*modrat2,rhoa2,taua2,tsol2,tsen2,tau_pred_min2,tau_pred_max2,1);

        for (iw=0; iw<nwave; iw++) {
            rhoa[iw] = rhoa1[iw]*(1-wt) + rhoa2[iw]*wt;
            taua[iw] = taua1[iw]*(1-wt) + taua2[iw]*wt;
            tsol[iw] = tsol1[iw]*(1-wt) + tsol2[iw]*wt;
            tsen[iw] = tsen1[iw]*(1-wt) + tsen2[iw]*wt;
        }
        *eps = eps1*(1-wt) + eps2*wt;

    } else {

        for (iw=0; iw<nwave; iw++) {
            rhoa[iw] = rhoa1[iw];
            taua[iw] = taua1[iw];
            tsol[iw] = tsol1[iw];
            tsen[iw] = tsen1[iw];
        }
        *eps = eps1;
    }

    // store model info for the dominant RH only

    //if (wt < 0.5) {
    //    *modmin = modmin1;
    //    *modmax = modmax1;
    //    *modrat = modrat1;
    //} else {
    //    *modmin = modmin2;
    //    *modmax = modmax2;
    //    *modrat = modrat2;
    //}

    free(taua1);
    free(taua2);
    free(tsol1);
    free(tsol2);
    free(tsen1);
    free(tsen2);
    free(rhoa1);
    free(rhoa2);
    free(tau_pred_min1);
    free(tau_pred_max1);
    free(tau_pred_min2);
    free(tau_pred_max2);

    return(0);
}


/* ---------------------------------------------------------------------------------------- */
/* fixedaer() - compute aerosol reflectance for fixed aerosol model                         */
/*                                                                                          */
/* B. Franz, August 2004.                                                                   */
/* ---------------------------------------------------------------------------------------- */
int fixedaer(int32_t sensorID, int32_t modnum, float wave[], int32_t nwave, int32_t iwnir_s, int32_t iwnir_l, 
             char models[MAXAERMOD][32], int32_t nmodels,
             float solz, float senz, float phi, float wv, float rhoa[],float *epsnir)
{
    float *eps;
    float *rhoas;
    int    iw;

    if ((rhoas = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for rhoas.\n");
        exit(1);
    }

    if (rhoa[iwnir_l] < 0.0) {
        *epsnir = BAD_FLT;
//        for (iw=0; iw<nwave; iw++)
//            rhoas[iw] = BAD_FLT;
        free(rhoas);
        return(1);
    }

    /* get model epsilon for all wavelengths at this geometry */
    eps = model_epsilon(modnum,iwnir_l,wave,nwave,solz,senz,phi);

    /* get SS aerosol reflectance at longest wavelength */
    if (rhoa_to_rhoas(sensorID,modnum,solz,senz,phi,wv,rhoa,wave,nwave,iwnir_l,iwnir_l,rhoas) != 0) {
        printf("Error getting rhoas\n");
        free(rhoas);
        return(1);
    }

    /* compute SS aerosol reflectance in visible bands */
    for (iw=0; iw<nwave; iw++) {
        rhoas[iw] = rhoas[iwnir_l]*eps[iw];
    }

    /* compute MS aerosol reflectance in visible bands */
    rhoas_to_rhoa(sensorID,modnum,solz,senz,phi,wv,rhoas,wave,nwave,0,nwave-1,rhoa);

    if (iwnir_s == iwnir_l) 
        *epsnir = eps[iwnir_l-1];
    else
        *epsnir = eps[iwnir_s];

    free(rhoas);

    return(0);
}


/* ---------------------------------------------------------------------------------------- */
/* fixedmodpair() - compute aerosol reflectance for fixed model pair                        */
/*                                                                                          */
/* B. Franz, August 2004.                                                                   */
/* ---------------------------------------------------------------------------------------- */
int fixedmodpair(
             int32_t sensorID, instr *input, float wave[], int32_t nwave, int32_t iwnir_s, int32_t iwnir_l, 
             float solz, float senz, float phi, float wv, int32_t modmin, int32_t modmax, float modrat, 
             float rhoa[],float *eps)
{
    float *rhoa1;
    float *rhoa2;
    float eps1;
    float eps2;
    int   iw;
    int   status;

    if ((rhoa1 = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for rhoa1.\n");
        exit(1);
    }
    if ((rhoa2 = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for rhoa2.\n");
        exit(1);
    }

    if (modmin < 0 || modmin >= input->naermodels ||
        modmax < 0 || modmax >= input->naermodels ||
        modrat < 0.0 || modrat > 1.0) {
        printf("Bogus input for fixed model pair: %d %d %f\n",modmin+1,modmax+1,modrat);
        exit(1);
    } 


    if  (rhoa[iwnir_l] > input->rhoamin) {

        rhoa2[iwnir_l] = rhoa1[iwnir_l] = rhoa[iwnir_l];

        /* get aerosol reflectance in visible for fixed model, and epsilon */
        status = fixedaer(sensorID,modmin,wave,nwave,iwnir_s,iwnir_l,
                      input->aermodels,input->naermodels,
                       solz,senz,phi,wv,rhoa1,&eps1);
        status = fixedaer(sensorID,modmax,wave,nwave,iwnir_s,iwnir_l,
                      input->aermodels,input->naermodels,
                      solz,senz,phi,wv,rhoa2,&eps2);

        /* convert aerosol relectance to radiance */
        if (status == 0) {
            for (iw=0; iw<nwave; iw++) {
	      if (iw != iwnir_l) // without this check tLw-La may go slight negative              
	        rhoa[iw] = MAX((1.0-modrat)*rhoa1[iw] + modrat*rhoa2[iw],0.0); 
  	    }
            *eps  = (1.0-modrat)*eps1 + modrat*eps2; 
	}

    } else if (rhoa[iwnir_l] > -(input->rhoamin)) {

	/* if input NIR is near zero, assume a small white aerosol signal */
        *eps    = 1.0;
        for (iw=0; iw<nwave; iw++) {
	    rhoa[iw] = MAX(rhoa[iwnir_l],1e-6);
	}

  	status = 0;

    } else {

	/* if input NIR is very negative, fail the case */
  	status = 1;
    }

    free(rhoa1);
    free(rhoa2);

    return(status);
}


/* ---------------------------------------------------------------------------------------- */
/* fixedaot() - compute aerosol reflectance for fixed aot(lambda)                           */
/*                                                                                          */
/* B. Franz, August 2004.                                                                   */
/* ---------------------------------------------------------------------------------------- */
int fixedaot(int32_t sensorID, float aot[], float wave[], int32_t nwave, int32_t iwnir_s, int32_t iwnir_l, 
             float solz, float senz, float phi, float wv, int32_t *modmin, int32_t *modmax, float *modrat, 
             float rhoa[],float *epsnir)
{
    static int firstCall = 1;
    static int angst_band1 = -1;
    static int angst_band2 = -1;

    float *phase1;
    float *phase2;
    float *f1   ;
    float *f2   ;
    float *lnf1 ;
    float *lnf2 ;
    float *rhoas1;
    float *rhoas2;
    float *rhoa1;
    float *rhoa2;
    float eps1;
    float eps2;
    float angstrom;
    int   iw, iwtab;
    float maxwave;

    maxwave=MAX(aertab->nwave,nwave);

    if ((rhoa1 = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for rhoa1.\n");
        exit(1);
    }
    if ((rhoa2 = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for rhoa2.\n");
        exit(1);
    }
    if ((rhoas1 = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for rhoas1.\n");
        exit(1);
    }
    if ((rhoas2 = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for rhoas2.\n");
        exit(1);
    }
    if ((f1 = (float *)calloc(maxwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for f1.\n");
        exit(1);
    }
    if ((f2 = (float *)calloc(maxwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for f2.\n");
        exit(1);
    }
    if ((lnf1 = (float *)calloc(maxwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for lnf1.\n");
        exit(1);
    }
    if ((lnf2 = (float *)calloc(maxwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for lnf2.\n");
        exit(1);
    }

    if (firstCall) {
        angst_band1 = windex(520,wave,nwave);
        angst_band2 = windex(865,wave,nwave);
        firstCall = 0;
    }

    /* bail on negative input aot */
    for (iw=0; iw<nwave; iw++) {
        if (aot[iw] < 0.0) {
            free(rhoa1);
            free(rhoa2);
            free(rhoas1);
            free(rhoas2);
            free(f1);
            free(f2);
            free(lnf1);
            free(lnf2);
            return(1);
        }
    }

    /* compute angstrom and use to select bounding models */
    if (aot[iwnir_l] > 0.0)
        angstrom = -log(aot[angst_band1]/aot[angst_band2])/
                    log(wave[angst_band1]/wave[angst_band2]);
    else
        angstrom = 0.0;

    model_select_angstrom(angstrom,modmin,modmax,modrat);


    /* get model phase function for all wavelengths at this geometry for the two models */
    phase1 = model_phase(*modmin,solz,senz,phi);
    phase2 = model_phase(*modmax,solz,senz,phi);

    /* compute factor for SS approximation, set-up for interpolation */
    for (iw=0; iw<aertab->nwave; iw++) {
        f1[iw] = aertab->model[*modmin]->albedo[iw]*phase1[iw]/4.0/mu0/mu;
        f2[iw] = aertab->model[*modmax]->albedo[iw]*phase2[iw]/4.0/mu0/mu;
        if (interp) {
	    lnf1[iw] = log(f1[iw]);
	    lnf2[iw] = log(f2[iw]);
	}
    }

    /* compute SS aerosol reflectance */
    if (interp) {
        for (iw=0; iw<nwave; iw++) {
	    iwtab = iwatab[iw];
	    if (aertab->wave[iwtab] != wave[iw] && wave[iw] > 0) {
	        rhoas1[iw] = aot[iw]*exp(linterp(aertab->wave,lnf1,aertab->nwave,wave[iw]));
	        rhoas2[iw] = aot[iw]*exp(linterp(aertab->wave,lnf2,aertab->nwave,wave[iw]));
	    } else {
                rhoas1[iw] = aot[iw]*f1[iwtab];
                rhoas2[iw] = aot[iw]*f2[iwtab];
	    }
	} 
    } else {
        for (iw=0; iw<nwave; iw++) {
	    iwtab = iwatab[iw];
            rhoas1[iw] = aot[iw]*f1[iwtab];
            rhoas2[iw] = aot[iw]*f2[iwtab];
	}
    }
    eps1 = rhoas1[iwnir_s]/rhoas1[iwnir_l];
    eps2 = rhoas2[iwnir_s]/rhoas2[iwnir_l];

    /* compute MS aerosol reflectance */
    rhoas_to_rhoa(sensorID,*modmin,solz,senz,phi,wv,rhoas1,wave,nwave,0,nwave-1,rhoa1);
    rhoas_to_rhoa(sensorID,*modmax,solz,senz,phi,wv,rhoas2,wave,nwave,0,nwave-1,rhoa2);

    for (iw=0; iw<nwave; iw++) {
        rhoa[iw] = (1.0-*modrat) * rhoa1[iw] + *modrat * rhoa2[iw];
    }
    *epsnir = (1.0-*modrat) * eps1 + *modrat * eps2;

    free(rhoa1);
    free(rhoa2);
    free(rhoas1);
    free(rhoas2);
    free(f1);
    free(f2);
    free(lnf1);
    free(lnf2);

    return(0);
}


/* ---------------------------------------------------------------------------------------- */
/* aerosol() - compute aerosol reflectance using specified algorithm                        */
/*                                                                                          */
/* B. Franz, 1 June 2004.                                                                   */
/* ---------------------------------------------------------------------------------------- */
int aerosol(l2str *l2rec, int32_t aer_opt_in, aestr *aerec, int32_t ip,
            float wave[], int32_t nwave, int32_t iwnir_s_in, int32_t iwnir_l_in, 
            float F0_in[], float La1_in[], float La2_out[], 
            float t_sol_out[], float t_sen_out[], float *eps, float taua_out[],
            int32_t *modmin, int32_t *modmax, float *modrat,
            int32_t *modmin2, int32_t *modmax2, float *modrat2)
{
    static int firstCall = 1;
    static int32_t mindx[MAXAERMOD];

    int   status = 1;
    float *rhoa;
    float *radref;
    float temp;
    float *aot;
    float angstrom;

    int iw;

    float *F0   ;
    float *taur ;
    float *La1  ;
    float *La2  ;
    float *t_sol;
    float *t_sen;
    float *taua ;
    float *taua_pred_min;
    float *taua_pred_max;

    int32_t  sensorID = l2rec->sensorID;
    float solz     = l2rec->solz  [ip];
    float senz     = l2rec->senz  [ip];
    float phi      = l2rec->delphi[ip];
    float wv       = l2rec->wv    [ip];
    float rh       = l2rec->rh    [ip];
    float pr       = l2rec->pr    [ip];

    instr *input = l2rec->input;

    if (firstCall == 1) {
        Nbands = nwave;
        Maxband = nwave + 1;    /* Must be >= Nbands */
    }

    if ((rhoa = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for rhoa.\n");
        exit(1);
    }
    if ((radref = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for raderef.\n");
        exit(1);
    }
    if ((F0 = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for F0.\n");
        exit(1);
    }
    if ((taur = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for taur.\n");
        exit(1);
    }
    if ((La1 = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for La1.\n");
        exit(1);
    }
    if ((La2 = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for La2.\n");
        exit(1);
    }
    if ((t_sol = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for t_sol.\n");
        exit(1);
    }
    if ((t_sen = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for t_sen.\n");
        exit(1);
    }
    if ((taua = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for taua.\n");
        exit(1);
    }
    if ((taua_pred_min = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for taua_pred_min.\n");
        exit(1);
    }
    if ((taua_pred_max = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("Unable to allocate space for taua_pred_max.\n");
        exit(1);
    }

    /* set static global evaluation level */  
    evalmask = input->evalmask;
    aer_opt  = aer_opt_in;

    /* transfer inputs per band to inputs per wavelength */
    for (iw=0; iw<nwave; iw++) {
 	F0  [iw] = F0_in [iw];
        taur[iw] = l2rec->Tau_r[iw];
        La1 [iw] = La1_in[iw];
        if (iw == iwnir_s_in) iwnir_s = iw;
        if (iw == iwnir_l_in) iwnir_l = iw;
    }

    /* compute total airmass (static global) */
    mu0     = cos(solz/radeg);
    mu      = cos(senz/radeg);
    airmass = 1.0/mu0 + 1.0/mu;
    if ((evalmask & TRANSSPHER) != 0) {
        airmass_plp = pp_airmass(solz) + pp_airmass(senz);
        airmass_sph = ky_airmass(solz) + ky_airmass(senz);
    }

    /* initialize epsilon and diffuse transmittance to defaults */
    *eps = 1.0;
    *modmin = BAD_INT;
    *modmax = BAD_INT;
    *modrat = BAD_FLT;
    *modmin2 = BAD_INT;
    *modmax2 = BAD_INT;
    *modrat2 = BAD_FLT;
    for (iw=0; iw<nwave; iw++) {
        taua  [iw] = 0.0;
        t_sol [iw] = 1.0;
        t_sen [iw] = 1.0;
        La2   [iw] = 0.0;
        radref[iw] = pi / F0[iw] / mu0;
    }
    
    /* load aerosol model tables */
    if (!loaded) {
        int32_t im;
        load_aermod(sensorID,wave,nwave,input->aermodfile,input->aermodels,input->naermodels);
        for (im=0; im<aertab->nmodel; im++) mindx[im] = im;
        if (have_rh && aertab->nmodel >= 30) {
	  printf("\nLimiting aerosol models based on RH.\n");
          use_rh = 1;
	}
    }

    /* Change the aerosol option if need be */
    if (use_rh)
      switch (aer_opt) {
        case AERWANG:     aer_opt = AERRH;  break;
        case AERWANGNIR:  aer_opt = AERRHNIR;  break;
        case AERWANGSWIR: aer_opt = AERRHSWIR; break;
        case AERMUMM:     aer_opt = AERRHMUMM; break;
      }


    if (firstCall) {
        if (aer_opt > 0 && aer_opt <= MAXAERMOD) {
            printf("\nUsing fixed aerosol model #%d (%s)\n",aer_opt,input->aermodels[aer_opt-1]);
	    printf("Extrapolating from %4.1f nm\n",wave[iwnir_l]);
        } else {
            switch (aer_opt) {
              case AERWHITE: 
                printf("\nUsing white-aerosol approximation\n");
   	        printf("Extrapolating from %4.1f nm\n",wave[iwnir_l]);
                break;
              case AERWANG: 
              case AERRH: 
                printf("\nUsing Gordon & Wang aerosol model selection\n");
   	        printf("Using bands at %4.1f and %4.1f nm for model selection\n",wave[iwnir_s],wave[iwnir_l]);
   	        printf("Extrapolating from %4.1f nm\n",wave[iwnir_l]);
                break;
              case AERWANGNIR: 
              case AERRHNIR: 
                printf("\nUsing Gordon & Wang aerosol model selection\n");
                printf("  and NIR correction with up to %d iterations\n",input->aer_iter_max);
   	        printf("Using bands at %4.1f and %4.1f nm for model selection\n",wave[iwnir_s],wave[iwnir_l]);
   	        printf("Extrapolating from %4.1f nm\n",wave[iwnir_l]);
                break;
              case AERWANGSWIR: 
              case AERRHSWIR: 
                printf("\nUsing Gordon & Wang aerosol model selection with NIR/SWIR switching.\n");
                printf("NIR correction with up to %d iterations\n",input->aer_iter_max);
   	        printf("NIR bands at %d and %d nm\n",input->aer_wave_short,input->aer_wave_long);
   	        printf("SWIR bands at %d and %d nm\n\n",input->aer_swir_short, input->aer_swir_long);
                break;
              case AERMUMM: 
              case AERRHMUMM: 
                printf("\nUsing Gordon & Wang aerosol model selection\n");
                printf("  and MUMM correction\n");
   	        printf("Using bands at %4.1f and %4.1f nm for model selection\n",wave[iwnir_s],wave[iwnir_l]);
   	        printf("Extrapolating from %4.1f nm\n",wave[iwnir_l]);
                break;
              case AERRHFRNIR: 
                printf("\nUsing Ahmad & Franz aerosol model selection.\n");
                printf("  and NIR correction with up to %d iterations\n",input->aer_iter_max);
   	        printf("Using bands at %4.1f and %4.1f nm for model selection\n",wave[iwnir_s],wave[iwnir_l]);
   	        printf("Extrapolating from %4.1f nm\n",wave[iwnir_l]);
                break;
              case AERRHMSEPS: 
                printf("\nUsing multi-scattering aerosol model selection.\n");
                printf("  and NIR correction with up to %d iterations\n",input->aer_iter_max);
   	        printf("Using bands at %4.1f and %4.1f nm for model selection\n",wave[iwnir_s],wave[iwnir_l]);
   	        printf("Extrapolating from %4.1f nm\n",wave[iwnir_l]);
                break;
              case FIXAOT: 
                printf("\nUsing fixed, input aerosol optical thickenesses for aerosol selection.\n");
                break;
              case FIXMODPAIR: 
                printf("\nUsing multiple scattering aerosols with fixed model pair\n");
                break;
              case FIXMODPAIRNIR: 
                printf("\nUsing multiple scattering aerosols with fixed model pair\n");
                printf("  and NIR iteration with up to %d iterations\n",input->aer_iter_max);
   	        printf("Extrapolating from %4.1f nm\n",wave[iwnir_l]);
                break;
              case FIXANGSTROM: 
                printf("\nUsing fixed aerosol model based on predefined Angstrom exponent\n");
   	        printf("Extrapolating from %4.1f nm\n",wave[iwnir_l]);
                break;
              case FIXANGSTROMNIR: 
                printf("\nUsing fixed aerosol model based on predefined Angstrom exponent\n");
                printf("  and NIR iteration with up to %d iterations\n",input->aer_iter_max);
   	        printf("Extrapolating from %4.1f nm\n",wave[iwnir_l]);
                break;
              default:
                printf("\nErroneous aerosol option, %d\n",aer_opt);
                exit(FATAL_ERROR);
                break;
	    }
	}
        firstCall = 0;
    }


    /* convert input aerosol radiances to relectance */
    for (iw=iwnir_s; iw<=iwnir_l; iw+=MAX(iwnir_l-iwnir_s,1))
        rhoa[iw] = La1[iw] * radref[iw];


    /* compute aerosol using requested method */
    /* -------------------------------------- */

    switch (aer_opt) {

      case AERWANG: case AERWANGNIR: case AERWANGSWIR: case AERMUMM:

        if (iwnir_l <= iwnir_s || wave[iwnir_s] < 600) {
  	    printf("Aerosol selection bands must be greater than 600nm with short wave less than long wave (%d,%d)\n",iwnir_l,iwnir_s);
            exit(1);
	}

        /* Multi-Scattering with Gordon & Wang Aerosol Selection */
        /* ----------------------------------------------------- */

        /* convert input NIR aerosol radiances to reflectance */
        for (iw=iwnir_s; iw<=iwnir_l; iw++) {
	    rhoa[iw] = La1[iw] * radref[iw];
        }

        /* require sufficient signal in two NIR bands */
	if  (rhoa[iwnir_s] > input->rhoamin && rhoa[iwnir_l] > input->rhoamin) {

	    /* require MS epsilon to be reasonable */
	    if (La1[iwnir_s]/La1[iwnir_l] > 0.1) {

	      status = wangaer(sensorID,wave,nwave,iwnir_s,iwnir_l,
                               aertab->nmodel,mindx,
                               solz,senz,phi,wv,rhoa,modmin,modmax,modrat,eps,taua,taua);  // this taua is not used //

              if (status == 0)
		  for (iw=0; iw<nwave; iw++)
    	    	        La2[iw] = rhoa[iw] / radref[iw];
	    }

	} else if (rhoa[iwnir_s] > -(input->rhoamin) && rhoa[iwnir_l] > -(input->rhoamin)) {

	    /* if input NIR is near zero, assume a small white aerosol signal */
            *eps    = 1.0;
            *modmin = aertab->nmodel;
            *modmax = aertab->nmodel;
            *modrat = 0.0;
            temp    = MAX(rhoa[iwnir_l],1e-6);
            for (iw=0; iw<nwave; iw++) {
	        rhoa[iw] = temp;
 	        La2 [iw] = rhoa[iw] / radref[iw];
	    }

  	    status = 0;

        } else {

	    /* if input NIR is very negative, fail the case */
  	    status = 1;
	}

	break;

      case AERRH:
      case AERRHNIR:
      case AERRHSWIR:
      case AERRHMUMM:
      case AERRHFRNIR:
      case AERRHMSEPS:

        if (iwnir_l <= iwnir_s || wave[iwnir_s] < 600) {
            printf("Aerosol selection bands must be greater than 600nm with short wave less than long wave");
            exit(1);
	}

        /* Multi-Scattering with RH-based selection              */
        /* ----------------------------------------------------- */

        /* convert input NIR aerosol radiances to relectance */
        for (iw=iwnir_s; iw<=iwnir_l; iw++) {
	    rhoa[iw] = La1[iw] * radref[iw];
        }

        /* require sufficient signal in two NIR bands */
	if  (rhoa[iwnir_s] > input->rhoamin && rhoa[iwnir_l] > input->rhoamin) {

	    /* require MS epsilon to be reasonable */
       	    if (rhoa[iwnir_s]/rhoa[iwnir_l] > 0.1 && rhoa[iwnir_s]/rhoa[iwnir_l] < 10.0) {

                  status = rhaer(sensorID,wave,nwave,iwnir_s,iwnir_l,
			     solz,senz,phi,wv,rh,pr,taur,rhoa,
				 modmin,modmax,modrat,modmin2,modmax2,modrat2,eps,taua,t_sol,t_sen);

              if (status == 0)
		  for (iw=0; iw<nwave; iw++)
    	    	        La2[iw] = rhoa[iw] / radref[iw];
	    }

	} else if (rhoa[iwnir_s] > -(input->rhoamin) && rhoa[iwnir_l] > -(input->rhoamin)) {

            /* if input NIR is near zero, assume a small white aerosol signal */
            *eps    = 1.0;
            *modmin = aertab->nmodel;
            *modmax = aertab->nmodel;
            *modrat = 0.0;
            *modmin2 = aertab->nmodel;
            *modmax2 = aertab->nmodel;
            *modrat2 = 0.0;
            temp    = MAX(rhoa[iwnir_l],1e-6);
            for (iw=0; iw<nwave; iw++) {
                rhoa[iw] = temp;
                La2 [iw] = rhoa[iw] / radref[iw];
            }

            diff_tran(sensorID,wave,nwave,iwnir_l,solz,senz,phi,wv,pr,taur,
		      *modmin,*modmax,*modrat,rhoa,taua,t_sol,t_sen,taua_pred_min,taua_pred_max,0);

  	    status = 0;

        } else {

	    /* if input NIR is very negative, fail the case */
  	    status = 1;
	}

	break;

     case AERWHITE:

        /* White Aerosol */
        /* ------------- */

        if ( La1[iwnir_l] > 0.0) {

            *eps    = 1.0;
            *modmin = 0;
            *modmax = 0;
            *modrat = 0.0;

            for (iw=0; iw<nwave; iw++) {
 	        La2 [iw] = *eps * F0[iw]/F0[iwnir_l]*La1[iwnir_l];
	        rhoa[iw] = La2[iw] * radref[iw];
	    }


            status = 0;
	}
        break;

      case FIXMODPAIR: case FIXMODPAIRNIR:

        /* Multi-Scattering with Fixed Model Pair */
        /* -------------------------------------- */

        if (aerec != NULL && aerec->mode == ON) {
            *modmin = aerec->mod_min[ip] - 1;
            *modmax = aerec->mod_max[ip] - 1;
            *modrat = aerec->mod_rat[ip];
	} else {
	    *modmin = input->aermodmin - 1;
	    *modmax = input->aermodmax - 1;
	    *modrat = input->aermodrat;
	}

        status = fixedmodpair(sensorID,input,wave,nwave,iwnir_s,iwnir_l,solz,senz,phi,wv, 
                          *modmin,*modmax,*modrat,rhoa,eps);
        if (status == 0) {
	    for (iw=0; iw<nwave; iw++)
    	        La2[iw] = rhoa[iw] / radref[iw];
        }

        break;

      case FIXANGSTROM: case FIXANGSTROMNIR:

        if (input->aer_angstrom > -2) {
	    angstrom = input->aer_angstrom;
	} else {
            angstrom = bin_climatology(input->aerbinfile,*(l2rec->day),l2rec->lon[ip],l2rec->lat[ip],"angstrom");
	}

        if (angstrom > -2) {

            model_select_angstrom(angstrom,modmin,modmax,modrat);

            status = fixedmodpair(sensorID,input,wave,nwave,iwnir_s,iwnir_l,solz,senz,phi,wv, 
                          *modmin,*modmax,*modrat,rhoa,eps);
	} else 
            status = 1;        

        if (status == 0) {
	    for (iw=0; iw<nwave; iw++)
    	        La2[iw] = rhoa[iw] / radref[iw];
        }

        break;

      case FIXAOT:

        /* Multi-Scattering with fixed AOT */
        /* ------------------------------- */
        if (aerec != NULL && aerec->mode == ON)
	    aot = &aerec->taua[ip*Nbands];
        else
	    aot = input->taua;

        status = fixedaot(sensorID,aot,wave,nwave,iwnir_s,iwnir_l,solz,senz,phi,wv, 
                          modmin,modmax,modrat,rhoa,eps);
        if (status == 0) {
            for (iw=0; iw<nwave; iw++)
    	        La2[iw] = rhoa[iw] / radref[iw];
        }

        break;

      default:

        /* Multi-Scattering with Fixed Model */
        /* --------------------------------- */

        *modmin = aer_opt-1;
        *modmax = aer_opt-1;
        *modrat = 0.0;

        if (*modmin < 0 || *modmin > input->naermodels-1) {
	    printf("Invalid aerosol option: %d\n",*modmin);
	    exit(1);
	}

        /* convert input NIR aerosol radiance to relectance */
	rhoa[iwnir_l] = La1[iwnir_l]*radref[iwnir_l];

        /* get aerosol reflectance in visible for fixed model, and epsilon */
        status = fixedaer(sensorID,*modmin,wave,nwave,iwnir_s,iwnir_l,
                          input->aermodels,input->naermodels,
                          solz,senz,phi,wv,rhoa,eps);

        /* convert aerosol relectance to radiance */
        if (status == 0)
 	    for (iw=0; iw<nwave; iw++)
 	        La2[iw]  = rhoa[iw] / radref[iw];
	break;
    }

    /* compute diffuse transmittance through aerosol and Rayleigh, if not yet computed */
    if (status == 0 && aer_opt != AERRHNIR && aer_opt != AERRHFRNIR && aer_opt != AERRHSWIR && aer_opt != AERRH && aer_opt != AERRHMSEPS) {

        diff_tran(sensorID,wave,nwave,iwnir_l,solz,senz,phi,wv,pr,taur,
                  *modmin,*modmax,*modrat,rhoa,taua,t_sol,t_sen,taua_pred_min,taua_pred_max,0);
    }

    /* transfer outputs per wavelength to outputs per band */
    for (iw=0; iw<nwave; iw++) {
        La2_out  [iw] = La2  [iw];
        taua_out [iw] = taua [iw];
        t_sol_out[iw] = t_sol[iw];
        t_sen_out[iw] = t_sen[iw];
    }

    /* model numbers are reported as 1-based */
    *modmin  = *modmin+1;
    *modmax  = *modmax+1;
    *modmin2 = *modmin2+1;
    *modmax2 = *modmax2+1;

    free(rhoa);
    free(radref);
    free(F0);
    free(taur);
    free(La1);
    free(La2);
    free(t_sol);
    free(t_sen);
    free(taua);
    free(taua_pred_min);
    free(taua_pred_max);

    return (status);
}



/* --------------------------------------------------------------- */
/* get_angstrom.c - compute angstrom coefficient                   */
/*                                                                 */
/* Inputs:                                                         */
/*     l2rec - level-2 structure containing one complete scan      */
/*             after atmospheric correction.                       */
/*     band  = band number 0-7                                     */
/* Outputs:                                                        */
/*     angst - angstrom coefficient                                */
/*                                                                 */
/* Written By: B. A. Franz, SAIC, August 2004                      */
/*                                                                 */
/* --------------------------------------------------------------- */
void get_angstrom(l2str *l2rec, int band, float angst[])
{
    static int   firstCall = 1;
    static int32_t  ib2;
    static float wave2;

    int32_t  ip;
    int32_t  ib1;
    float wave1;
    float aot1;
    float aot2;

    if (firstCall) {
        ib2 = windex(865.0,l2rec->fwave,l2rec->nbands);
        wave2 = l2rec->fwave[ib2];
        firstCall = 0;
    }

    if (band < 0)
        ib1 = windex(443.0,l2rec->fwave,l2rec->nbands);
    else
        ib1 = band;
    wave1 = l2rec->fwave[ib1];

    for (ip=0; ip<l2rec->npix; ip++) {
        aot1 = l2rec->taua[ip*l2rec->nbands+ib1];
        aot2 = l2rec->taua[ip*l2rec->nbands+ib2];
        if (aot1 > 0.0 && aot2 > 0.0) 
            angst[ip] = -log(aot1/aot2)/log(wave1/wave2);
        else if (aot1 == 0.0 && aot2 == 0.0) 
	    angst[ip] = 0.0;
        else {
	    angst[ip] = BAD_FLT;
            l2rec->flags[ip] |= PRODFAIL;
        }
    }

    return;
}


/* --------------------------------------------------------------- */
/* get_ms_epsilon.c -                                              */
/* --------------------------------------------------------------- */
void get_ms_epsilon(l2str *l2rec, float eps[])
{
  int32_t  ip, ipb1, ipb2;

    for (ip=0; ip<l2rec->npix; ip++) {
      ipb1 = ip*l2rec->nbands+iwnir_s;
      ipb2 = ip*l2rec->nbands+iwnir_l;
      if (l2rec->La[ipb2] > 0.0) {
          eps[ip] = l2rec->La[ipb1]/l2rec->La[ipb2]
	          * l2rec->Fo[iwnir_l]/l2rec->Fo[iwnir_s];
      } else {
	  /* "should" indicate ATMFAIL, but just to be safe */
          eps[ip] = BAD_FLT;
          l2rec->flags[ip] |= PRODFAIL;
       }
    }

    return;
}
