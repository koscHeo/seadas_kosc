/* --------------------------------------------------------------------------------------- */
/* Loisel & Stramski IOP model                                                             */
/*                                                                                         */
/* References:                                                                             */
/*                                                                                         */
/*    H. Loisel, J.M. Nicolas, A. Sciandra, D. Stramski, and A. Poteau. 2006. The          */
/*    spectral dependency of optical backscattering by marine particles from satellite     */
/*    remote sensing of the global ocean. Journal of Geophysical Research, 111, C09024,    */
/*    doi: 10.1029/2005JC003367.                                                           */
/*                                                                                         */
/*    H. Loisel, D. Stramski, B. G. Mitchell, F. Fell, V. Fournier-Sicre, B. Lemasle       */
/*    and M. Babin. 2001. Comparison of the ocean inherent optical properties obtained     */
/*    from measurements and inverse modeling. Applied Optics.40: 2384-2397.                */
/*                                                                                         */
/*    H. Loisel, E. Bosc D. Stramski, K. Oubelker and P-Y. Deschamps. 2001. Seasonal       */
/*    variability of the backscattering coefficients in the Mediterranean Sea based        */
/*    on Satellite SeaWIFS imagery. Geophysical Research letters. 28: 4203-4206.           */
/*                                                                                         */
/*    H. Loisel, J.M. Nicolas, P.Y.Deschamps, and R. Frouin. 2002. Seasonal and            */
/*    inter-annual variability of the particulate matter in the global ocean.              */
/*    Geophysical Research letters 29: 2996.                                               */
/*                                                                                         */
/*    H. Loisel and D. Stramski. 2000. Estimation of the inherent optical properties       */
/*    of natural waters from irradiance attenuation coefficient and reflectance in         */
/*    the presence of Raman scattering. Applied Optics. 39: 3001-3011.                     */
/*                                                                                         */
/*    Fortran code provided by H. Loisel, May 2008.                                        */
/*                                                                                         */
/* Implementation:                                                                         */
/*                                                                                         */
/*    B. Franz, NASA/OBPG, June 2008.                                                      */
/* --------------------------------------------------------------------------------------- */

#include <stdlib.h>
#include <math.h>
#include "l12_proto.h"
#include "l2prod.h"
#include "amoeba.h"

#define LASNSOL 3
#define LASNKC  3
#define LASNRC  2

typedef float m_array[LASNSOL][LASNRC];

typedef struct las_table_struc {
  int nwav;
  int nsol;
  int nrc;
  int nkc;
  float *wav;
  float sol[LASNSOL];
  m_array *rc ;
  m_array *kc ;
} lastabstr;

static lastabstr tab;        // R(0-) and <Kd>1 coefficient table 

static int LastRecNum = -1;
static float badData  = BAD_FLT;

static int nwave      = -1;  // number of wavelengths to fit
static float *wave;          // sensor wavelengths to fit
static float *a;             // total absorption coefficient per band and pixel           
static float *b;             // total scattering coefficient  per band and pixel              
static float *c;             // beam attenuation coefficient per band and pixel               
static float *bb;            // backscatter coefficient per band and pixel                  
static float *bbp;           // particulate backscatter coefficient per band and pixel
static float *aw;            // pure-water total absorption  per band           
static float *bbw;           // pure-water backscattering per band              



/* ------------------------------------------------------------------------------- */
/* read_las_tables() - load Loisel & Stramski R0 and Kd coefficient tables         */ 
/* ------------------------------------------------------------------------------- */
void read_las_tables(int sensorID)
{
    FILE *fp;
    char *tmp_str;
    char fname[FILENAME_MAX];
    char sensor_name[FILENAME_MAX];
    int iwav, isol;
    static int firstCall = 1;


    tab.nsol = LASNSOL;
    tab.nkc  = LASNKC;
    tab.nrc  = LASNRC;

    if ((tmp_str = getenv("OCDATAROOT")) == NULL) {
       printf("OCDATAROOT environment variable is not defined.\n");
       exit(1);
    }

    /* sensor specific LUTs used - only available SeaWiFS */
    switch (sensorID){
      case SEAWIFS:
         sprintf(sensor_name,"seawifs");
         tab.nwav = 5;
         break;
      default: /* interpolate from seawifs wavelengths */
         sprintf(sensor_name,"seawifs");
         tab.nwav = 5;
         break;
    }

    if (tab.rc == NULL) {
        if ( (tab.rc = (m_array *)calloc(tab.nwav,sizeof(m_array))) == NULL) {
            printf("-E- : Error allocating memory to read_las_tables\n");
            exit(FATAL_ERROR);
        }
        if ( (tab.kc = (m_array *)calloc(tab.nwav,sizeof(m_array))) == NULL) {
            printf("-E- : Error allocating memory to read_las_tables\n");
            exit(FATAL_ERROR);
        }
        if ( (tab.wav = (float *)calloc(tab.nwav,sizeof(float))) == NULL) {
            printf("-E- : Error allocating memory to read_las_tables\n");
            exit(FATAL_ERROR);
        }
    }
    printf("Loading Loisel & Stramski IOP model tables:\n");

    sprintf(fname,"%s/%s/iop/las/%s",tmp_str,sensor_name,"LUT_RRS");
    printf("  %s\n",fname);
    if((fp = fopen(fname,"r")) == NULL){
        printf("Error opening %s\n", fname);
        exit(1);
    } 
    for (isol=0; isol<tab.nsol; isol++) for (iwav=0; iwav<tab.nwav; iwav++)  
        fscanf(fp,"%f %f %f %f",&tab.wav[iwav],&tab.sol[isol],&tab.rc[iwav][isol][0],&tab.rc[iwav][isol][1]);
    fclose(fp);

    sprintf(fname,"%s/%s/iop/las/%s",tmp_str,sensor_name,"LUT_KD");
    printf("  %s\n",fname);
    if((fp = fopen(fname,"r")) == NULL){
        printf("Error opening %s\n", fname);
        exit(1);
    }
    for (isol=0; isol<tab.nsol; isol++) for (iwav=0; iwav<tab.nwav; iwav++)  
        fscanf(fp,"%f %f %f %f %f",&tab.wav[iwav],
            &tab.sol[isol],&tab.kc[iwav][isol][0],&tab.kc[iwav][isol][1],&tab.kc[iwav][isol][2]);
    fclose(fp);
}


/* ------------------------------------------------------------------------------- */
/* computes R(0-) and <Kd>1 at one wavelength from Rrs and Rrs(443)/Rrs(555)       */ 
/* ------------------------------------------------------------------------------- */
void get_R0_Kd(float rat, float Rrs, float wave, float solz, float *R0, float *Kd)
{
    float log_q = log10(rat);
    int   iwav  = windex(wave,tab.wav,tab.nwav);
    float dwav  = wave-tab.wav[iwav];

    float R01, R02;
    float Kd1, Kd2;

    if (fabs(dwav) > 1.0) {

        // if the input wavelength is not one of the table wavelengths, call the 
        // function recursively for bounding table wavelengths and interpolate

        int iwav1, iwav2;
        if (dwav > 0.0) {
            iwav2 = MIN(iwav+1,tab.nwav-1);
	    iwav1 = iwav2-1;
	} else {
	    iwav1 = MAX(iwav-1,0);
            iwav2 = iwav1+1;
	}
        get_R0_Kd(rat,Rrs,tab.wav[iwav1],solz,&R01,&Kd1);
        get_R0_Kd(rat,Rrs,tab.wav[iwav2],solz,&R02,&Kd2);
        *R0 = R01 + (wave-tab.wav[iwav1])/(tab.wav[iwav2]-tab.wav[iwav1])*(R02-R01);
        *Kd = Kd1 + (wave-tab.wav[iwav1])/(tab.wav[iwav2]-tab.wav[iwav1])*(Kd2-Kd1);

    } else {

        // compute R(0-) and <Kd>1 for bounding solar zenith angles of the 
        // coefficient table and interpolate to input solar zenith angle

        int   isol1, isol2;
        float solz1, solz2;

        isol2 = 0;
        while (tab.sol[isol2] < solz) isol2++;
        isol2 = MAX(MIN(isol2,tab.nsol-1),1);
        isol1 = isol2-1;

        solz1 = tab.sol[isol1];
        solz2 = tab.sol[isol2];

        R01 = tab.rc[iwav][isol1][0] * pow(Rrs,tab.rc[iwav][isol1][1]);
        R02 = tab.rc[iwav][isol2][0] * pow(Rrs,tab.rc[iwav][isol2][1]);

        Kd1 = pow(10.0,(tab.kc[iwav][isol1][0]*log_q+tab.kc[iwav][isol1][1])/(tab.kc[iwav][isol1][2]+log_q));    
        Kd2 = pow(10.0,(tab.kc[iwav][isol2][0]*log_q+tab.kc[iwav][isol2][1])/(tab.kc[iwav][isol2][2]+log_q));

        *R0 = R01 + (solz-solz1)*(R02-R01)/(solz2-solz1);
        *Kd = Kd1 + (solz-solz1)*(Kd2-Kd1)/(solz2-solz1);
    }

} 


/* ------------------------------------------------------------------------------- */
/* have we run for this scan line?                                                 */
/* ------------------------------------------------------------------------------- */
int las_ran(int recnum)
{                                                                                
    if ( recnum == LastRecNum )
        return 1;
    else
        return 0;                                                                              
}


/* ------------------------------------------------------------------------------- */
/* allocates private arrays                                                        */
/* ------------------------------------------------------------------------------- */
void alloc_las(int npix, int nbands)
{
    a   = (float*) calloc(npix*nbands,sizeof(float));
    b   = (float*) calloc(npix*nbands,sizeof(float));
    c   = (float*) calloc(npix*nbands,sizeof(float));
    bb  = (float*) calloc(npix*nbands,sizeof(float));
    bbp = (float*) calloc(npix*nbands,sizeof(float));
    aw  = (float*) calloc(nbands,sizeof(float));
    bbw = (float*) calloc(nbands,sizeof(float));
}


/* ------------------------------------------------------------------------------- */
/* computes eta = bw/b                                                             */
/* ------------------------------------------------------------------------------- */
float eta_func(float bbw, float bb)
{
    float bw = bbw*2.0;
    float bbp, b;
    
    if (bb <= bbw)
        bbp = bbw/10.0; 
    else 
        bbp = bb-bbw;

    b = bw + bbp/0.0183;

    return(bw/b);
}


/* ------------------------------------------------------------------------------- */
/* computes alpha exponent of bb function                                          */
/* ------------------------------------------------------------------------------- */
float alpha_func(float solz, float eta)
{
    float a = -0.0007*solz+0.1024;
    float b = -0.0042*solz+0.3594;
    return(a*log10(eta)+b);
}      


/* ------------------------------------------------------------------------------- */
/* computes delta exponent of bb function                                          */
/* ------------------------------------------------------------------------------- */
float delta_func(float eta)
{
    return(0.871+0.40*eta-1.83*pow(eta,2.0));   
}


/* ------------------------------------------------------------------------------- */
/* run the model, collect the output                                               */
/* ------------------------------------------------------------------------------- */
void run_las(l2str *l2rec)
{
    static int firstCall = 1;
    static int ib443 = -1;
    static int ib555 = -1;
    static int maxiter = 4;
    static float nw  = 1.334;
    float *Rrs, ratio;
    float R0, Kd;
    float solz, mu;
    float alpha, delta, eta;
    float h, m, minh, maxh, ec50;
    int status;
    int ip,ib,ipb,it;

    if (firstCall) {

        firstCall = 0;

        wave  = l2rec->fwave;
        ib443 = windex(443,wave,l2rec->nbands);
        ib555 = windex(555,wave,l2rec->nbands);
        nwave = ib555+1; 

        alloc_las(l2rec->npix,nwave);
        read_las_tables(l2rec->sensorID);
        get_aw_bbw(l2rec,wave,nwave,aw,bbw);
    }

    if ((Rrs = (float *)calloc(nwave,sizeof(float))) == NULL) {
        printf("-E- : Error allocating memory to run_las\n");
        exit(FATAL_ERROR);
    }

    for (ip=0; ip<l2rec->npix; ip++) {

        status = 0;

        // check for negative radiances

        for (ib=0; ib<nwave; ib++) {
            ipb = ip*nwave+ib;
            Rrs[ib] = l2rec->Rrs[ipb];
            if (Rrs[ib] <= 0.0) status=1;
	}

        // if no negatives and pixel not already masked, run model   
     
        if ( status == 0 && !l2rec->mask[ip] ) {

	    if ((l2rec->input->brdf_opt & FOQMOREL) != 0) 
              solz = 0.0;
            else
	      solz  = MIN(l2rec->solz[ip],60.0);
            ratio = Rrs[ib443]/Rrs[ib555];
            mu    = cos(asin(sin(solz/RADEG)/nw));

            for (ib=0; ib<nwave; ib++) {

  	        ipb = ip*nwave+ib;

                get_R0_Kd(ratio, Rrs[ib], wave[ib], solz, &R0, &Kd);

                eta = 0.05;  
                for (it=0; it<maxiter; it++) {
                    alpha   = alpha_func(solz,eta);      
                    delta   = delta_func(eta);   
                    bb[ipb] = Kd*pow(10,alpha)*pow(R0,delta);
		    eta     = eta_func(bbw[ib],bb[ipb]);
		}
        
                minh = -0.0034*solz+1.248;
                maxh = -0.0084*solz+1.7223;
                ec50 =  0.0033*solz-1.9837;

                m = 0.0034*pow(solz,2)-0.0674*solz-9.34;
                h = pow(10.0,minh+(maxh-minh)/(1+pow(10,(ec50-log10(eta))*m)));

                a[ipb] = mu*Kd/sqrt(1+h*R0/(1-R0));
                b[ipb] = bbw[ib]*2.0/eta;
                c[ipb] = a[ipb] + b[ipb];
                bbp[ipb] = bb[ipb] - bbw[ib];

                if ( !finite(a[ipb]) || !finite(bb[ipb]) ) {
		    a  [ipb] = badData;
		    b  [ipb] = badData;
		    c  [ipb] = badData;
		    bb [ipb] = badData;
		    bbp[ipb] = badData;
		}
	    }

        } else {

            for (ib=0; ib<nwave; ib++) {
                ipb = ip*nwave+ib;
                a  [ipb]  = badData; 
                b  [ipb]  = badData; 
                c  [ipb]  = badData; 
                bb [ipb]  = badData;
                bbp[ipb]  = badData;
            }
	}        
    }

    LastRecNum = l2rec->iscan; 
    free(Rrs);
    return;
}


/* ------------------------------------------------------------------- */
/*                                                                     */
/* ------------------------------------------------------------------- */
double las_eta_amb(FITSTRUCT *ambdata, double par[])
{
  int iw;

  ambdata->merit = 0.0;
  for (iw=0; iw < nwave; iw++) {
      ambdata->yfit[iw] = par[0]*pow(wave[0]/wave[iw],par[1]);
      ambdata->merit += pow((ambdata->y[iw] - ambdata->yfit[iw]), 2) * ambdata->wgt[iw];
  }

  return(ambdata->merit);
}


/* ------------------------------------------------------------------- */
/*                                                                     */
/* ------------------------------------------------------------------- */
#define NPARETALAS 2
float fit_las_eta_amb(float *bbp)
{
  static int firstCall = 1;
  static float tol = 1.e-6;             /* fractional change in chisqr */
  static FITSTRUCT ambdata;             /* amoeba interface structure  */
  static double init[NPARETALAS*(NPARETALAS+1)]; /* initial simplex    */
  static double par [NPARETALAS];
  static double len [NPARETALAS];
  static double *wts ;
  static double *y   ;
  static double *yfit;
  static int npar = NPARETALAS;   
  static int maxiter = 100;

  int   i, j;
  short isml;
  int status = 0;

  if (firstCall) {
      firstCall = 0;
      if ( (wts = (double *)calloc(nwave,sizeof(double))) == NULL) {
          printf("-E- : Error allocating memory to fit_las_eta_amb\n");
          exit(FATAL_ERROR);
      }
      if ( (y = (double *)calloc(nwave,sizeof(double))) == NULL) {
          printf("-E- : Error allocating memory to fit_las_eta_amb\n");
          exit(FATAL_ERROR);
      }
      if ( (yfit = (double *)calloc(nwave,sizeof(double))) == NULL) {
          printf("-E- : Error allocating memory to fit_las_eta_amb\n");
          exit(FATAL_ERROR);
      }
      for (i=0; i<nwave; i++) {
	  wts[i] = 1.0;
      }
      ambdata.nfunc =  npar;      /* number of model parameters */
      ambdata.npnts = nwave;      /* number of wavelengths      */
      ambdata.wgt   =   wts;      /* Input weights on values    */
      ambdata.meta  =  NULL;

  }

  for (i=0; i<nwave; i++) {
      y[i] = (double) bbp[i];
  }

  ambdata.y     =     y;      /* Input values                  */
  ambdata.yfit  =  yfit;      /* Output model predicted values */

  par[0] = y[0]; len[0] = 0.01;
  par[1] = 0.0;  len[1] = 2.0;

  /* initialize simplex with first guess model parameters */
  for (j=0; j<npar+1; j++) 
      for (i=0; i<npar; i++) 
          init[j*npar+i] = par[i];

  for (i=0; i<npar; i++) {
      init[npar+i*(npar+1)] += len[i];
      par[i] = 0.0;
  }

  /* run optimization */
  isml = amoeba(init, &ambdata, las_eta_amb, tol); 

  for (i = 0; i < npar; i++) {
      par[i] = init[npar * isml + i];
  }

  /* check convergence and record parameter results */
  if (ambdata.niter >= maxiter || par[0] == BAD_FLT || !finite(par[1]))
      return(BAD_FLT);
  else
      return(par[1]);
}




/* ------------------------------------------------------------------------------- */
/* interface to l2_hdf_generic()                                                   */
/* ------------------------------------------------------------------------------- */
void get_las(l2str *l2rec, l2prodstr *p, float prod[])
{
    int   prodID = p->cat_ix;
    int   ib  = p->prod_ix;
    int   ip, ipb;

    if ( !las_ran(l2rec->iscan) )
        run_las(l2rec);

    for (ip=0; ip<l2rec->npix; ip++) {

        if (ib >= nwave) { 
	    prod[ip] = p->badData;
            l2rec->flags[ip] |= PRODFAIL;
            continue;
	}

        ipb = ip*nwave+ib;

        switch (prodID) {

        case CAT_a_las : 
            if ( a[ipb] > badData )
                prod[ip] = a[ipb];
            else {
                prod[ip] = p->badData;
                l2rec->flags[ip] |= PRODFAIL;
	    }
            break;

        case CAT_b_las : 
            if ( b[ipb] > badData )
                prod[ip] = b[ipb];
            else {
                prod[ip] = p->badData;
                l2rec->flags[ip] |= PRODFAIL;
	    }
            break;

        case CAT_c_las : 
            if ( c[ipb] > badData )
                prod[ip] = c[ipb];
            else {
                prod[ip] = p->badData;
                l2rec->flags[ip] |= PRODFAIL;
	    }
            break;

        case CAT_bb_las : 
	    if ( bb[ipb] > badData )
                prod[ip] = bb[ipb];
	    else {
                prod[ip] = p->badData;
                l2rec->flags[ip] |= PRODFAIL;
	    }
            break;

        case CAT_bbp_las : 
	    if ( bbp[ipb] > badData )
                prod[ip] = bbp[ipb];
	    else {
                prod[ip] = p->badData;
                l2rec->flags[ip] |= PRODFAIL;
	    }
            break;

        case CAT_bbps_las : 
	    prod[ip] = get_bbp_las_eta(l2rec,ip);
	    if (prod[ip] == BAD_FLT) {
                prod[ip] = p->badData;
                l2rec->flags[ip] |= PRODFAIL;
	    }
            break;

        default:
            printf("-E- %s line %d : erroneous product ID %d passed to Loisel & Stramski model().\n",
                __FILE__,__LINE__,prodID);
            exit(1);
        }
    }

    return;
}



/* ------------------------------------------------------------------------------- */
/* interface to convl12() to return iops                                           */
/* ------------------------------------------------------------------------------- */
void iops_las(l2str *l2rec)
{
    int ib, ip;

    if ( !las_ran(l2rec->iscan) )
        run_las(l2rec);

    for (ip=0; ip<l2rec->npix; ip++) {
        for (ib=0; ib<nwave; ib++) {
            l2rec->a [ip*nwave+ib] = a [ip*nwave+ib];
            l2rec->bb[ip*nwave+ib] = bb[ip*nwave+ib];
        }
        for (ib=nwave; ib<l2rec->nbands; ib++) {
            l2rec->a [ip*nwave+ib] = badData;
            l2rec->bb[ip*nwave+ib] = badData;
        }
    }

    return;
}


/* ------------------------------------------------------------------------------- */
/* interface to giop()                                                             */
/* ------------------------------------------------------------------------------- */
int get_bbp_las(l2str *l2rec, int ip, float tab_wave[], float tab_bbp[], int tab_nwave)
{
    int   ipb, iw;

    if ( !las_ran(l2rec->iscan) )
        run_las(l2rec);

    for (iw=0; iw<nwave; iw++) {
        ipb = ip*nwave+iw;
        if (bbp[ipb] < 0)
	    return(0);
    }

    ipb = ip*nwave;
    for (iw=0; iw<tab_nwave; iw++) {
        tab_bbp[iw] = linterp(wave,&bbp[ipb],nwave,tab_wave[iw]);
    }

    return(1);
}


/* ------------------------------------------------------------------------------- */
/* interface to giop()                                                             */
/* ------------------------------------------------------------------------------- */
float get_bbp_las_eta(l2str *l2rec, int ip)
{
    int   ipb, iw;

    if ( !las_ran(l2rec->iscan) )
        run_las(l2rec);

    for (iw=0; iw<nwave; iw++) {
        ipb = ip*nwave+iw;
        if (bbp[ipb] < 0)
	  return(BAD_FLT);
    }

    ipb = ip*nwave;
    return(fit_las_eta_amb(&bbp[ipb]));
}


