#include <stdlib.h>
#include <math.h>
#include "l12_proto.h"
#include "opp.h"

/**
 * Calculate the Primary Productivity using 1 of 3 algorithms
 *  1) Behrenfeld-Falkowski: (BeFa)
 *  2) Behrenfeld-Falkowski  algorithm, but
             modifies the pb_opt function after Eppley (as
             implemented by Antoine and Morel): (Eppley)
    3) Primary Productivity using a chl:Carbon ratio.
                     This is a spectrally resolved version of the cbpm, using nine separate
                     wavelengths: (Updated CbPM)

      Adapted from Oregon State U.
      http://www.science.oregonstate.edu/ocean.productivity/

      by R. Healy at NASA
      January 2015
 */
void get_opp(l2str *l2rec, int prodnum, float prod[])
{
    static int32_t  ib440;
    int32_t ip;
    float *kd, *par;
    static float32 *parin,*lat,*lon;
    static int nlat,nlon;
    float sst, chl, trise, tset, hrl, mld, bbp,k490,zno3, irr;
    static int firstCall = 1, havefile;
    char *parfile = l2rec->input->parfile;


    char name[H4_MAX_NC_NAME];
    char sdsname[H4_MAX_NC_NAME];
    int ncid, grpid, ndims, nvars, ngatts, unlimdimid;
    int32 sd_id;
    int32 sds_id;
    int32 rank;
    int32 sds_index;
    int32 nt;
    int32 dims[H4_MAX_VAR_DIMS];
    int32 nattrs;
    int status;
    nc_type rh_type;                 /* variable type */
    int dimids[H4_MAX_VAR_DIMS];    /* dimension IDs */
    int natts;                      /* number of attributes */
    int nsz;                 /* number of dims */
    size_t length;
    float kdtmp;

    if (firstCall && strcmp(parfile,"") != 0) {
            /* try netCDF first */
            if (nc_open(parfile, NC_NOWRITE, &ncid) == 0) {


                strcpy(sdsname,"lat");

                status = nc_inq_varid(ncid, sdsname, &sds_id);
                if (status != NC_NOERR) {
                    fprintf(stderr,"-E- %s line %d:  Error reading %s from %s.\n",
                            __FILE__,__LINE__,sdsname,parfile);
                    exit(1);
                }

                status = nc_inq_var (ncid, sds_id, 0, &rh_type, &ndims, dimids,
                                          &natts);

                if (nc_inq_dimlen(ncid, dimids[0], &length) != NC_NOERR) {
                    char name[H4_MAX_NC_NAME];
                    nc_inq_dim(ncid, dimids[0], name, &length);
                    fprintf(stderr,
                            "-E- %s line %d: could not get size of demension \"%s\" in netCDF File.\n",
                            __FILE__, __LINE__, name);
                    exit(1);
                }

                nlat = length;

                if ( (lat = (float *)calloc(nlat,sizeof(float))) == NULL) {
                    printf("-E- : Error allocating memory to tindx\n");
                    exit(FATAL_ERROR);
                }

                if (nc_get_var(ncid, sds_id, lat) !=0){
                    fprintf(stderr,"-E- %s line %d:  Error reading %s from %s.\n",
                            __FILE__,__LINE__,sdsname,parfile);
                    exit(1);
                }

                strcpy(sdsname,"lon");

                status = nc_inq_varid(ncid, sdsname, &sds_id);
                if (status != NC_NOERR) {
                    fprintf(stderr,"-E- %s line %d:  Error reading %s from %s.\n",
                            __FILE__,__LINE__,sdsname,parfile);
                    exit(1);
                }

                status = nc_inq_var (ncid, sds_id, 0, &rh_type, &ndims, dimids,
                                          &natts);
                if (status != NC_NOERR) {
                    fprintf(stderr,"-E- %s line %d:  Error reading %s from %s.\n",
                            __FILE__,__LINE__,sdsname,parfile);
                    exit(1);
                }

                if (nc_inq_dimlen(ncid, dimids[0], &length) != NC_NOERR) {
                    char name[H4_MAX_NC_NAME];
                    nc_inq_dim(ncid, dimids[0], name, &length);
                    fprintf(stderr,
                            "-E- %s line %d: could not get size of demension \"%s\" in netCDF File.\n",
                            __FILE__, __LINE__, name);
                    exit(1);
                }

                nlon = length;

                if ( (lon = (float *)calloc(nlon,sizeof(float))) == NULL) {
                    printf("-E- : Error allocating memory to tindx\n");
                    exit(FATAL_ERROR);
                }

                if (nc_get_var(ncid, sds_id, lon) !=0){
                    fprintf(stderr,"-E- %s line %d:  Error reading %s from %s.\n",
                            __FILE__,__LINE__,sdsname,parfile);
                    exit(1);
                }

                strcpy(sdsname,"par");

                status = nc_inq_varid(ncid, sdsname, &sds_id);
                if (status != NC_NOERR) {
                    fprintf(stderr,"-E- %s line %d:  Error reading %s from %s.\n",
                            __FILE__,__LINE__,sdsname,parfile);
                    exit(1);
                }

                status = nc_inq_var (ncid, sds_id, 0, &rh_type, &ndims, dimids,
                                          &natts);
                if (status != NC_NOERR) {
                    fprintf(stderr,"-E- %s line %d:  Error reading %s from %s.\n",
                            __FILE__,__LINE__,sdsname,parfile);
                    exit(1);
                }
                if (ndims !=2 ) {
                    fprintf(stderr,"-E- %s line %d:  Wrong number of dimensions for %s.  Need 2 got %d.\n",
                            __FILE__,__LINE__,sdsname,ndims);
                    exit(1);

                }
                printf("PARFILE: %s nlat=%d nlon=%d\n",parfile,nlat,nlon);

                if ( (parin = (float *)calloc(nlat*nlon,sizeof(float))) == NULL) {
                    printf("-E- : Error allocating memory to tindx\n");
                    exit(FATAL_ERROR);
                }

                if (nc_get_var(ncid, sds_id, parin) !=0){
                    fprintf(stderr,"-E- %s line %d:  Error reading %s from %s.\n",
                            __FILE__,__LINE__,sdsname,parfile);
                    exit(1);
                }

                havefile = 1;

        } else {
            fprintf(stderr,"-E- %s line %d:  Error opening parfile = %s.\n",
                    __FILE__,__LINE__,parfile);
            exit(1);
        }
    }

    if ( (kd = (float *)calloc(l2rec->npix,sizeof(float))) == NULL) {
        printf("-E- : Error allocating memory to par in get_opp\n");
        exit(FATAL_ERROR);
    }
    if ( (par = (float *)calloc(l2rec->npix,sizeof(float))) == NULL) {
        printf("-E- : Error allocating memory to par in get_opp\n");
        exit(FATAL_ERROR);
    }

    if (prodnum == CAT_opp_cbpm2) {
        if (l2rec->input->iop_opt == IOPNONE) {
            printf("IOP-based Kd_lee product requires iop model selection (iop_opt).  ");
            printf("Using default model.\n");
            l2rec->input->iop_opt = IOPDEFAULT;
            get_iops(l2rec,l2rec->input->iop_opt);
        }
        Kd490_obpg(l2rec,kd);
    }

    /* get irradiance E/D/m^2*/

    if (havefile)
        get_par_clim(parin,lat,lon,nlat,nlon,l2rec->lat,l2rec->lon,l2rec->npix,par);
    else
        get_par(l2rec, par);

    for (ip=0; ip<l2rec->npix; ip++) {
        sst = l2rec->sstref[ip];
        chl = l2rec->chl[ip];
        /*  Get the rise and set time for this day to get the number of hours of daylight:*/
        triseset(l2rec->day, l2rec->lon[ip], l2rec->lat[ip], &trise, &tset);

        switch (prodnum) {
        case CAT_opp_mld:
            prod[ip]  = get_mld (l2rec->lon[ip],l2rec->lat[ip],*l2rec->day);
            break;
        case CAT_opp_zno3:
            prod[ip]  = get_zno3(l2rec->lon[ip],l2rec->lat[ip],*l2rec->day);
            break;
        case CAT_opp_bbp:
            prod[ip]  = l2rec->bb[l2rec->nbands*ip+ib440];
            break;
        case CAT_opp_par:
            prod[ip]  = par[ip];
            break;
        default:

            if (par[ip] != BAD_FLT && chl > 0) {

                switch (prodnum) {
                case CAT_opp_befa:
                    if (sst > -2)
                        prod[ip] = opp_befa(chl, par[ip], sst, tset-trise );
                    else
                        l2rec->flags[ip] |= PRODFAIL;

                    break;
                case CAT_opp_eppley:
                    if (sst > -2)
                        prod[ip] = opp_eppley(chl, par[ip], sst, tset-trise );
                    else
                        l2rec->flags[ip] |= PRODFAIL;
                    //if (ip % 100 == 0)
                    //    printf("RJH: %f %f %f %f %f\n",prod[ip],chl, par[ip], sst, tset-trise);
                    break;
                case CAT_opp_cbpm2:
                    if (firstCall) {
                        ib440 = bindex_get(440);
                        if (ib440 < 0) {
                           printf("opp_cbpm2: incompatible sensor wavelengths (no 440 for backscatter coefficient).\n");
                           exit(1);
                        }
                    }
                    mld  = get_mld (l2rec->lon[ip],l2rec->lat[ip],*l2rec->day);
                    zno3 = get_zno3(l2rec->lon[ip],l2rec->lat[ip],*l2rec->day);
                    bbp = l2rec->bb[l2rec->nbands*ip+ib440];
                    irr = par[ip];
//                    bbp = 0.005;
//                    mld = 50;
//                    zno3=125;
//                    irr=53;
//                    tset=19;
//                    trise=5;
//                    kdtmp=0.03; //kd[ip]
//                    chl = 0.1;
//                    if (mld == BAD_FLT) mld = 50;
//                    if (zno3 == BAD_FLT) zno3 = 50;
                    if (bbp > 0 && mld != BAD_FLT && zno3 != BAD_FLT) {
                            prod[ip] = opp_cbpm2(chl, bbp, irr, kd[ip], mld, zno3, tset-trise);
                     //if (ip % 1 == 0)
                       // printf("RJH: %d %d %f %f %f %f %f %f %f %f %d ",ip,l2rec->iscan,chl, irr, sst, bbp, kd[ip], mld, zno3, tset-trise, *l2rec->day);
                        //if (prod[ip] == 0) prod[ip] = BAD_FLT;
                       // if (ip % 1 == 0)
                        //    printf("%f\n",prod[ip]);
                    }else{
                        prod[ip] = BAD_FLT;
                        //l2rec->flags[ip] |= PRODFAIL;
                    }
                    break;
                default:
                    printf("Error: %s : Unknown product specifier: %d\n",__FILE__, prodnum);
                    exit(1);
                    break;
                }
            } else {
                l2rec->flags[ip] |= PRODFAIL;
                prod[ip] = BAD_FLT;
            }
            break;
        }


    }

    free(kd);
    free(par);

    firstCall = 0;
}

void get_par_clim(float *parin,float *lat,float *lon,int nlat,int nlon,float *latp,float *lonp,int32_t npix,float *par) {

    int32_t i,j,n,incx,incy,startx,starty;
    float dx,dy;

    dy = *(lat+1) - *(lat);
    dx = *(lon+1) - *(lon);

    if (dx > 0) {
        incx=-1;
        startx = nlon-1;
    } else {
        incx = 1;
        startx = 0;
    }
    if (dy > 0) {
        incy=-1;
        starty = nlat-1;
    } else {
        incy = 1;
        starty = 0;
    }

    for (n=0; n<npix; n++) {
        for( i = startx; i >=0 && i<nlon && (fabs(lonp[n] - lon[i]) > fabs(dx)); i+=incx);
        for( j = starty; j >=0 && j<nlat && (fabs(latp[n] - lat[j]) > fabs(dy)); j+=incy);
//            printf("RJH: lat[%d]=%f latp[%d]=%f diff=%f dy=%f\n",j,lat[j],n,latp[n],fabs(latp[n] - lat[j]),dy);
        if (j<nlat && i < nlon && i>0 && j>0)
           par[n] = parin[j*nlon + i];
        else
           par[n] = BAD_FLT;

//        printf("RJH: n=%d par=%f latp=%f lonp=%f lat[%d]=%f lon[%d]=%f\n",n,par[n], latp[n],lonp[n],j,lat[j],i,lon[i]);
    }
}
/**
!C--------------------------------------------------------------------------*\

   !Description:     opp_befa - computes daily primary productivity using
                     the Behrenfeld-Falkowski (BeFa) algorithm.  The BeFa
                     algorithm estimates productivity using surface chl
                     (mg m-3), surface irradiance (Einsteins m-2 d-1),
                     sea surface temperature (C).
             Pb_opt is modelled as a polynomial function of SST.

   !Input Parameters:  
      @param[in] chl            Chlorophyll_a surface concentration in milligrams
                                chlorophyl per cubic meter
      @param[in] irr            Photosynthetically available radiation in Einsteins per
                                day per square meter
      @param[in] sst            Sea surface temperature in degrees Centigrade
      @param[in] dayL           Length day in decimal hours.

   !Output Parameters: 
      @param[out]       Primary productivity in milligrams Carbon per square meter
                        per hour

   !Revision History:  

      First programmed up by Monica Chen at Rutgers
      (1996)
      
      Revised by K. Turpie at NASA 
      (August 1997)
      
      Maintained by Don Shea at NASA
      
      Now maintained by Robert O'Malley at Oregon State University  
      (April, 2005 - present)

      Modified for inclusion in l2gen for OBPG by R. Healy at NASA
      (January 2015)

   !References and Credits
   
      Behrenfeld,M.J; Falkowski,P.G.; 1997. Photosynthetic Rates Derived
      from Satellite-Based Chlorophyll Concentration.  Limnology and 
      Oceanography, Volume 42, Number 1
      
!END------------------------------------------------------------------------*\
**/


float opp_befa( float chl, float irr, float sst, float dayL ) {

   float  chl_tot,
           z_eu   ,
           pb_opt ,
           irrFunc,
           npp;


   /* Calculate euphotic depth (z_eu) with Morel's Case I model.            */
   /* Calculate chl_tot from Satellite Surface Chlorophyll Data.            */

   if (chl <  1.0)
     chl_tot = 38.0 * pow( chl, 0.425 );
   else
     chl_tot = 40.2 * pow( chl, 0.507 );


   z_eu = 200.0 * pow( chl_tot, (-.293) );

   if (z_eu <= 102.0)
     z_eu = 568.2 * pow( chl_tot, (-.746) );


   /* Calculate the Pb_opt from satellite sea surface temperature (sst).    */
   
   if (sst < -10.0)
     pb_opt = 0.00;
   else if (sst <  -1.0)
     pb_opt = 1.13;
   else if (sst >  28.5)
     pb_opt = 4.00;
   else {
     pb_opt = 1.2956 + 2.749e-1*sst + 6.17e-2*pow(sst,2) - 2.05e-2*pow(sst, 3)
       + 2.462e-3*pow(sst,4) - 1.348e-4*pow(sst,5) + 3.4132e-6*pow(sst,6) 
       - 3.27e-8*pow(sst,7);
   }


   /* calculate the irradiance function */
   
   irrFunc = 0.66125 * irr / ( irr + 4.1 );


   /* Return the primary production calculation.                            */
   
   npp = pb_opt * chl * dayL * irrFunc * z_eu;

   return npp;
}

/*
!C--------------------------------------------------------------------------*\

   !Description:     opp_eppley - computes daily primary productivity using
                     the Behrenfeld-Falkowski (BeFa) algorithm, but
             modifies the pb_opt function after Eppley (as
             implemented by Antoine and Morel).  The BeFa
                     algorithm estimates productivity using surface chl
                     (mg m-3), surface irradiance (Einsteins m-2 d-1),
                     sea surface temperature (C), and day length (hours).
             Pb_opt is modelled as an exponential function of SST.

   !Input Parameters:
      @param[in] chl            Chlorophyll_a surface concentration in milligrams
                                chlorophyl per cubic meter
      @param[in] irr            Photosynthetically available radiation in Einsteins per
                                day per square meter
      @param[in] sst            Sea surface temperature in degrees Centigrade
      @param[in] dayL           Length day in decimal hours.

   !Output Parameters:
      @param[out]       Primary productivity in milligrams Carbon per square meter
                        per hour

   !Revision History:

      First programmed up by Monica Chen at Rutgers
      (1996)

      Revised by K. Turpie at NASA
      (August 1997)

      Maintained by Don Shea at NASA

      Now maintained by Robert O'Malley at Oregon State University
      (April, 2005 - present)

      Modified for inclusion in l2gen for OBPG by R. Healy at NASA
      (January 2015)

   !References and Credits

      Behrenfeld,M.J; Falkowski,P.G.; 1997. Photosynthetic Rates Derived
      from Satellite-Based Chlorophyll Concentration.  Limnology and
      Oceanography, Volume 42, Number 1

      Eppley, R.W.; 1972.  Temperature and Phytoplankton Growth in the Sea.
      Fishery Bulletin, Volume 79, Number 4

      Antoine, D.; Morel, A.; 1996.  Oceanic Primary Production
      1.  Adatptation of a Spectral Light-Photosynthesis Model
      in view of Application to Satellite Chlorophyll Observations

!END------------------------------------------------------------------------*\
*/

float opp_eppley( float chl,
           float irr,
           float sst,
           float dayL ) {

   float chl_tot,
          z_eu,
          pb_opt,
          irrFunc,
          npp;


   /* Calculate euphotic depth (z_eu) with Morel's Case I model.            */
   /* Calculate chl_tot from Satellite Surface Chlorophyll Data.            */

   if (chl <  1.0)
     chl_tot = 38.0 * pow( chl, 0.425 );
   else
     chl_tot = 40.2 * pow( chl, 0.507 );


   z_eu = 200.0 * pow( chl_tot, (-.293) );

   if (z_eu <= 102.0)
     z_eu = 568.2 * pow( chl_tot, (-.746) );


   /* Calculate the Pb_opt from satellite sea surface temperature (sst).    */

   pb_opt = 1.54 * pow(10, 0.0275 * sst - 0.07);


   /* calculate the irradiance function */

   irrFunc = 0.66125 * irr / ( irr + 4.1 );


   /* Return the primary production calculation.                            */

   npp = pb_opt * chl * dayL * irrFunc * z_eu;

   return npp;
}
/**

   !Description:     opp_cbpm2 - computes daily primary productivity using a chl:Carbon ratio.
                     This is a spectrally resolved version of the cbpm, using nine separate
                     wavelengths.  It is also depth resolved, integrating the effects from
                     the surface down to a fixed depth of 200 m.

                     The cbpm2 algorithm estimates productivity using chl (m-1), bbp (m-1),
             surface irradiance (Einsteins m-2 d-1), k490 (m-1), mld (m), zno3 (m)
                     and day length (hours).

Net primary productivity is carbon * growth rate, where carbon is proportional to particulate
backscatter

    carbon = 13000 * (bbp - 0.00035)

and growth rate is a function of nutrient and temperature stress (f(nut,T) and photoacclimation
(f(Ig))

    growth rate (u) = umax * f(nut,T) * f(Ig)

where:

    umax = 2

    f(nut,T) = ((Chl/C)sat - y0) / ((Chl/C)max - y0)

    f(Ig) = 1 - exp (-5 * Ig)


and:

    (Chl/C)sat = ratio of satellite observed chl and carbon (carbon from bbp)

    (Chl/C)max = 0.022 + (0.045-0.022) * exp (-3 * Ig)

    Ig = median mixed layer light level
       = surface irradiance * exp (-k(lambda) * MLD/2)

The above items are analyzed for nine separate wavelengths, and is vertically resolved to a depth
of 200 m.

For more details, please see the paper by Westberry, et al (2008)


   !Input Parameters:
      @param[in] chl            chlorophyll concentration
      @param[in] bbp            backscatter
      @param[in] irr            Photosynthetically available radiation in Einsteins per
                     day per square meter
      @param[in] k490           absorbence at 490nm
      @param[in] mld            mixing layer depth in meters
      @param[in] zno3           depth of the nitrocline
      @param[in] daylength      length of the day in decimal hours.

   !Output Parameters:
      @param[out]        Primary productivity in milligrams Carbon per square meter
                     per day

   !Dependencies:
      function austinPetzold_1986 ( float lambda, float K490 )

         given a reference k490 vlaue, determine k(lambda) for a specified lambda

         ref:
            Austin, R. W., and T. J. Petzold (1986), Spectral dependence of the diffuse
            attenuation coefficient of light in ocean waters, Opt. Eng., 25, 473 – 479

   !Revision History:

   08-16-2010 first release version (Robert O'Malley)
      [original code written in matlab by T. Westberry]

   01-05-2011   O'Malley
      add uMax trap on mu[m]
      correct z_eu determination

   !References and Credits

      Westberry, T. Behrenfeld, M.J., Siegel, D.A., and Boss, E.; 2008.  Carbon-based
      primary productivity modeling with vertically resolved photoacclimation.  Global
      Biogeochemical Cycles, Vol. 22, GB2024, doi:10.1029/2007GB003078

*/

double opp_cbpm2( double chl,
		  double bbp,
		  double irr,
		  double k490,
		  double mld,
		  double zno3,
		  double daylength) {

  double austinPetzold_1986( double, double );

  double uMax;			/* max growth rate */
  double chlCarbonMax;		/* max chl:carbon ration */
  double nutTempFunc;		/* f(nut,T) */
  double chlCarbonSat;		/* satalite chl:carbon ratio */
  double carbon;		/* bbp converted to carbon */
  double IgFunc;		/* f(Ig) */
  double IgFuncz;               /* f(Ig) below the mixed layer depth */
  double z_eu;			/* euphotic depth at 1% light level */
  double npp;                   /* net primary production */

/* --------------------- */
/*   spectral variables  */
/* --------------------- */

  double lambda[] = { 400, 412, 443, 490, 510, 555, 625, 670, 700 };
  double parFraction[] = { 0.0029, 0.0032, 0.0035, 0.0037, 0.0037, 0.0036, 0.0032, 0.0030, 0.0024 };
  double X[] = { .11748, .122858, .107212, .07242, .05943, .03996, .04000, .05150, .03000 };
  double e[] = { .64358, .653270, .673358, .68955, .68567, .64204, .64700, .69500, .60000 }; 
  double Kw[]= { .01042, .007932, .009480, .01660, .03385, .06053, .28400, .43946, .62438 };
  double Kd[9];
  double Kbio;
  double Kdif[9];

  double Klambda[9];
  double Eo[9];
  double Ez_mld[9];
  double par_mld;
  double delChlC;

  double y0;

/* --------------------------- */
/*   depth resolved variables  */
/* --------------------------- */

  double z[200];               /* depths */
  double chl_C[200];           /* chl:c ratio */
  double chlz[200];            /* chl */
  double mu[200];              /* growth */
  double Ezlambda[9][200];     /* fraction of light at nine wavelengths */
  double parz[200];            /* total light */
  double prcnt[200];           /* percent light */
  double Cz[200];              /* carbon */
  double ppz[200];             /* npp */

  int i;
  int m;
  int mzeu;
  double r;
  double prcnt0;
  double prcnt1;
  double z0;
  double z1;
  double numerator;
  double denominator;
  double fraction;
  double deltaZ;

  if(irr <= 0.0){
    return 0.0;
  }

  /* --------------------- */
  /*   initialize values   */
  /* --------------------- */

  z_eu = -9999;     //  1.05.2011
  y0 = 0.0003;                     /* min  chl:c  when  mu = 0 */
  for (i = 0; i<200; i++){
    z[i]= (float)(i+1);
  }
  r = 0.1;
  
  uMax = 2.0;                      /* after Banse (1991) */
  npp  = 0.0;
  mzeu = 0.0;

  for (i=0; i<9; i++) {
    Klambda[i]= austinPetzold_1986(lambda[i],k490);
    Eo[i] = irr * parFraction[i];
    Ez_mld[i] = Eo[i] * 0.975 * exp(-Klambda[i] * mld / 2.0);
  }

  /* ----------------------------- */
  /*   reintegrate to get par at   */
  /*   depth ...                   */
  /*   do trapezoidal integration  */
  /* ----------------------------- */

  par_mld = 0.0;
  for (i = 0; i < 8; i++ ){
    par_mld += (lambda[i+1]-lambda[i])*(Ez_mld[i+1]+Ez_mld[i])/2;
  }

  par_mld /= daylength;

  IgFunc = 1 - exp(-5.0 * par_mld);

  if(bbp < 0.00035)
    bbp = 0.00036;
  carbon = 13000.0 * (bbp - 0.00035);
  
  chlCarbonSat = chl / carbon;

  if ( chlCarbonSat < y0 ) {
    chlCarbonSat = y0;
  }

  chlCarbonMax = 0.022 + (0.045-0.022) * exp(-3.0 * par_mld);
  delChlC = chlCarbonMax - chlCarbonSat;

  nutTempFunc = (chlCarbonSat - y0) / (chlCarbonMax - y0);

  /* ''''''''''''''''''''''''' */
  /*   calculate Kd offset     */
  /*   carry through to depth  */
  /*   non-chl attenuation     */
  /* ------------------------- */

  for (i=0; i<9; i++) {
    Kbio = X[i] * pow(chl,e[i]);
    Kd[i] = Kw[i] + Kbio;
    Kdif[i] = Klambda[i] - Kd[i];
  }

  /* ''''''''''''''''''''''''''''''''''' */
  /*   integrate down the water column   */
  /*   in one-meter steps                */
  /* ----------------------------------- */

  for (m=0; m<200; m++) {

    /* ---------------------------------------------- */
    /*   if you are in the mixed layer, do this way   */
    /* ---------------------------------------------- */

    if ( z[m] < mld ) {
      chl_C[m] = chlCarbonSat;
      chlz[m] = chl_C[m] * carbon;
      mu[m] = uMax * nutTempFunc * IgFunc;

      if ( mu[m] > uMax ) {     //  1.05.2011
	mu[m] = uMax;           //  1.05.2011
      }                         //  1.05.2011
  
      for ( i=0; i<9; i++) {
	Ezlambda[i][m] = Eo[i]*0.975*exp(-Klambda[i]*z[m]);
      }

      parz[m] = 0.0;
      for (i = 0; i < 8; i++ ){
	parz[m] += (lambda[i+1]-lambda[i])*(Ezlambda[i+1][m]+Ezlambda[i][m])/2;
      }

      Cz[m] = carbon;

    } else {

      /* '''''''''''''''''''''''''''''''''''''''''''''''''''''''''' */
      /*   if below mixed layer must treat properties differently   */
      /* ---------------------------------------------------------- */

      for (i=0; i<9; i++) {
	Kbio = X[i] * pow(chlz[m-1],e[i]);     /*  after Morel & Maritorena (2001)  */
	Kd[i] = Kw[i] + Kbio;
	Kd[i] += Kdif[i];
	Ezlambda[i][m] = Ezlambda[i][m-1]*exp(-Kd[i]*1.0);
      }

      parz[m] = 0.0;
      for (i = 0; i < 8; i++ ){
	parz[m] += (lambda[i+1]-lambda[i])*(Ezlambda[i+1][m]+Ezlambda[i][m])/2;
      }

      deltaZ = zno3 - z[m];
      if ( deltaZ < 0 ) {
	deltaZ = 0;
      }

      chl_C[m] = (0.022 + (0.045-0.022) * exp(-3.0 * parz[m] / daylength));
      chl_C[m] -= delChlC * (1-exp(-0.075*deltaZ));

      IgFuncz = 1 - exp(-5.0 * parz[m]/daylength);
      mu[m] = uMax * nutTempFunc * IgFuncz;

      if ( mu[m] > uMax ) {     //  1.05.2011
	mu[m] = uMax;           //  1.05.2011
      }                         //  1.05.2011
  
      if (mu[m-1] >= r ) {
	Cz[m] = carbon;
      } else {
	Cz[m] = carbon * mu[m-1] / r;
      }

      chlz[m] = chl_C[m] * Cz[m];

    }
	   
    prcnt[m] = parz[m] / (irr * 0.975);

    /*  track this to get to the euphotic depth  */

    if ( prcnt[m] >= 0.01 ) {
      mzeu = m;
    } else {

      /* ''''''''''''''''''''''''''' */
      /*   now find 1% light depth   */
      /*   in case the user wants    */
      /*   to use this information   */
      /* --------------------------- */

      if (z_eu == -9999 ) {     // 01.05.11
        prcnt0 = prcnt[mzeu];
        prcnt1 = prcnt[mzeu+1];
        z0 = z[mzeu];
        z1 = z[mzeu+1];
        numerator = prcnt0 - 0.01;
        denominator = prcnt0 - prcnt1;
        fraction = numerator / denominator;
        z_eu = z0 + (z1-z0)*fraction;
      }
    }

    ppz[m] = mu[m] * Cz[m];

  }

  /* ------------------------------- */
  /*   do trapezoidal integration    */
  /*   from m = 0 to m = 200         */
  /* ------------------------------- */

  //  note:  186 m is the euphotic depth for pure water

  if ( mzeu < 186 ) {
    npp = 0;
    for (i = 0; i < 199; i++ ){
      npp += (z[i+1]-z[i])*(ppz[i+1]+ppz[i])/2;
    }
  } else {
    npp = BAD_FLT;
  }

  return npp;
}

/* =================================================================  */

double austinPetzold_1986 ( double lambda,
                             double K490 ) {

  double wave[] = { 350, 360, 370, 380, 390, 400,
                    410, 420, 430, 440, 450, 460, 470, 480, 490, 500,
                    510, 520, 530, 540, 550, 560, 570, 580, 590, 600,
                    610, 620, 630, 640, 650, 660, 670, 680, 690, 700 };

  double M[] = { 2.1442, 2.0504, 1.9610, 1.8772, 1.8009, 1.7383,
		 1.7591, 1.6974, 1.6108, 1.5169, 1.4158, 1.3077, 1.1982, 1.0955, 1.0000, 0.9118, 
		 0.8310, 0.7578, 0.6924, 0.6350, 0.5860, 0.5457, 0.5146, 0.4935, 0.4840, 0.4903, 
		 0.5090, 0.5380, 0.6231, 0.7001, 0.7300, 0.7301, 0.7008, 0.6245, 0.4901, 0.2891 };

  double Kdw[] = { 0.0510, 0.0405, 0.0331, 0.0278, 0.0242, 0.0217, 
		   0.0200, 0.0189, 0.0182, 0.0178, 0.0176, 0.0176, 0.0179, 0.0193, 0.0224, 0.0280, 
		   0.0369, 0.0498, 0.0526, 0.0577, 0.0640, 0.0723, 0.0842, 0.1065, 0.1578, 0.2409, 
		   0.2892, 0.3124, 0.3296, 0.3290, 0.3559, 0.4105, 0.4278, 0.4521, 0.5116, 0.6514 };

  double l0;
  double l1;
  double k0;
  double k1;
  double m0;
  double m1;
  double kdiff;
  double mdiff;
  double num;
  double den;
  double frac;
  double Kdw_l;
  double M_l;
  double Kd;

  int ref;
  int i;

  // -- INTERPOLATE TO WAVELENGTH OF INTEREST --  //

  for (i = 1; i < 36; i++) {
    if ( wave[i] >= lambda ) {
      l1 = wave[i];
      k1 = Kdw[i];
      m1 = M[i];
      l0 = wave[i-1];
      k0 = Kdw[i-1];
      m0 = M[i-1];
      break;
    }
  }

  num = lambda - l0;
  den = l1 - l0;
  frac = num / den;

  kdiff = k1 - k0;
  Kdw_l = k0 + frac*kdiff;

  mdiff = m1 - m0;
  M_l = m0 + frac*mdiff;
  

  // -- GET REFERENCE WAVELENGTH (=490 FOR NOW) AND APPLY MODEL -- //

  ref = 14;

  Kd = (M_l/M[ref]) * (K490 - Kdw[ref]) + Kdw_l;

  return Kd;

}
