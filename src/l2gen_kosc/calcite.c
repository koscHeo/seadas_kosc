/*---------------------------------------------------------------------*/
/* calcite.c -  get calcium carbonate concentration.                   */
/*                                                                     */
/* Inputs:                                                             */
/*     l2rec - level-2 structure containing one complete scan after    */
/*             atmospheric correction.                                 */
/* Outputs:                                                            */
/*     caco3 - calcium carbonate concentration, per pixel        .     */
/*                                                                     */
/* Written by: W. Robinson, GSC, 7 Jun 2000.                           */
/*             S. Bailey, OCDPG, July 2004, conversion to C.           */
/*             B. Franz, OCDPG, Sep 2004, sensor generalization and    */
/*                 implementation of 2-Band algorithm.                 */
/*                                                                     */
/*             2014: Standardized to use common table for 2-band alg   */
/*             and adjust green nLw as needed for sensor.              */
/*                                                                     */
/*             2014: Changed bbstar from 4 to 1.628.                   */
/*                                                                     */
/*---------------------------------------------------------------------*/

#include <stdlib.h>
#include <math.h>
#include "l12_proto.h"
#include "l2_flags.h"

#define BAD_CACO3 BAD_FLT

static int32_t  caco3_msk = LAND | HIGLINT | CLOUD;
static float pi        = PI;
static float radeg     = RADEG;
static float bbstr     = 1.628;
static float caco3min  = 1.18e-5;
static float caco3hi   = 0.003;

/* --------------------------------------------------------------------- */
/* calcite_3b() - calcium carbonate concentration from 3-Band algorith.. */
/*                                                                       */
/* Gordon, H.R. Boynton, G.C., Balch, W.M., Groom, S.B., Harbour, D.S.,  */
/* Smyth, T.J., Retrieval of Coccolithophore Calcite Concentration from  */
/* SeaWiFS Imagery, GRL, 28, 8, 1587-1590.                               */
/*                                                                       */
/* --------------------------------------------------------------------- */

float calcite_3b(l2str *l2rec, int32_t ip)
{
    static int   firstCall = 1;
    static int   maxiter   = 10;
    static float ftrans    = 6.179;               /* (1/.298)*(1/.543)   */

    static float wave[3]   = {670.,760.,870.};    /* approx. wavelengths */
    static int   bx  [3];
    static float aw  [3];
    static float bbw [3];
    static float bbc [3];
    static float t   [3];
    static float b68diff;
    static float b78diff;
    static float fw1,fw2;
    
    static float oobswf[3][8] = { 
    {0.000313529, 0.000770558, 0.00152194,  0.000155573, 
     0.00116455,  0.0,         0.000445433, 0.000124172},
    {0.000201709, 6.96143e-05, 7.00147e-06, 2.28957e-07, 
     4.17788e-05, 0.00159814,  0.0,         0.00536827},
    {0.000463807, 8.54003e-05, 2.47401e-05, 0.000755890, 
     0.00587073,  0.00021686,  0.0111331,   0.0} };

    int32_t ipb, ib, i;
    float  *rhoaw;
    float  rho[3];
    float  bbc_cclth, r8_cclth, aeps_cclth;
    float  mu0;
    float  bbcinit, bbctol;
    int    numiter;
    int32_t   nwave, status = 0;
    float *awptr, *bbwptr;
    float caco3;
    float bbw546;

    nwave = l2rec->nbands;

    if (firstCall) {

        /* save coeffs for the three bands, resolve actual sensor wave */
        for (i=0; i<3; i++) {
            bx  [i] = windex(wave[i],l2rec->fwave,nwave);
            wave[i] = l2rec->fwave[bx[i]];

            bbc [i] = 0.0;
	}

        b68diff = wave[2]-wave[0];
        b78diff = wave[1]-wave[0];

        /* spectral dependence of bbc */
        fw1 = pow(wave[0]/wave[1],1.35);
        fw2 = pow(wave[0]/wave[2],1.35);

        firstCall = 0;
    }

    /* set aw & bbw */

    ipb = ip*nwave;
    awptr  = &l2rec->sw_a_avg[ipb];
    bbwptr = &l2rec->sw_bb_avg[ipb];

    for (i=0; i<3; i++) {
        aw  [i] = awptr [bx[i]];
        bbw [i] = bbwptr[bx[i]];
    }
    bbw546 = seawater_bb(546.0,l2rec->sstref[ip],l2rec->sssref[ip]);

    status    = 0;
    numiter   = 0;
    bbctol    = 100.;
    bbcinit   = 0.00085;
    caco3     = BAD_CACO3;
    bbc[0]    = 0.000;

    /* skip pixel if already masked (this should not include ATMFAIL) */
    if ((l2rec->flags[ip] & caco3_msk) != 0) {
        return(caco3);
    }

    mu0 = cos( l2rec->solz[ip]/radeg );

    if ( (rhoaw = (float *)calloc( l2rec->nbands,sizeof(float))) == NULL) {
        printf("-E- : Error allocating memory to rhoaw\n");
        exit(FATAL_ERROR);
    }

    /* compute the aerosol/water reflectance (include out-of-band correction for SeaWiFS) */
    if (l2rec->sensorID == SEAWIFS) {
        for( ib = 0; ib < l2rec->nbands; ib++ ) {
            ipb = nwave*ip + ib;
            rhoaw[ib] = ( ( l2rec->Lt[ipb]/l2rec->tg_sol[ipb]/l2rec->tg_sen[ipb] - l2rec->tLf[ipb]
                      - l2rec->Lr[ipb] ) / l2rec->t_o2[ipb] - l2rec->TLg[ipb]) * pi/l2rec->Fo[ib]/mu0;
        }
        for (i = 0; i < 3; i++) {
	    rho[i] = rhoaw[bx[i]];
            for( ib = 0; ib < l2rec->nbands; ib++ ) {
                rho[i] -= rhoaw[ib]*oobswf[i][ib]; 
                if (rho[i] <= 0.0) status = 1;
            }
        }
    } else {
        for( i = 0; i < 3; i++ ) {
            ib  = bx[i];
            ipb = nwave*ip + ib;
            rho[i] = ( (l2rec->Lt[ipb]/l2rec->tg_sol[ipb]/l2rec->tg_sen[ipb] - l2rec->tLf[ipb]
                   - l2rec->Lr[ipb] ) / l2rec->t_o2[ipb] - l2rec->TLg[ipb]) * pi/l2rec->Fo[ib]/mu0;
            if (rho[i] <= 0.0) status = 1;
        }
    }

    /* skip pixel on negative surface reflectance */
    if (status != 0) {
        free(rhoaw);
        return (caco3);
    }

    /* compute total transmittance */
    for (i=0; i<3; i++) {        
        ipb  = nwave*ip + bx[i];
        t[i] = l2rec->tg_sol[ipb]*l2rec->tg_sen[ipb]*l2rec->t_sol[ipb]*l2rec->t_sen[ipb];
    }

    /* compute backscatter at 546 nm */
    while (bbctol > 5. && numiter < maxiter) {

        numiter++;

        bbc[1] = bbc[0]*fw1;
        bbc[2] = bbc[0]*fw2;

        /* reflectance at longest wavelength */        
        r8_cclth = rho[2] - (bbw[2]+bbc[2])/(aw[2]+bbw[2]+bbc[2])/ftrans*t[2];

        if ((r8_cclth > 0.09) || (r8_cclth < 0)){
            status = 1;
            bbc[0] = 0;
            break;
        }

        /* atmospheric epsilon at two longest wavelengths */        
        aeps_cclth = log((rho[1] - (bbw[1]+bbc[1])/(aw[1]+bbw[1]+bbc[1])/ftrans*t[1])/r8_cclth)/b78diff;

        if (aeps_cclth > 0.4){
            status = 1;
            bbc[0] = 0;
            break;
        }

        /* --------------- */
        bbc[0] = (rho[0] - r8_cclth * exp(aeps_cclth*b68diff))/t[0] * (aw[0]+bbw[0]+bbc[0])*ftrans - bbw[0];

        if ((bbc[0] <= 0) || isnan(bbc[0])) {
            status = 1;
            bbc[0] = 0;
            break;
        }
       
        bbctol  = fabs((bbcinit - bbc[0])/bbcinit)*100.;
        bbcinit = bbc[0];                        
    }


    if (status == 0) {

        bbc_cclth = bbc[0]/pow((546./wave[0]),1.35) - bbw546;

        // convert to calcite in moles/m^3 (Balch 2005)
        caco3 = bbc_cclth / bbstr;                   
        if (caco3 < caco3min ) {
	    caco3 = caco3min;                        
        }
    } 

    free(rhoaw);

    return(caco3);
}


/* --------------------------------------------------------------------- */
/* calcite_2b() - calcium carbonate concentration from 2-Band algorith.. */
/*                                                                       */
/* Gordon, H.R. and Balch, W.M., MODIS Detached Coccolith Concentration  */
/* Algorithm Theoretical Basis Document,  April 30, 1999                 */
/*                                                                       */
/* --------------------------------------------------------------------- */

#define N443 490
#define N550 456

float calcite_2b(l2str *l2rec, int32_t ip)
{
    static int   firstCall = 1;
    static int   bandShift = 0;
    static float t443[N443]; 
    static float t550[N550]; 
    static float tbb [N443][N550];
    static int32_t  ib443;
    static int32_t  ib550;

    float caco3;
    int   i443 = 0, i550 = 0;
    float x,y,a,b;
    int32_t  nwave, i,nc;
    float x443, x550;
    float bb1, bb2, bb;

    nwave = l2rec->nbands;

    if (firstCall) {

        FILE *fp;
        char *filedir;
        char filename[FILENAME_MAX];
        char line [80];

        strcpy(filename,l2rec->input->picfile);
        if (strlen(filename) == 0) {
            printf("-E- %s line %d: No picfile specified.\n",__FILE__,__LINE__);
            exit(1);
        }
        printf("Loading PIC 2-band algorithm table %s\n",filename);

        ib443 = windex(443.,l2rec->fwave,nwave);
        ib550 = windex(550.,l2rec->fwave,nwave);

        if (strstr(filename,"common") != NULL) {
            printf("Assuming PIC table is for 443nm and 555nm.\n");
            bandShift = 1;
            ib443 = bindex_get(443);
            ib550 = bindex_get_555();
            if (ib443 < 0 || ib550 < 0) {
                printf("-E- %s line %d: required bands not available PIC\n",
                __FILE__,__LINE__);
                exit(1);
	    }            
	}

        if ( (fp = fopen(filename,"r")) == NULL ) {
            fprintf(stderr,"-E- %s line %d: unable to open %s for reading\n",
                    __FILE__,__LINE__,filename);
            exit(1);
        }

        // Skip comment lines
        nc = 0;
        while ( fgets(line,80,fp ) ) {
	    if (line[0] != '#' && line[0] != '\n' ) {
                break;
	    }
	    nc++;
	}
        rewind(fp);
        for (i=0; i<nc; i++)
            fgets(line,80,fp);

        // Load table
        
        for (i443=0; i443<N443; i443++)
          for (i550=0; i550<N550; i550++) {
            fscanf(fp,"%f %f %f %f\n",&x,&y,&a,&b);
            t443[i443] = x;
            t550[i550] = y;
            tbb [i443][i550] = a;
	  }

        firstCall = 0;
    }

    caco3 = BAD_CACO3;

    // skip pixel if already masked (this includes ATMFAIL)
    if (l2rec->mask[ip]) {
      return(caco3);
    }

    x443 = l2rec->nLw[ip*nwave+ib443];
    x550 = l2rec->nLw[ip*nwave+ib550];

    // if required radiances are negative, set to background level
    if (x550 <= 0.0) {
        return(caco3);
    }
    if (x443 <= 0.0) {
        return(caco3);
    }

    // adjust nLw to 555 based on Rrs555/Rrs5xx ratio
    if (bandShift) {
        float Rrs555 = conv_rrs_to_555(l2rec->Rrs[ip*nwave+ib550],l2rec->fwave[ib550]);
        x550 *= Rrs555/l2rec->Rrs[ip*nwave+ib550];
    }

    // locate bounding table indices
    for (i=0; i<N443; i++) {
        if (x443 < t443[i]) {
            i443 = i;
            break;
	}
    }
    for (i=0; i<N550; i++) {
        if (x550 < t550[i]) {
	    i550 = i;
	    break;
	}
    }

    // radiances less than table entries, fail
    if (i443 <=0) {
        return(caco3);
    }
    if (i550 <=0) {
        return(caco3);
    }

    // radiances greater than table entries, fail
    if (i443 >= N443 || i550 >= N550) {
        return(caco3);
    }

    // radiances associated with missing table entries, fail
    if (tbb[i443-1][i550-1] > 998.9 || tbb[i443  ][i550-1] > 998.9 || 
        tbb[i443-1][i550  ] > 998.9 || tbb[i443  ][i550  ] > 998.9) {
        return(caco3);
    }

    // interpolate to get bb(546) 
    bb1 = tbb[i443-1][i550-1] + (x443 - t443[i443-1])*
      (tbb[i443][i550-1] - tbb[i443-1][i550-1])/(t443[i443] - t443[i443-1]);

    bb2 = tbb[i443-1][i550  ] + (x443 - t443[i443-1])*
      (tbb[i443][i550  ] - tbb[i443-1][i550  ])/(t443[i443] - t443[i443-1]);

    bb = bb1 + (x550 - t550[i550-1]) * (bb2 - bb1)/(t550[i550] - t550[i550-1]);

    // convert to calcite in moles/m^3 (Balch 2005)
    caco3 = bb / bbstr;   
     
    if (caco3 < caco3min ) {   
        caco3 = caco3min;     
    }

    return(caco3);
}



/* --------------------------------------------------------------------- */
/* calcite_c() - calcium carbonate concentration (combined algorithm)    */

/* --------------------------------------------------------------------- */
float calcite_c(l2str *l2rec, int32_t ip) {
    float caco3 = BAD_CACO3;
    int32_t npix = l2rec->npix;

    // if the 2-band algorithm fails for any reason, try the 3-band 
    // but only accept the high-value retrievals from the 3-band

    int nwave = l2rec->nbands;
    int32_t ib443 = windex(443.,l2rec->fwave,nwave);
    int32_t ib550 = windex(550.,l2rec->fwave,nwave);
    float nLw443 = l2rec->nLw[ip * nwave + ib443];
    float nLw550 = l2rec->nLw[ip * nwave + ib550];

    // if required radiances are negative, set to background level
    if (nLw550 <= 0.0 || nLw443 <= 0) {
        caco3 = calcite_3b(l2rec, ip);
        if (caco3 < caco3hi)
            caco3 = BAD_CACO3;
    } else {
        caco3 = calcite_2b(l2rec, ip);
    }
    // if valid value
    if (caco3 > BAD_CACO3) {
        caco3 = MAX(caco3, caco3min);
    }

    return (caco3);
}

/* ------------------------------------------------------------------- */
/* calcite() - l2_hdf_generic interface for calcite (pic)              */
/* ------------------------------------------------------------------- */
void calcite(l2str *l2rec, l2prodstr *p, float prod[])
{
    int32_t ip;

    for (ip=0; ip<l2rec->npix; ip++) {
        switch (p->cat_ix) {
          case CAT_calcite:
            prod[ip] = calcite_c(l2rec,ip);
            break;
          case CAT_calcite_2b:
            prod[ip] = calcite_2b(l2rec,ip);
            break;
          case CAT_calcite_3b:
            prod[ip] = calcite_3b(l2rec,ip);
            break;
          default:
            printf("Error: %s : Unknown product specifier: %d\n",__FILE__,p->cat_ix);
            exit(1);
            break;
        }

        if (prod[ip] == BAD_CACO3)
            l2rec->flags[ip] |= PRODFAIL;
    }

    return;
}

// Latitude-dependent bbstar, this not currently in use
float get_bbstar(float lat) {

    static float c[] = {8.701E-01,1.200E-01,-5.999E-04,-6.606E-04,
            4.202E-05,-1.150E-06,1.614E-08,-1.138E-10,3.196E-13};

    float bbstr;
    float alat  = fabs(lat);

    if (alat > 55) 
        bbstr = 1.429;
    else if (alat < 1) 
        bbstr = 0.877;
    else
        bbstr = c[0]+alat*(c[1]+alat*(c[2]+alat*(c[3]+alat*(c[4]+alat*(c[5]+alat*(c[6]+alat*(c[7]+alat*c[8])))))));

    return (bbstr);
}

