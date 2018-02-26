/* =================================================================== */
/* Algorithms for computing diffuse attenuation coefficient for MSl12. */
/*                                                                     */
/* B. Franz, NASA Ocean Biology Processing Group, SAIC, March 2005.    */
/* =================================================================== */

#include <stdlib.h>
#include <math.h>
#include "l12_proto.h"

#define KD_MAX  6.4
#define KD_MIN  0.016

static float kdbad = BAD_FLT;

/* ------------------------------------------------------------------- */
/* Kd490_KD2 -  diffuse attenuation at 490nm (2-band polynomial).      */
/*                                                                     */
/* Inputs:                                                             */
/*     l2rec - level-2 structure containing one complete scan after    */
/*             atmospheric correction.                                 */
/* Outputs:                                                            */
/*     k490 - diffuse attenuation coefficient, 1 value per pixel.      */
/*                                                                     */
/* Algorithm: P.J.Werdell, June 2009                                   */
/* ------------------------------------------------------------------- */
void Kd490_KD2(l2str *l2rec, float *Kd)
{
    static int32_t  *w  = NULL;
    static float *a  = NULL;
    static int   ib1 = -1;
    static int   ib2 = -1;

    int32_t   ip, ipb;
    float R;

    if (w == NULL) {
        w = l2rec->input->kd2w;
        a = l2rec->input->kd2c;
        if (w[0] < 0 || w[1] < 0) {
            printf("Kd490_KD2: algorithm coefficients not provided for this sensor.\n");
            exit(1);
	}
        ib1 = bindex_get(w[0]);
        ib2 = bindex_get(w[1]);
        if (ib1 < 0 || ib2 < 0) {
            printf("Kd490_KD2: incompatible sensor wavelengths for this algorithm\n");
            exit(1);
        }
    }

    for (ip=0; ip<l2rec->npix; ip++) {

        ipb = l2rec->nbands*ip;

        if (l2rec->mask[ip] || 
            l2rec->Rrs[ipb+ib1] <= 0.0 || l2rec->Rrs[ipb+ib2] <= 0.0) {
            Kd[ip] = kdbad;
            l2rec->flags[ip] |= PRODFAIL;
        } else {
	    R = log10(l2rec->Rrs[ipb+ib1]/l2rec->Rrs[ipb+ib2]);
            if (isnan(R)) {
                Kd[ip] = kdbad;
                l2rec->flags[ip] |= PRODFAIL;
	    } else {
                Kd[ip] = a[0]+pow(10.0,a[1] + R*(a[2]+R*(a[3]+R*(a[4]+R*a[5]))));
                if (Kd[ip] > KD_MAX) {
                    Kd[ip] = KD_MAX;
                    l2rec->flags[ip] |= PRODWARN;
	        }
	    }
	}
    }
}


/* ------------------------------------------------------------------- */
/* Kd490_mueller() -  diffuse attenuation at 490nm (J. Mueller).       */
/*                                                                     */
/* Inputs:                                                             */
/*     l2rec - level-2 structure containing one complete scan after    */
/*             atmospheric correction.                                 */
/* Outputs:                                                            */
/*     k490 - diffuse attenuation coefficient, 1 value per pixel.      */
/*                                                                     */
/* Algorithm provided by: J. Mueller                                   */
/* OCTS coefficents: S. Bailey, 16 July 2001                           */
/* ------------------------------------------------------------------- */
void Kd490_mueller(l2str *l2rec, float k490[])
{
    int32_t   ip, ipb;

    static float a  =  0.15645;
    static float b  = -1.5401;
    static float Kw =  0.016;

    static float badval = 0.0;
    static float maxval = 6.4;

    if (l2rec->sensorID == OCTS) {
        /* Coefs for 490/565 band combination */ 
        a =  0.2166;
        b = -1.6355;
    }

    for (ip=0; ip<l2rec->npix; ip++) {

        ipb = l2rec->nbands*ip;

        if (l2rec->mask[ip]) {                   /* pixel was masked          */
            k490[ip] = kdbad;
            l2rec->flags[ip] |= PRODFAIL;
        } else if (l2rec->nLw[ipb+2] <= 0.0 ) {  /* unknown, high attenuation */
            k490[ip] = maxval;
        } else if (l2rec->nLw[ipb+4] <= 0.0 ) {  /* unknown, low  attenuation */
            k490[ip] = Kw;
        } else {
            k490[ip] = a*pow(l2rec->nLw[ipb+2]/l2rec->nLw[ipb+4],b) + Kw;
            if (k490[ip] > maxval) {
                k490[ip] = maxval;
	    }
	}
    }
}


/* ------------------------------------------------------------------- */
/* Kd490_obpg -  diffuse attenuation at 490nm (Mueller & Werdell).     */
/*                                                                     */
/* Inputs:                                                             */
/*     l2rec - level-2 structure containing one complete scan after    */
/*             atmospheric correction.                                 */
/* Outputs:                                                            */
/*     k490 - diffuse attenuation coefficient, 1 value per pixel.      */
/*                                                                     */
/* Algorithm provided by: J. Mueller                                   */
/* New coefficents and updated form: P.J.Werdell, February 2005.       */
/* ------------------------------------------------------------------- */
void Kd490_obpg(l2str *l2rec, float k490[])
{
    int32_t   ip, ipb;

    static float a;
    static float b;
    static int32_t ib490;
    static int32_t ib555;

    static int firstCall = 1;

    if (firstCall) {

        /* select 490/555 or 490/565 fit, depending on proximity of */
        /* sensor bands to fitted bands                             */

        if ((ib555 = bindex_get(551)) > 0) {
            a =  0.1853;
            b = -1.349;
        } else if ((ib555 = bindex_get(565)) > 0) {
            a =  0.1787;
            b = -1.122;
	} else {
	  printf("Kd_obpg: incompatible sensor wavelengths (no 555 or 565).\n");
          exit(1);
	}

        if ((ib490 = bindex_get(490)) < 0) {
	  printf("Kd_obpg: incompatible sensor wavelengths (no 490).\n");
          exit(1);
	}

        firstCall = 0;
    }

    for (ip=0; ip<l2rec->npix; ip++) {

        ipb = l2rec->nbands*ip;

        if (l2rec->mask[ip] || 
            l2rec->nLw[ipb+ib490] <= 0.0 || l2rec->nLw[ipb+ib555] <= 0.0) {
            k490[ip] = kdbad;
            l2rec->flags[ip] |= CHLFAIL;
            l2rec->flags[ip] |= PRODFAIL;
        } else {
            k490[ip] = a*pow(l2rec->nLw[ipb+ib490]/l2rec->nLw[ipb+ib555],b);
            if (k490[ip] > KD_MAX) {
                k490[ip] = KD_MAX;
                l2rec->flags[ip] |= PRODWARN;
	    }
	    /* not until reprocessing
            if (k490[ip] < KD_MIN) {
                k490[ip] = KD_MIN;
                l2rec->flags[ip] |= PRODWARN;
	    }
	    */
	}

    }
}


/* ------------------------------------------------------------------- */
/* Kd490_morel() - diffuse attenuation at 490 using Morel (2007)       */
/*                                                                     */
/* Inputs:                                                             */
/*     l2rec - level-2 structure containing one complete scan after    */
/*             atmospheric correction.                                 */
/*     band  - waveband number (0 - nbands-1) at which Kd computed.    */
/*                                                                     */
/* Outputs:                                                            */
/*     Kd - diffuse attenuation at 490nm, 1 value per pixel.           */
/*                                                                     */
/* Description:                                                        */
/*  This produces the estimate of diffuse attenation at 490            */
/*  using the satellite derived chlorophyll.                           */
/*                                                                     */
/* Reference:                                                          */
/*                                                                     */
/* Morel, A., Y. Huot, B. Gentili, P.J. Werdell, S.B. Hooker (2007).   */
/* Consistency of products derived from various ocean color sensors:   */
/* An examination before merging these products and extending their    */ 
/* applications, Remote Sensing of Environment, to be submitted.       */
/*                                                                     */
/* New equation 8 (combined LOV and NOMAD derivation)                  */
/*                                                                     */
/* Original Implementation: B. Franz, February 2007                    */
/*---------------------------------------------------------------------*/
void Kd490_morel(l2str *l2rec, float *Kd)
{
    static float Kw = 0.01660;
    static float X  = 0.077746;
    static float e  = 0.672846;

    float chl;
    int32_t  ip;

    for (ip=0; ip<l2rec->npix; ip++) {

        chl = l2rec->chl[ip];

        if (l2rec->mask[ip] || chl <= 0.0) {
	    Kd[ip] = kdbad;
            l2rec->flags[ip] |= PRODFAIL;
	} else {
            Kd[ip] = Kw + X * pow(chl,e);
            if (Kd[ip] > KD_MAX) {
                Kd[ip] = KD_MAX;
                l2rec->flags[ip] |= PRODWARN;
	    } else
            if (Kd[ip] < KD_MIN) {
                Kd[ip] = KD_MIN;
                l2rec->flags[ip] |= PRODWARN;
	    }
	}
    }

    return;
}


/* ------------------------------------------------------------------- */
/* Kd_PAR_morel() - spectrally integrated attenuation using Morel      */
/*                                                                     */
/* Inputs:                                                             */
/*     l2rec - level-2 structure containing one complete scan after    */
/*             atmospheric correction.                                 */
/*     depth  - layer depth 1=1/k490, 2=2/k490                         */
/*                                                                     */
/* Outputs:                                                            */
/*     Kd - Kd(PAR), 1 value per pixel.                                */
/*                                                                     */
/* Description:                                                        */
/*                                                                     */
/* Reference:                                                          */
/*                                                                     */
/* Morel, A., Y. Huot, B. Gentili, P.J. Werdell, S.B. Hooker (2007).   */
/* Consistency of products derived from various ocean color sensors:   */
/* An examination before merging these products and extending their    */ 
/* applications, Remote Sensing of Environment, to be submitted.       */
/*                                                                     */
/* Original Implementation: B. Franz, November 2006                    */
/*---------------------------------------------------------------------*/
void Kd_PAR_morel(l2str *l2rec, int depth, float *Kd)
{
    static float *Kd490 = NULL;

    int32_t  ip;

    if (Kd490 == NULL) {
        if ((Kd490 = malloc(l2rec->npix*sizeof(float))) == NULL) {
            printf("-E- %s:  Error allocating space for %d records.\n",__FILE__,l2rec->npix);
            exit(1);
	}
    }

    Kd490_morel(l2rec,Kd490);

    for (ip=0; ip<l2rec->npix; ip++) {

        if (l2rec->mask[ip] || Kd490[ip] <= 0.0) {
	    Kd[ip] = kdbad;
            l2rec->flags[ip] |= PRODFAIL;
	} else {
	    switch (depth) {
              case 1: 
                Kd[ip] = 0.0864 + 0.884 * Kd490[ip] - 0.00137/Kd490[ip];
	        break;
              case 2: 
                Kd[ip] = 0.0665 + 0.874 * Kd490[ip] - 0.00121/Kd490[ip];
	        break;
	      default:
                printf("-E- %s: Invalid depth for Kd(PAR) (1 or 2).\n",__FILE__);
                exit(1);	    
  	        break;
	    }
	}
    }

    return;
}

/* ------------------------------------------------------------------- */
/* Zhl_morel() - heated layer depth using Morel (2007)                 */
/*                                                                     */
/* Inputs:                                                             */
/*     l2rec - level-2 structure containing one complete scan after    */
/*             atmospheric correction.                                 */
/*                                                                     */
/* Outputs:                                                            */
/*     Zhl - heated layer depth (m)                                    */
/*                                                                     */
/* Description:                                                        */
/*                                                                     */
/* Reference:                                                          */
/*                                                                     */
/* Morel, A., Y. Huot, B. Gentili, P.J. Werdell, S.B. Hooker (2007).   */
/* Consistency of products derived from various ocean color sensors:   */
/* An examination before merging these products and extending their    */ 
/* applications, Remote Sensing of Environment, to be submitted.       */
/*                                                                     */
/* Original Implementation: B. Franz, November 2006                    */
/*---------------------------------------------------------------------*/
void Zhl_morel(l2str *l2rec, float *Zhl)
{
    static float *KdPAR = NULL;

    int32_t  ip;

    /* need Kd(PAR) at 2nd optical depth */
    if (KdPAR == NULL) {
        if ((KdPAR = malloc(l2rec->npix*sizeof(float))) == NULL) {
            printf("-E- %s:  Error allocating space for %d records.\n",__FILE__,l2rec->npix);
            exit(1);
	}
    }
    Kd_PAR_morel(l2rec,2,KdPAR);

    for (ip=0; ip<l2rec->npix; ip++) {

        if (l2rec->mask[ip] || KdPAR[ip] <= 1e-5) {
	    Zhl[ip] = -999;
            l2rec->flags[ip] |= PRODFAIL;
	} else
            Zhl[ip] = 2.0/KdPAR[ip];
    }

    return;
}


/* ------------------------------------------------------------------- */
/* Kd532() - spectral diffuse attenuation using Mueller/Austin&Petzold */
/*                                                                     */
/* Inputs:                                                             */
/*     l2rec - level-2 structure containing one complete scan after    */
/*             atmospheric correction.                                 */
/*     flag  - return Kd in meters (1) or 1/meters                     */
/*     Kd - diffuse attenuation at 490 nm.                             */
/* Outputs:                                                            */
/*     Kd - diffuse attenuation at 532 nm.                             */
/*                                                                     */
/* Description:                                                        */
/*  This produces the estimate of diffuse attenation at 532 nm from    */
/*  the estimate of diffuse attenuation at 490 nm using the spectral   */
/*  K algorithm of Austin and Petzold.                                 */
/*                                                                     */
/* Reference:                                                          */
/*  Austin, R. W and T. J. Petzold, "Spectral Dependence of the        */
/*  Diffuse Attenuation Coefficient of Light in Ocean Waters",         */
/*  SPIE Vol 489 Ocean Optics (1984) pp 168-178                        */
/*                                                                     */
/* Original Implementation: P. Martinolich, NRL/Stennis, May 2005      */
/*---------------------------------------------------------------------*/

void Kd532(l2str *l2rec, int flag, float k532[])
{
    const float M532  = 0.68052;
    const float KW490 = 0.0224;
    const float KW532 = 0.05356;

    static float maxval = 6.4;

    int   ip;

    float temp;

    for (ip=0; ip<l2rec->npix; ip++) {

        if ( k532[ip] < 0.0 ) 
            k532[ip] = kdbad;
        else if ( k532[ip] >= maxval )
            k532[ip] = maxval;
        else if (l2rec->mask[ip])
            k532[ip] = kdbad;
	else {
            temp = M532*(k532[ip]-KW490) + KW532;
            if ( flag > 0 )
                k532[ip] = 1.0 / temp;
            else
		k532[ip] = temp;
            if ( k532[ip] < 0.0 )
                k532[ip] = kdbad;
            if ( k532[ip] > maxval )
                k532[ip] = maxval;
        }

    }
}

/* ------------------------------------------------------------------- */
/* Kd_lee() - spectral diffuse attenuation using Lee, et. (2005)       */
/*                                                                     */
/* Inputs:                                                             */
/*     l2rec - level-2 structure containing one complete scan after    */
/*             atmospheric correction.                                 */
/*     band  - waveband number (0 - nbands-1) at which Kd computed.    */
/*                                                                     */
/* Outputs:                                                            */
/*     Kd - diffuse attenuation at specified band, 1 value per pixel.  */
/*                                                                     */
/* Description:                                                        */
/*  This produces the estimate of diffuse attenation at the given band */
/*  using the satellite derived total absorption and backscattering.   */
/*                                                                     */
/* Reference:                                                          */
/*  Lee, Z.P., M. Darecki, K.L. Carder, C.O.Davis,                     */
/*  D. Stramski, W.J. Rhea, "Diffuse Attenuation coefficient of        */
/*  downwelling irradiance:  An evalution of remote sensing methods".  */
/*                                                                     */
/* Original Implementation: P. Martinolich, NRL/Stennis, March 2005    */
/*---------------------------------------------------------------------*/
void Kd_lee(l2str *l2rec, int band, float *Kd)
{

    const float m1 =   4.18;
    const float m2 =   0.52;
    const float m3 = -10.80;

    float m0;
    int ip, ipb;

    if (l2rec->input->iop_opt == IOPNONE) {
        printf("IOP-based Kd_lee product requires iop model selection (iop_opt).  ");
        printf("Using default model.\n");
        l2rec->input->iop_opt = IOPDEFAULT;
        get_iops(l2rec,l2rec->input->iop_opt);
    }

    for (ip=0; ip<l2rec->npix; ip++) {

        ipb = ip*l2rec->nbands + band;

        if (l2rec->mask[ip] || 
            l2rec->a[ipb] <= 0.0 || l2rec->bb[ipb] <= 0.0 ) {
	    Kd[ip] = kdbad;
            l2rec->flags[ip] |= PRODFAIL;
	} else {
	    m0 = 1.0 + 0.005 * l2rec->solz[ip];
	    Kd[ip] = m0 * l2rec->a[ipb]
		       + m1 * (1.0 - m2 * exp( m3 * l2rec->a[ipb])) * l2rec->bb[ipb];
            if (Kd[ip] > KD_MAX) {
                Kd[ip] = KD_MAX;
                l2rec->flags[ip] |= PRODWARN;
	    } else
            if (Kd[ip] < KD_MIN) {
                Kd[ip] = KD_MIN;
                l2rec->flags[ip] |= PRODWARN;
	    }
        }
    }

    return;
}



/* ------------------------------------------------------------------- */
/* Kd_PAR_lee() - spectrally integrated attenuation using ZP Lee       */
/*                                                                     */
/* Inputs:                                                             */
/*     l2rec - level-2 structure containing one complete scan after    */
/*             atmospheric correction.                                 */
/*                                                                     */
/* Outputs:                                                            */
/*     Kd - Kd(PAR), 1 value per pixel.                                */
/*                                                                     */
/* Description:                                                        */
/*                                                                     */
/* Reference:                                                          */
/*                                                                     */
/* Original Implementation: B. Franz, September 2007                   */
/*---------------------------------------------------------------------*/
void Kd_PAR_lee(l2str *l2rec, float *Kd)
{
    void Zphotic_lee(l2str *l2rec, l2prodstr *p, float Zp[]);

    float *Zp = Kd;
    l2prodstr p;
    int32_t ip;

    // get first optical depth from Lee
    p.cat_ix   = CAT_Zphotic_lee;
    p.prod_ix  = -3;
    Zphotic_lee(l2rec, &p, Zp);
    
    // invert to get Kd(PAR) at 1st optical depth
    for (ip=0; ip<l2rec->npix; ip++) {
        if (Zp[ip] > 0.0)
	    Kd[ip] = 1.0/Zp[ip];
        else {
	    Kd[ip] = kdbad;
	    l2rec->flags[ip] |= PRODFAIL;
	}
    }
}

/* ------------------------------------------------------------------- */
/* Kd_jamet() - spectral diffuse attenuation using Jamet et al, (2012) */
/*                                                                     */
/* Inputs:                                                             */
/*     l2rec - level-2 structure containing one complete scan after    */
/*             atmospheric correction.                                 */
/*     band  - waveband number (0 - nbands-1) at which Kd computed.    */
/*                                                                     */
/* Outputs:                                                            */
/*     Kd - diffuse attenuation at specified band, 1 value per pixel.  */
/*                                                                     */
/* Description:                                                        */
/*  This derives the Kd using the full set of Rrs available and only   */
/*  for the SeaWiFS instrument (at this time).  A neural network is    */
/*  used to perform this.                                              */
/*                                                                     */
/* Reference:                                                          */
/*  Jamet, C., H. Loisel, and D. Dessailly (2012), Retrieval of the    */
/*  spectral diffuse attenuation coefficient Kd(lambda) in open and    */
/*  coastal ocean waters using a neural network inversion, J Geophys.  */
/*  Res., 117, C10023, doi: 10.1029/2012JC008076.                      */
/*                                                                     */
/* Original Implementation: C. Jamet                                   */
/*---------------------------------------------------------------------*/
void Kd_jamet(l2str *l2rec, int band, float *Kd)
  {
  static int   firstCall = 1;
  int ip, ipb;
  static float alg_lambda[16], alg_inp[16];
  int32 i, j, bad_prd;
  static float *b1, *w1, *w2, *mean, *sd, b2;
  static float *xnorm, *alayer;
  static int32 nrrs, n_input, n_norm, n_nodes;
  static int32 *bnd_ptr;
  float y;
 /*
  *  Initialize the coefficients for the computation
  */
  if( firstCall )
    {
    FILE *fp;
    char *filedir;
    char filename[FILENAME_MAX];
    char line [801];
    int com_ln;
    int poub;

    firstCall = 0;
    if ( (bnd_ptr = (int32 *)calloc(l2rec->nbands,sizeof(int32))) == NULL) {
        printf("-E- : Error allocating memory to bnd_ptr in get_Kd\n");
        exit(FATAL_ERROR);
    }
    if( ( filedir = getenv("OCDATAROOT")) == NULL)
      {
      printf("-E- %s, %d: OCDATAROOT env variable undefined.\n", 
        __FILE__, __LINE__ );
      exit(1);
      }
    strcpy(filename, filedir);
    switch (l2rec->sensorID)
      {
      case SEAWIFS:
        strcat(filename, "/seawifs/iop/kd_jamet/seawifs_kd_jamet.dat");
        break;
      case MERIS:
        strcat(filename, "/meris/iop/kd_jamet/meris_kd_jamet.dat");
        break;
      default:
        printf( 
      "-E- %s, %d: Kd_jamet can only be generated for the SEAWIFS instrument\n",
           __FILE__, __LINE__ );
        exit(1);
        break;
      }
    printf("Loading Kd coefficient table %s\n", filename );

    if( ( fp = fopen(filename,"r")) == NULL )
      {
      fprintf(stderr,"-E- %s, %d: unable to open %s for reading\n",
        __FILE__,__LINE__,filename);
      exit(1);
      }
   /*
    *  read the # bands for the Rrs and their values
    */
    com_ln = 1;
    while( com_ln == 1 )
      {
      if( fgets(line,800,fp) == NULL ) {
        fprintf( stderr,"-E- %s, %d:  Error reading %s\n",
          __FILE__, __LINE__, filename );
        exit(1);
        }
      if (line[0] != ';' ) com_ln = 0;
      }
    if( sscanf( line, "%d", &nrrs ) != 1 ) {
      fprintf( stderr,"-E- %s, %d:  Error reading %s\n",
        __FILE__, __LINE__, filename );
      exit(1);
      }
    for( i = 0; i < nrrs; i++ )
      {
      fgets(line,800,fp);
      sscanf( line, "%f", &alg_lambda[i] );
      }
    n_input = nrrs + 1;
    n_norm = n_input + 1;
   /*
    * get the # nodes in the hidden layer
    */
    com_ln = 1;
    while( com_ln == 1 )
      {
      fgets(line,800,fp);
      if (line[0] != ';' ) com_ln = 0;
      }
    sscanf( line, "%d", &n_nodes );

    b1 = malloc( n_nodes * sizeof( float ) );
    w1 = malloc( n_nodes * n_norm * sizeof( float ) );
    w2 = malloc( n_nodes * sizeof( float ) );
    mean = malloc( n_norm * sizeof( float ) );
    sd = malloc( n_norm * sizeof( float ) );
    xnorm =  malloc( n_input *sizeof( float ) );
    alayer = malloc( n_nodes * sizeof( float ) );

    if( ( b1 == NULL ) || ( w1 == NULL ) || ( w2 == NULL ) ||
        ( mean == NULL ) || ( sd == NULL ) ||
        ( xnorm == NULL ) || ( alayer == NULL ) )
      {
      fprintf(stderr,"-E- %s, %d: allocation of Kd weights failed\n",
        __FILE__,__LINE__ );
      exit(1);
      }
   /*
    *  get the weights for the hidden and output layers
    */
    com_ln = 1;
    while( com_ln == 1 )
      {
      fgets(line,800,fp);
      if (line[0] != ';' ) com_ln = 0;
      }
   /* ( this line is a throw-away line) */
    for(i=0; i<n_nodes; i++)
      {
      fgets(line,800,fp);
      sscanf(line,"%d %d %f",&poub,&poub,&b1[i]);
      }
    fgets(line,800,fp);
    sscanf(line,"%d %d %f",&poub,&poub,&b2);
    for(j=0; j<n_input; j++)
      for(i=0; i<n_nodes; i++)
        {
        fgets(line,800,fp);
        sscanf(line,"%d %d %f",&poub,&poub, ( w1 + i + n_nodes * j ) );
        }
    for(j=0; j<n_nodes; j++)
      {
      fgets(line,800,fp);
      sscanf(line,"%d %d %f",&poub,&poub,&w2[j]);
      }
   /*
    *  lastly, read the mean and standard deviation data
    */
    com_ln = 1;
    while( com_ln == 1 )
      {
      fgets(line,800,fp);
      if (line[0] != ';' ) com_ln = 0;
      }
    sscanf(line, "%f",&mean[0]);
    for(i=1; i< n_norm; i++) {
      fgets(line,800,fp);
      sscanf(line, "%f",&mean[i]);
      }
    com_ln = 1;
    while( com_ln == 1 )
      {
      fgets(line,800,fp);
      if (line[0] != ';' ) com_ln = 0;
      }
    sscanf( line, "%f",&sd[0]);
    for(i=1; i< n_norm; i++) {
      fgets(line,800,fp);
      sscanf( line, "%f",&sd[i]);
      }
    fclose( fp );
   /*
    *  set up correct band locations
    */
    printf( "Kd_jamet band assignment:\n" );
    printf( "Alg wave   inst wave   inst band\n" );
    for( i = 0; i < nrrs; i++ ) {
      bnd_ptr[i] = windex( alg_lambda[i], l2rec->fwave, l2rec->nbands );
      printf( "%7.2f    %7.2f    %3d\n", alg_lambda[i], 
        l2rec->fwave[bnd_ptr[i]], bnd_ptr[i] );
      }
    printf( "\n" );
    }
 /*
  *  derive the Kd
  */
  for (ip=0; ip<l2rec->npix; ip++) {
    ipb = ip * l2rec->nbands;
   /*
    *  Normalize the Rrs, lambda
    */
    bad_prd = 0;
    if( l2rec->mask[ip] ) bad_prd = 1;
/*
    for( i = 0; i < nrrs; i++ )
      {
      alg_inp[i] = l2rec->Rrs[ ipb + bnd_ptr[i] ];
      if( i <= 1 ) {
        if( alg_inp[i] <= -0.01 ) bad_prd = 1;
        }
      else {
        if( alg_inp[i] <= 0.0 ) bad_prd = 1;
        }
      }
*/
    for( i = 0; i < nrrs; i++ )
      {
      alg_inp[i] = l2rec->Rrs[ ipb + bnd_ptr[i] ];
      if( alg_inp[i] < BAD_FLT + 1 )  bad_prd = 1;
      }
    alg_inp[nrrs] = l2rec->fwave[band];

    if( bad_prd == 1 ){
      Kd[ip] = kdbad;
      l2rec->flags[ip] |= PRODFAIL;
      }
    else {
      for( i = 0; i < n_input; i++ )
        xnorm[i] = ( ( 2./3. ) *( alg_inp[i] - mean[i] ) ) / sd[i];
     /*
      *  Apply the nn algorithm
      */
      for( i = 0; i < n_nodes; i++ ){
        alayer[i] = 0.0;
        for( j = 0; j < n_input; j++ ){
          alayer[i] += ( xnorm[j] * *( w1 + i + n_nodes * j ) );
        }
        alayer[i] = 1.715905 * (float)tanh( ( 2./3.) * 
          ( double)( alayer[i] + b1[i] ) );
      }
    
      for( j = 0, y = 0.; j < n_nodes; j++ )
        y += ( alayer[j] * w2[j] );
     /*
      *  De-normalize the log(Kd) and make Kd
      */
      y = 1.5 * ( y + b2) * sd[n_input] + mean[n_input];
      *( Kd + ip )= (float)pow( 10., (double)y );
      }
    }
  return;
  }

/* ------------------------------------------------------------------- */
/* Kd_rhos() - spectral diffuse attenuation using Stumpf, et. (2013) */
/*                                                                     */
/* Inputs:                                                             */
/*     l2rec - level-2 structure containing one complete scan after    */
/*             atmospheric correction.                                 */
/*     band  - waveband number (0 - nbands-1) at which Kd computed.    */
/*                                                                     */
/* Outputs:                                                            */
/*     Kd - diffuse attenuation at specified band, 1 value per pixel.  */
/*                                                                     */
/* Description:                                                        */
/*  This produces the estimate of diffuse attenation at the given band */
/*  using the satellite derived total absorption and backscattering.   */
/*                                                                     */
/* Reference:                                                          */
/* Original Implementation: R. Healy, NASA-GSFC, August 2015           */
/*---------------------------------------------------------------------*/
void Kd_rhos(l2str *l2rec,float *Kd)
{

    const float factor =   0.7;

    static int ib443,ib490,ib620,ib665,ib645,ib469,ib859, ib865,firstCall=1;
    int ip, ipb;

    if (firstCall == 1) {

        ib645 = bindex_get(645);
        ib469 = bindex_get(469);
        ib859 = bindex_get(859);
        firstCall = 0;

            ib443 = -1;
        if (ib645 < 0 || ib469 < 0 || ib859 < 0) {
                // try the Meris wavelengths

                ib443 = bindex_get(443);
                ib490 = bindex_get(490);
                ib620 = bindex_get(620);
                ib665 = bindex_get(665);
                ib865 = bindex_get(865);
                firstCall = 0;


                if (ib443 < 0 || ib490 < 0 || ib620 < 0 || ib665 < 0 || ib865 < 0) {
            printf("Kd_rhos: incompatible sensor wavelengths for this algorithm\n");
            exit(1);
                } else {
                    printf("Kd_rhos: Using the Meris band averaging for this product\n");
        }

    }

    }

    if (ib443 < 0) {
        for (ip=0; ip<l2rec->npix; ip++) {

            ipb = l2rec->nbands*ip;

       // if (l2rec->mask[ip] ||
         if(       l2rec->Rrs[ipb+ib645] <= 0.0 || l2rec->Rrs[ipb+ib469] <= 0.0 || l2rec->Rrs[ipb+ib859] <= 0.0) {
            Kd[ip] = kdbad;
            l2rec->flags[ip] |= PRODFAIL;
        } else {
                Kd[ip] = factor * (l2rec->rhos[ipb+ib645] - l2rec->rhos[ipb+ib859])
                        /(l2rec->rhos[ipb+ib469] - l2rec->rhos[ipb+ib859]);
            }
        }
    } else {
        for (ip=0; ip<l2rec->npix; ip++) {

            ipb = l2rec->nbands*ip;

           // if (l2rec->mask[ip] ||
             if(l2rec->Rrs[ipb+ib665] <= 0.0 || l2rec->Rrs[ipb+ib620] <= 0.0 || l2rec->Rrs[ipb+ib443] <= 0.0 || l2rec->Rrs[ipb+ib490] <= 0.0 || l2rec->Rrs[ipb+ib865] <= 0.0) {
                Kd[ip] = kdbad;
                l2rec->flags[ip] |= PRODFAIL;
            } else {
                Kd[ip] = factor * ((l2rec->rhos[ipb+ib620]+l2rec->rhos[ipb+ib665])/2 - l2rec->rhos[ipb+ib865])
                        /((l2rec->rhos[ipb+ib443]+l2rec->rhos[ipb+ib490])/2 - l2rec->rhos[ipb+ib865]);
            }

        }
    }
    return;
}

/* ------------------------------------------------------------------- */
/* get_Kd() - l2_hdf_generic interface for Kd                          */
/* ------------------------------------------------------------------- */
void get_Kd(l2str *l2rec, l2prodstr *p, float prod[])
{
    switch (p->cat_ix) {
    case CAT_Kd_mueller:
        Kd490_mueller(l2rec,prod);
        break;
    case CAT_Kd_obpg:
        Kd490_obpg(l2rec,prod);
        break;
    case CAT_Kd_KD2:
        Kd490_KD2(l2rec,prod);
        break;
    case CAT_Kd_lee:
        Kd_lee(l2rec,p->prod_ix,prod);
        break;
    case CAT_Kd_532:
        Kd490_mueller(l2rec,prod);
        Kd532(l2rec,p->prod_ix,prod);
        break;
    case CAT_Kd_morel:
        Kd490_morel(l2rec,prod);
        break;
    case CAT_Kd_jamet:
        Kd_jamet(l2rec,p->prod_ix,prod);
        break;
    case CAT_Kd_rhos:
        Kd_rhos(l2rec,prod);
        break;
    case CAT_KPAR_morel:
        Kd_PAR_morel(l2rec,p->prod_ix,prod);
        break;
    case CAT_KPAR_lee:
        Kd_PAR_lee(l2rec,prod);
        break;
    case CAT_Zhl_morel:
        Zhl_morel(l2rec,prod);
        break;
    default:
        printf("Error: %s : Unknown product specifier: %d\n",__FILE__,p->cat_ix);
        exit(1);
        break;
    }

    return;
}





/* =================================================================== */
/* defunct versions                                                    */
/* =================================================================== */


/* ------------------------------------------------------------------- */
/* Kd_morel() - spectral diffuse attenuation using Morel (2007)        */
/*                                                                     */
/* Inputs:                                                             */
/*     l2rec - level-2 structure containing one complete scan after    */
/*             atmospheric correction.                                 */
/*     band  - waveband number (0 - nbands-1) at which Kd computed.    */
/*                                                                     */
/* Outputs:                                                            */
/*     Kd - diffuse attenuation at specified band, 1 value per pixel.  */
/*                                                                     */
/* Description:                                                        */
/*  This produces the estimate of diffuse attenation at the given band */
/*  using the satellite derived chlorophyll.                           */
/*                                                                     */
/* Reference:                                                          */
/*                                                                     */
/* Morel, A., Y. Huot, B. Gentili, P.J. Werdell, S.B. Hooker (2007).   */
/* Consistency of products derived from various ocean color sensors:   */
/* An examination before merging these products and extending their    */ 
/* applications, Remote Sensing of Environment, to be submitted.       */
/*                                                                     */
/* Equation 8                                                          */
/*                                                                     */
/* Original Implementation: B. Franz, August 2006                      */
/*---------------------------------------------------------------------*/
/*
void Kd_morel(l2str *l2rec, int band, float *Kd)
{
    typedef struct kd_morel_table {
      float wave;
      float Kw;
      float X;
      float e;
    } kdtabstr;

    static kdtabstr *kdtab = NULL;
    static int ntab = 0;
                
    float chl;
    float Kd1, Kd2;
    float wt;             
    float wave;
    int32_t  ip, itab, i1, i2;

    if (kdtab == NULL) {
        FILE *fp = NULL;
        char *tmp_str;
        char  file   [FILENAME_MAX] = "";
        if ((tmp_str = getenv("OCDATAROOT")) == NULL) {
            printf("OCDATAROOT environment variable is not defined.\n");
            exit(1);
        }
        strcpy(file,tmp_str); strcat(file,"/common/kd_morel.dat");
        if ( (fp = fopen(file,"r")) == NULL ) {
            printf("-E- %s:  Error opening file %s.\n",__FILE__,file);
            exit(1);
        }
        printf("\nLoading Morel spectral Kd table\n    %s\n",file);
        fscanf(fp,"%d",&ntab);
        if (ntab <= 0) {
            printf("-E- %s:  Error reading %s.\n",__FILE__,file);
            exit(1);
        }
        if ((kdtab = malloc(ntab*sizeof(kdtabstr))) == NULL) {
            printf("-E- %s:  Error allocating space for %d records.\n",__FILE__,ntab);
            exit(1);
	} 
        for (itab=0; itab<ntab; itab++) {
            fscanf(fp,"%f %f %f %f",
                &kdtab[itab].wave,
                &kdtab[itab].Kw,
                &kdtab[itab].X,
		&kdtab[itab].e);
            printf("%f %f %f %f\n",
                kdtab[itab].wave,
                kdtab[itab].Kw,
                kdtab[itab].X,
		kdtab[itab].e);
        }
        fclose(fp);
        printf("\n");
    }

    if (band >= 0) 
        wave = l2rec->fwave[band];
    else
        wave = 490.;

    for (itab=0; itab<ntab; itab++)
        if (wave <= kdtab[itab].wave)
	    break;

    if (wave == kdtab[itab].wave) {
        i1 = itab;
        i2 = itab;
    } else if (itab == ntab) {
        i1 = itab-2; 
        i2 = itab-1;
        wt = (wave-kdtab[i1].wave)/(kdtab[i2].wave - kdtab[i1].wave); 
    } else {
        i1 = itab; 
        i2 = itab+1;
        wt = (wave-kdtab[i1].wave)/(kdtab[i2].wave - kdtab[i1].wave); 
    }             

    for (ip=0; ip<l2rec->npix; ip++) {

        chl = l2rec->chl[ip];

        if (l2rec->mask[ip] || chl <= 0.0) {
	    Kd[ip] = kdbad;
            l2rec->flags[ip] |= PRODFAIL;
	} else {
	    if (i1 == i2) {
	        Kd[ip] = kdtab[i1].Kw + kdtab[i1].X * pow(chl,kdtab[i1].e);
	    } else {
	        Kd1 = kdtab[i1].Kw + kdtab[i1].X * pow(chl,kdtab[i1].e);
	        Kd2 = kdtab[i2].Kw + kdtab[i2].X * pow(chl,kdtab[i2].e);
                Kd[ip] = Kd1 + (Kd2-Kd1)*wt;
	    }
            if (Kd[ip] > KD_MAX) {
                Kd[ip] = KD_MAX;
                l2rec->flags[ip] |= PRODWARN;
	    } else
            if (Kd[ip] < KD_MIN) {
                Kd[ip] = KD_MIN;
                l2rec->flags[ip] |= PRODWARN;
	    }
	}
    }

    return;
}
*/
