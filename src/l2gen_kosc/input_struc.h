#ifndef _INPUT_STR_H
#define _INPUT_STR_H

#include <stdio.h>
#include "hdf.h"
#include "l12_parms.h"
#include "l2prod.h"
#include "filter.h"

#ifdef __cplusplus
extern "C" {
#endif
    
typedef struct input_struct {

  int32_t   sensorID;
  int32_t   subsensorID;
  int32_t   format;
  int32_t   nbands;

  int       modis_subset_compat; /* force modis file to be read as subsetted */
  char      oformat[20];  // output file type
  char      oformat_depth[20];  // output file color depth l1brsgen only

  char      ifile   [MAX_IFILES][FILENAME_MAX];
  char      ofile   [MAX_OFILES][FILENAME_MAX];
  char      l2prod  [MAX_OFILES][PRODSTRLEN];
  char      suite   [32];
  char      ilist   [FILENAME_MAX];
  char      calfile [FILENAME_MAX];
  char      xcalfile[FILENAME_MAX];
  char      polfile [FILENAME_MAX];
  char      cldfile [FILENAME_MAX];  
  char      geofile [FILENAME_MAX];
  char      metafile[FILENAME_MAX];
  char      btfile  [FILENAME_MAX];
  char      fqfile  [FILENAME_MAX];
  char      parfile  [FILENAME_MAX];
  char      sstcoeffile [FILENAME_MAX];
  char      sst4coeffile[FILENAME_MAX];
  char      sst3coeffile[FILENAME_MAX];
  char      sstssesfile [FILENAME_MAX];
  char      sst4ssesfile[FILENAME_MAX];
  char      sst3ssesfile[FILENAME_MAX];
  char      picfile[FILENAME_MAX];
  char      owtfile[FILENAME_MAX];
  char      owtchlerrfile[FILENAME_MAX];
  char      vcnnfile[FILENAME_MAX];
  char      def_l2prod[MAX_OFILES][PRODSTRLEN];

  char      aermodfile [FILENAME_MAX];
  char      aermodels [MAXAERMOD][32];
  int32_t   naermodels;
  int32_t   aermodmin;
  int32_t   aermodmax;
  float     aermodrat;

  int32_t   evalmask;
  int32_t   mode;           /* 0: Forward processing only               */
                            /* 1: Inverse calibration mode, given Lw    */
                            /* 2: Inverse calibration mode, given nLw   */
                            /* 3: Inverse calibration mode, assuming nLw=0 */
  int32_t   spixl;          /* starting pixel no. of the input (1-rel)  */
  int32_t   epixl;          /* ending pixel no. of the input (1-rel)    */
  int32_t   dpixl;          /* pixel subsampling increment              */
  int32_t   sline;          /* starting line no. of the input (1-rel)   */
  int32_t   eline;          /* ending line no. of the input (1-rel)     */
  int32_t   dline;          /* line subsampling increment               */
  int32_t   ctl_pt_incr;    /* control-point reduction factor           */

  int32_t   proc_ocean;     /* 1=perform ocean-specific processing      */
  int32_t   proc_land;      /* 1=perform land-specific processing       */
  int32_t   proc_sst;       /* 1=perform sst-specific processing        */
  int32_t   atmocor;        /* 1=perform atmospheric correction         */
  int32_t   seawater_opt;   /* pure seawater IOP options                */
  int32_t   aer_opt;        /* aerosol model option                     */
  int32_t   aer_wave_short; /* shortest wavelength for model selection  */
  int32_t   aer_wave_long;  /* longest wavelength for model selection   */
  int32_t   aer_swir_short; /* shortest wavelength foe swir to nir corr */
  int32_t   aer_swir_long;  /* longest wavelength for swir to nir corr  */
  float     aer_rrs_short;  /* Rrs at shortest aerosol wavelength       */
  float     aer_rrs_long;   /* Rrs at longest aerosol wavelength        */
  float     aer_angstrom;   /* angstrom for aerosol model selection     */
  int32_t   gas_opt;        /* 1=ozone, 2=co2, 4=no2, 8=h2o             */
  int32_t   atrem_opt;      /* 1=ozone, 2=co2, 4=no2, 8=co, 16=ch4,
                               32=o2, 64=n2o                            */
  int32_t   atrem_full;     /* 1=do full, explicit atrem calc, 0=off    */
  int32_t   atrem_geom;     /* 1=do atrem geometry every pixel,
                               0=only calculate when error > 0.05 %     */
  int32_t   atrem_model;    /* 0=determine model from latitude and date
                               1-6 use one of the standard models       */
  int32_t   brdf_opt;       /* 1=Fresnel, 2=Rgoth, 3=f/Q*Rgoth          */
  int32_t   iop_opt;        /* Base IOP model for downstream products   */
  int32_t   pol_opt;        /* polarization correction option           */
  int32_t   rad_opt;        /* radcor switch for MERIS smile correction */
  int32_t   ocrvc_opt;      /* OCR-VC switch                            */
  int32_t   absaer_opt;     /* absorbing aerosol test option            */
  int32_t   aer_iter_max;   /* aerosol_iteration_limit                  */
  int32_t   glint_opt;      /* 1=apply glint correction                 */
  int32_t   outband_opt;    /* 1=apply seawifs out-of-band correction   */
  int32_t   oxaband_opt;    /* 1=apply seawifs 765 nm Oxygen corr       */
  int32_t   filter_opt;     /* 1=apply filtering in filter_file         */
  int32_t   cirrus_opt;     /* 1=apply cirrus cloud filtering           */
  int32_t   resolution;     /* process at this nadir pixel res (meters) */
                         /* 250, 500, 1000, -1=native (modis only)   */
  float  *taua;   /* Input aerosol optical thickness          */

  int32_t xcal_nwave;    /* number of wavelengths to which xcal applied */
  int32_t *xcal_opt; /* xcal option per band              */
  float   *xcal_wave; /* sensor wavelengths to which xcal applied */
  int32_t band_shift_opt; /* band fill, 0=lin.interp. 1=bio.opt.band shift */

  float add_noise_sigma; /* uncertainties, 0=none, !=0 => std dev for gaussian noise */
  float noise_scale; /* scaling factor for gain perturbation by scaling SNR*/ 
  int32_t wiggle_band; /*  band to be perturbed for analytical uncertainty calculation */
  float wiggle_by;/* noise level designated band perturbation */

  int32_t   sl_pixl;        /* seawifs straylight pixel limit           */
  float  sl_frac;        /* seawifs straylight Ltyp fraction         */

  char   filter_file [FILENAME_MAX];    /* filter specification file */
  fctlstr fctl;

  char   aerfile[FILENAME_MAX];         /* input aerosol spec file   */
  char   tgtfile[FILENAME_MAX];         /* input cal target file     */

  char   met1   [FILENAME_MAX];  /* Meteorological ancillary file    */
  char   met2   [FILENAME_MAX];  /* Meteorological ancillary file    */
  char   met3   [FILENAME_MAX];  /* Meteorological ancillary file    */
  char   ozone1 [FILENAME_MAX];  /* Ozone ancillary file             */
  char   ozone2 [FILENAME_MAX];  /* Ozone ancillary file             */
  char   ozone3 [FILENAME_MAX];  /* Ozone ancillary file             */
  char   anc_cor_file[FILENAME_MAX];  /* ancillary correction file   */
  char   land   [FILENAME_MAX];  /* Land mask file                   */
  char   water  [FILENAME_MAX];  /* Bathymetry mask file             */
  char   demfile[FILENAME_MAX];  /* Digital elevation map file       */
  char   elevfile[FILENAME_MAX]; /* Elevation file                   */
  char   elev_auxfile[FILENAME_MAX];  /* Auxiliary elevation file    */
  char   icefile[FILENAME_MAX];  /* Ice mask/fraction file           */
  char   sstfile[FILENAME_MAX];  /* SST reference file               */
  char   sssfile[FILENAME_MAX];  /* SSS file                         */
  char   no2file[FILENAME_MAX];  /* NO2 file                         */
  char   alphafile[FILENAME_MAX];/* angstrom climatology file        */
  char   tauafile[FILENAME_MAX]; /* AOT climatology file             */
  char   aerbinfile[FILENAME_MAX];  /* Bin file for aerosol inputs   */
  char   owmcfile[FILENAME_MAX]; /* Ocean water classification       */
  char   prodXMLfile[FILENAME_MAX]; /* product XML output file name  */
  char   breflectfile[FILENAME_MAX]; /* bottom reflectance input file*/
  int    sstreftype;

  float  *gain    ;       /* Vicarious calibration gain       */
  float  *gain_unc;       /* Vicarious gain uncertainty       */
  float  *offset  ;       /* Vicarious calibration offset     */
  float  albedo;                 /* cloud reflectance threshold      */
  float  cloud_wave;             /* cloud test wavelength            */
  float  cloud_eps;              /* cloud reflectance ratio          */
  float  glint;                  /* glint threshold                  */
  float  sunzen;                 /* solar zenith angle threshold     */
  float  satzen;                 /* sensor zenith angle threshold    */
  float  epsmin;                 /* min epsilon for atm corr failure */
  float  epsmax;                 /* max epsilon for atm corr failure */
  float  tauamax;                /* max tau 865 for hi-taua flagging */
  float  nlwmin;                 /* min nlw 555 for low lw flagging  */
  float  hipol;                  /* high polarization threshold      */
  float  wsmax;                  /* max windspeed for whitecap corr  */
  float  coccolith[8];           /* coccolithophore algorithm coefs. */
  float  absaer;                 /* threshold for abs aerosol index  */
  float  rhoamin;                /* low aerosol threshold            */
  float  cirrus_thresh[2];       /* cirrus reflectance thresholds    */

  float  windspeed;              /* use fixed windspeed as specified */
  float  windangle;              /* use fixed wind dir as specified  */
  float  pressure;               /* use fixed pressure as specified  */
  float  ozone;                  /* use fixed ozone as specified     */
  float  watervapor;             /* use fixed pr. water as specified */
  float  relhumid;               /* use fixed rh as specified        */
  float  ice_threshold;          /* fraction above which is flag ice */

  int32_t   landmask;               /* 0=off, 1=on */
  int32_t   bathmask;               /* 0=off, 1=on */
  int32_t   cloudmask;              /* 0=off, 1=on */
  int32_t   glintmask;              /* 0=off, 1=on */
  int32_t   sunzenmask;             /* 0=off, 1=on */
  int32_t   satzenmask;             /* 0=off, 1=on */
  int32_t   hiltmask;               /* 0=off, 1=on */
  int32_t   stlightmask;            /* 0=off, 1=on */

  char   program_name[128];
  char   pro_control[4096];
  char   input_parms[32768];
  char   input_files[6144];
  char   mask_names[1024];
  char   pversion[1024];
  char   rflag[1024];

  /* Vicarious calibration */
  float *vcal_nLw;
  float *vcal_Lw ;
  float vcal_chl;
  float vcal_solz;
  int   vcal_opt;
  float   vcal_depth;  /*  vcaltarget depth mask value */
  int32_t vcal_min_nbin;  /* min # samples in bin to accept */
  int32_t vcal_min_nscene;  /* min # scenes in bin to accept */

  /* MUMM control */
  float   mumm_alpha;             
  float   mumm_gamma;             
  float   mumm_epsilon;            

  /* QAA IOP model control */
  float   qaa_adg_s;
  int     qaa_wave[5];

  /* GSM IOP model control */
  int32_t gsm_opt;
  float   gsm_adg_s;
  float   gsm_bbp_s;
  float   *gsm_aphw;
  float   *gsm_aphs;
  int32_t gsm_fit;

  /* GIOP IOP model control */
  char   giop_aph_file[FILENAME_MAX];
  char   giop_adg_file[FILENAME_MAX];
  char   giop_bbp_file[FILENAME_MAX];
  char   giop_acdom_file[FILENAME_MAX];
  char   giop_anap_file[FILENAME_MAX];
  char   giop_bbph_file[FILENAME_MAX];
  char   giop_bbnap_file[FILENAME_MAX];
  int    giop_maxiter;
  int    giop_fit_opt;
  int    giop_aph_opt;
  int    giop_adg_opt;
  int    giop_acdom_opt;
  int    giop_anap_opt;
  int    giop_bbp_opt;
  int    giop_bbph_opt;
  int    giop_bbnap_opt;
  int    giop_rrs_opt;
  int    giop_iterate;
  float  giop_aph_s;
  float  giop_adg_s;
  float  giop_bbp_s;
  float  giop_aph_w;
  float  giop_adg_w;
  float  giop_bbp_w;
  float  giop_grd[2];
  float  giop_rrs_diff;
  float  *giop_wave;
  float  *giop_rrs_unc;

  /* empirical chlorophyll algorithm coeffs */
  int32_t chloc2w[2];
  float   chloc2c[5];
  int32_t chloc3w[3];
  float   chloc3c[5];
  int32_t chloc4w[4];
  float   chloc4c[5];

  int32_t kd2w[2];
  float   kd2c[6];

  float   flh_offset;

  /* sst stuff */
  int32_t viirsnv7; /* =1 to use the VIIRSN V7 high satz latband equation and coeffs */
  int32_t viirsnosisaf; /* =1 to use the VIIRSN OSI-SAF equation and coeffs (sort of v5 like) */
  float sstrefdif;   /* tighter threshold to match sst with reference */
  int32_t      newavhrrcal;  /* new avhrr calibration equation */
  float  ch22detcor[10];	 /* channel 22 detector corrections  */
  float  ch23detcor[10];	 /* channel 23 detector corrections  */

  /* the following fields support inverse (calibration) processing */
  char    il2file  [MAX_OFILES][FILENAME_MAX];
  char    flaguse[1024];
  float   chlthreshold;
  float   aotthreshold;
  float   maxpointdist;     /* Provide max distance between L1 and L2 pixels 
  			    (-1. - use average resolution of L1 data
			    default=max{L1 resolution, L2 resolution} */
  
  int32_t xcalbox;          /* Pixel size of the central box in the L1 scene (e.g. 5 pixels around MOBY) to be extracted into xcalfile, default=0-whole L1 */
  int32_t xcalboxcenter[2]; /* Centeral [ipix, iscan] of the box in the L1 scene, default =[0,0] - center of the L1 scene */
  int32_t xcalpervalid;     /* Minimum percent of valid cross-calibration pixels within the box or the L1 scene */
  int32_t xcalsubsmpl;      /* Subsampling rate for the data to be used for the cross-calibration  */

  /* the following fields support l1mapgen */
  int32_t    stype; /* scaling type 0=log, 1=linear*/
  float   datamin;
  float   datamax;
  float   west;
  float   east;
  float   north;
  float   south;
  int32_t   width;
  float	  threshold;
  int32_t   rgb[3];
  int	subsamp;
  int32_t  xbox;    /* number of pixels to retrieve around a point */
  int32_t  ybox;

  int32_t  deflate;
   
  int32_t raman_opt;  /*RAMAN Rrs correction model*/ 
  
  char   viirscalparfile[FILENAME_MAX];  /* VIIRS calibration parfile      */
} instr;

#ifdef __cplusplus
}
#endif

#endif



