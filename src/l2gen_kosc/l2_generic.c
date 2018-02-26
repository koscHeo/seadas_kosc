/* =========================================================== */
/* Module l2_generic.c                                         */
/*                                                             */
/* Functions to open and write a multi-sensor (generic) l2     */
/* file in HDF/NCDF format.                                    */
/*                                                             */ 
/* Written By:                                                 */
/*     Bryan A. Franz, SAIC GSC, March 1998.                   */
/*     Gary Fu,        SAIC GSC, March 1999.                   */
/*     Joel M. Gales,  Futuretech, Sept. 1999.                 */
/*     Gene Eplee, SAIC GSC, SeaWiFS Project, December 2000.   */
/*           Update time correction and mirror side            */
/*                       factors.                              */
/*     Gene Eplee, SAIC, SeaWiFS Project, March 2004.          */
/*           Convert time correction and mirror side           */
/*           factors to simultaneous exponentials.             */
/*     Joel Gales, Futuretech, OBPG Project, Sept 2012.        */
/*           Add support for L2 NETCDF4 output                 */
/*     Joel Gales, Futuretech, OBPG Project, Feb 2013.         */
/*           Add support for INT8,UINT8 datatypes for NETCDF4  */
/*     Joel Gales, Futuretech, OBPG Project, Nov 2013.         */
/*           Add support for CF-compliant metadata             */
/* =========================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h> 
#include <time.h>
#include <math.h>
#include "l12_proto.h"
#include "l2_generic.h"
#include <timeutils.h>
#include "l2prod.h"
#include "l1_aci_hdf.h"
#include "flags_sst.h"
#include "flags_iop.h"
#include "mph_flags.h"
#include "version.h"
#include <dfutils.h>


#define ROUND(x) ((((x) >= 0)?0.5:-0.5) + (x))

/* Global variables to facilitate communication */
static float bad_float  = BAD_FLT;
//static float fill_float = FIL_FLT;
//static float fill_int   = FIL_INT;
//static float fill_byte  = FIL_BYT;
static int32_t  *Lambda    , *Lambda_p;
static int32 numScans;
static int32 numPixels;
static int32 numBands;
static int32 numBandsIR;
static int32 spix; 
static int32 cpix; 
static int32 epix;
static int32 cscan; 
static int32 nctl; 
static int32 *ictl;
static int32 *jctl;

static FILE *fp_meta = NULL;

// MOVED TO l12_proto.h
//#define   INSTITUTION "NASA Goddard Space Flight Center, Ocean Ecology Laboratory, Ocean Biology Processing Group"
//
//#define   LICENSE "http://science.nasa.gov/earth-science/earth-science-data/data-information-policy/"
//#define   NAMING_AUTHORITY "gov.nasa.gsfc.sci.oceandata"
//#define   KEYWORDS_VOCABULARY "NASA Global Change Master Directory (GCMD) Science Keywords"
//#define   KEYWORDS_OC "Oceans > Ocean Chemistry > Chlorophyll; Oceans > Ocean Optics > Ocean Color"
//#define   KEYWORDS_SST "Oceans > Ocean Temperature > Sea Surface Temperature"
//#define   STDNAME_VOCABULARY "NetCDF Climate and Forecast (CF) Metadata Convenention"
//#define   CREATOR_NAME "NASA/GSFC/OBPG"
//#define   CREATOR_EMAIL "data@oceancolor.gsfc.nasa.gov"
//#define   CREATOR_URL "http://oceandata.sci.gsfc.nasa.gov"
//#define   PROJECT "Ocean Biology Processing Group (NASA/GSFC/OBPG)"
//#define   PUBLISHER_NAME "NASA/GSFC/OBPG"
//#define   PUBLISHER_EMAIL "data@oceancolor.gsfc.nasa.gov"
//#define   PUBLISHER_URL "http://oceandata.sci.gsfc.nasa.gov"

#define   GEOBOX_INC 20.0

/* -------------------------------------------------------- */
/* -------------------------------------------------------- */
#define RFACTOR 100
int32 get_ctl(int32_t ctl_pt_fact, int32 ictl[], int32 jctl[])
{
    int32 nctl;
    int32 dctl;
    int32 i;

    if (ctl_pt_fact <= 0) {

        if (numPixels < RFACTOR) {
            nctl = numPixels;
            for (i = 0; i < numPixels; i++) {
                ictl[i] = i+1;
            }

        } else {
            dctl = (int32_t) numPixels/RFACTOR;
            nctl = MIN((int32_t) numPixels/dctl + 1,numPixels);
            for (i = 0; i < (nctl-1); i++) {
                ictl[i] = i * dctl + 1;
            }
            ictl[nctl-1] = numPixels;
        }

    } else {
        dctl = ctl_pt_fact;
        nctl = MIN((int32_t) numPixels/dctl + 1,numPixels);

        for (i = 0; i < (nctl-1); i++) {
            ictl[i] = i * dctl + 1;
       }
        ictl[nctl-1] = numPixels;
    }

    for (i = 0; i < numScans; i++)
        jctl[i] = i+1;

    return(nctl);
}


/* -------------------------------------------------------- */
/* Assign SDSes to Vgroups                                  */
/* -------------------------------------------------------- */
int MakeVgroups(filehandle *l2file){

    int32	h_id;
    int32 v_id;
    int32 sd_id = l2file->sd_id;
    int   i;

    /* Do we have the extra meta-data for SeaWiFS */
    int seawifs_meta = 0;
    if (l2file->sensorID == SEAWIFS) {
        int32 sds_id;
        if ( sd_select(sd_id,"scan_ell",&sds_id) == 0 )
            seawifs_meta = 1;
    }

    h_id = Hopen(l2file->name, DFACC_RDWR, 0);
    if(h_id == FAIL){
        fprintf(stderr,"-E- %s line %d: Hopen() failed for file, %s.\n",
                __FILE__,__LINE__,l2file->name);
        return(HDF_FUNCTION_ERROR);
    }
    Vstart(h_id);

    /* Scan-Line Attributes */
    DPTB( v_attach(h_id, &v_id)			);
    Vsetclass(v_id, "Per File Data");
    Vsetname(v_id, "Sensor Band Parameters");
    DPTB( AddSdsToVgroup(sd_id, v_id, "wavelength")	);
    DPTB( AddSdsToVgroup(sd_id, v_id, "vcal_gain")	);
    DPTB( AddSdsToVgroup(sd_id, v_id, "vcal_offset")	);
    DPTB( AddSdsToVgroup(sd_id, v_id, "F0")	  	);
    DPTB( AddSdsToVgroup(sd_id, v_id, "k_oz")		);
    DPTB( AddSdsToVgroup(sd_id, v_id, "Tau_r")		);
    Vdetach(v_id);

    if (seawifs_meta) {

        /* Sensor Tilt */
        DPTB( v_attach(h_id, &v_id)			);
        Vsetclass(v_id, "Per File Data");
        Vsetname(v_id, "Sensor Tilt");
        DPTB( AddSdsToVgroup(sd_id, v_id, "ntilts")		);
        DPTB( AddSdsToVgroup(sd_id, v_id, "tilt_flags")	);
        DPTB( AddSdsToVgroup(sd_id, v_id, "tilt_ranges")	);
        Vdetach(v_id);

    }

    /* Scan-Line Attributes */
    DPTB( v_attach(h_id, &v_id)			);
    Vsetclass(v_id, "Per Scan Data");
    Vsetname(v_id, "Scan-Line Attributes");
    DPTB( AddSdsToVgroup(sd_id, v_id, "year")		);
    DPTB( AddSdsToVgroup(sd_id, v_id, "day")		);
    DPTB( AddSdsToVgroup(sd_id, v_id, "msec")		);
    DPTB( AddSdsToVgroup(sd_id, v_id, "slon")		);
    DPTB( AddSdsToVgroup(sd_id, v_id, "clon")		);
    DPTB( AddSdsToVgroup(sd_id, v_id, "elon")		);
    DPTB( AddSdsToVgroup(sd_id, v_id, "slat")		);
    DPTB( AddSdsToVgroup(sd_id, v_id, "clat")		);
    DPTB( AddSdsToVgroup(sd_id, v_id, "elat")		);
    DPTB( AddSdsToVgroup(sd_id, v_id, "csol_z")		);
    Vdetach(v_id);

    /* Geophysical Data */
    DPTB( v_attach(h_id, &v_id));
    Vsetclass(v_id, "Per Scan Data");
    Vsetname(v_id, "Geophysical Data");
    for (i=0; i<l2file->tot_prod; i++) {
        DPTB( AddSdsToVgroup(sd_id, v_id, l2file->l2_prod_names[i])	);
    }
    Vdetach(v_id);

    /* Navigation */
    DPTB( v_attach(h_id, &v_id)			);
    Vsetclass(v_id, "Per Scan Data");
    Vsetname(v_id, "Navigation Data");
    DPTB( AddSdsToVgroup(sd_id, v_id, "longitude")	);
    DPTB( AddSdsToVgroup(sd_id, v_id, "latitude")		);
    DPTB( AddSdsToVgroup(sd_id, v_id, "cntl_pt_cols")	);
    DPTB( AddSdsToVgroup(sd_id, v_id, "cntl_pt_rows")	);
    DPTB( AddSdsToVgroup(sd_id, v_id, "tilt")		);
    if (seawifs_meta) {
        DPTB( AddSdsToVgroup(sd_id, v_id, "orb_vec")	);
        DPTB( AddSdsToVgroup(sd_id, v_id, "sun_ref")	);
        DPTB( AddSdsToVgroup(sd_id, v_id, "att_ang")	);
        DPTB( AddSdsToVgroup(sd_id, v_id, "sen_mat")	);
        DPTB( AddSdsToVgroup(sd_id, v_id, "scan_ell")	);
        DPTB( AddSdsToVgroup(sd_id, v_id, "nflag")	);
    }
    Vdetach(v_id);

    Vend(h_id);

    if(Hclose(h_id) != SUCCEED){
        fprintf(stderr,"-E- %s line %d: Hclose(%d) failed for file, %s .\n",
                __FILE__,__LINE__,h_id,l2file->name);
        return(HDF_FUNCTION_ERROR);
    }
    return(LIFE_IS_GOOD);
}


/*----------------------------------------------------------------- */
/* Open an L2 file and store global attributes in it.               */
/* ---------------------------------------------------------------- */
int openl2( filehandle *l2file)   
{
    float  *Gain  ;
    float  *Offset;
    float  *Fonom ;
    float  *Fobar , *Fobar_p;
    float  *Tau_r , *Tau_r_p;
    float  *k_oz  , *k_oz_p;
    float  *k_no2 , *k_no2_p;
    float  *aw    ;
    float  *bbw   ;

    char   title[256];
    char   soft_id[200];     /* software version info */
    int32_t   sd_id;
    int32_t   sds_id;
    int32_t   tot_prod;
    VOIDP   pbuf;

    int    i, k, ii;
    static int firstCall=1;
    l2prodstr *p;
    char   *tmp_ptr, tmp_str[2048];
    int32_t    n;
    char    wavestr[5];
    char    prefix[256];
    char    suffix[256];
    char    avhrrbird[10];
    int     flagbits[32];
    int16_t sstflagbits[16];
    char    buf1[1024], buf2[32];
    char   *p1, *p2;
    int32_t    nbvis;
    int32_t year, day, hr, mn, sc;
    int16_t month, dom;

    idDS ds_id;

    int32 dm[3];
    const char dm_name[3][80];
    char *end_str;


    // strings for the different named global attributes
    char* calibrationDataStr;
    char* equatorCrossingLonStr;
    char* historyStr;
    char* inputFilesStr;
    char* maskNamesStr;
    char* orbitNumberStr;
    char* pixelControlPointsStr;
    char* processingVersionStr;
    char* productNameStr;
    char* scanControlPointsStr;
    char* softwareNameStr;
    char* softwareVersionStr;
    char* titleStr;

    char* numberOfScanLinesStr;
    char* pixelsPerScanLineStr;
    char* totalBandNumberStr;
    char* bandNumberStr;

    if (strcmp(l2file->input->metafile,"") != 0) {
        fp_meta = fopen(l2file->input->metafile,"w");
        if (fp_meta == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: Unable to open specified meta-data file, %s .\n",
                    __FILE__,__LINE__,l2file->input->metafile);
            return(HDF_FUNCTION_ERROR);
        }
    }   

    /* Create the L2 file */
    ds_format_t fileFormat;
    int32 nt_chr, nt_i16, nt_i32, nt_f32;
    if ( l2file->format == FMT_L2HDF) {
        fileFormat = DS_HDF;
        nt_chr = DFNT_CHAR;
        nt_i16 = DFNT_INT16;
        nt_i32 = DFNT_INT32;
        nt_f32 = DFNT_FLOAT32;

        calibrationDataStr = "Calibration Data";
        equatorCrossingLonStr = "Orbit Node Longitude";
        historyStr = "Processing Control";
        inputFilesStr = "Input Files";
        maskNamesStr = "Mask Names";
        orbitNumberStr = "Orbit Number";
        pixelControlPointsStr = "Number of Pixel Control Points";
        processingVersionStr = "Processing Version";
        productNameStr = "Product Name";
        scanControlPointsStr = "Number of Scan Control Points";
        softwareNameStr = "Software Name";
        softwareVersionStr = "Software Version";
        titleStr = "Title";
        numberOfScanLinesStr = "Number of Scan Lines";
        pixelsPerScanLineStr = "Pixels per Scan Line";
        totalBandNumberStr = "total band number";
        bandNumberStr = "band number";


    } else if ( l2file->format == FMT_L2NCDF) {
        fileFormat = DS_NCDF;
        nt_chr = NC_CHAR;
        nt_i16 = NC_SHORT;
        nt_i32 = NC_INT;
        nt_f32 = NC_FLOAT;

        calibrationDataStr = "calibration_data";
        equatorCrossingLonStr = "equatorCrossingLongitude";
        historyStr = "history";
        inputFilesStr = "source";
        maskNamesStr = "mask_names";
        orbitNumberStr = "orbit_number";
        pixelControlPointsStr = "pixel_control_points";
        processingVersionStr = "processing_version";
        productNameStr = "product_name";
        scanControlPointsStr = "scan_control_points";
        softwareNameStr = "software_name";
        softwareVersionStr = "software_version";
        titleStr = "title";
        numberOfScanLinesStr = "number_of_lines";
        pixelsPerScanLineStr = "pixels_per_line";
        totalBandNumberStr = "number_of_bands";
        bandNumberStr = "number_of_reflective_bands";
    }
    ds_id = startDS(l2file->name, fileFormat, DS_WRITE, l2file->deflate);
    if ( ds_id.fid == FAIL) return(HDF_FUNCTION_ERROR); 
    l2file->sd_id = ds_id.fid;

    /* Get number of bands, pixels, scans */
    numPixels = l2file->npix;
    numScans = l2file->nscan;
    numBands = l2file->nbands;
    numBandsIR = l2file->nbandsir;
    spix     = 0;
    cpix     = numPixels/2;
    epix     = numPixels-1;
    cscan    = numScans/2;

    if (firstCall == 1) {
        firstCall = 0;
        if ( (Lambda = calloc(numBands+NBANDSIR,sizeof(int32_t))) == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: Unable to allocate control Lambda array.\n",
                    __FILE__,__LINE__);
            return(MEMORY_ALLOCATION_ERROR);
        }

    }


    if ((Gain = (float *) calloc(numBands,sizeof(float ))) == NULL) {
        printf("-E- %s line %d : error allocating memory for l2_generic:openl2.\n",
                __FILE__,__LINE__);
        exit(1);
    }
    if ((Offset = (float *) calloc(numBands,sizeof(float ))) == NULL) {
        printf("-E- %s line %d : error allocating memory for l2_generic:openl2.\n",
                __FILE__,__LINE__);
        exit(1);
    }
    if ((Fonom = (float *) calloc(numBands,sizeof(float ))) == NULL) {
        printf("-E- %s line %d : error allocating memory for l2_generic:openl2.\n",
                __FILE__,__LINE__);
        exit(1);
    }
    if ((Fobar = (float *) calloc(numBands,sizeof(float ))) == NULL) {
        printf("-E- %s line %d : error allocating memory for l2_generic:openl2.\n",
                __FILE__,__LINE__);
        exit(1);
    }
    if ((Tau_r = (float *) calloc(numBands,sizeof(float ))) == NULL) {
        printf("-E- %s line %d : error allocating memory for l2_generic:openl2.\n",
                __FILE__,__LINE__);
        exit(1);
    }
    if ((k_oz = (float *) calloc(numBands,sizeof(float ))) == NULL) {
        printf("-E- %s line %d : error allocating memory for l2_generic:openl2.\n",
                __FILE__,__LINE__);
        exit(1);
    }
    if ((k_no2 = (float *) calloc(numBands,sizeof(float ))) == NULL) {
        printf("-E- %s line %d : error allocating memory for l2_generic:openl2.\n",
                __FILE__,__LINE__);
        exit(1);
    }
    if ((aw = (float *) calloc(numBands,sizeof(float ))) == NULL) {
        printf("-E- %s line %d : error allocating memory for l2_generic:openl2.\n",
                __FILE__,__LINE__);
        exit(1);
    }
    if ((bbw = (float *) calloc(numBands,sizeof(float ))) == NULL) {
        printf("-E- %s line %d : error allocating memory for l2_generic:openl2.\n",
                __FILE__,__LINE__);
        exit(1);
    }

    /* Get control-point array */
    if ( (ictl = calloc(numPixels,sizeof(int32))) == NULL) {
        fprintf(stderr,
                "-E- %s line %d: Unable to allocate control-point array.\n",
                __FILE__,__LINE__);
        return(MEMORY_ALLOCATION_ERROR);
    }
    if ( (jctl = calloc(numScans,sizeof(int32))) == NULL) {
        fprintf(stderr,
                "-E- %s line %d: Unable to allocate %d scan control-point array.\n",
                __FILE__,__LINE__,numScans);
        return(MEMORY_ALLOCATION_ERROR);
    }
    nctl = get_ctl(l2file->ctl_pt_incr,ictl,jctl);

    /* Get sensor-specific attributes */
    if ((n = rdsensorinfo(l2file->sensorID,l2file->input->evalmask,NULL,NULL)) != numBands) {
        fprintf(stderr,"-E- %s Line %d:  Error reading sensor table. %d %d\n",
                __FILE__,__LINE__,numBands,n);
        return(-1);
    }
    rdsensorinfo(l2file->sensorID,l2file->input->evalmask,"Lambda",(void **) &Lambda_p);
    rdsensorinfo(l2file->sensorID,l2file->input->evalmask,"Fobar", (void **) &Fobar_p );
    rdsensorinfo(l2file->sensorID,l2file->input->evalmask,"Tau_r", (void **) &Tau_r_p );
    rdsensorinfo(l2file->sensorID,l2file->input->evalmask,"k_oz",  (void **) &k_oz_p  );
    rdsensorinfo(l2file->sensorID,l2file->input->evalmask,"k_no2",  (void **) &k_no2_p  );


    for (i=0; i<numBands; i++) {
        Lambda [i] = Lambda_p[i];
        // multiply by 10 to put into W/m2/nm, since internally all radiances are mW/cm2/um
        Fobar  [i] = Fobar_p [i] * 10.0;
        Tau_r  [i] = Tau_r_p [i];   
        k_oz   [i] = k_oz_p  [i];
        k_no2  [i] = k_no2_p [i];
        aw     [i] = aw_spectra (Lambda[i],BANDW);
        bbw    [i] = bbw_spectra(Lambda[i],BANDW);
        Gain   [i] = l2file->input->gain[i];
        Offset [i] = l2file->input->offset[i];
        if (l2file->input->outband_opt >= 2) {
            get_f0_thuillier_ext(Lambda[i],BANDW,&Fonom[i]);
            Fonom[i] *= 10.0;
        } else {
            Fonom[i] = Fobar[i];
        }
    }
    for (i=numBands; i<numBands+numBandsIR; i++)
        Lambda [i] = Lambda_p[i];    

    /*                                                                  */
    /* Build list of dataset names from input parameters                */
    /* ---------------------------------------------------------------- */
    /*                                                                  */
    tot_prod = prodlist(l2file->sensorID,l2file->input->evalmask,l2file->l2prod,
            l2file->def_l2prod,l2file->l2_prod_names);

    strcpy(l2file->l2_prod_names[tot_prod++],"l2_flags"); 

    printf("\n\nThe following products will be included in %s.\n",l2file->name);
    for (i=0; i<tot_prod; i++)
        printf("%d %s\n",i,l2file->l2_prod_names[i]);

    l2file->tot_prod = tot_prod;

    // cache the product structures
    l2file->prodptr = (l2prodstr*) allocateMemory(tot_prod * sizeof(l2prodstr), "l2file->prodptr");
    for (i=0; i<tot_prod; i++) {
        if ((p = get_l2prod_index(l2file->l2_prod_names[i],l2file->sensorID,
                numBands+numBandsIR,numPixels,numScans,Lambda)) == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: product index failure.\n",
                    __FILE__,__LINE__);
            return(1);
        }
        l2file->prodptr[i] = *p; // need to actually copy all the memory
    }

    if ( l2file->format == FMT_L2NCDF) {
        int dumdim;
        if ( nc_def_dim(ds_id.fid, numberOfScanLinesStr, numScans, &dumdim)
                != NC_NOERR) exit(1);
        if ( nc_def_dim(ds_id.fid, pixelsPerScanLineStr, numPixels, &dumdim)
                != NC_NOERR) exit(1);
        if ( nc_def_dim(ds_id.fid, pixelControlPointsStr, nctl, &dumdim)
                != NC_NOERR) exit(1);
        if ( nc_def_dim(ds_id.fid, totalBandNumberStr, numBands+numBandsIR, &dumdim)
                != NC_NOERR) exit(1);
        if ( nc_def_dim(ds_id.fid, bandNumberStr, numBands, &dumdim)
                != NC_NOERR) exit(1);

        nc_def_grp( ds_id.fid, "sensor_band_parameters", &l2file->grp_id[0]);
        nc_def_grp( ds_id.fid, "scan_line_attributes", &l2file->grp_id[2]);
        nc_def_grp( ds_id.fid, "geophysical_data", &l2file->grp_id[3]);
        nc_def_grp( ds_id.fid, "navigation_data", &l2file->grp_id[4]);
        nc_def_grp( ds_id.fid, "processing_control", &l2file->grp_id[5]);
    }


    /*                                                                  */
    /* Create the scan-line datasets                                    */
    /* ---------------------------------------------------------------- */
    /*                                                                  */
    if ( l2file->format == FMT_L2NCDF) ds_id.fid = l2file->grp_id[2];

    dm[0] = numScans;
    strcpy((char *) dm_name[0], numberOfScanLinesStr);

    PTB( createDS(ds_id, (int) l2file->sensorID, "year", dm, dm_name));
    PTB( createDS(ds_id, (int) l2file->sensorID, "day", dm, dm_name));
    PTB( createDS(ds_id, (int) l2file->sensorID, "msec", dm, dm_name));
    PTB( createDS(ds_id, (int) l2file->sensorID, "detnum", dm, dm_name));
    PTB( createDS(ds_id, (int) l2file->sensorID, "mside", dm, dm_name));
    PTB( createDS(ds_id, (int) l2file->sensorID, "slon", dm, dm_name));
    PTB( createDS(ds_id, (int) l2file->sensorID, "clon", dm, dm_name));
    PTB( createDS(ds_id, (int) l2file->sensorID, "elon", dm, dm_name));
    PTB( createDS(ds_id, (int) l2file->sensorID, "slat", dm, dm_name));
    PTB( createDS(ds_id, (int) l2file->sensorID, "clat", dm, dm_name));
    PTB( createDS(ds_id, (int) l2file->sensorID, "elat", dm, dm_name));
    PTB( createDS(ds_id, (int) l2file->sensorID, "csol_z", dm, dm_name));

    if ( l2file->format == FMT_L2NCDF) ds_id.fid = l2file->grp_id[4];

    dm[1] = nctl;
    strcpy((char *) dm_name[1], pixelControlPointsStr);

    PTB( createDS(ds_id, (int) l2file->sensorID, "longitude", dm, dm_name));
    PTB( createDS(ds_id, (int) l2file->sensorID, "latitude", dm, dm_name));

    dm[0] = nctl;
    strcpy((char *) dm_name[0], pixelControlPointsStr);
    PTB( createDS(ds_id, (int) l2file->sensorID, "cntl_pt_cols", dm, dm_name));

    dm[0] = numScans;
    strcpy((char *) dm_name[0], numberOfScanLinesStr);
    PTB( createDS(ds_id, (int) l2file->sensorID, "cntl_pt_rows", dm, dm_name));
    PTB( createDS(ds_id, (int) l2file->sensorID, "tilt", dm, dm_name));

    /*                                                                  */
    /* Create the geophysical datasets                                  */
    /* ---------------------------------------------------------------- */
    /*                                                                  */
    if ( l2file->format == FMT_L2NCDF) ds_id.fid = l2file->grp_id[3];

    dm[0] = numScans;
    strcpy((char *) dm_name[0], numberOfScanLinesStr);
    dm[1] = numPixels;
    strcpy((char *) dm_name[1], pixelsPerScanLineStr);

    for (i=0; i<tot_prod; i++) {
        // Skip parameters already included if user requested them specifically
        if (!strcmp(l2file->l2_prod_names[i],"detnum") ||
            !strcmp(l2file->l2_prod_names[i],"mside")  ||
            !strcmp(l2file->l2_prod_names[i],"year")  ||
            !strcmp(l2file->l2_prod_names[i],"day")  ||
            !strcmp(l2file->l2_prod_names[i],"msec")  ||
            !strcmp(l2file->l2_prod_names[i],"slon")  ||
            !strcmp(l2file->l2_prod_names[i],"clon")  ||
            !strcmp(l2file->l2_prod_names[i],"elon")  ||
            !strcmp(l2file->l2_prod_names[i],"slat")  ||
            !strcmp(l2file->l2_prod_names[i],"clat")  ||
            !strcmp(l2file->l2_prod_names[i],"elat")  ||
            !strcmp(l2file->l2_prod_names[i],"csol_z") ||
            !strcmp(l2file->l2_prod_names[i],"latitude")  ||
            !strcmp(l2file->l2_prod_names[i],"longitude")  ||
            !strcmp(l2file->l2_prod_names[i],"tilt")
        )
            continue;

        PTB( createDS(ds_id, (int) l2file->sensorID, l2file->l2_prod_names[i],
                dm, dm_name));

        ds_id.sid = selectDS(ds_id, l2file->l2_prod_names[i]);

        /*                                                              */
        /* Add flag-name attributes if this is a flag product           */
        /*                                                              */
        if ((strcmp(l2file->l2_prod_names[i],"l2_flags") == 0) ||
                (strcmp(l2file->l2_prod_names[i],"flags_sst") == 0) ||
                (strcmp(l2file->l2_prod_names[i],"flags_sst_triple") == 0) ||
                (strcmp(l2file->l2_prod_names[i],"flags_sst4") == 0) ||
                (strcmp(l2file->l2_prod_names[i],"qual_sst") == 0) ||
                (strcmp(l2file->l2_prod_names[i],"qual_sst_triple") == 0) ||
                (strcmp(l2file->l2_prod_names[i],"qual_sst4") == 0) ||
                (strcmp(l2file->l2_prod_names[i],"MPH_flags") == 0)) {

            if ((strcmp(l2file->l2_prod_names[i],"l2_flags") == 0)){
                tmp_str[0] = '\0';
                for (k=0; k<NFLAGS; k++) {
                    strcat(tmp_str,l2_flag_lname[k]);
                    flagbits[k] = pow(2,k);
                    if (k < NFLAGS-1)
                        strcat(tmp_str," ");
                    /*
                     * Keep the old flag attribute set for HDF4 files
                     */
                    if (l2file->format == FMT_L2HDF){
                        PTB( setAttr( ds_id, l2_flag_sname[k], nt_chr,
                                strlen(l2_flag_lname[k])+1, (VOIDP)l2_flag_lname[k]) );
                    }
                }
                PTB( setAttr( ds_id, "flag_masks", nt_i32, NFLAGS, (VOIDP)flagbits));

            } else if ((strcmp(l2file->l2_prod_names[i],"flags_sst") == 0) ||
                    (strcmp(l2file->l2_prod_names[i],"flags_sst_triple") == 0) ||
                    (strcmp(l2file->l2_prod_names[i],"flags_sst4") == 0) ) {
                tmp_str[0] = '\0';
                for (k=0; k<NSSTFLAGS; k++) {
                    sstflagbits[k] =  pow(2,k);
                    if (l2file->sensorID == AVHRR) {
                        strcat(tmp_str,avhrr_sst_flag_lname[k]);
                    } else if (l2file->sensorID == VIIRS) {
                        strcat(tmp_str,viirs_sst_flag_lname[k]);
                    } else {
                        strcat(tmp_str,sst_flag_lname[k]);
                    }
                    if (k < NSSTFLAGS-1)
                        strcat(tmp_str," ");
                }
                PTB( setAttr( ds_id, "flag_masks", nt_i16, NSSTFLAGS, (VOIDP)sstflagbits));

            } else if ((strcmp(l2file->l2_prod_names[i],"qual_sst") == 0) ||
 		    (strcmp(l2file->l2_prod_names[i],"qual_sst_triple") == 0)  ||
                    (strcmp(l2file->l2_prod_names[i],"qual_sst4") == 0) ) {
                tmp_str[0] = '\0';
                for (k=0; k<NQSSTFLAGS; k++) {
                    sstflagbits[k] =  k;
                    strcat(tmp_str,qual_sst_flag_lname[k]);
                    if (k < NQSSTFLAGS-1)
                        strcat(tmp_str," ");
                }
                PTB( setAttr( ds_id, "flag_masks", nt_i16, NQSSTFLAGS, (VOIDP)sstflagbits));

            } else if ((strcmp(l2file->l2_prod_names[i],"MPH_flags") == 0) ) {
                tmp_str[0] = '\0';
                for (k=0; k<NMPHFLAGS; k++) {
                    sstflagbits[k] =  pow(2,k);
                    strcat(tmp_str,mph_flag_lname[k]);
                    if (k < NMPHFLAGS-1)
                        strcat(tmp_str," ");
                }
                PTB( setAttr( ds_id, "mph_flags", nt_i16, NMPHFLAGS, (VOIDP)sstflagbits));
            }
            PTB( setAttr( ds_id, "flag_meanings", nt_chr, strlen(tmp_str)+1, (VOIDP)tmp_str));
            
        }
        /*                                                              */
        /* Add solar irradiance meta-data if this is the nLw or Rrs     */
        /*                                                              */
        p = l2file->prodptr + i; // get product structure from cache
        switch (p->cat_ix) {
          case CAT_nLw:
          case CAT_Rrs:
            PTB( setAttr( ds_id, "solar_irradiance", nt_f32, 1, &Fonom[p->prod_ix]));
        }
        /*                                                            */
        /* Add bad value attributes to HDF4 files                     */
        /*                                                            */
        if ((l2file->format == FMT_L2HDF) && ((strcmp(l2file->l2_prod_names[i],"l2_flags") != 0) ||
                (strcmp(l2file->l2_prod_names[i],"flags_sst") != 0) ||
                (strcmp(l2file->l2_prod_names[i],"flags_sst_triple") != 0) ||
                (strcmp(l2file->l2_prod_names[i],"flags_sst4") != 0))) {

            switch (p->datatype) {
              case DFNT_UINT8:
                pbuf = (VOIDP) float2uint8((VOIDP)&bad_float,0,1,1,p->slope,p->offset);
                setAttr( ds_id,"bad_value_scaled", nt_chr,    1, pbuf);
                pbuf = (VOIDP) unscale_sds(pbuf,p,0,1,1);
                setAttr( ds_id,"bad_value_unscaled", nt_f32, 1, pbuf);
                break;
              case DFNT_INT16:
                pbuf = (VOIDP) float2int16((VOIDP)&bad_float,0,1,1,p->slope,p->offset);
                setAttr( ds_id,"bad_value_scaled", nt_i16,   1, pbuf);
                pbuf = (VOIDP) unscale_sds(pbuf,p,0,1,1);
                setAttr( ds_id,"bad_value_unscaled", nt_f32, 1, pbuf);
                break;
              case DFNT_INT32:
                break;
              case DFNT_FLOAT32:
                setAttr( ds_id,"bad_value_scaled", nt_f32, 1, &bad_float);
                setAttr( ds_id,"bad_value_unscaled", nt_f32, 1, &bad_float);
                break;
              default:
                break;
            }
        }
        
        endaccessDS(ds_id);        
    }

    /*                                                                  */
    /* Create the Sensor Band Parameters datasets                       */
    /* ---------------------------------------------------------------- */
    /*                                                                  */
    if ( l2file->format == FMT_L2NCDF) ds_id.fid = l2file->grp_id[0];

    dm[0] = numBands+numBandsIR;
    strcpy((char *) dm_name[0], totalBandNumberStr);
    PTB( createDS(ds_id, l2file->sensorID, "wavelength", dm, dm_name));

    dm[0] = numBands;
    strcpy((char *) dm_name[0], bandNumberStr);

    if (numBands > 0) {
        PTB( createDS(ds_id, l2file->sensorID, "vcal_gain", dm, dm_name));
        PTB( createDS(ds_id, l2file->sensorID, "vcal_offset", dm, dm_name));
        PTB( createDS(ds_id, l2file->sensorID, "F0", dm, dm_name));
        PTB( createDS(ds_id, l2file->sensorID, "aw", dm, dm_name));
        PTB( createDS(ds_id, l2file->sensorID, "bbw", dm, dm_name));
        PTB( createDS(ds_id, l2file->sensorID, "k_oz", dm, dm_name));
        PTB( createDS(ds_id, l2file->sensorID, "k_no2", dm, dm_name));
        PTB( createDS(ds_id, l2file->sensorID, "Tau_r", dm, dm_name));
    }

    /*                                                                  */
    /* Write out some global attributes                                 */
    /* ---------------------------------------------------------------- */
    /*                                                                  */
    sprintf(title,"%s Level-2 Data",sensorName[l2file->sensorID]);
    //strcpy(soft_id, VERSION);
    sprintf(soft_id, "%d.%d.%d-r%d",L2GEN_VERSION_MAJOR,L2GEN_VERSION_MINOR,L2GEN_VERSION_PATCH_LEVEL,SVN_REVISION);

    if ( l2file->format == FMT_L2NCDF) ds_id.fid = l2file->sd_id;

    PTB( SetChrGA (ds_id, titleStr, title)                               );
    PTB( SetChrGA (ds_id, productNameStr, basename(l2file->name))       );
    PTB( SetChrGA (ds_id, processingVersionStr, l2file->input->pversion));
    if (l2file->orbit_node_lon > -180.0 && l2file->orbit_node_lon < 180.0)
        PTB( SetF32GA (ds_id, equatorCrossingLonStr, l2file->orbit_node_lon)  );
    if (l2file->orbit_number > 0)
        PTB( SetI32GA (ds_id, orbitNumberStr, l2file->orbit_number)  );
    PTB( SetChrGA (ds_id, historyStr, l2file->pro_control)    );
    if ( l2file->format == FMT_L2HDF) {
        PTB( SetI32GA (ds_id, scanControlPointsStr, numScans)     );
        PTB( SetI32GA (ds_id, pixelControlPointsStr, nctl)        );
    }

    // Processing Control attibutes
    if ( l2file->format == FMT_L2NCDF) ds_id.fid = l2file->grp_id[5];
    PTB( SetChrGA (ds_id, softwareNameStr, PROGRAM)                     );
    PTB( SetChrGA (ds_id, softwareVersionStr, soft_id)                  );
    PTB( SetChrGA (ds_id, inputFilesStr, l2file->input_files)           );

    // cleanup and populate calfile metadata
    char *calfile_str;
    if ((calfile_str = malloc(strlen(l2file->calfile) + 1)) == NULL) {
        fprintf(stderr, "-E- %s line %d: Unable to copy calfile string.\n",
        __FILE__, __LINE__);
        exit(1);
    }
    strcpy(calfile_str,"\0");
    char *tmp_calfile;
    if ((tmp_calfile = strdup(l2file->calfile)) == NULL) {
        fprintf(stderr, "-E- %s line %d: Unable to copy calfile string.\n",
        __FILE__, __LINE__);
        exit(1);
    }

    char *token = strtok_r(tmp_calfile, ",", &end_str);
    while (token != NULL) {
        strcpy(tmp_str, token);
        trimBlanks(tmp_str);
        strcat(calfile_str, basename(tmp_str));
        token = strtok_r(NULL, ",", &end_str);
        if (token != NULL)
            strcat(calfile_str, ", ");
    }
    free(tmp_calfile);

    PTB( SetChrGA (ds_id, calibrationDataStr, calfile_str));

    PTB( SetChrGA (ds_id, maskNamesStr, l2file->mask_names));

    if(l2file->format == FMT_L2NCDF) {
        // write out netCDF specific attributes

        // global attr
        ds_id.fid = l2file->sd_id;
        PTB( SetChrGA (ds_id, "instrument" , instrumentName[l2file->sensorID]) );

        if (l2file->sensorID == AVHRR) {
            strcpy(avhrrbird,"NOAA-");
            strncat(avhrrbird,xsatid2name(l2file->subsensorID)+2,2);
            PTB( SetChrGA (ds_id, "platform" , avhrrbird) );
        } else {
            PTB( SetChrGA (ds_id, "platform" , platformName[l2file->sensorID]) );
        }
        PTB( SetChrGA (ds_id, "Conventions", "CF-1.6")                      );
        PTB( SetChrGA (ds_id, "Metadata_Conventions","Unidata Dataset Discovery v1.0"));
        PTB( SetChrGA (ds_id, "license", LICENSE)                           );
        PTB( SetChrGA (ds_id, "naming_authority", NAMING_AUTHORITY)         );
        // create the id -
        if (strcmp(l2file->input->pversion,"Unspecified") != 0){
            strcpy(buf1,l2file->input->pversion);
            strcat(buf1,"/L2/");
        } else {
            strcpy(buf1,"L2/");
        }
        strcat(buf1,basename(l2file->name));
        PTB( SetChrGA (ds_id, "id", buf1 ));


        time_t tnow;
        time(&tnow);
        strcpy(buf1, unix2isodate(tnow, 'G'));
        PTB( SetChrGA (ds_id, "date_created", buf1)                         );

        PTB( SetChrGA (ds_id, "keywords_vocabulary", KEYWORDS_VOCABULARY)   );
        if (strstr(l2file->l2prod,"sst")){
            PTB( SetChrGA (ds_id, "keywords", KEYWORDS_SST)                 );
        } else {
            PTB( SetChrGA (ds_id, "keywords", KEYWORDS_OC)                  );
        }
        PTB( SetChrGA (ds_id, "standard_name_vocabulary", STDNAME_VOCABULARY)     );
        PTB( SetChrGA (ds_id, "institution", INSTITUTION)                   );
        PTB( SetChrGA (ds_id, "creator_name", CREATOR_NAME)                 );
        PTB( SetChrGA (ds_id, "creator_email", CREATOR_EMAIL)               );
        PTB( SetChrGA (ds_id, "creator_url", CREATOR_URL)                   );
        PTB( SetChrGA (ds_id, "project", PROJECT)                           );
        PTB( SetChrGA (ds_id, "publisher_name", PUBLISHER_NAME)             );
        PTB( SetChrGA (ds_id, "publisher_url", PUBLISHER_URL)               );
        PTB( SetChrGA (ds_id, "publisher_email", PUBLISHER_EMAIL)           );
        //Some missions have DOIs
        if (l2file->sensorID == SEAWIFS
                || l2file->sensorID == OCTS
                || l2file->sensorID == CZCS
                || l2file->sensorID == HMODISA
                || l2file->sensorID == HMODIST) {
            PTB( SetChrGA (ds_id, "identifier_product_doi_authority","http://dx.doi.org"));

            switch (l2file->sensorID) {
            case SEAWIFS:
                PTB( SetChrGA (ds_id, "identifier_product_doi","10.5067/ORBVIEW-2/SEAWIFS_OC.2014.0"));
                break;
            case OCTS:
                PTB( SetChrGA (ds_id, "identifier_product_doi","10.5067/ADEOS/OCTS_OC.2014.0"));
                break;
            case CZCS:
                PTB( SetChrGA (ds_id, "identifier_product_doi","10.5067/NIMBUS-7/CZCS_OC.2014.0"));
                break;
            case HMODISA:
                PTB( SetChrGA (ds_id, "identifier_product_doi","10.5067/AQUA/MODIS_OC.2014.0"));
                break;
            case HMODIST:
                PTB( SetChrGA (ds_id, "identifier_product_doi","10.5067/TERRA/MODIS_OC.2014.0"));
                break;
            }
        }
        PTB( SetChrGA (ds_id, "processing_level", "L2")                     );
        PTB( SetChrGA (ds_id, "cdm_data_type", "swath")                     );
        if (strlen(l2file->node_crossing_time) > 1) {
            double eqcrosstm = zulu2unix(l2file->node_crossing_time);
            strcpy(buf1, unix2isodate(eqcrosstm, 'G'));
            PTB( SetChrGA (ds_id, "equatorCrossingDateTime", buf1)  );
        }

        PTB( SetChrGA (ds_id, "spatialResolution",l2file->spatialResolution));

        // Write input parameters metadata
        ds_id.fid = l2file->grp_id[5];
        int32_t grp_id_input_parms;
        nc_def_grp( l2file->grp_id[5], "input_parameters", &grp_id_input_parms);
        ds_id.fid = grp_id_input_parms;

        char *tmp_parms;
        if ((tmp_parms = strdup(l2file->input_parms)) == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: Unable to copy input_parms string.\n",
                    __FILE__, __LINE__);
            exit(1);
        }

        char *token = strtok_r(tmp_parms, "\n", &end_str);
        while (token != NULL) {
            char *end_token;
            strcpy(tmp_str, token);
            char *name = strtok_r(token, "=", &end_token);
            for (i = 0; i < strlen(name); i++) {
                if (name[i] == ' ') {
                    name[i] = 0;
                    break;
                }
            }
            strcpy(tmp_str, strtok_r(NULL, ";", &end_token));
            trimBlanks(tmp_str);
            if (name[0] != '#')
                PTB(SetChrGA(ds_id, name, tmp_str));
            token = strtok_r(NULL, "\n", &end_str);
        }
        free(tmp_parms);


    } else {

        // write out HDF4 specific attributes
        PTB( SetChrGA (ds_id, "Sensor Name", sensorName[l2file->sensorID])  );
        PTB( SetChrGA (ds_id, "Mission", platformName[l2file->sensorID])  );
        PTB( SetI32GA (ds_id, "Number of Bands", numBands)                  );
        PTB( SetI32GA (ds_id, "Number of Scan Lines", numScans)             );
        PTB( SetI32GA (ds_id, "Pixels per Scan Line", numPixels)            );
        PTB( SetI32GA (ds_id, "Scene Center Scan Line", cscan)              );
        PTB( SetChrGA (ds_id, "Node Crossing Time", l2file->node_crossing_time)  );
        PTB( SetChrGA (ds_id, "Processing Time", ydhmsf(now(),'G'))         );
        PTB( SetChrGA (ds_id, "Input Parameters", l2file->input_parms)      );

    }

    /*                                                                  */
    /* Write out some global datasets                                   */
    /* ---------------------------------------------------------------- */
    /*                                                                  */
    if ( l2file->format == FMT_L2NCDF) ds_id.fid = l2file->grp_id[0];
    PTB( writeDS(ds_id, "wavelength",   Lambda, 0, 0, 0, numBands+numBandsIR, 0,0) );

    if (numBands > 0) {
        PTB( writeDS(ds_id, "vcal_gain",    Gain,   0, 0, 0, numBands, 0,0) );
        PTB( writeDS(ds_id, "vcal_offset",  Offset, 0, 0, 0, numBands, 0,0) );
        PTB( writeDS(ds_id, "F0",           Fobar,  0, 0, 0, numBands, 0,0) );
        PTB( writeDS(ds_id, "aw",           aw,     0, 0, 0, numBands, 0,0) );
        PTB( writeDS(ds_id, "bbw",          bbw,    0, 0, 0, numBands, 0,0) );
        PTB( writeDS(ds_id, "k_oz",         k_oz,   0, 0, 0, numBands, 0,0) );
        PTB( writeDS(ds_id, "k_no2",        k_no2,  0, 0, 0, numBands, 0,0) );
        PTB( writeDS(ds_id, "Tau_r",        Tau_r,  0, 0, 0, numBands, 0,0) );

        if ( l2file->format == FMT_L2NCDF) ds_id.fid = l2file->grp_id[4];
        PTB( writeDS(ds_id, "cntl_pt_cols", ictl,   0, 0, 0, nctl,     0,0) );
        PTB( writeDS(ds_id, "cntl_pt_rows", jctl,   0, 0, 0, numScans, 0,0) );
    }

    free(Gain);
    free(Offset);
    free(Fonom);
    free(Fobar);
    free(Tau_r);
    free(k_oz);
    free(k_no2);
    free(aw);
    free(bbw);
    free(calfile_str);

    return(LIFE_IS_GOOD);
}


/*----------------------------------------------------------------- */
/* Update scan-line datasets for the specified scan.                */
/* ---------------------------------------------------------------- */

int writel2( filehandle *l2file, int32_t recnum, l2str *l2rec, int outfile_number)
{
    static float32 *lon  = NULL;
    static float32 *lat  = NULL;
    static int32   *buf  = NULL;
    static float    fsol = -1.0;
    VOIDP   pbuf;
    int32   i,j,ip;
    l2prodstr *p;
    uint32_t mask;
    idDS ds_id;
    static int32_t sst_flag_cnt[NSSTFLAGS];
    static int32_t sst_triple_flag_cnt[NSSTFLAGS];
    static int32_t sst4_flag_cnt[NSSTFLAGS];
    static int32_t giop_flag_cnt[NSSTFLAGS];
    static int32_t qualsst_flag_cnt[NQSSTFLAGS];
    static int32_t qualsst_triple_flag_cnt[NQSSTFLAGS];
    static int32_t qualsst4_flag_cnt[NQSSTFLAGS];
    static int32_t qaa_flag_cnt[1]={0};
    static int32_t carder_flag_cnt[1]={0};
    static int32_t niwa_flag_cnt[1]={0};
    static const char *flag_lname[1] = {"PRODFAIL"};

    static uint8_t first=1;
    static float last_lat;
    static float geobox[4][100];
    static int32 geobox_cnt=0;
    float gring_fval[100];
    int32 gring_ival[100];

    ds_id.deflate = 0;
    ds_id.fid = l2file->sd_id;
    if(l2file->format == FMT_L2NCDF)
        ds_id.fftype = DS_NCDF;
    else
        ds_id.fftype = DS_HDF;

    if (recnum >= numScans) {
        fprintf(stderr,"-W- %s line %d: ", __FILE__,__LINE__);
        fprintf(stderr,"attempt to write rec %d of %d\n",recnum,numScans);
        return(1);
    } 

    if (recnum >= cscan && fsol < 0.0) {
        fsol = l2rec->fsol;
    }

    /* allocate buffer space */
    if (buf == NULL) {
        if ( (buf = calloc(numPixels,sizeof(int32))) == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: Unable to allocate buffer space.\n",
                    __FILE__,__LINE__);
            exit(1);
        }
    }
    if (lon == NULL) {
        if ( (lon = calloc(numPixels,sizeof(float32))) == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: Unable to allocate buffer space.\n",
                    __FILE__,__LINE__);
            exit(1);
        }
    }
    if (lat == NULL) {
        if ( (lat = calloc(numPixels,sizeof(float32))) == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: Unable to allocate buffer space.\n",
                    __FILE__,__LINE__);
            exit(1);
        }
    }


    /* Set arrays that depend on the number of control points */
    for (i = 0; i < nctl; i++) {
       lon [i] = l2rec->lon  [ictl[i]-1];
       lat [i] = l2rec->lat  [ictl[i]-1];
    }

    /* Write the scan-line data */
    if ( l2file->format == FMT_L2NCDF) ds_id.fid = l2file->grp_id[2];
    PTB( writeDS(ds_id, "year", l2rec->year, recnum, 0,0,1,1,1) );
    PTB( writeDS(ds_id, "day" , l2rec->day , recnum, 0,0,1,1,1) );
    PTB( writeDS(ds_id, "msec", l2rec->msec, recnum, 0,0,1,1,1) );
    PTB( writeDS(ds_id, "mside",&l2rec->mside,recnum, 0,0,1,1,1) );
    PTB( writeDS(ds_id, "detnum",&l2rec->detnum,recnum, 0,0,1,1,1) );
    PTB( writeDS(ds_id, "slon",&(l2rec->lon[spix]),recnum, 0,0,1,1,1) );
    PTB( writeDS(ds_id, "clon",&(l2rec->lon[cpix]),recnum, 0,0,1,1,1) );
    PTB( writeDS(ds_id, "elon",&(l2rec->lon[epix]),recnum, 0,0,1,1,1) );
    PTB( writeDS(ds_id, "slat",&(l2rec->lat[spix]),recnum, 0,0,1,1,1) );
    PTB( writeDS(ds_id, "clat",&(l2rec->lat[cpix]),recnum, 0,0,1,1,1) );
    PTB( writeDS(ds_id, "elat",&(l2rec->lat[epix]),recnum, 0,0,1,1,1) );
    PTB( writeDS(ds_id, "csol_z",&(l2rec->solz[cpix]),recnum, 0,0,1,1,1) );

    if ( l2file->format == FMT_L2NCDF) ds_id.fid = l2file->grp_id[4];
    PTB( writeDS(ds_id, "longitude", lon, recnum, 0,0,1,nctl,1) );
    PTB( writeDS(ds_id, "latitude" , lat, recnum, 0,0,1,nctl,1) );
    PTB( writeDS(ds_id, "tilt",&(l2rec->tilt), recnum, 0, 0, 1, 1, 1) );

    if (outfile_number == 0){
        if ((first == 1) || (fabs(last_lat - l2rec->lat[cpix]) > GEOBOX_INC)
                || (recnum == (numScans - 1))) {
            // make sure the points used in the gring are valid
            if (!(l2rec->flags[cpix] & NAVFAIL)) {

                if ((!(l2rec->flags[spix] & NAVFAIL)
                        && !(l2rec->flags[epix] & NAVFAIL))) {

                    first = 0;
                    geobox[0][geobox_cnt] = l2rec->lon[spix];
                    geobox[1][geobox_cnt] = l2rec->lat[spix];
                    geobox[2][geobox_cnt] = l2rec->lon[epix];
                    geobox[3][geobox_cnt] = l2rec->lat[epix];
                    last_lat = l2rec->lat[cpix];
                    geobox_cnt++;
                }
            }
        }
        // just in case the last line is buggered...
        if ((first == 0) && (recnum < (numScans - 1))
                && (l2rec->lat[cpix] != last_lat)) {   // && (geobox_cnt < 2)) {
            if (!(l2rec->flags[cpix] & NAVFAIL)) {

                if ((!(l2rec->flags[spix] & NAVFAIL)
                        && !(l2rec->flags[epix] & NAVFAIL))) {
                    geobox[0][geobox_cnt] = l2rec->lon[spix];
                    geobox[1][geobox_cnt] = l2rec->lat[spix];
                    geobox[2][geobox_cnt] = l2rec->lon[epix];
                    geobox[3][geobox_cnt] = l2rec->lat[epix];
                }
            }
        }
        if (recnum == (numScans - 1) && geobox_cnt == 1) geobox_cnt++;
    }
    /*                                                                  */
    /* Write the geophysical data                                       */
    /* ---------------------------------------------------------------- */
    /*                                                                  */
    if ( l2file->format == FMT_L2NCDF) ds_id.fid = l2file->grp_id[3];
    for (i=0; i<l2file->tot_prod; i++) {
        // Skip parameters already included if user requested them specifically
        if (!strcmp(l2file->l2_prod_names[i],"detnum") ||
            !strcmp(l2file->l2_prod_names[i],"mside")  ||
            !strcmp(l2file->l2_prod_names[i],"year")  ||
            !strcmp(l2file->l2_prod_names[i],"day")  ||
            !strcmp(l2file->l2_prod_names[i],"msec")  ||
            !strcmp(l2file->l2_prod_names[i],"slon")  ||
            !strcmp(l2file->l2_prod_names[i],"clon")  ||
            !strcmp(l2file->l2_prod_names[i],"elon")  ||
            !strcmp(l2file->l2_prod_names[i],"slat")  ||
            !strcmp(l2file->l2_prod_names[i],"clat")  ||
            !strcmp(l2file->l2_prod_names[i],"elat")  ||
            !strcmp(l2file->l2_prod_names[i],"csol_z") ||
            !strcmp(l2file->l2_prod_names[i],"latitude")  ||
            !strcmp(l2file->l2_prod_names[i],"longitude")  ||
            !strcmp(l2file->l2_prod_names[i],"tilt")
        )
            continue;

        // get the product index record
        p = l2file->prodptr + i; // get product structure from cache

        // extract the product and scale (if needed)
        if (p->slope == 1.0 && p->offset == 0.0) {
            pbuf = prodgen(p,l2rec);
        }
        else {
           pbuf = scale_sds( (float *) prodgen(p,l2rec), p, 1);
        }
        /* update flag counters when appropriate */
        if ( (strcmp(l2file->l2_prod_names[i],"flags_sst") == 0)) {
             update_flag_cnts16(sst_flag_cnt,pbuf, NSSTFLAGS, l2rec->npix, 1L);
         } else if   ( (strcmp(l2file->l2_prod_names[i],"flags_sst4") == 0) ) {
             update_flag_cnts16(sst4_flag_cnt,pbuf, NSSTFLAGS, l2rec->npix, 1L);
         } else if   ( (strcmp(l2file->l2_prod_names[i],"flags_sst_triple") == 0) ) {
             update_flag_cnts16(sst_triple_flag_cnt,pbuf, NSSTFLAGS, l2rec->npix, 1L);
         } else if   ( (strcmp(l2file->l2_prod_names[i],"flags_giop") == 0) ) {
             update_flag_cnts(giop_flag_cnt,pbuf, NGIOPFLAGS, l2rec->npix, 1L);
         } else if   ( (strcmp(l2file->l2_prod_names[i],"flags_qaa") == 0) ) {
             update_flag_cnts(qaa_flag_cnt,pbuf, 1, l2rec->npix, PRODFAIL);
         } else if   ( (strcmp(l2file->l2_prod_names[i],"flags_carder") == 0) ) {
             update_flag_cnts(carder_flag_cnt,pbuf, 1, l2rec->npix, PRODFAIL);
         } else if   ( (strcmp(l2file->l2_prod_names[i],"flags_niwa") == 0) ) {
             update_flag_cnts(niwa_flag_cnt,pbuf, 1, l2rec->npix, PRODFAIL);
         } else if   ( (strcmp(l2file->l2_prod_names[i],"qual_sst") == 0) ) {
             update_qual_cnts(qualsst_flag_cnt,pbuf, NQSSTFLAGS, l2rec->npix);
         } else if   ( (strcmp(l2file->l2_prod_names[i],"qual_sst_triple") == 0) ) {
             update_qual_cnts(qualsst_triple_flag_cnt,pbuf, NQSSTFLAGS, l2rec->npix);
         } else if   ( (strcmp(l2file->l2_prod_names[i],"qual_sst4") == 0) ) {
             update_qual_cnts(qualsst4_flag_cnt,pbuf, NQSSTFLAGS, l2rec->npix);
         }

         // write to L2 file
         if (p->rank == 2) {
           PTB( writeDS(ds_id,l2file->l2_prod_names[i],pbuf,
                    recnum,0,0,1,numPixels,1) );
        } else {
            PTB( writeDS(ds_id,l2file->l2_prod_names[i],pbuf,
                    recnum,0,0,1,1,1) );
        }
    }

    /* Update global flag counter */
    update_flag_cnts(l2file->flag_cnt,l2rec->flags, NFLAGS, l2rec->npix, 1L);

    if ( (l2file->format == FMT_L2NCDF) && !(recnum % 200)){
        nc_sync(l2file->sd_id);
    }

    /* Write global attributes */
    if (recnum == (numScans- 1)) {
        float flag_perc[NFLAGS];
        if ( l2file->format == FMT_L2NCDF) {

            // write out to netCDF file
            ds_id.fid = l2file->sd_id;

            scene_meta_write(ds_id);

            PTB( SetF64GA(ds_id, "earth_sun_distance_correction",fsol));

            // Write flag percentages metadata
            int32_t grp_id_flag_percentages;
            ds_id.sid = NC_GLOBAL;

            /* Report flag percentages */
            /* determine if there are any flag products */
            for (i=0; i<l2file->tot_prod; i++) {
                if ( (strcmp(l2file->l2_prod_names[i],"flags_sst") == 0)) {
                    nc_def_grp( l2file->grp_id[5], "sst_flag_percentages", &grp_id_flag_percentages);
                    ds_id.fid = grp_id_flag_percentages;
                    printf("\nSST: Percentage of pixels flagged:\n");
                    if (l2file->sensorID == AVHRR)
                        write_flag_pcnts(ds_id, fp_meta, sst_flag_cnt,NSSTFLAGS,avhrr_sst_flag_lname, numScans, numPixels);
                    else if (l2file->sensorID == VIIRS)
                        write_flag_pcnts(ds_id, fp_meta, sst_flag_cnt,NSSTFLAGS,viirs_sst_flag_lname, numScans, numPixels);
                    else
                        write_flag_pcnts(ds_id, fp_meta, sst_flag_cnt,NSSTFLAGS,sst_flag_lname, numScans, numPixels);
                }else if   ( (strcmp(l2file->l2_prod_names[i],"flags_sst4") == 0) ) {
                    nc_def_grp( l2file->grp_id[5], "sst4_flag_percentages", &grp_id_flag_percentages);
                    ds_id.fid = grp_id_flag_percentages;
                    printf("\nSST4: Percentage of pixels flagged:\n");
                        write_flag_pcnts(ds_id, fp_meta, sst4_flag_cnt,NSSTFLAGS,sst_flag_lname, numScans, numPixels);
                }else if   ( (strcmp(l2file->l2_prod_names[i],"flags_sst_triple") == 0) ) {
                    nc_def_grp( l2file->grp_id[5], "sst_triple_flag_percentages", &grp_id_flag_percentages);
                    ds_id.fid = grp_id_flag_percentages;
                    printf("\nsst_triple: Percentage of pixels flagged:\n");
		    write_flag_pcnts(ds_id, fp_meta, sst_triple_flag_cnt,NSSTFLAGS,viirs_sst_flag_lname, numScans, numPixels);
                }else if   ( (strcmp(l2file->l2_prod_names[i],"flags_giop") == 0) ) {
                    nc_def_grp( l2file->grp_id[5], "giop_flag_percentages", &grp_id_flag_percentages);
                    ds_id.fid = grp_id_flag_percentages;
                    printf("\nGIOP: Percentage of pixels flagged:\n");
                    write_flag_pcnts(ds_id, fp_meta, giop_flag_cnt,NGIOPFLAGS,giop_flag_lname, numScans, numPixels);
                }else if   ( (strcmp(l2file->l2_prod_names[i],"flags_qaa") == 0) ) {
                    nc_def_grp( l2file->grp_id[5], "qaa_flag_percentages", &grp_id_flag_percentages);
                    ds_id.fid = grp_id_flag_percentages;
                    printf("\nQAA: Percentage of pixels flagged:\n");
                    write_flag_pcnts(ds_id, fp_meta, qaa_flag_cnt,1,flag_lname, numScans, numPixels);
                }else if   ( (strcmp(l2file->l2_prod_names[i],"flags_carder") == 0) ) {
                    nc_def_grp( l2file->grp_id[5], "carder_flag_percentages", &grp_id_flag_percentages);
                    ds_id.fid = grp_id_flag_percentages;
                    printf("\nCARDER: Percentage of pixels flagged:\n");
                    write_flag_pcnts(ds_id, fp_meta, carder_flag_cnt,1,flag_lname, numScans, numPixels);
                }else if   ( (strcmp(l2file->l2_prod_names[i],"flags_niwa") == 0) ) {
                    nc_def_grp( l2file->grp_id[5], "niwa_flag_percentages", &grp_id_flag_percentages);
                    ds_id.fid = grp_id_flag_percentages;
                    printf("\nNIWA: Percentage of pixels flagged:\n");
                    write_flag_pcnts(ds_id, fp_meta, niwa_flag_cnt,1,flag_lname, numScans, numPixels);
                }else if   ( (strcmp(l2file->l2_prod_names[i],"qual_sst") == 0) ) {
                    nc_def_grp( l2file->grp_id[5], "qual_sst_percentages", &grp_id_flag_percentages);
                    ds_id.fid = grp_id_flag_percentages;
                    printf("\nQUAL_SST: Percentage of pixels flagged:\n");
                    write_qual_flag_pcnts(ds_id, fp_meta, qualsst_flag_cnt,NQSSTFLAGS,qual_sst_flag_lname);
                }else if   ( (strcmp(l2file->l2_prod_names[i],"qual_sst4") == 0) ) {
                    nc_def_grp( l2file->grp_id[5], "qual_sst4_percentages", &grp_id_flag_percentages);
                    ds_id.fid = grp_id_flag_percentages;
                    printf("\nQUAL_SST4: Percentage of pixels flagged:\n");
                    write_qual_flag_pcnts(ds_id, fp_meta, qualsst4_flag_cnt,NQSSTFLAGS,qual_sst_flag_lname);
                }else if   ( (strcmp(l2file->l2_prod_names[i],"qual_sst_triple") == 0) ) {
                    nc_def_grp( l2file->grp_id[5], "qual_sst_triple_percentages", &grp_id_flag_percentages);
                    ds_id.fid = grp_id_flag_percentages;
                    printf("\nQUAL_sst_triple: Percentage of pixels flagged:\n");
                    write_qual_flag_pcnts(ds_id, fp_meta, qualsst_triple_flag_cnt,NQSSTFLAGS,qual_sst_flag_lname);
                }

            }
            nc_def_grp( l2file->grp_id[5], "flag_percentages", &grp_id_flag_percentages);
            ds_id.fid = grp_id_flag_percentages;
            printf("\nPercentage of pixels flagged:\n");
            write_flag_pcnts(ds_id, fp_meta, l2file->flag_cnt,NFLAGS,l2_flag_lname, numScans, numPixels);


            // Geobox attributes
            if ( l2file->format == FMT_L2NCDF) ds_id.fid = l2file->sd_id;
            j = 1;
            gring_fval[0] = geobox[0][0];
            for (i=0; i<geobox_cnt; i++) {
                gring_fval[j++] = geobox[2][i];
            }
            for (i=0; i<geobox_cnt-1; i++) {
                gring_fval[j++] = geobox[0][geobox_cnt-1-i];
            }
            if ( l2file->format == FMT_L2NCDF) ds_id.fid = l2file->grp_id[4];
            PTB(setAttr(ds_id,"gringpointlongitude",NC_FLOAT,j,(VOIDP)gring_fval));

            j = 1;
            gring_fval[0] = geobox[1][0];
            gring_ival[0] = j;
            for (i=0; i<geobox_cnt; i++) {
                gring_ival[j] = j+1;
                gring_fval[j++] = geobox[3][i];
            }
            for (i=0; i<geobox_cnt-1; i++) {
                gring_ival[j] = j+1;
                gring_fval[j++] = geobox[1][geobox_cnt-1-i];
            }
            if ( l2file->format == FMT_L2NCDF) ds_id.fid = l2file->grp_id[4];
            PTB(setAttr(ds_id,"gringpointlatitude",NC_FLOAT,j,(VOIDP)gring_fval));
            PTB(setAttr(ds_id,"gringpointsequence",NC_INT,j,(VOIDP)gring_ival));

        } else {

            // write out to HDF4
            scene_meta_write(ds_id);

            PTB( SetF64GA(ds_id, "Earth-Sun Distance Correction",fsol));

            /* Report flag percentages */
            printf("\nPercentage of pixels flagged:\n");
            for (i=0; i<NFLAGS; i++) {
                flag_perc[i] = ((float) l2file->flag_cnt[i])/numScans/numPixels*100.0;
                printf("Flag #%2d: %16s %10d %8.4f\n",
                        i+1,l2_flag_lname[i],l2file->flag_cnt[i],flag_perc[i]);
                if (fp_meta != NULL)
                    fprintf(fp_meta,"Flag #%2d: %16s %10d %8.4f\n",
                            i+1,l2_flag_lname[i],l2file->flag_cnt[i],flag_perc[i]);
            }
            PTB(sd_setattr(ds_id.fid, "Flag Percentages", DFNT_FLOAT32,
                    NFLAGS,(VOIDP)flag_perc));

        }

    }   

    return(LIFE_IS_GOOD);
}

/* -------------------------------------------------------- */
/* Finish access for the current file.                      */
/* -------------------------------------------------------- */
int closel2(filehandle *l2file)
{
    idDS ds_id;
    ds_id.deflate = 0;
    ds_id.fid = l2file->sd_id;
    if(l2file->format == FMT_L2NCDF)
        ds_id.fftype = DS_NCDF;
    else
        ds_id.fftype = DS_HDF;

    if ( l2file->format == FMT_L2HDF) {
        PTB( MakeVgroups(l2file) );
    }

    if (endDS(ds_id)) {
        fprintf(stderr,"-E- %s line %d: endDS(%d) failed for file, %s.\n",
                __FILE__,__LINE__,ds_id.fid,l2file->name);
        return(HDF_FUNCTION_ERROR);
    }

    if (fp_meta != NULL)
        fclose(fp_meta);

    if(l2file->prodptr)
        free(l2file->prodptr);
    l2file->prodptr = NULL;

    return(LIFE_IS_GOOD);
}

void  update_flag_cnts(int32_t *flag_cnt, int32_t*flags, int32_t nflags, int32_t npix, uint32_t init_mask) {
    int32_t i, ip;
    uint32_t mask;

    /* Update flag counter */
    for (ip=0; ip<npix; ip++) {
        mask = init_mask;
        for (i=0; i<nflags; i++) {
            flag_cnt[i] += ((flags[ip] & mask) > 0);
            mask *= 2L;
        }
    }
    return;
}

void  update_flag_cnts16(int32_t *flag_cnt, int16_t *flags, int32_t nflags, int32_t npix, uint32_t init_mask) {
    int32_t i, ip;
    uint32_t mask;

    /* Update flag counter */
    for (ip=0; ip<npix; ip++) {
        mask = init_mask;
        for (i=0; i<nflags; i++) {
            flag_cnt[i] += ((flags[ip] & mask) > 0);
            mask *= 2L;
        }
    }
    return;
}

void  update_qual_cnts(int32_t *flag_cnt, int8_t* flags, int32_t nflags, int32_t npix) {
    int32_t i, ip;
    uint32_t mask;

    /* Update flag counter */
    for (ip=0; ip<npix; ip++) {
        mask = 0;
        for (i=0; i<nflags; i++) {
            flag_cnt[i] += ((flags[ip] == mask));
            mask ++;
        }
    }
    return;
}


int write_flag_pcnts(idDS ds_id, FILE *fpmeta, int32_t *flag_cnt,int32_t nflags, const char *flag_lname[], int32_t numScans, int32_t numPixels)
{
    int32_t i, status;
    float *flag_perc;

    flag_perc = (float *)calloc(nflags,sizeof(float));

    /* Report flag percentages */
    for (i=0; i<nflags; i++) {
        flag_perc[i] = ((float) flag_cnt[i])/numScans/numPixels*100.0;
        printf("Flag #%2d: %16s %10d %8.4f\n",i+1,flag_lname[i],flag_cnt[i],flag_perc[i]);
        if (fp_meta != NULL)
            fprintf(fpmeta,"Flag #%2d: %16s %10d %8.4f\n",i+1,flag_lname[i],flag_cnt[i],flag_perc[i]);

        PTB(setAttr(ds_id, flag_lname[i], NC_FLOAT,1,(VOIDP) (flag_perc+i)));
    }

    return 0;
}

int write_qual_flag_pcnts(idDS ds_id, FILE *fpmeta, int32_t *flag_cnt,int32_t nflags, const char *flag_lname[])
{
    int32_t i, status,sumflags=0;
    float *flag_perc;

    flag_perc = (float *)calloc(nflags,sizeof(float));

    for (i=0; i<nflags; i++) sumflags+=flag_cnt[i];

    /* Report flag percentages */
    for (i=0; i<nflags; i++) {
        flag_perc[i] = ((float) flag_cnt[i])/sumflags*100.0;
        printf("Flag #%2d: %16s %10d %8.4f\n",i+1,flag_lname[i],flag_cnt[i],flag_perc[i]);
        if (fp_meta != NULL)
            fprintf(fpmeta,"Flag #%2d: %16s %10d %8.4f\n",i+1,flag_lname[i],flag_cnt[i],flag_perc[i]);

        PTB(setAttr(ds_id, flag_lname[i], NC_FLOAT,1,(VOIDP) (flag_perc+i)));
    }

    return 0;
}
