/**
 *  @file l1_meris_N1.h
 *  @brief MERIS reader
 *  @author Paul Lyon
 *  @author Naval Research Laboratory, Stennis Space Center, MS
 */

/* ============================================================================ */
/* module l1_meris_N1.c - functions to read MERIS   Reduced RES L2 for MSL12    */
/* NOTE!! THIS IS A TEMPORARY READER FOR L2 FILES, INTO L1 SEADAS RECORD        */
/* Written By:  Paul E. Lyon NRL, Oct. 2006.                                    */
/*                                                                              */
/* ============================================================================ */

#include <stdbool.h>
#include "l1_meris_N1.h"
#include "l12_proto.h"
#include "epr_api.h"
#include "epr_field.h"
#include <math.h>

#define MERIS_NBANDS 15

#define MERIS_BANDINFO_FILENAME      "band_info.txt"
#define MERIS_WAVELENGTH_FR_FILENAME "central_wavelen_fr.txt"
#define MERIS_WAVELENGTH_RR_FILENAME "central_wavelen_rr.txt"
#define MERIS_SUN_FLUX_FR_FILENAME   "sun_spectral_flux_fr.txt"
#define MERIS_SUN_FLUX_RR_FILENAME   "sun_spectral_flux_rr.txt"

// max # of pixels in a reduced resolution scan line
#define MERIS_REDUCED_RESOLUTION_MAX 1121
// # of detectors in a full resolution scan line
#define MERIS_FR_DETECTORS 3700
// # of detectors in a full resolution scan line
#define MERIS_RR_DETECTORS 925

// max line length in a config table file
#define MERIS_LINE_MAX 255

/* ------------------------------------------------------------------------ */
/* These L1 flag values were derrived from the MERIS specification          */
/*                                                                          */
/*  http://earth.esa.int/pub/ESA_DOC/ENVISAT/Vol11_Meris_5b.pdf             */
/*  page 64,  Table 11.4.1.7.4.2-3 MER_RR__1P - Quality Flag Coding         */
/*                                                                          */
/* dshea                                                                    */
/* ------------------------------------------------------------------------ */
#define MERIS_L1FLAG_COSMETIC   0x01
#define MERIS_L1FLAG_DUPLICATED 0x02
#define MERIS_L1FLAG_GLINT      0x04
#define MERIS_L1FLAG_SUSPECT    0x08
#define MERIS_L1FLAG_LAND       0x10
#define MERIS_L1FLAG_BRIGHT     0x20
#define MERIS_L1FLAG_COASTLINE  0x40
#define MERIS_L1FLAG_INVALID    0x80

/* ------------------------------------------------------------------------ */
/* These L2 flag values were derrived from the MERIS specification          */
/*                                                                          */
/*  http://earth.esa.int/pub/ESA_DOC/ENVISAT/Vol11_Meris_5b.pdf             */
/*  page 105,  Table 11.5.1.7.4.8-2 Description of the Flag Coding          */
/*                                                                          */
/* dshea                                                                    */
/* ------------------------------------------------------------------------ */
#define MERIS_L2FLAG_WHITE_SCATTER   0x000001
#define MERIS_L2FLAG_PRESSURE_CONF   0x000002
#define MERIS_L2FLAG_HIGH_GLINT      0x000004
#define MERIS_L2FLAG_DDV             0x000008
#define MERIS_L2FLAG_MEDIUM_GLINT    0x000010
#define MERIS_L2FLAG_ICE_HAZE        0x000020
#define MERIS_L2FLAG_CASE2_Y         0x000040
#define MERIS_L2FLAG_CASE2_ANOM      0x000080
#define MERIS_L2FLAG_CASE2_S         0x000100
#define MERIS_L2FLAG_ABSOA_DUST      0x000200
#define MERIS_L2FLAG_OOADB           0x000400
#define MERIS_L2FLAG_SUSPECT         0x000800
#define MERIS_L2FLAG_COSMETIC        0x001000
#define MERIS_L2FLAG_COASTLINE       0x002000
#define MERIS_L2FLAG_PCD_19          0x004000
#define MERIS_L2FLAG_PCD_18          0x008000
#define MERIS_L2FLAG_PCD_17          0x010000
#define MERIS_L2FLAG_PCD_16          0x020000
#define MERIS_L2FLAG_PCD_15          0x040000
#define MERIS_L2FLAG_PCD_14          0x080000
#define MERIS_L2FLAG_PCD_1_13        0x100000
#define MERIS_L2FLAG_WATER           0x200000
#define MERIS_L2FLAG_CLOUD           0x400000
#define MERIS_L2FLAG_LAND            0x800000

static EPR_SProductId *fileID = NULL;
static int spix = 0;
static int file_npix;
static int         year, day, msec;
static double      time_interval;
static bool         isLevel1;

/* ------------------------------------------------------------------------ */
/*                                                                          */
/* calc_smile_delta()                                                       */
/*   calculates the smile delta for all of the bands of a pixel             */
/*                                                                          */
/* int *shouldCorrect  true if the reflectance correction should be made    */
/* int *indexes1       index of lower band to use for interpolation         */
/* int *indexes2       upper band to use for intrepolation                  */
/* float *radiances    original measurment                                  */
/* float *theoretWLs   theoretical wavelengths for each band                */
/* float *theoretE0s   theoretical sun spectral fluxes for each band        */
/* float *detectorWLs  detector wavelengths for each band                   */
/* float *detectorE0s  sun spectral fluxes for the detector wavelengths     */
/* float *delta        resulting smile delta correction                     */
/*                                                                          */
/* D. Shea, SAIC, Jan 2009.                                                 */
/* ------------------------------------------------------------------------ */
void calc_smile_delta(int   *shouldCorrect, 
        int   *indexes1,
        int   *indexes2,
        float *radiances,
        float *theoretWLs,
        float *theoretE0s,
        float *detectorWLs,
        float *detectorE0s,
        float *delta)
{
    double r0, r1, r2, rc, dl, dr;
    int i0, i1, i2;

    for(i0=0; i0<MERIS_NBANDS; i0++) {

        // perform irradiance correction
        r0 = radiances[i0] / detectorE0s[i0];
        rc = r0 * theoretE0s[i0];
        delta[i0] = rc - radiances[i0];
        if (shouldCorrect[i0]) {
            // perform reflectance correction
            i1 = indexes1[i0];
            i2 = indexes2[i0];
            r1 = radiances[i1] / detectorE0s[i1];
            r2 = radiances[i2] / detectorE0s[i2];
            dl = (theoretWLs[i0] - detectorWLs[i0]) / (detectorWLs[i2] - detectorWLs[i1]);
            dr = (r2 - r1) * dl * theoretE0s[i0];
            delta[i0] += dr;
        }
    } // for bands

}


/* ------------------------------------------------------------------------ */
/* radcor()                                                                 */
/*     loads smile calibration deltas into the radcor array of l1rec        */
/*                                                                          */
/* l1rec  level 1 record to set radcor in                                   */
/* ip     pixel number used to calculate the smile correction               */
/* land   0=ocean, 1=land                                                   */
/*                                                                          */
/*                                                                          */
/* D. Shea, SAIC, Jan 2009.                                                 */
/* ------------------------------------------------------------------------ */
void radcor(l1str *l1rec, int32_t ip, int land)
{
    static int   firstCall = 1;
    static int   switch_land[MERIS_NBANDS];
    static int   lower_land[MERIS_NBANDS];
    static int   upper_land[MERIS_NBANDS];
    static int   switch_water[MERIS_NBANDS];
    static int   lower_water[MERIS_NBANDS];
    static int   upper_water[MERIS_NBANDS];
    static float theoWL[MERIS_NBANDS];
    static float theoE0[MERIS_NBANDS];
    static float *detectorWL;
    static float *detectorE0;
    static int   fullResolution = 0;
    static int   numDetectors;

    int ib;
    int ipb;
    int *shouldCorrect;
    int correctionPossible;
    int *index1;
    int *index2;
    int detectorIndex;
    float *Ltemp;

    // bail if it is not the MERIS sensor
    if(l1rec->sensorID != MERIS)
        return;

    // bail if the rad_opt is set to 0
    // if set to -1 (default) then do the correction for MERIS
    if(l1rec->input->rad_opt == 0)
        return;
    if ((Ltemp = (float *) calloc(l1rec->nbands,sizeof(float ))) == NULL) {
        printf("-E- %s line %d : error allocating memory for l1_generic_write:open1_write.\n",
                __FILE__,__LINE__);
        exit(1);
    }

    if (firstCall) {
        char *tmp_str;
        char  line   [MERIS_LINE_MAX];
        char  path   [FILENAME_MAX] = "";
        char  file   [FILENAME_MAX] = "";
        FILE *fp = NULL;
        int count;
        int dummy;
        int band;
        int detector;
        float *floatp;

        if ((tmp_str = getenv("OCDATAROOT")) == NULL) {
            printf("OCDATAROOT environment variable is not defined.\n");
            exit(1);
        }

        strcpy(path,tmp_str); strcat(path,"/"); strcat(path,sensorDir[l1rec->sensorID]); strcat(path,"/cal/");

        /*-------------------------------------------------------------*/
        // read the band_info file
        strcpy(file, path); strcat(file, MERIS_BANDINFO_FILENAME);
        if ( (fp = fopen(file,"r")) == NULL ) {
            fprintf(stderr,
                    "-E- %s line %d: unable to open %s for reading\n",__FILE__,__LINE__,file);
            exit(1);
        }

        // disgard the first line of column labels
        fgets( line, MERIS_LINE_MAX, fp );
        for(band=0; band<MERIS_NBANDS; band++) {
            if(!fgets( line, MERIS_LINE_MAX, fp ) ) {
                fprintf(stderr,
                        "-E- %s line %d: unable to read band %d from file %s\n",__FILE__,__LINE__,band,file);
                exit(1);
            }
            count = sscanf(line,"%d %d %d %d %d %d %d %f %f", &dummy, &switch_land[band], &lower_land[band],
                    &upper_land[band], &switch_water[band], &lower_water[band], &upper_water[band],
                    &theoWL[band], &theoE0[band]);
            if(count != 9) {
                fprintf(stderr,
                        "-E- %s line %d: unable to read band %d line from file %s, count = %d\n",
                        __FILE__,__LINE__,band,file,count);
                exit(1);
            }

            // now change the indexes to 0 based array indexes
            lower_land[band]--;
            upper_land[band]--;
            lower_water[band]--;
            upper_water[band]--;

        } // for band
        fclose(fp);


        /*-------------------------------------------------------------*/
        // determine if we need the full resolution or reduced resolution tables
        if(file_npix > MERIS_REDUCED_RESOLUTION_MAX) {
            fullResolution = 1;
            detectorWL = (float*) malloc(sizeof(float)*MERIS_NBANDS*MERIS_FR_DETECTORS);
            detectorE0 = (float*) malloc(sizeof(float)*MERIS_NBANDS*MERIS_FR_DETECTORS);
        } else {
            fullResolution = 0;
            detectorWL = (float*) malloc(sizeof(float)*MERIS_NBANDS*MERIS_RR_DETECTORS);
            detectorE0 = (float*) malloc(sizeof(float)*MERIS_NBANDS*MERIS_RR_DETECTORS);
        }

        /*-------------------------------------------------------------*/
        // read the central wavelength table
        if(fullResolution) {
            strcpy(file, path); strcat(file, MERIS_WAVELENGTH_FR_FILENAME);
            numDetectors = MERIS_FR_DETECTORS;
        } else {
            strcpy(file, path); strcat(file, MERIS_WAVELENGTH_RR_FILENAME);
            numDetectors = MERIS_RR_DETECTORS;
        }

        if ( (fp = fopen(file,"r")) == NULL ) {
            fprintf(stderr,
                    "-E- %s line %d: unable to open %s for reading\n",__FILE__,__LINE__,file);
            exit(1);
        }

        // disgard the first line of column labels
        fgets( line, MERIS_LINE_MAX, fp );

        floatp = detectorWL;
        for(detector=0; detector<numDetectors; detector++) {
            if(!fgets( line, MERIS_LINE_MAX, fp ) ) {
                fprintf(stderr,
                        "-E- %s line %d: unable to read detector %d from file %s\n",__FILE__,__LINE__,detector,file);
                exit(1);
            }

            /* I really should read all of the bands defined by MERIS_NBANDS, but I am just reading
	 15 bands for now.  dshea */
            count = sscanf(line,"%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f", &dummy, floatp, floatp+1,
                    floatp+2, floatp+3, floatp+4, floatp+5, floatp+6, floatp+7, floatp+8, floatp+9,
                    floatp+10, floatp+11, floatp+12, floatp+13, floatp+14);
            if(count != 16) {
                fprintf(stderr,
                        "-E- %s line %d: unable to read whole detector %d line from file %s, count = %d\n",
                        __FILE__,__LINE__,detector,file,count);
                exit(1);
            }

            floatp += MERIS_NBANDS;
        } // for detector
        fclose(fp);


        /*-------------------------------------------------------------*/
        // read the sun spectral flux table
        if(fullResolution) {
            strcpy(file, path); strcat(file, MERIS_SUN_FLUX_FR_FILENAME);
            numDetectors = MERIS_FR_DETECTORS;
        } else {
            strcpy(file, path); strcat(file, MERIS_SUN_FLUX_RR_FILENAME);
            numDetectors = MERIS_RR_DETECTORS;
        }

        if ( (fp = fopen(file,"r")) == NULL ) {
            fprintf(stderr,
                    "-E- %s line %d: unable to open %s for reading\n",__FILE__,__LINE__,file);
            exit(1);
        }

        // disgard the first line of column labels
        fgets( line, MERIS_LINE_MAX, fp );

        floatp = detectorE0;
        for(detector=0; detector<numDetectors; detector++) {
            if(!fgets( line, MERIS_LINE_MAX, fp ) ) {
                fprintf(stderr,
                        "-E- %s line %d: unable to read detector %d from file %s\n",__FILE__,__LINE__,detector,file);
                exit(1);
            }

            /* I really should read all of the bands defined by MERIS_NBANDS, but I am just reading
	 15 bands for now.  dshea */
            count = sscanf(line,"%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f", &dummy, floatp, floatp+1,
                    floatp+2, floatp+3, floatp+4, floatp+5, floatp+6, floatp+7, floatp+8, floatp+9,
                    floatp+10, floatp+11, floatp+12, floatp+13, floatp+14);
            if(count != 16) {
                fprintf(stderr,
                        "-E- %s line %d: unable to read whole detector %d line from file %s, count = %d\n",
                        __FILE__,__LINE__,detector,file,count);
                exit(1);
            }

            floatp += MERIS_NBANDS;
        } // for detector
        fclose(fp);

        if(want_verbose)
            printf("\nSmile corrections enabled.\n\n");

        firstCall = 0;

    } // if firstCall


    /*-------------------------------------------------------------*/
    // do the correction

    detectorIndex = l1rec->pixdet[ip];
    correctionPossible = ((detectorIndex >= 0) && (detectorIndex < numDetectors));

    if(correctionPossible) {
        if(land) {
            shouldCorrect = switch_land;
            index1 = lower_land;
            index2 = upper_land;
        } else {
            shouldCorrect = switch_water;
            index1 = lower_water;
            index2 = upper_water;
        }

        // correct all bands for ozone
        for (ib=0; ib<l1rec->nbands; ib++) {
            ipb = ip*l1rec->nbands + ib;

            /* Correct for ozone absorption.  We correct for inbound and outbound here, */
            /* then we put the inbound back when computing Lw.                          */
            Ltemp[ib] = l1rec->Lt[ipb]/l1rec->tg_sol[ipb]/l1rec->tg_sen[ipb];
        } // for ib

        calc_smile_delta(shouldCorrect,
                index1,
                index2,
                Ltemp,
                theoWL,
                theoE0,
                &(detectorWL[detectorIndex*MERIS_NBANDS]),
                &(detectorE0[detectorIndex*MERIS_NBANDS]),
                &(l1rec->radcor[ip*l1rec->nbands]));

    } // correction possible

    free(Ltemp);

}



/**
 * @brief opens a MERIS file for reading to load into L1 record
 * @param[in]   file   file handle to MERIS file
 *
 * @author Paul E. Lyon NRL, Oct. 2006.
 */

int
openl1_meris_N1(filehandle * file)
{
    const char         *fltime;
    const char         *fname;
    char                monthstr[10];
    static char         months_list[12][4] =
    { "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL",
            "AUG", "SEP", "OCT", "NOV", "DEC" };
    unsigned int        source_w, source_h;
    int                 minute, hour, month, sec, i;
    float               fsec, milisec;
    const EPR_SRecord  *record;
    const EPR_SField   *field;
    EPR_SBandId        *band_id = NULL;
    char               *sph_names[] =
    { "NUM_BANDS", "LINE_TIME_INTERVAL", "NUM_SLICES" };
    int                 nf = 3;

    // Initialize the API
    epr_init_api(e_log_debug, NULL, NULL);

    // Open the N1 input file 

    fileID = epr_open_product(file->name);
    if (fileID == NULL) {
        fprintf(stderr, "-E- %s line %d: epr_open_product(%s) failed.\n",
                __FILE__, __LINE__, file->name);
        return (1);
    }

    band_id  = epr_get_band_id(fileID, "l1_flags");
    isLevel1 = (band_id != NULL) ? true : false;

    // Get pixel and scan dimensions 

    source_h = fileID->scene_height;
    source_w = fileID->scene_width;
    file_npix = source_w;
    if(want_verbose) {
        if (isLevel1)
            printf("MERIS Level1B Npix  :%u Nscans:%u\n", source_w, source_h);
        else
            printf("MERIS Level2  Npix  :%u Nscans:%u\n", source_w, source_h);
    } // want_verbose

    // get specific product header (SPH)

    record = epr_get_sph(fileID);
    if (record == NULL) {
        fprintf(stderr, "-E- %s line %d: epr_get_sph(fileID) failed.\n",
                __FILE__, __LINE__);
        return (1);
    }


    // get year, day, msecs 

    field = epr_get_field(record, "FIRST_LINE_TIME");
    if (field == NULL) {
        fprintf(stderr,
                "-E- %s line %d: epr_get_field(record,FIRST_LINE_TIME) failed.\n",
                __FILE__, __LINE__);
        return (1);
    }

    fltime = epr_get_field_elem_as_str(field);
    sscanf(fltime, "%02d-%3s-%04d %02d:%02d:%f", &day, monthstr, &year,
            &hour, &minute, &fsec);
    monthstr[4] = '\0';
    month = 1;
    for (i = 0; i < 12; i++)
        if (strncmp(monthstr, months_list[i], 3) == 0)
            month = i + 1;
    sec = (int) trunc((double) fsec);
    milisec = (float) (fsec - sec) * 1e3;
    ymdhms2ydmsec(year, month, day, hour, minute, sec, (int32_t *) &year, (int32_t *) &day, (int32_t *) &msec);
    msec += (int) milisec;

    field = epr_get_field(record, "LINE_TIME_INTERVAL");
    if (field == NULL) {
        fprintf(stderr,
                "-E- %s line %d: epr_get_field(record,FIRST_LINE_TIME) failed.\n",
                __FILE__, __LINE__);
        return (1);
    }
    time_interval = (double) epr_get_field_elem_as_uint(field, 0);        // microseconds;
    if(want_verbose)
        printf("MERIS time interval %lf\n", time_interval);

    // dump a few more items for infomational purposes

    for (i = 0; i < nf; i++) {
        field = epr_get_field(record, sph_names[i]);
        if (field == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: epr_get_field(record,%s) failed.\n",
                    __FILE__, __LINE__, sph_names[i]);
            return (1);
        }
        if(want_verbose)
            printf("\t%s = %d\n", sph_names[i],
                    (int)epr_get_field_elem_as_uint(field, 0));
    }

    // get the START_PIXEL and END_PIXEL that were injected into the MPH
    // by the extractor.  START_PIXEL and END_PIXEL are 1 based.  The epr
    // api reverses the pixels on a line, so spix = source_w - epix.
    record = epr_get_mph(fileID);
    if (record) {
        field = epr_get_field(record,"PRODUCT");
        if (field) {
            fname = epr_get_field_elem_as_str(field);
            if (strncmp(fname,"MER_RR",6) == 0)
                strcpy(file->spatialResolution,"1.2 km");
            else
                strcpy(file->spatialResolution,"300 m");
        }

        field = epr_get_field(record, "END_PIXEL");
        if (field) {
            int tmp_epix = epr_get_field_elem_as_uint(field, 0);
            field = epr_get_field(record, "START_PIXEL");
            if (field) {
                int tmp_spix = epr_get_field_elem_as_uint(field, 0);
                if((tmp_spix < source_w) && (tmp_spix < tmp_epix)) {
                    spix = source_w - tmp_epix;
                    source_w = tmp_epix - tmp_spix + 1; // spix and epix are inclusive
                    if(want_verbose)
                        printf("OBPG Extract - spix=%d, npix=%d\n",
                                spix+1, source_w);
                } // spix, epix reasonable
            } // spix good
        } // epix good
    } // found the MPH

    // define number of input products

    if (isLevel1) {
        file->nbands = 15;
    } else {
        file->n_inprods = 7;
        file->nbands = 15;
    }

    file->npix =  (int) source_w;
    file->nscan = (int) source_h;


    return (LIFE_IS_GOOD);
}

/**
 * @brief reads 1 scan line from MERIS file, loads l1rec
 * @param[in]   file   file handle to MERIS file
 * @param[in]   scan   scan number to read
 * @apram[out]  l1rec  output l1rec
 *
 * @author Paul E. Lyon NRL, Oct. 2006.
 * W. Robinson, SAIC, 22 May 2012  account for msec going to next day
 */

int readl1_meris_N1(filehandle *file, int32 scan, l1str *l1rec)
{
    static int          firstCall = 1;
    int                 err_code;

    // only look at SeaWiFS bands until we have rayleigh/aerosol tables
    // NOTE: the following names are part of the API and NOT the file.  See libmeris.
    static char        *l2_names[] = { 
            "reflec_1",  "reflec_2",  "reflec_3",  "reflec_4",    // 412.7, 442.6, 489.9, 509.8
            "reflec_5",  "reflec_6",  "reflec_7",  "reflec_8",    // 559.7, 619.6, 664.5, 680.8
            "reflec_9",  "reflec_10", "reflec_11", "reflec_12",   // 708.3, 753.4, 761.5, 778.4
            "reflec_13", "reflec_14", "reflec_15"                 // 864.9, 884.9, 900.0
    };
    static char        *l1_names[] = { 
            "radiance_1",  "radiance_2",  "radiance_3",  "radiance_4",  // 412.7, 442.6, 489.9, 509.8
            "radiance_5",  "radiance_6",  "radiance_7",  "radiance_8",  // 559.7, 619.6, 664.5, 680.8
            "radiance_9",  "radiance_10", "radiance_11", "radiance_12", // 708.3, 753.4, 761.5, 778.4
            "radiance_13", "radiance_14", "radiance_15"                 // 864.9, 884.9, 900.0
    };
    static char       **names;
    int32               npix;
    int32               nbands;
    int32               ip, ib, ipb;
    instr              *input;
    EPR_SBandId        *band_id = NULL;
    EPR_SRaster        *temp_raster, *bm_raster;
    epr_uint            flag;
    static char        *invalid_flag;
    float64 recsec, sec70;
    int16 sh_year, sh_day;

    input = file->input;

    if (firstCall) {
        // signal that we don't do any atmospheric correction using atmocor flag
        if (!isLevel1) {
            input->atmocor = 0;
            names = l2_names;
        } else {
            //            input->atmocor = 1;
            names = l1_names;
            invalid_flag = malloc(sizeof(char) * file->npix);
        }
        if(want_verbose)
            printf("file->nbands = %d, l1rec->nbands = %d\n", (int)file->nbands,
                    (int)l1rec->nbands );
        firstCall = 0;
    }

    nbands = file->nbands;
    npix = file->npix;

    /*  set time for this scan and account for overflow of msec to next day  */
    recsec = (float64)( msec + time_interval * scan / 1e3 ) / 1.e3;
    sec70 = ymds2unix( (int16)year, 1, (int16)day, recsec );
    unix2yds( sec70, &sh_year, &sh_day, &recsec );
    *(l1rec->year) = (int32) sh_year;
    *(l1rec->day)  = (int32) sh_day;
    *(l1rec->msec) = (int32) ( recsec * 1e3 );

    // read L1 flags and set invalid_flag if suspect or invalid set

    if(isLevel1) {
        band_id = epr_get_band_id(fileID, "l1_flags");
        if (band_id == NULL) {
            printf("-E- %s line %d: Error finding band_id:l1_flags\n", __FILE__,
                    __LINE__);
            return (1);
        }
        temp_raster = epr_create_compatible_raster(band_id, file_npix, 1, 1, 1);
        if (temp_raster == NULL) {
            printf("-E- %s line %d: Error allocating raster space.\n", __FILE__,
                    __LINE__);
            return (1);
        }
        err_code = epr_read_band_raster(band_id, 0, scan, temp_raster);
        for (ip = 0; ip < npix; ip++) {
            flag = epr_get_pixel_as_uint(temp_raster, spix+ip, 0);
            if((flag & MERIS_L1FLAG_SUSPECT) || (flag & MERIS_L1FLAG_INVALID)) {
                invalid_flag[ip] = 1;
            } else {
                invalid_flag[ip] = 0;
            }
        }

        epr_free_raster(temp_raster);
    } 


    // read latitude

    band_id = epr_get_band_id(fileID, "latitude");
    if (band_id == NULL) {
        printf("-E- %s line %d: Error finding band_id:latitude\n", __FILE__,
                __LINE__);
        return (1);
    }
    temp_raster = epr_create_compatible_raster(band_id, file_npix, 1, 1, 1);
    if (temp_raster == NULL) {
        printf("-E- %s line %d: Error allocating raster space.\n", __FILE__,
                __LINE__);
        return (1);
    }
    err_code = epr_read_band_raster(band_id, 0, scan, temp_raster);
    for (ip = 0; ip < npix; ip++)
        l1rec->lat[ip] = epr_get_pixel_as_float(temp_raster, spix+ip, 0);

    // read longitude

    band_id = epr_get_band_id(fileID, "longitude");
    if (band_id == NULL) {
        printf("-E- %s line %d: Error finding band_id:longitude\n", __FILE__,
                __LINE__);
        return (1);
    }
    err_code = epr_read_band_raster(band_id, 0, scan, temp_raster);
    for (ip = 0; ip < npix; ip++)
        l1rec->lon[ip] = epr_get_pixel_as_float(temp_raster, spix+ip, 0);
    epr_free_raster(temp_raster);

    // read sun zenith

    band_id = epr_get_band_id(fileID, "sun_zenith");
    if (band_id == NULL) {
        printf("-E- %s line %d: Error finding band_id:sun_zenith", __FILE__,
                __LINE__);
        return (1);
    }
    temp_raster = epr_create_compatible_raster(band_id, file_npix, 1, 1, 1);
    if (temp_raster == NULL) {
        printf("-E- %s line %d: Error allocating raster space.\n", __FILE__,
                __LINE__);
        return (1);
    }
    err_code = epr_read_band_raster(band_id, 0, scan, temp_raster);
    for (ip = 0; ip < npix; ip++)
        l1rec->solz[ip] = epr_get_pixel_as_float(temp_raster, spix+ip, 0);
    epr_free_raster(temp_raster);

    // read sun azimuth

    band_id = epr_get_band_id(fileID, "sun_azimuth");
    if (band_id == NULL) {
        printf("-E- %s line %d: Error finding band_id:sun_azimuth", __FILE__,
                __LINE__);
        return (1);
    }
    temp_raster = epr_create_compatible_raster(band_id, file_npix, 1, 1, 1);
    if (temp_raster == NULL) {
        printf("-E- %s line %d: Error allocating raster space.\n", __FILE__,
                __LINE__);
        return (1);
    }
    err_code = epr_read_band_raster(band_id, 0, scan, temp_raster);
    for (ip = 0; ip < npix; ip++)
        l1rec->sola[ip] = epr_get_pixel_as_float(temp_raster, spix+ip, 0);
    epr_free_raster(temp_raster);

    // read view zenith

    band_id = epr_get_band_id(fileID, "view_zenith");
    if (band_id == NULL) {
        printf("-E- %s line %d: Error finding band_id:view_zenith", __FILE__,
                __LINE__);
        return (1);
    }
    temp_raster = epr_create_compatible_raster(band_id, file_npix, 1, 1, 1);
    if (temp_raster == NULL) {
        printf("-E- %s line %d: Error allocating raster space.\n", __FILE__,
                __LINE__);
        return (1);
    }
    err_code = epr_read_band_raster(band_id, 0, scan, temp_raster);
    for (ip = 0; ip < npix; ip++)
        l1rec->senz[ip] = epr_get_pixel_as_float(temp_raster, spix+ip, 0);
    epr_free_raster(temp_raster);

    // read view azimuth

    band_id = epr_get_band_id(fileID, "view_azimuth");
    if (band_id == NULL) {
        printf("-E- %s line %d: Error finding band_id:view_azimuth", __FILE__,
                __LINE__);
        return (1);
    }
    temp_raster = epr_create_compatible_raster(band_id, file_npix, 1, 1, 1);
    if (temp_raster == NULL) {
        printf("-E- %s line %d: Error allocating raster space.\n", __FILE__,
                __LINE__);
        return (1);
    }
    err_code = epr_read_band_raster(band_id, 0, scan, temp_raster);
    for (ip = 0; ip < npix; ip++)
        l1rec->sena[ip] = epr_get_pixel_as_float(temp_raster, spix+ip, 0);
    epr_free_raster(temp_raster);

    // set pixnum and check for navigation failure

    for (ip = 0; ip < npix; ip++) {
        l1rec->pixnum[ip] = spix + ip;
        if (l1rec->lon[ip] < -181.0 || l1rec->lon[ip] > 181.0 ||
                l1rec->lat[ip] < -91.0 || l1rec->lat[ip] > 91.0)
            l1rec->navfail[ip] = 1;
    }

    if ( isLevel1 ) {

        // read detector_index
        band_id  = epr_get_band_id(fileID, "detector_index");
        if (band_id == NULL) {
            printf("-E- %s line %d: Error finding band_id:detector_index\n", __FILE__,
                    __LINE__);
            return (1);
        }
        temp_raster = epr_create_compatible_raster(band_id, file_npix, 1, 1, 1);
        if (temp_raster == NULL) {
            printf("-E- %s line %d: Error allocating raster space.\n", __FILE__,
                    __LINE__);
            return (1);
        }
        err_code = epr_read_band_raster(band_id, 0, scan, temp_raster);
        for (ip = 0; ip < npix; ip++) {
            l1rec->pixdet[ip] = epr_get_pixel_as_uint(temp_raster, spix+ip, 0);
        }
        epr_free_raster(temp_raster);


        // if Level-2 file only get data processed as "water"

    } else {

        band_id = epr_get_band_id(fileID, "l2_flags");
        if (band_id == NULL) {
            printf("-E- %s line %d: Error finding band_id:l2_flags", __FILE__,
                    __LINE__);
            return (1);
        }
        bm_raster = epr_create_compatible_raster(band_id, file_npix, 1, 1, 1);
        if (bm_raster == NULL) {
            printf("-E- %s line %d: Error allocating bitmap raster space.\n",
                    __FILE__, __LINE__);
            return (1);
        }
        err_code = epr_read_band_raster(band_id, 0, scan, bm_raster);
        if (err_code) {
            printf("-E- %s line %d: Error reading bitmap raster space.\n",
                    __FILE__, __LINE__);
            return (1);
        }
        for (ip = 0; ip < npix; ip++)
            l1rec->in_flags[ip] = epr_get_pixel_as_uint(bm_raster, spix+ip, 0);
        epr_free_raster(bm_raster);

        if (file->n_inprods > 0) {

            // read in algal_1 and store in l1rec->in_prods[0]

            band_id = epr_get_band_id(fileID, "algal_1");
            if (band_id == NULL) {
                l1rec->in_prods[0][0] = BAD_FLT; /* flag as not found */
            } else {
                temp_raster = 
                        epr_create_compatible_raster(band_id, file_npix, 1, 1, 1);
                if (temp_raster == NULL) {
                    printf("-E- %s line %d: Error allocating raster space.\n",
                            __FILE__, __LINE__);
                    return (1);
                }
                err_code = epr_read_band_raster(band_id, 0, scan, temp_raster);
                for (ip = 0; ip < npix; ip++)
                    l1rec->in_prods[0][ip] =
                            epr_get_pixel_as_float(temp_raster, spix+ip, 0);
                epr_free_raster(temp_raster);
            }

            // read in algal_2 and store in l1rec->in_prods[1]

            band_id = epr_get_band_id(fileID, "algal_2");
            if (band_id == NULL) {
                l1rec->in_prods[1][0] = BAD_FLT; /* flag as not found */
            } else {
                temp_raster =
                        epr_create_compatible_raster(band_id, file_npix, 1, 1, 1);
                if (temp_raster == NULL) {
                    printf("-E- %s line %d: Error allocating raster space.\n",
                            __FILE__, __LINE__);
                    return (1);
                }
                err_code = epr_read_band_raster(band_id, 0, scan, temp_raster);
                for (ip = 0; ip < npix; ip++)
                    l1rec->in_prods[1][ip] =
                            epr_get_pixel_as_float(temp_raster, spix+ip, 0);
                epr_free_raster(temp_raster);
            }

            // read in yellow_subs and store in l1rec->in_prods[2]

            band_id = epr_get_band_id(fileID, "yellow_subs");
            if (band_id == NULL) {
                l1rec->in_prods[2][0] = BAD_FLT; /* flag as not found */
            } else {
                temp_raster =
                        epr_create_compatible_raster(band_id, file_npix, 1, 1, 1);
                if (temp_raster == NULL) {
                    printf("-E- %s line %d: Error allocating raster space.\n",
                            __FILE__, __LINE__);
                    return (1);
                }
                err_code = epr_read_band_raster(band_id, 0, scan, temp_raster);
                for (ip = 0; ip < npix; ip++)
                    l1rec->in_prods[2][ip] =
                            epr_get_pixel_as_float(temp_raster, spix+ip, 0);
                epr_free_raster(temp_raster);
            }

            // read in total_susp and store in l1rec->in_prods[3]

            band_id = epr_get_band_id(fileID, "total_susp");
            if (band_id == NULL) {
                l1rec->in_prods[3][0] = BAD_FLT; /* flag as not found */
            } else {
                temp_raster =
                        epr_create_compatible_raster(band_id, file_npix, 1, 1, 1);
                if (temp_raster == NULL) {
                    printf("-E- %s line %d: Error allocating raster space.\n",
                            __FILE__, __LINE__);
                    return (1);
                }
                err_code = epr_read_band_raster(band_id, 0, scan, temp_raster);
                for (ip = 0; ip < npix; ip++)
                    l1rec->in_prods[3][ip] =
                            epr_get_pixel_as_float(temp_raster, spix+ip, 0);
                epr_free_raster(temp_raster);
            }

            // read in photosyn_rad and store in l1rec->in_prods[4]

            band_id = epr_get_band_id(fileID, "photosyn_rad");
            if (band_id == NULL) {
                l1rec->in_prods[4][0] = BAD_FLT; /* flag as not found */
            } else {
                temp_raster =
                        epr_create_compatible_raster(band_id, file_npix, 1, 1, 1);
                if (temp_raster == NULL) {
                    printf("-E- %s line %d: Error allocating raster space.\n",
                            __FILE__, __LINE__);
                    return (1);
                }
                err_code = epr_read_band_raster(band_id, 0, scan, temp_raster);
                for (ip = 0; ip < npix; ip++)
                    l1rec->in_prods[4][ip] =
                            epr_get_pixel_as_float(temp_raster, spix+ip, 0);
                epr_free_raster(temp_raster);
            }

            // read in aero_alpha and store in l1rec->in_prods[5]

            band_id = epr_get_band_id(fileID, "aero_alpha");
            if (band_id == NULL) {
                l1rec->in_prods[5][0] = BAD_FLT; /* flag as not found */
            } else {
                temp_raster =
                        epr_create_compatible_raster(band_id, file_npix, 1, 1, 1);
                if (temp_raster == NULL) {
                    printf("-E- %s line %d: Error allocating raster space.\n",
                            __FILE__, __LINE__);
                    return (1);
                }
                err_code = epr_read_band_raster(band_id, 0, scan, temp_raster);
                for (ip = 0; ip < npix; ip++)
                    l1rec->in_prods[5][ip] =
                            epr_get_pixel_as_float(temp_raster, spix+ip, 0);
                epr_free_raster(temp_raster);
            }

            // read in aero_opt_thick and store in l1rec->in_prods[6]

            band_id = epr_get_band_id(fileID, "aero_opt_thick");
            if (band_id == NULL) {
                l1rec->in_prods[6][0] = BAD_FLT; /* flag as not found */
            } else {
                temp_raster =
                        epr_create_compatible_raster(band_id, file_npix, 1, 1, 1);
                if (temp_raster == NULL) {
                    printf("-E- %s line %d: Error allocating raster space.\n",
                            __FILE__, __LINE__);
                    return (1);
                }
                err_code = epr_read_band_raster(band_id, 0, scan, temp_raster);
                for (ip = 0; ip < npix; ip++)
                    l1rec->in_prods[6][ip] =
                            epr_get_pixel_as_float(temp_raster, spix+ip, 0);
                epr_free_raster(temp_raster);
            }
        }
    }

    // read in data

    for (ib = 0; ib < nbands; ib++) {

        band_id = epr_get_band_id(fileID, names[ib]);
        if (band_id == NULL) {
            if(isLevel1) {
                printf("-E- %s line %d: Error finding band_id:%s\n", __FILE__,
                        __LINE__, names[ib]);
                return (1);
            }

            // fill missing L2 bands with 0

            if((ib == 10) || (ib == 14)) {
                for (ip = 0; ip < npix; ip++) {
                    ipb = ip*nbands + ib;
                    l1rec->Lt[ipb] = 0.0;
                }
            } else {
                printf("-E- %s line %d: Error finding band_id:%s\n", __FILE__,
                        __LINE__, names[ib]);
                return (1);
            }

        } else {
            temp_raster = epr_create_compatible_raster(band_id, file_npix, 1, 1, 1);
            if (temp_raster == NULL) {
                printf("-E- %s line %d: Error allocating raster space.\n",
                        __FILE__, __LINE__);
                return (1);
            }
            err_code = epr_read_band_raster(band_id, 0, scan, temp_raster);
            if (!err_code) {

                // copy to Lt record.  Note that this might actually be Rrs
                // if the input file is a MERIS L2 file

                for (ip = 0; ip < npix; ip++) {

                    ipb = ip*nbands + ib;
                    l1rec->Lt[ipb] = epr_get_pixel_as_float(temp_raster, spix+ip, 0);
                    l1rec->detnum = spix+ip;

                    if (isLevel1) {

                        // set the value of all of the bands to BAD_FLT if
                        // navfail has been flagged.
                        if(invalid_flag[ip])
                            l1rec->Lt[ipb] = BAD_FLT;
                        else
                            l1rec->Lt[ipb] /= 10.0;     // units conversion

                        // mark negative input data as HILT
                        if (l1rec->Lt[ipb] < 0.0)
                            l1rec->hilt[ip] = 1;

                    } else {

                        // MERIS L2 divides each pixel into three classes and 
                        // we only want the water pixels

                        if (l1rec->in_flags[ip] & MERIS_L2FLAG_LAND) { /* bit 23 = land */
                            l1rec->Lt[ipb]  = 0.0;
                        } else if (l1rec->in_flags[ip] & MERIS_L2FLAG_CLOUD) { /* bit 22 = cloud */
                            l1rec->Lt[ipb]  = 10.0;
                            l1rec->hilt[ip] = 1;
                        }
                    }
                }
            }
            epr_free_raster(temp_raster);

        } // if band_id found

    } // for ib

    l1rec->sensorID = file->sensorID;
    l1rec->npix     = file->npix;

    return (LIFE_IS_GOOD);
}


/**
 * @brief reads 1 scan line from MERIS file, loads l1rec
 * @param[in]   file   file handle to MERIS file
 * @param[in]   scan   scan number to read
 * @apram[out]  l1rec  output l1rec
 *
 * @author Paul E. Lyon NRL, Oct. 2006.
 */

int readl1_lonlat_meris_N1(filehandle *file, int32 scan, l1str *l1rec)
{
    static int          firstCall = 1;

    static char       **names;
    int32               npix;
    int32               ip, ib, ipb;
    static EPR_SBandId        *lon_band_id;
    static EPR_SBandId        *lat_band_id;
    static EPR_SRaster        *temp_raster;


    npix = file->npix;


    if (firstCall) {

        lon_band_id = epr_get_band_id(fileID, "longitude");
        if (lon_band_id == NULL) {
            printf("-E- %s line %d: Error finding band_id:longitude\n", 
                    __FILE__, __LINE__);
            return (1);
        }

        lat_band_id = epr_get_band_id(fileID, "latitude");
        if (lat_band_id == NULL) {
            printf("-E- %s line %d: Error finding band_id:latitude\n", 
                    __FILE__, __LINE__);
            return (1);
        }

        temp_raster = epr_create_compatible_raster(lat_band_id, file_npix, 1, 1, 1);
        if (temp_raster == NULL) {
            printf("-E- %s line %d: Error allocating raster space.\n", 
                    __FILE__, __LINE__);
            return (1);
        }

        firstCall = 0;
    }


    // read latitude

    epr_read_band_raster(lat_band_id, 0, scan, temp_raster);
    for (ip = 0; ip < npix; ip++)
        l1rec->lat[ip] = epr_get_pixel_as_float(temp_raster, spix+ip, 0);

    // read longitude

    epr_read_band_raster(lon_band_id, 0, scan, temp_raster);
    for (ip = 0; ip < npix; ip++)
        l1rec->lon[ip] = epr_get_pixel_as_float(temp_raster, spix+ip, 0);


    return (LIFE_IS_GOOD);
}


/**
 * @brief closes MERIS file, loads l1rec
 * @param[in]   file   file handle to MERIS file
 *
 * @author Paul E. Lyon NRL, Oct. 2006.
 */

int
closel1_meris_N1(filehandle *file)
{
    if (epr_close_product(fileID)) {
        fprintf(stderr,
                "-E- %s line %d: epr_close_product failed for file, %s.\n",
                __FILE__, __LINE__, file->name);
        return (1);
    }

    return (LIFE_IS_GOOD);
}

