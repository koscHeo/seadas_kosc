/* ============================================================================ */
/* module l1_viirs_nc.c - functions to read VIIRS L1A for MSL12                 */
/* Written By:  Joel M. Gales GSFC Futuretech, Sep. 2011.                       */
/*                                                                              */
/* ============================================================================ */

// #include <stdbool.h>
#include "l12_proto.h"
#include "calibrate_viirs.h"
#include <Calibrate_Viirs_Connector.h>

#define MBAND_NUM_DETECTORS 16

static int32 lastScan = -1;

static int scanLineGrp;
static int scanStartTimeId;
static int HAMSideId;

static int geoFileId;
static int geoNavigationGrp;
static int geoGeolocationGrp;
static int geoScanLineGrp;
static int lonId, latId, senzId, senaId, solzId, solaId, esdistId, scanQualityId, pixelQualityId;

static double ***f_cal_corr = NULL; /* f table correction [band][det][ms] */

static float *Fobar; // reflectance to radiance conversion factors
static int extract_pixel_start = 0;
static int extract_pixel_stop = 0;

static short *tmpShort;
static unsigned char *tmpByte;
static int nscan, nline, npix;

static double starttime;
static double lastvalidtime;

static double scan_start_tai = -999;

int openl1_viirs_nc(filehandle * file) {
    char *fltime;

    size_t tmpSize;
    int fileID, dimid, status;
    size_t att_len; // Change to size_t JMG 08/12/13
    int orbit_number;

    // Open the netcdf4 input file
    status = nc_open(file->name, NC_NOWRITE, &fileID);
    if (status == FAIL) {
        fprintf(stderr, "-E- %s line %d: nc_open(%s) failed.\n",
        __FILE__, __LINE__, file->name);
        return (1);
    }

    // Get pixel and scan dimensions
    status = nc_get_att_int(fileID, NC_GLOBAL, "number_of_filled_scans", &nscan);
    check_err(status, __LINE__, __FILE__);
    nline = nscan*MBAND_NUM_DETECTORS;
    
    status = nc_inq_dimid(fileID, "Mband_pixels", &dimid);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_dimlen(fileID, dimid, &tmpSize);
    check_err(status, __LINE__, __FILE__);
    npix = tmpSize;

    nc_type vr_type;  /* attribute type */
    size_t  vr_len;   /* attribute length */
    if ((nc_inq_att (fileID, NC_GLOBAL, "extract_pixel_start", &vr_type, &vr_len) == 0)) {
        status = nc_get_att_int(fileID, NC_GLOBAL, "extract_pixel_start", &extract_pixel_start);
        check_err(status, __LINE__, __FILE__);
        extract_pixel_start--; // Attribute is one-based
        status = nc_get_att_int(fileID, NC_GLOBAL, "extract_pixel_stop", &extract_pixel_stop);
        check_err(status, __LINE__, __FILE__);
        extract_pixel_stop--; // Attribute is one-based
        if(npix != (extract_pixel_stop - extract_pixel_start + 1)) {
            printf("-E- Problem with the extracted L1A file pixel dimension.\n");
            printf("    npix(%d), extract_pixel_stop(%d), extract_pixel_start(%d) do not work together.\n",
                    npix, extract_pixel_stop, extract_pixel_start);
            exit(EXIT_FAILURE);
        }
    }

    if (want_verbose) {
        printf("VIIRS L1A Npix  :%d Nlines:%d\n", npix, nline);
    } // want_verbose

    // get start and end time
    status = nc_inq_attlen(fileID, NC_GLOBAL, "time_coverage_start", &att_len);
    check_err(status, __LINE__, __FILE__);

    /* allocate required space before retrieving values */
    fltime = (char *) malloc(att_len + 1); /* + 1 for trailing null */

    /* get attribute values */
    status = nc_get_att_text(fileID, NC_GLOBAL, "time_coverage_start", fltime);
    check_err(status, __LINE__, __FILE__);
    fltime[att_len] = '\0';
//    isodate2ydmsec(fltime, (int32_t *) &year,(int32_t *) &day, (int32_t *) &msec);

    starttime = lastvalidtime = isodate2unix(fltime);
    free(fltime);

    status = nc_inq_attlen(fileID, NC_GLOBAL, "time_coverage_end", &att_len);
    check_err(status, __LINE__, __FILE__);

    /* allocate required space before retrieving values */
    fltime = (char *) malloc(att_len + 1); /* + 1 for trailing null */

    /* get attribute values */
    status = nc_get_att_text(fileID, NC_GLOBAL, "time_coverage_end", fltime);
    check_err(status, __LINE__, __FILE__);
    fltime[att_len] = '\0';

    double stoptime = isodate2unix(fltime);
    free(fltime);

    if ((nc_inq_att(fileID, NC_GLOBAL, "orbit_number", &vr_type, &vr_len) == 0)) {
        status = nc_get_att_int(fileID, NC_GLOBAL, "orbit_number", &orbit_number);
        check_err(status, __LINE__, __FILE__);
    } else {
        status = nc_get_att_int(fileID, NC_GLOBAL, "OrbitNumber", &orbit_number);
        check_err(status, __LINE__, __FILE__);
    }

    file->sd_id = fileID;
    file->nbands = 10;
    file->npix = npix;
    file->nscan = nline;
    file->ndets = MBAND_NUM_DETECTORS;
    file->terrain_corrected = 1; // presumed.
    file->orbit_number = orbit_number;
    strcpy(file->spatialResolution, "750 m");

    rdsensorinfo(file->sensorID, file->input->evalmask,
                 "Fobar", (void **) &Fobar);

    if (want_verbose)
        printf("file->nbands = %d\n", (int) file->nbands);

    status = nc_inq_ncid(file->sd_id, "scan_line_attributes", &scanLineGrp);
    check_err(status, __LINE__, __FILE__);

    // get tai93
    status = nc_inq_varid(scanLineGrp, "scan_start_time", &scanStartTimeId);
    check_err(status, __LINE__, __FILE__);

    // get mirror side
    status = nc_inq_varid(scanLineGrp, "HAM_side", &HAMSideId);
    check_err(status, __LINE__, __FILE__);

    // Setup geofile pointers
    if(file->geofile[0] != '\0') {
        
        status = nc_open(file->geofile, NC_NOWRITE, &geoFileId);
        check_err(status, __LINE__, __FILE__);

        status = nc_inq_grp_ncid(geoFileId, "geolocation_data", &geoGeolocationGrp);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(geoGeolocationGrp, "longitude", &lonId);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(geoGeolocationGrp, "latitude", &latId);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(geoGeolocationGrp, "sensor_zenith", &senzId);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(geoGeolocationGrp, "sensor_azimuth", &senaId);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(geoGeolocationGrp, "solar_zenith", &solzId);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(geoGeolocationGrp, "solar_azimuth", &solaId);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(geoGeolocationGrp, "quality_flag", &pixelQualityId);
        check_err(status, __LINE__, __FILE__);

        status = nc_inq_grp_ncid(geoFileId, "navigation_data", &geoNavigationGrp);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(geoNavigationGrp, "earth_sun_distance", &esdistId);
        check_err(status, __LINE__, __FILE__);

        status = nc_inq_grp_ncid(geoFileId, "scan_line_attributes", &geoScanLineGrp);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(geoScanLineGrp, "scan_quality", &scanQualityId);
        check_err(status, __LINE__, __FILE__);
    } // geofile

    return (LIFE_IS_GOOD);
}

int readl1_viirs_nc(filehandle *file, int32 line, l1str *l1rec)
{
    static int firstCall = 1;
   
    int32 ip, ib, ipb;
    int i;
    double scan_sec;
    int16 scan_year, scan_day;
    double f_corr;

    int status, ncid, varid;
    size_t start[] = { 0,0,0 };
    size_t count[] = { 1,1,1 };

    int32 scan = line / MBAND_NUM_DETECTORS;
    
    for (ip = 0; ip < npix; ip++) {
      l1rec->pixnum[ip] = ip + extract_pixel_start;
      flag_bowtie_deleted(l1rec, ip, extract_pixel_start);
    }

    if (firstCall) {
        firstCall = 0;

        // init calibration only if a geo file was given
        if(file->geofile[0] != '\0') {

            /*----- Calibration LUT -----*/
            double tai2unixoffset =  725846400.0;
            if (file->input->calfile[0] != '\0') {
                double grantime= starttime-tai2unixoffset;
                // convert TAI93 to "UTC93"
                grantime -= leapseconds_since_1993(grantime); 
                grantime += 1104537600.0; // convert "UTC93" to UTC (1958) (verify)
                grantime *= 1000000.0; // convert to IET
                load_fcal_lut(file->input->calfile, (int64_t) grantime, &f_cal_corr);
            }

            // Initialize L1A calibration
            VcstViirsCal_initialize(file->name, file->viirscalparfile);

            tmpShort = (short *) malloc(npix * sizeof(short));
            tmpByte = (unsigned char *) malloc(npix);
        }       
    }

    l1rec->sensorID = file->sensorID;
    l1rec->npix = file->npix;

    // set time for this scan - if scan_start_time value not properly set, estimate from scene start time.
    if(scan != lastScan) {
        start[0] = scan;
        status = nc_get_var1_double(scanLineGrp, scanStartTimeId, start, &scan_start_tai);
        check_err(status, __LINE__, __FILE__);
    }
    
    double tai2unixoffset =  725846400.0;
    if (scan_start_tai > 0){
        lastvalidtime = scan_start_tai+tai2unixoffset;
        unix2yds(lastvalidtime, &scan_year, &scan_day, &scan_sec);
        *(l1rec->year) = (int32_t) scan_year;
        *(l1rec->day) = (int32_t) scan_day;
        *(l1rec->msec) = (int32_t) (scan_sec*1e3);
        *(l1rec->msec) -= leapseconds_since_1993(lastvalidtime-tai2unixoffset)*1e3;
    } else {
        for (ip = 0; ip < npix; ip++) 
            l1rec->flags[ip] |= NAVFAIL;
        *(l1rec->year) = -999;
        *(l1rec->day) = -999;
        *(l1rec->msec) = -999;
    }

    //------------------------------------
    // if there is no geo file just return
    // This is used for l1info with only a L1A file and no GEO
    //-------------------------------------
    if(file->geofile[0] == '\0') {
        return 0;
    }

    // first check the scan quality flag
    // 1   SCE_side_A_B 
    // 2   SCE_side_invalid 
    // 4   Sector_rotation 
    // 8   Encoder_degraded 
    // 16  SAA 
    // 32  Solar_eclipse 
    // 64  Lunar_eclipse 
    // 128 HAM_side
    //
    // Sector_rotation
    short scanQualityWarnMask = 2 | 8 | 128;
    short scanQualityFailMask = 4;
    static short scanQualityFlag = 0;
    
    if(scan != lastScan) {
        start[0] = scan;
        status = nc_get_var1_short(geoScanLineGrp, scanQualityId, start, &scanQualityFlag);
        check_err(status, __LINE__, __FILE__);
    }
    if(scanQualityFlag & scanQualityFailMask) {
        for (ip = 0; ip < npix; ip++) 
            l1rec->flags[ip] |= NAVFAIL;
        return 0;
    }
    if(scanQualityFlag & scanQualityWarnMask) {
        for (ip = 0; ip < npix; ip++) 
            l1rec->flags[ip] |= NAVWARN;
    }
   
    static unsigned char HAMSideVal = 0;
    if(scan != lastScan) {
        start[0] = scan;
        status = nc_get_var1_uchar(scanLineGrp, HAMSideId, start, &HAMSideVal);
        check_err(status, __LINE__, __FILE__);
    }
    l1rec->mside = HAMSideVal;
    
    // set up to read all pixels of the line.
    start[0] = line;
    start[1] = 0;
    count[0] = 1;
    count[1] = npix; // read all pixels
    
    status = nc_get_vara_float(geoGeolocationGrp, latId, start, count, l1rec->lat);
    check_err(status, __LINE__, __FILE__);

    status = nc_get_vara_float(geoGeolocationGrp, lonId, start, count, l1rec->lon);
    check_err(status, __LINE__, __FILE__);

    status = nc_get_vara_short(geoGeolocationGrp, solzId, start, count, tmpShort);
    check_err(status, __LINE__, __FILE__);
    for (i = 0; i < npix; i++)
      if (tmpShort[i] != -32767) l1rec->solz[i] = tmpShort[i]*0.01;

    status = nc_get_vara_short(geoGeolocationGrp, solaId, start, count, tmpShort);
    check_err(status, __LINE__, __FILE__);
    for (i = 0; i < npix; i++)
      if (tmpShort[i] != -32767) l1rec->sola[i] = tmpShort[i]*0.01;

    status = nc_get_vara_short(geoGeolocationGrp, senzId, start, count, tmpShort);
    check_err(status, __LINE__, __FILE__);
    for (i = 0; i < npix; i++)
      if (tmpShort[i] != -32767) l1rec->senz[i] = tmpShort[i]*0.01;

    status = nc_get_vara_short(geoGeolocationGrp, senaId, start, count, tmpShort);
    check_err(status, __LINE__, __FILE__);
    for (i = 0; i < npix; i++)
      if (tmpShort[i] != -32767) l1rec->sena[i] = tmpShort[i]*0.01;
    
    status = nc_get_vara_uchar(geoGeolocationGrp, pixelQualityId, start, count, tmpByte);
    check_err(status, __LINE__, __FILE__);
    // 1 Input_invalid
    // 2 Pointing_bad
    // 4 Terrain_bad
    unsigned char qualityFailMask = 1 || 2; 
    unsigned char qualityWarnMask = 4;
    for (i = 0; i < npix; i++) {
      if (tmpByte[i] & qualityFailMask) 
          l1rec->flags[i] |= NAVFAIL;
      if (tmpByte[i] & qualityWarnMask) 
          l1rec->flags[i] |= NAVWARN;
    }   

    static float esdist = -999.9;
    if(scan != lastScan) {
        start[0] = scan;
        status = nc_get_var1_float(geoNavigationGrp, esdistId, start, &esdist);
        check_err(status, __LINE__, __FILE__);
    }
    
    /* Earth-sun distance correction for this scan */
    if(esdist > -999.0)
        l1rec->fsol = pow(1.0 / esdist, 2);

    int nbands=16, nRSBbands=10, nCIRbands=1, nTEBbands=5;

    // read in calibrated L1B data
    static float *l1bptrs[16]={0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0,
                               0x0,0x0,0x0,0x0,0x0,0x0,0x0,0x0};

    static char *bandType[16]={"RSB","RSB","RSB","RSB","RSB","RSB",
                               "RSB","RSB","CIR","RSB","RSB","TEB",
                               "TEB","TEB","TEB","TEB"};

    // Note: l1bptrs arrays are 3200 pixels wide
    int oldVerbose = want_verbose;
    want_verbose = 0;
    VcstViirsCal_calibrateMOD(line, nbands, l1bptrs);
    want_verbose = oldVerbose;
    
    l1rec->detnum = line % l1rec->ndets;

    int irsb=0, iteb=0;
    int l1bptrs_scan_idx;
    for (ib = 0; ib < nbands; ib++) {

        /* get specific f table cal correction  */
        f_corr = (f_cal_corr == NULL) ? 1.0
                : f_cal_corr[ib][l1rec->detnum][l1rec->mside];

        if (strcmp(bandType[ib], "TEB") == 0) {

            for (ip = 0; ip < npix; ip++) {
                ipb = ip * NBANDSIR + iteb;
                l1rec->Ltir[ipb] = 0;

                if (l1bptrs[ib][l1rec->detnum * npix + ip] != -32767) {
                    l1rec->Ltir[ipb] = l1bptrs[ib][l1rec->detnum * npix + ip] / 10.0;

                    /* Apply F-factor */
                    l1rec->Ltir[ipb] *= f_corr;
                }

            }
            iteb++;

        } else if (strcmp(bandType[ib], "CIR") == 0) {

            for (ip = 0; ip < npix; ip++) {

                l1bptrs_scan_idx = l1rec->detnum * 3200 + ip + extract_pixel_start;

                if (l1bptrs[ib][l1bptrs_scan_idx] != -32767) {
                    l1rec->rho_cirrus[ip] = l1bptrs[ib][l1bptrs_scan_idx];

                    /* Apply F-factor */
                    l1rec->rho_cirrus[ip] *= f_corr;
                }

            }

        } else if (strcmp(bandType[ib], "RSB") == 0) {

            l1rec->Fo[irsb] = Fobar[irsb] * l1rec->fsol;

            // copy to Lt record.
            for (ip = 0; ip < npix; ip++) {
                ipb = ip * file->nbands + irsb;

                l1bptrs_scan_idx = l1rec->detnum * 3200 + ip + extract_pixel_start;

                if (l1bptrs[ib][l1bptrs_scan_idx] != -32767) {
                    l1rec->Lt[ipb] = l1bptrs[ib][l1bptrs_scan_idx];

                    /* convert from reflectance to radiance */
                    l1rec->Lt[ipb] *= l1rec->Fo[irsb] / PI;

                    /* Apply F-factor */
                    l1rec->Lt[ipb] *= f_corr;
                }

            }
            irsb++;
        } // if RSB

    } // for ib

    radiance2bt(l1rec, -1); // calculate brightness temperature
    
    lastScan = scan;
    return (LIFE_IS_GOOD);
}

int readl1_lonlat_viirs_nc(filehandle *file, int32 line, l1str *l1rec)
{
    double scan_sec;
    int16 scan_year, scan_day;
    int32 ip;
    int status;
    size_t start[] = { 0,0,0 };
    size_t count[] = { 1,1,1 };

    int32 scan = line / MBAND_NUM_DETECTORS;
    
    if(file->geofile[0] == '\0') {
        printf("-E- Geolocation file needs to be set\n");
        exit(EXIT_FAILURE);
    }
    
    l1rec->sensorID = file->sensorID;
    l1rec->npix = file->npix;

    // set time for this scan - if scan_start_time value not properly set, estimate from scene start time.
    if(scan != lastScan) {
        start[0] = scan;
        status = nc_get_var1_double(scanLineGrp, scanStartTimeId, start, &scan_start_tai);
        check_err(status, __LINE__, __FILE__);
    }
        
    double tai2unixoffset =  725846400.0;
    if (scan_start_tai > 0){
        lastvalidtime = scan_start_tai+tai2unixoffset;
        unix2yds(lastvalidtime, &scan_year, &scan_day, &scan_sec);
        *(l1rec->year) = (int32_t) scan_year;
        *(l1rec->day) = (int32_t) scan_day;
        *(l1rec->msec) = (int32_t) (scan_sec*1e3);
        *(l1rec->msec) -= leapseconds_since_1993(lastvalidtime-tai2unixoffset)*1e3;
    } else {
        for (ip = 0; ip < npix; ip++) 
            l1rec->flags[ip] |= NAVFAIL;
        *(l1rec->year) = -999;
        *(l1rec->day) = -999;
        *(l1rec->msec) = -999;
    }

    // set up to read all pixels of the line.
    start[0] = line;
    start[1] = 0;
    count[0] = 1;
    count[1] = npix; // read all pixels
    
    status = nc_get_vara_float(geoGeolocationGrp, latId, start, count, l1rec->lat);
    check_err(status, __LINE__, __FILE__);

    status = nc_get_vara_float(geoGeolocationGrp, lonId, start, count, l1rec->lon);
    check_err(status, __LINE__, __FILE__);

    lastScan = scan;
    return (LIFE_IS_GOOD);
}

int closel1_viirs_nc(filehandle *file) {
    int status;

    status = nc_close(file->sd_id);
    check_err(status, __LINE__, __FILE__);
    
    if(file->geofile[0] == '\0') {
        status = nc_close(geoFileId);
        check_err(status, __LINE__, __FILE__);

        free(tmpShort);
        free(tmpByte);
    }
    
    return (LIFE_IS_GOOD);
}

