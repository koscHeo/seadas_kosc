/* ============================================================================ */
/* module l1_meris_CC.c - functions to read MERIS CC (coastal color) for MSL12  */
/* Written By:  Joel M. Gales GSFC Futuretech, Sep. 2011.                       */
/*                                                                              */
/*  W. Robinson, SAIC, 22 May 2012 make repair to start, end and scan times     */
/* ============================================================================ */

// #include <stdbool.h>
#include "l12_proto.h"

static int ncol;
static int16 nline, npix;

static int32 spix = 0;
static double starttime;
static double *scan_start_tai, *scan_stop_tai;
static double lastvalidtime;
static int lastvalidscan = 0;
static double time_interval;

int openl1_orca(filehandle * file) {
    char *fltime;

    size_t source_w, source_h;
    int32 nscan;
    int fileID, xid, yid, geo_id, retval, grp_id, sds_id;
    int minute, hour, month, sec, i;
    int year2, day2, msec2;
    float fsec, milisec, fnpix;
    size_t att_len; // Change to size_t JMG 08/12/13
    double dblsec, stoptime;
    int16 y16, m16, d16;
    size_t nctlpix = 0;
    size_t nctlscan = 0;
    size_t start[3], count[3];
    int16_t year=2014, day=127, secs=43200; // 05/07/2014

    // Open the netcdf4 input file
    retval = nc_open(file->name, NC_NOWRITE, &fileID);
    if (retval == FAIL) {
        fprintf(stderr, "-E- %s line %d: nc_open(%s) failed.\n",
        __FILE__, __LINE__, file->name);
        return (1);
    }

    // Get pixel and scan dimensions
    retval = nc_inq_dimid(fileID, "number_of_lines", &yid);
    retval = nc_inq_dimid(fileID, "pixels_per_line", &xid);
    retval = nc_inq_dimlen(fileID, xid, &source_w);
    retval = nc_inq_dimlen(fileID, yid, &source_h);

    if (want_verbose) {
        printf("ORCA L1B Npix  :%d Nscans:%d\n", (int) source_w,
                (int) source_h);
    } // want_verbose
    npix = (int32) source_w;
    nscan = (int32) source_h;
    nline = nscan;

    file->sd_id = fileID;
    file->nbands = 91;
    file->npix = npix;
    file->nscan = nscan;
    strcpy(file->spatialResolution, "750 m");

    retval = nc_inq_dimid(file->sd_id, "tp_x", &xid);
    retval = nc_inq_dimlen(file->sd_id, xid, &nctlpix);
    ncol = nctlpix;

    retval = nc_inq_dimid(file->sd_id, "tp_y", &yid);
    retval = nc_inq_dimlen(file->sd_id, yid, &nctlscan);
    scan_start_tai = (double *) calloc(nscan,sizeof(double));
    scan_stop_tai  = (double *) calloc(nscan,sizeof(double));

    start[0] = 0;
    count[0] = nscan;

    nc_inq_ncid(file->sd_id,"scan_line_attributes" , &grp_id);

    nc_inq_varid(grp_id,"scan_start_time",&sds_id);
    nc_get_vara_double(grp_id, sds_id, start, count, scan_stop_tai);

    /*
     *  Temporary fix for scan_start and stop
     */

    for (i=0; i<nscan; i++) {
        scan_start_tai[i] = yds2unix(year, day, secs+i);
    }

    starttime = scan_start_tai[0];

//    nc_inq_varid(grp_id,"scan_end_time",&sds_id);
//    nc_get_vara_double(grp_id, sds_id, start, count, scan_stop_tai);

    stoptime = scan_start_tai[nscan];

    time_interval = ( stoptime - starttime ) / ( nscan - 1 ); /* in sec */

/*
    retval = navigation_meris(fileID);
*/
    return (LIFE_IS_GOOD);
}

int readl1_orca(filehandle *file, int32 scan, l1str *l1rec)
{
    static int firstCall = 1;
    int err_code, retval;
    float *lon, *lat, *senz, *sena, *solz, *sola;

    int32_t ip, ib, ipb, Ibp;
    int32_t nbands = l1rec->nbands;
    epr_uint flag;
    size_t start[3], count[3];
    int i;
    int xid, band_id, sds_id, grp_id;
    float *data_ctl1, *data_ctl2, *ctlfpix;
    float *rad_data;
    static float gain[15];
    double tai2unixoffset =  725846400.0;
    double msec;// recsec, sec70;
    int16 scan_year, scan_day;

    if (firstCall) {
        if (want_verbose)
            printf("file->nbands = %d, l1rec->nbands = %d\n",
                    (int) file->nbands, (int) l1rec->nbands);
        firstCall = 0;

        for (ip = 0; ip < npix; ip++) {
            l1rec->pixnum[ip] = ip;
            l1rec->flags[ip] = 0;
        }

    }


//      set time for this scan - if scan_start_time value not properly set, estimate from scene start time.
  // printf("RJH: scan_start_tai=%f\n",scan_start_tai[scan]);
    if (scan_start_tai[scan] > 0){
        lastvalidtime = scan_start_tai[scan]+tai2unixoffset;
        lastvalidscan = scan;
        unix2yds(lastvalidtime, &scan_year, &scan_day, &msec);
    } else {
        unix2yds(lastvalidtime+(time_interval * (scan-lastvalidscan)),&scan_year,&scan_day, &msec);
    }
    *(l1rec->year)  = (int32_t) scan_year;
    *(l1rec->day)   = (int32_t) scan_day;
    *(l1rec->msec)  = (int32_t) msec*1000;
    *(l1rec->msec) -= leapseconds_since_1993(lastvalidtime-tai2unixoffset)*1000;

    retval = nc_inq_ncid(file->sd_id,"navigation_data" , &grp_id);
    if (retval == FAIL) {
        fprintf(stderr,
                "-E- %s line %d: nc_inq_ncid failed for file, %s  group, %s.\n",
                __FILE__, __LINE__, file->name, "navigation_data");
        return (1);
    }

    start[0] = scan;
    start[1] = 0;
    start[2] = 0;
    count[0] = 1;
    count[1] = npix;
    count[2] = 0;

    lon  = (float *) calloc(npix, sizeof(float));
    lat  = (float *) calloc(npix, sizeof(float));
    senz = (float *) calloc(npix, sizeof(float));
    sena = (float *) calloc(npix, sizeof(float));
    solz = (float *) calloc(npix, sizeof(float));
    sola = (float *) calloc(npix, sizeof(float));

    retval = nc_inq_varid(grp_id, "lon", &band_id);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, file->name, "lon");
        return (1);
    }
    retval = nc_get_vara_float(grp_id, band_id, start, count,lon);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, file->name, "lon");
        return (1);
    }
    retval = nc_inq_varid(grp_id, "lat", &band_id);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, file->name, "lon");
        return (1);
    }
    retval = nc_get_vara_float(grp_id, band_id, start, count,lat);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, file->name, "lat");
        return (1);
    }
    //until geolocation is read, set fill values -
    for (ip = spix; ip < npix; ip++) {
        l1rec->lon[ip] = lon[ip];
        l1rec->lat[ip] = lat[ip];
        l1rec->solz[ip] = -999;//solz[scan * npix + ip];
        l1rec->sola[ip] = -999;//sola[scan * npix + ip];
        l1rec->senz[ip] = -999;//senz[scan * npix + ip];
        l1rec->sena[ip] = -999;//sena[scan * npix + ip];
    }

    retval = nc_inq_varid(grp_id, "sena", &band_id);
    retval = nc_get_vara_float(grp_id, band_id, start, count,sena);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, file->name, "sena");
        return (1);
    }
    retval = nc_inq_varid(grp_id, "senz", &band_id);
    retval = nc_get_vara_float(grp_id, band_id, start, count,senz);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, file->name, "senz");
        return (1);
    }
    retval = nc_inq_varid(grp_id, "sola", &band_id);
    retval = nc_get_vara_float(grp_id, band_id, start, count,sola);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, file->name, "sola");
        return (1);
    }
    retval = nc_inq_varid(grp_id, "solz", &band_id);
    retval = nc_get_vara_float(grp_id, band_id, start, count,solz);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, file->name, "solz");
        return (1);
    }

    for (ip = spix; ip < npix; ip++) {
      l1rec->lon[ip] = lon[ip];
      l1rec->lat[ip] = lat[ip];
      l1rec->solz[ip] = solz[ip];
      l1rec->sola[ip] = sola[ip];
      l1rec->senz[ip] = senz[ip];
      l1rec->sena[ip] = sena[ip];
      //printf("RJH: ip=%d lon=%lf lat=%lf solz=%lf sola=%lf senz=%lf sena=%lf\n",ip,lon[ip],lat[ip], solz[ip],sola[ip],senz[ip],sena[ip] );
    }

    // read in radiance data
    rad_data = (float *) calloc(npix*nbands, sizeof(float));

    start[0] = 0;
    start[1] = scan;
    start[2] = 0;
    count[0] = nbands;
    count[1] = 1;
    count[2] = npix;

    nc_inq_ncid(file->sd_id,"earth_view_data" , &grp_id);

    retval = nc_inq_varid(grp_id, "Lt_visnir", &band_id);
    retval = nc_get_vara_float(grp_id, band_id, start, count,rad_data);
    if (retval == FAIL) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, file->name, "Lt_visnir");
        return (1);
    }

    for (ib = 0; ib < nbands; ib++) {

        // copy to Lt record.
        for (ip = spix; ip < npix; ip++) {
            ipb = ip * nbands + ib;
            Ibp = ib * npix + ip;
            l1rec->Lt[ipb] = rad_data[Ibp];
            l1rec->detnum = ip;

            // mark negative input data as HILT
            // navfail commented out for Lt < 0 - too limiting for hyperspectral - RJH
            if (l1rec->Lt[ipb] < 0.0)
                l1rec->Lt[ipb] = 0.0001;
//                l1rec->navfail[ip] = 1;
        }
    } // for ib
    free(rad_data);
    free(lat);
    free(lon);
    free(solz);
    free(sola);
    free(senz);
    free(sena);

    l1rec->sensorID = file->sensorID;
    l1rec->npix = file->npix;

    return (LIFE_IS_GOOD);
}

int closel1_orca(filehandle *file) {
    int retval;

    retval = nc_close(file->sd_id);
    if (retval == FAIL) {
        fprintf(stderr, "-E- %s line %d: nc_close failed for file, %s.\n",
        __FILE__, __LINE__, file->name);
        return (1);
    }
    free(scan_start_tai);

    return (LIFE_IS_GOOD);
}

