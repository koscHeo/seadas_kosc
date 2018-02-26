/* ============================================================================ */
/* module l1_olci.c - functions to read OLCI (coastal color) for MSL12          */
/* Written By:  Richard Healy (SAIC) July 29, 2015.                             */
/*                                                                              */
/* ============================================================================ */

// #include <stdbool.h>
#include "l12_proto.h"
#include "olci.h"
#include "math.h"

#define NUMRADFILES 21
#define TIE_ROW_PTS      1
#define TIE_COL_PTS     64
#define NEW_CACHE_SIZE 32000000
#define NEW_CACHE_NELEMS 1217
#define NEW_CACHE_PREEMPTION .75

static int16 nline, npix;

static int32 spix = 0;
static int64_t *scan_start_tai;
static double lastvalidtime;
static int lastvalidscan = 0;
static int32_t olci_sd[MAXOLCI_RADFILES],geoFileID,coordFileID,tcoordFileID,instrumentFileID;
//static int bindx[15] = { 1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 15, 16, 17, 18 };
static double pi = PI;

double deg2rad(double);
double angular_distance(double lat1, double lon1, double lat2, double lon2);

olci_t* createPrivateData_olci (int numBands) {

    olci_t* data = (olci_t*)calloc(1, sizeof(olci_t));
    if(data == NULL) {
        fprintf(stderr,"-E- %s line %d: unable to allocate private data for olci\n",
                __FILE__,__LINE__);
        exit(1);
    }

    data->wave = (double *) malloc(numBands*sizeof(double) );
    data->fwhm = (double *) malloc(numBands*sizeof(double) );
    if(data->wave==NULL || data->fwhm==NULL) {
        fprintf(stderr,"-E- %s line %d: unable to allocate scale/offset data for olci\n",
                __FILE__,__LINE__);
        exit(1);
    }

    data->numRadFiles = NUMRADFILES;
    return data;
}
/*
 *  W. Robinson, SAIC, 18 Feb 2016  create correct name for radiance dataset 
 *  for I/O
 *
 *  Put olci_varname into olci struct instead of using static global. - rjh 2/29/2016
 */
int openl1_olci(filehandle * l1file) {
    char *fltime;
    olci_t *data = l1file->private_data;
    size_t source_w, source_h, source_b;
    int32 nscan;
    int timefileID, xid, yid, bid,geo_id, retval, grp_id, sds_id;
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
    unsigned short orbit_number;

    if (!data) {
    	data = l1file->private_data = createPrivateData_olci(MAXOLCI_RADFILES);
        if ((data->instrumentFile = (char *) malloc((sizeof("instrument_data.nc"))*sizeof(char))) == NULL) {
            printf("%s, %d - E - unable to allocate instrument data filename \n",
                    __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        if ((data->time_coordinatesFile = (char *) malloc((sizeof("time_coordinates.nc"))*sizeof(char))) == NULL) {
            printf("%s, %d - E - unable to allocate time coordinates filename \n",
                    __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        if ((data->geoCoordinatesFile = (char *) malloc((sizeof("geo_coordinates.nc"))*sizeof(char))) == NULL) {
            printf("%s, %d - E - unable to allocate geo coordinates filename \n",
                    __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        if ((data->tieGeometriesFile = (char *) malloc((sizeof("tie_geometries.nc"))*sizeof(char))) == NULL) {
            printf("%s, %d - E - unable to allocate tie geometries filename \n",
                    __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        if ((data->tieGeoCoordinatesFile = (char *) malloc((sizeof("tie_geo_coordinates.nc"))*sizeof(char))) == NULL) {
            printf("%s, %d - E - unable to allocate tie geo coordinates filename \n",
                    __FILE__, __LINE__);
            exit(EXIT_FAILURE);
        }
        sprintf(data->instrumentFile,"%s","instrument_data.nc");
        sprintf(data->time_coordinatesFile,"%s","time_coordinates.nc");
        sprintf(data->tieGeometriesFile,"%s","tie_geometries.nc");
        sprintf(data->geoCoordinatesFile,"%s","geo_coordinates.nc");
        sprintf(data->tieGeoCoordinatesFile,"%s","tie_geo_coordinates.nc");
        for (i=0;i<data->numRadFiles;i++){
            if ((data->olci_radfiles[i] = (char *) malloc((sizeof("Oa01_radiance.nc"))*sizeof(char))) == NULL) {
                printf("%s, %d - E - unable to allocate radiance filename \n",
                        __FILE__, __LINE__);
                exit(EXIT_FAILURE);
            }
           sprintf(data->olci_radfiles[i],"Oa%2.2d_radiance.nc",i+1);
           //printf("olci rad=%s\n",data->olci_radfiles[i]);
        }
    } else {
    }
    // Open the netcdf4 input file
    if (nc_set_chunk_cache(NEW_CACHE_SIZE, NEW_CACHE_NELEMS,
              NEW_CACHE_PREEMPTION)) {
        fprintf(stderr, "-E- %s line %d: nc_set_chunk_cache (%s) failed.\n",
        __FILE__, __LINE__, l1file->name);
        return (1);
    }
    retval = nc_open(data->instrumentFile, NC_NOWRITE, &instrumentFileID);
    if (retval == FAIL) {
        fprintf(stderr, "-E- %s line %d: nc_open(%s) failed.\n",
        __FILE__, __LINE__, l1file->name);
        return (1);
    }

    // Get pixel and scan dimensions
    retval = nc_inq_dimid(instrumentFileID, "rows", &yid);
    retval = nc_inq_dimid(instrumentFileID, "columns", &xid);
    retval = nc_inq_dimid(instrumentFileID, "bands", &bid);
    retval = nc_inq_dimlen(instrumentFileID, xid, &source_w);
    retval = nc_inq_dimlen(instrumentFileID, yid, &source_h);
    retval = nc_inq_dimlen(instrumentFileID, bid, &source_b);
    retval = nc_get_att_ushort(instrumentFileID, NC_GLOBAL, "absolute_orbit_number", &orbit_number);

//    source_b = 15; // Temporary value for use with MERIS default files/data;
//    printf("RJH: Using Temporary nbands=%d for use with MERIS data copied into run/data/olci", (int)source_b);

    if (want_verbose) {
        printf("OLCI L1B Npix  :%d Nscans:%d\n", (int) source_w,
                (int) source_h);
    } // want_verbose
    npix = (int32) source_w;
    nscan = (int32) source_h;
    nline = nscan;

    l1file->orbit_number = orbit_number;
    l1file->nbands = source_b;
    l1file->npix = npix;
    l1file->nscan = nscan;
    strcpy(l1file->spatialResolution, "750 m");
    data->numBands = source_b;

    scan_start_tai = (int64_t *) calloc(nscan,sizeof(int64_t));

    start[0] = 0;
    count[0] = nscan;

    retval = nc_open(data->time_coordinatesFile, NC_NOWRITE, &timefileID);
    if (retval == FAIL) {
        fprintf(stderr, "-E- %s line %d: nc_open(%s) failed.\n",
        __FILE__, __LINE__, data->time_coordinatesFile);
        return (1);
    }

    nc_inq_varid(timefileID,"time_stamp",&sds_id); //microseconds since 1/1/2000
    nc_get_vara_long(timefileID, sds_id, start, count, scan_start_tai);

    retval = nc_open(data->tieGeometriesFile, NC_NOWRITE, &geoFileID);
    if (retval == FAIL) {
        fprintf(stderr, "-E- %s line %d: nc_open(%s) failed.\n",
        __FILE__, __LINE__, data->tieGeometriesFile);
        return (1);
    }

    retval = nc_open(data->geoCoordinatesFile, NC_NOWRITE, &coordFileID);
    if (retval == FAIL) {
        fprintf(stderr, "-E- %s line %d: nc_open(%s) failed.\n",
        __FILE__, __LINE__, data->geoCoordinatesFile);
        return (1);
    }
    retval = nc_open(data->tieGeoCoordinatesFile, NC_NOWRITE, &tcoordFileID);
    if (retval == FAIL) {
        fprintf(stderr, "-E- %s line %d: nc_open(%s) failed.\n",
        __FILE__, __LINE__, data->tieGeoCoordinatesFile);
        return (1);
    }

    for (i=0; i<data->numRadFiles; i++) {
          // Open each of the netcdf4 Lt files for each band
        retval = nc_open(data->olci_radfiles[i], NC_NOWRITE, &olci_sd[i]);
        if (retval == FAIL) {

              fprintf(stderr,
                      "-E- %s line %d: nc_open failed for file, %s.\n",
                      __FILE__, __LINE__, data->olci_radfiles[i]);
              return (1);
          }
          // NOT so good  strcpy(olci_varname[i],data->olci_radfiles[i]);
          sprintf( data->olci_varname[i], "Oa%02d_radiance", i+1 );
    }


//    for (i=0; i<l1file->nbands;i++)
//        printf("Band % 2d: %s\n", i, olci_varname[bindx[i]]);
//          printf("Band % 2d: %s\n", i, data->olci_varname[i]);
    return (LIFE_IS_GOOD);
}

int readl1_olci(filehandle *file, int32 scan, l1str *l1rec)
{
    static int firstCall = 1;
    static float64 utime;
    static double time_interval;
    static int64_t tai2unixoffset; //=  946684800;
    static float  *senz, *sena, *solz, *sola;
    static uint *tsolz,*tsenz;
    static int32_t *lon, *lat, *tsola, *tsena,  *tlon, *tlat;
    static double lon0,lat0,xdist, *dist,*sdist,*tielat,*tielon,*tiesolz,*tiesenz,*tiesola,*tiesena;
    static float scale_lon,scale_lat,scale_solz,scale_senz,scale_sola,scale_sena,scale_rad;
    static int16_t *detector_index;
    static float *lambda0,*fwhm;
    static size_t detectors;
    static unsigned int *rad_data;
    static size_t nctlpix,nctlscan;
    static unsigned short tie_row_pts, tie_col_pts;

    int err_code, retval;

    int32_t ip, ib, ipb, Ibp;
    int32_t nbands = l1rec->nbands;
    epr_uint flag;
    size_t start[3], count[3];
    int i;

    int xid, yid, band_id, sds_id;
    double msec;// recsec, sec70;
    int16 scan_year, scan_day;

    gsl_spline *spline[4];
    gsl_interp_accel *spl_acc;

    char varnam[14];
    olci_t *data = file->private_data;

    if (firstCall) {
        firstCall = 0;
        if (want_verbose)
            printf("file->nbands = %d, l1rec->nbands = %d\n",
                    (int) file->nbands, (int) l1rec->nbands);

        for (ip = 0; ip < npix; ip++) {
            l1rec->pixnum[ip] = ip;
            l1rec->flags[ip] = 0;
        }
        utime = yds2unix(2000,1,0);
        tai2unixoffset = utime; // This is the time offset between what the time utilities need (1/1/1970) and what
                                // the OLCI netcdf provides (1/1/2000)

        retval = nc_inq_dimid(tcoordFileID, "tie_columns", &xid);
        retval = nc_inq_dimlen(tcoordFileID, xid, &nctlpix);
        retval = nc_inq_dimid(tcoordFileID, "tie_rows", &yid);
        retval = nc_inq_dimlen(tcoordFileID, yid, &nctlscan);
        retval = nc_get_att_ushort(tcoordFileID, NC_GLOBAL, "ac_subsampling_factor", &tie_col_pts);
        if (retval) tie_col_pts = TIE_COL_PTS;
        retval = nc_get_att_ushort(tcoordFileID, NC_GLOBAL, "al_subsampling_factor", &tie_row_pts);
        if (retval) tie_row_pts = TIE_ROW_PTS;

        if (tie_row_pts != 1) {
            printf("-E- %s line %d: Sorry.  I can only handle tie_row_pts = 1 not tie_row_pts = %d\n",
                            __FILE__, __LINE__,tie_row_pts);
            exit(-1);
        }

        if ( ((nctlscan-1)*tie_row_pts+1) != file->nscan) {
            printf("-E- %s line %d: Sanity check failed - tie_rows (%d) x tie_row_pts (%d) != nscan (%d) in file %s\n",
                            __FILE__, __LINE__,(int)nctlscan,tie_row_pts, file->nscan,data->tieGeoCoordinatesFile);
            exit(-1);
        }
        if (((nctlpix-1)*tie_col_pts+1) != file->npix ) {
            printf("-E- %s line %d: Sanity check failed - tie_cols (%d) x tie_col_pts (%d) != npix (%d) in file %s\n",
                            __FILE__, __LINE__, (int)nctlpix,tie_col_pts, file->nscan,data->tieGeoCoordinatesFile);
            exit(-1);
        }


        tlon  = ( int32_t *) calloc(nctlpix, sizeof(int32_t));
        tlat  = ( int32_t *) calloc(nctlpix, sizeof(int32_t));
        tsenz = (uint32_t *) calloc(nctlpix, sizeof(uint32_t));
        tsena = ( int32_t *) calloc(nctlpix, sizeof(int32_t));
        tsolz = (uint32_t *) calloc(nctlpix, sizeof(uint32_t));
        tsola = ( int32_t *) calloc(nctlpix, sizeof(int32_t));

        tielon  = (double *) calloc(nctlpix, sizeof(double));
        tielat  = (double *) calloc(nctlpix, sizeof(double));
        tiesenz = (double *) calloc(nctlpix, sizeof(double));
        tiesena = (double *) calloc(nctlpix, sizeof(double));
        tiesolz = (double *) calloc(nctlpix, sizeof(double));
        tiesola = (double *) calloc(nctlpix, sizeof(double));
        dist    = (double *) calloc(nctlpix, sizeof(double));
        sdist   = (double *) calloc(nctlpix, sizeof(double));

        lon  = (int32_t *) calloc(npix, sizeof(int32_t));
        lat  = (int32_t *) calloc(npix, sizeof(int32_t));

        retval = nc_inq_dimid(instrumentFileID, "detectors", &xid);
        retval = nc_inq_dimlen(instrumentFileID, xid, &detectors);
        //printf("RJH: instrumentFile: detectors=%d \n",(int)detectors);

        detector_index = (int16_t*) calloc(npix*nline, sizeof(int16_t));
        lambda0 = (float*) calloc(nbands*detectors, sizeof(float));
        fwhm = (float*) calloc(nbands*detectors, sizeof(float));

        start[0] = 0;
        start[1] = 0;
        start[2] = 0;
        count[0] = nline;
        count[1] = npix;
        count[2] = 0;

        retval = nc_inq_varid(instrumentFileID, "detector_index", &xid);
        if (retval != NC_NOERR) {
            fprintf(stderr,
                    "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                    __FILE__, __LINE__, data->instrumentFile, "detector_index");
            exit (FATAL_ERROR);
        }
        retval = nc_get_vara_short(instrumentFileID, xid, start, count,detector_index);
        if (retval != NC_NOERR) {
            fprintf(stderr,
                    "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                    __FILE__, __LINE__, data->instrumentFile, "detector_index");
            exit (FATAL_ERROR);
        }

        start[0] = 0;
        start[1] = 0;
        start[2] = 0;
        count[0] = nbands;
        count[1] = detectors;
        count[2] = 0;

        retval = nc_inq_varid(instrumentFileID, "lambda0", &xid);
        if (retval != NC_NOERR) {
            fprintf(stderr,
                    "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                    __FILE__, __LINE__, data->instrumentFile, "lambda0");
            exit (FATAL_ERROR);
    }
        retval = nc_get_vara_float(instrumentFileID, xid, start, count,lambda0);
        if (retval != NC_NOERR) {
            fprintf(stderr,
                    "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                    __FILE__, __LINE__, data->instrumentFile, "lambda0");
            exit (FATAL_ERROR);
        }
        retval = nc_inq_varid(instrumentFileID, "FWHM", &xid);
        if (retval != NC_NOERR) {
            fprintf(stderr,
                    "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                    __FILE__, __LINE__, data->instrumentFile, "FWHM");
            exit (FATAL_ERROR);
        }
        retval = nc_get_vara_float(instrumentFileID, xid, start, count,fwhm);
        if (retval != NC_NOERR) {
            fprintf(stderr,
                    "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                    __FILE__, __LINE__, data->instrumentFile, "FWHM");
            exit (FATAL_ERROR);
        }

        rad_data = (unsigned int *) calloc(npix, sizeof(rad_data)); //BYSCAN

    }


//      set time for this scan - if scan_start_time value not properly set, estimate from scene start time.
    if (scan_start_tai[scan] > 0){
        lastvalidtime = scan_start_tai[scan]/1000000 + tai2unixoffset;
        if (scan > 0) time_interval = (scan_start_tai[scan] - lastvalidtime)/(scan-lastvalidscan);
        lastvalidscan = scan;
        unix2yds(lastvalidtime, &scan_year, &scan_day, &msec);
    } else {
        unix2yds(lastvalidtime+(time_interval * (scan-lastvalidscan)),&scan_year,&scan_day, &msec);
    }
    *(l1rec->year)  = (int32_t) scan_year;
    *(l1rec->day)   = (int32_t) scan_day;
    *(l1rec->msec)  = (int32_t) msec*1000;
 //   *(l1rec->msec) -= leapseconds_since_1993(lastvalidtime-tai2unixoffset)*1000;

    //printf("RJH: year=%d day=%d %lf %ld\n",scan_year, scan_day, msec,scan_start_tai[scan]);

    start[0] = scan;
    start[1] = 0;
    start[2] = 0;
    count[0] = 1;
    count[1] = nctlpix;
    count[2] = 0;

    retval = nc_inq_varid(tcoordFileID, "longitude", &xid);
    retval = nc_get_att_float(tcoordFileID, xid, "scale_factor", &scale_lon);
    retval = nc_get_vara_int(tcoordFileID, xid, start, count,tlon);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, data->tieGeoCoordinatesFile, "lon");
        exit (FATAL_ERROR);
    }
    retval = nc_inq_varid(tcoordFileID, "latitude", &yid);
    retval = nc_get_att_float(tcoordFileID, yid, "scale_factor", &scale_lat);
    retval = nc_get_vara_int(tcoordFileID, yid, start, count,tlat);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, data->tieGeoCoordinatesFile, "lat");
        exit (FATAL_ERROR);
    }

    retval = nc_inq_varid(geoFileID, "OAA", &band_id);
    retval = nc_get_att_float(geoFileID, band_id, "scale_factor", &scale_sena);
    retval = nc_get_vara_int(geoFileID, band_id, start, count,tsena);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, data->tieGeoCoordinatesFile, "sena");
        exit (FATAL_ERROR);
    }
    retval = nc_inq_varid(geoFileID, "OZA", &band_id);
    retval = nc_get_att_float(geoFileID, band_id, "scale_factor", &scale_senz);
    retval = nc_get_vara_uint(geoFileID, band_id, start, count,tsenz);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, data->tieGeoCoordinatesFile, "senz");
        exit (FATAL_ERROR);
    }
    retval = nc_inq_varid(geoFileID, "SAA", &band_id);
    retval = nc_get_att_float(geoFileID, band_id, "scale_factor", &scale_sola);
    retval = nc_get_vara_int(geoFileID, band_id, start, count,tsola);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, data->tieGeoCoordinatesFile, "sola");
        exit (FATAL_ERROR);
    }
    retval = nc_inq_varid(geoFileID, "SZA", &band_id);
    retval = nc_get_att_float(geoFileID, band_id, "scale_factor", &scale_solz);
    retval = nc_get_vara_uint(geoFileID, band_id, start, count,tsolz);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, data->tieGeoCoordinatesFile, "solz");
        exit (FATAL_ERROR);
    }

    for (ip = 0; ip < nctlpix; ip++) {
        tielat[ip]  = scale_lat*tlat[ip];
        tielon[ip]  = scale_lon*tlon[ip];
        tiesena[ip] = scale_sena*tsena[ip];
        tiesenz[ip] = scale_senz*tsenz[ip];
        tiesola[ip] = scale_sola*tsola[ip];
        tiesolz[ip] = scale_solz*tsolz[ip];

        dist[ip]    = angular_distance(tielat[ip],tielon[ip],tielat[0],tielon[0]);
        //printf("RJH: TIE: ip=%d dist=%f lon=%lf lat=%lf %lf %lf solz=%lf sola=%lf senz=%lf sena=%lf\n",ip,dist[ip],tielon[ip],tielat[ip],tielon[0],tielat[0], tiesolz[ip],tiesola[ip],tiesenz[ip],tiesena[ip] );

    }

    start[0] = scan;
    start[1] = 0;
    start[2] = 0;
    count[0] = 1;
    count[1] = npix;
    count[2] = 0;

    //until geolocation is read, set fill values -
    for (ip = spix; ip < npix; ip++) {
        l1rec->lon[ip] = -999;
        l1rec->lat[ip] = -999;

        l1rec->solz[ip] = -999;//solz[scan * npix + ip];
        l1rec->sola[ip] = -999;//sola[scan * npix + ip];
        l1rec->senz[ip] = -999;//senz[scan * npix + ip];
        l1rec->sena[ip] = -999;//sena[scan * npix + ip];
    }

    retval = nc_inq_varid(coordFileID, "longitude", &xid);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, data->geoCoordinatesFile, "lon");
        exit (FATAL_ERROR);
    }
    retval = nc_get_vara_int(coordFileID, xid, start, count,lon);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, data->geoCoordinatesFile, "lat");
        exit (FATAL_ERROR);
    }
    retval = nc_inq_varid(coordFileID, "latitude", &yid);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, data->geoCoordinatesFile, "lat");
        exit (FATAL_ERROR);
    }
    retval = nc_get_vara_int(coordFileID, yid, start, count,lat);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, data->geoCoordinatesFile, "lat");
        exit (FATAL_ERROR);
    }


     for (i=0; i<4; i++)
       spline[i] = gsl_spline_alloc(gsl_interp_cspline, nctlpix);

     spl_acc = gsl_interp_accel_alloc();

     // Sort distances and corresponding elevation values
     memcpy(sdist,dist,sizeof(double)*nctlpix);
     gsl_sort2( sdist, 1, tiesolz, 1, nctlpix);
     // Initiate spline
     gsl_spline_init(spline[0], sdist, tiesolz, nctlpix);

     memcpy(sdist,dist,sizeof(double)*nctlpix);
     gsl_sort2( sdist, 1, tiesenz, 1, nctlpix);
     gsl_spline_init(spline[1], sdist, tiesenz, nctlpix);

     memcpy(sdist,dist,sizeof(double)*nctlpix);
     gsl_sort2( sdist, 1, tiesola, 1, nctlpix);
     gsl_spline_init(spline[2], sdist, tiesola, nctlpix);

     memcpy(sdist,dist,sizeof(double)*nctlpix);
     gsl_sort2( sdist, 1, tiesena, 1, nctlpix);
     gsl_spline_init(spline[3], sdist, tiesena, nctlpix);

    for (ip = spix; ip < npix; ip++) {

      l1rec->lon[ip] = lon[ip]*scale_lon;
      l1rec->lat[ip] = lat[ip]*scale_lat;

      xdist    = angular_distance(l1rec->lat[ip],l1rec->lon[ip],tielat[0],tielon[0]);

      //printf("RJH:XDIST: ip=%d dist=%f lon=%lf lat=%lf %d %d \n",ip,xdist,l1rec->lon[ip],l1rec->lat[ip],lon[ip],lat[ip]);
      if (xdist <= dist[nctlpix-1]) {
          l1rec->solz[ip] = (float)gsl_spline_eval( spline[0], xdist, spl_acc);
          l1rec->senz[ip] = (float)gsl_spline_eval( spline[1], xdist, spl_acc);
          l1rec->sola[ip] = (float)gsl_spline_eval( spline[2], xdist, spl_acc);
          l1rec->sena[ip] = (float)gsl_spline_eval( spline[3], xdist, spl_acc);
          //printf("solz=%lf sola=%lf senz=%lf sena=%lf\n", l1rec->solz[ip],l1rec->sola[ip],l1rec->senz[ip],l1rec->sena[ip] );
      }
    }

    // read in radiance data

    for (ib = 0; ib < nbands; ib++) {
//        retval = nc_inq_varid(olci_sd[ib], olci_varname[bindx[ib]], &band_id);
        retval = nc_inq_varid(olci_sd[ib], data->olci_varname[ib], &band_id);
        if (retval < 0) {
          fprintf(stderr,
                    "-E- %s line %d: nc_inq_varid failed for file, %s  field, %s.\n",
                    __FILE__, __LINE__, data->olci_radfiles[ib], varnam);
            exit (FATAL_ERROR);
        }
        retval = nc_get_att_float(olci_sd[ib], band_id, "scale_factor", &scale_rad);
        if (retval == FAIL) {
            fprintf(stderr,
                    "-E- %s line %d: nc_get_att_float failed for file, %s  field, %s.\n",
                    __FILE__, __LINE__, data->olci_radfiles[ib], varnam);
            exit (FATAL_ERROR);
        }
        retval = nc_get_vara_uint(olci_sd[ib], band_id, start, count,rad_data); //BYSCAN
        if (retval == FAIL) {
            fprintf(stderr,
                    "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                    __FILE__, __LINE__, data->olci_radfiles[ib], varnam);
            exit (FATAL_ERROR);
        }
        // copy to Lt record.
        for (ip = spix; ip < npix; ip++) {
            ipb = ip * nbands + ib;
            l1rec->Lt[ipb] = rad_data[ip]*scale_rad/10.; //BYSCAN
            l1rec->detnum = ip;

            // mark negative input data as HILT
            // navfail commented out for Lt < 0 - too limiting for hyperspectral - RJH
            if (l1rec->Lt[ipb] < 0.0)
                l1rec->Lt[ipb] = 0.0001;
//                l1rec->navfail[ip] = 1;
        }
    } // for ib


//    for (ip=0;ip<npix;ip++) {
//        for (ib=0; ib<nbands; ib++){
//            if (detector_index[ip]>=0) {
//               ipb = ib*detectors+detector_index[ip];
//               l1rec->fwave[ib] = *(lambda0+ipb);
//               printf("%d %f/%f/%d \n",ib,*(lambda0+ipb),*(fwhm+ipb),detector_index[ip]);
//            }
//        }
//    }

    gsl_interp_accel_free(spl_acc);

    l1rec->sensorID = file->sensorID;
    l1rec->npix = file->npix;

    return (LIFE_IS_GOOD);
}

int closel1_olci(filehandle *file) {
    int retval;
    int i;
    olci_t *data = file->private_data;

    retval = nc_close(geoFileID);
    retval += nc_close(coordFileID);
    retval += nc_close(tcoordFileID);
    retval += nc_close(instrumentFileID);
    for (i=0; i< data->numRadFiles; i++)
        retval += nc_close(olci_sd[i]);

    if (retval != 0) {
        fprintf(stderr, "-E- %s line %d: nc_close failed for one or more OLCI files.\n",
        __FILE__, __LINE__);
        return (1);
    }
    free(scan_start_tai);

    return (LIFE_IS_GOOD);
}

/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
/*::  Compute the angular distance                                  :*/
/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/

double angular_distance(double lat1, double lon1, double lat2, double lon2) {
  double theta, dist;
  theta = lon1 - lon2;
  dist = sin(deg2rad(lat1)) * sin(deg2rad(lat2)) + cos(deg2rad(lat1)) * cos(deg2rad(lat2)) * cos(deg2rad(theta));
  dist = dist> 1? 1:dist;
  dist = dist<-1?-1:dist;
  return (acos(dist));
}

/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
/*::  This function converts decimal degrees to radians             :*/
/*:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
double deg2rad(double deg) {
  return (deg * pi / 180);
}

