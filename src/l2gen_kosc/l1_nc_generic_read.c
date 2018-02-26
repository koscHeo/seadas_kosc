#include "l12_proto.h"
#include "l1_nc_generic_read.h"

static int ncol;
static int16 nline, npix;
static int32 spix = 0;

int openl1_nc_generic(filehandle *file) {
    size_t source_w, source_h;
    int32 nscan;
    int fileID, xid, yid, retval, grp_id, sds_id;
    int i;
    size_t start[3], count[3];

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
        printf("L1B Npix  :%d Nscans:%d\n", (int) source_w,
                (int) source_h);
    } // want_verbose
    npix = (int32) source_w;
    nscan = (int32) source_h;
    nline = nscan;

    file->sd_id = fileID;
    file->npix = npix;
    file->nscan = nscan;
    file->terrain_corrected = 1; //assumed
    return (LIFE_IS_GOOD);
}

int readl1_nc_generic(filehandle *file, int32 scan, l1str *l1rec)
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
    float *rad_data;
    int msec;
    int scan_year, scan_day;
    static float *scale_factor, *add_offset;
    
    lon  = (float *) calloc(npix, sizeof(float));
    lat  = (float *) calloc(npix, sizeof(float));
    senz = (float *) calloc(npix, sizeof(float));
    sena = (float *) calloc(npix, sizeof(float));
    solz = (float *) calloc(npix, sizeof(float));
    sola = (float *) calloc(npix, sizeof(float));

    
    if (firstCall) {
        if (want_verbose)
            printf("file->nbands = %d, l1rec->nbands = %d\n",
                    (int) file->nbands, (int) l1rec->nbands);
        firstCall = 0;

        for (ip = 0; ip < npix; ip++) {
            l1rec->pixnum[ip] = ip;
            l1rec->flags[ip] = 0;
        }
        retval = nc_inq_ncid(file->sd_id,"navigation_data" , &grp_id);
        if (retval == FAIL) {
            fprintf(stderr,
                "-E- %s line %d: nc_inq_ncid failed for file, %s  group, %s.\n",
                __FILE__, __LINE__, file->name, "navigation_data");
            return (1);
        }
        scale_factor = (float *) calloc(4, sizeof(float));
        add_offset = (float *) calloc(4, sizeof(float));
        
        retval = nc_inq_varid(grp_id, "sena", &band_id);
        nc_get_att_float(grp_id, band_id, "scale_factor", &scale_factor[0]);
        nc_get_att_float(grp_id, band_id, "add_offset", &add_offset[0]);
        
        retval = nc_inq_varid(grp_id, "senz", &band_id);
        nc_get_att_float(grp_id, band_id, "scale_factor", &scale_factor[1]);
        nc_get_att_float(grp_id, band_id, "add_offset", &add_offset[1]);
        
        retval = nc_inq_varid(grp_id, "sola", &band_id);
        nc_get_att_float(grp_id, band_id, "scale_factor", &scale_factor[2]);
        nc_get_att_float(grp_id, band_id, "add_offset", &add_offset[2]); 
        
        retval = nc_inq_varid(grp_id, "solz", &band_id);
        nc_get_att_float(grp_id, band_id, "scale_factor", &scale_factor[3]);
        nc_get_att_float(grp_id, band_id, "add_offset", &add_offset[3]); 
    }

    start[0] = scan;
    count[0] = 1;

    nc_inq_ncid(file->sd_id,"scan_line_attributes" , &grp_id);

    nc_inq_varid(grp_id,"year",&sds_id);
    nc_get_vara_int(grp_id, sds_id, start, count, &scan_year);
    nc_inq_varid(grp_id,"day",&sds_id);
    nc_get_vara_int(grp_id, sds_id, start, count, &scan_day);
    nc_inq_varid(grp_id,"msec",&sds_id);
    nc_get_vara_int(grp_id, sds_id, start, count, &msec);
    
//    l1rec->scantime = yds2unix((int16) scan_year, (int16) scan_day, (double)(msec/1.e3));
    *l1rec->year = scan_year;
    *l1rec->day  = scan_day;
    *l1rec->msec = msec;
    
    start[0] = scan;
    start[1] = 0;
    start[2] = 0;
    count[0] = 1;
    count[1] = npix;
    count[2] = 1;
    
    retval = nc_inq_ncid(file->sd_id,"navigation_data" , &grp_id);
    if (retval == FAIL) {
        fprintf(stderr,
            "-E- %s line %d: nc_inq_ncid failed for file, %s  group, %s.\n",
            __FILE__, __LINE__, file->name, "navigation_data");
        return (1);
    }
    retval = nc_inq_varid(grp_id, "longitude", &band_id);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, file->name, "longitude");
        return (1);
    }
    retval = nc_get_vara_float(grp_id, band_id, start, count,lon);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, file->name, "longitude");
        return (1);
    }
    retval = nc_inq_varid(grp_id, "latitude", &band_id);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, file->name, "latitude");
        return (1);
    }
    retval = nc_get_vara_float(grp_id, band_id, start, count,lat);
    if (retval != NC_NOERR) {
        fprintf(stderr,
                "-E- %s line %d: nc_get_vara_float failed for file, %s  field, %s.\n",
                __FILE__, __LINE__, file->name, "latitude");
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
      l1rec->sena[ip] = sena[ip] * scale_factor[0] + add_offset[0];
      l1rec->senz[ip] = senz[ip] * scale_factor[1] + add_offset[1];
      l1rec->sola[ip] = sola[ip] * scale_factor[2] + add_offset[2];
      l1rec->solz[ip] = solz[ip] * scale_factor[3] + add_offset[3];
    }

    // read in radiance data
    rad_data = (float *) calloc(npix*nbands, sizeof(float));

    nc_inq_ncid(file->sd_id,"geophysical_data" , &grp_id);
    
    for (ib = 0; ib < nbands; ib++) {
        
        retval = nc_get_vara_float(grp_id, ib, start, count,rad_data);
        if (retval == FAIL) {
            fprintf(stderr,
                    "-E- %s line %d: nc_get_vara_float failed for file, %s  band %d\n",
                    __FILE__, __LINE__, file->name, ib);
            return (1);
        }

        // copy to Lt record.
        for (ip = spix; ip < npix; ip++) {
            ipb = ip * nbands + ib;
            l1rec->Lt[ipb] = rad_data[ip];
            l1rec->detnum = ip;

            if (l1rec->Lt[ipb] < 0.0)
                l1rec->Lt[ipb] = 0.0001;
        }
    } // for ib
    free(rad_data);
    free(lat);
    free(lon);
    free(solz);
    free(sola);
    free(senz);
    free(sena);
    if (scan == file->nscan) {
        free(scale_factor);
        free(add_offset);
    }

    l1rec->sensorID = file->sensorID;
    l1rec->npix = file->npix;

    return (LIFE_IS_GOOD);
}

int closel1_nc_generic(filehandle *file) {
    int retval;

    retval = nc_close(file->sd_id);
    if (retval == FAIL) {
        fprintf(stderr, "-E- %s line %d: nc_close failed for file, %s.\n",
        __FILE__, __LINE__, file->name);
        return (1);
    }

    return (LIFE_IS_GOOD);
}

