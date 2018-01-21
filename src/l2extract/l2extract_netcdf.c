/*
 * subroutine to extract a netCDF L2 file
 */

#include <netcdf.h>

#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <unistd.h>


#include <dfutils.h>
#include <genutils.h>

#include <l12_parms.h>
#include <l2_flags.h>

#include "l2extract.h"

#define MAX_VARIABLES 128

#define NCDIE(function){ \
  int status = function; \
  switch(status){ \
    case NC_NOERR: break; \
    default: printf("NetCDF error: file %s, line %d, %s\n", \
            __FILE__, __LINE__, nc_strerror(status)); \
    exit(1); \
  } \
}

void copyGlobalAttributes(int ncid_r, int ncid_w) {
    int numAtts;
    int i;
    char name[NC_MAX_NAME+1];

    NCDIE(nc_inq_natts(ncid_r, &numAtts));
    for(i=0; i<numAtts; i++) {
        NCDIE(nc_inq_attname(ncid_r, NC_GLOBAL, i, name));
        NCDIE(nc_copy_att(ncid_r, NC_GLOBAL, name, ncid_w, NC_GLOBAL));
    }
}

void copyVariableAttributes(int ncid_r, int varid_r, int ncid_w, int varid_w) {
    int numAtts;
    int i;
    char name[NC_MAX_NAME+1];
    int status;

    NCDIE(nc_inq_varnatts(ncid_r, varid_r, &numAtts));
    for(i=0; i<numAtts; i++) {
        NCDIE(nc_inq_attname(ncid_r, varid_r, i, name));
        NCDIE(nc_copy_att(ncid_r, varid_r, name, ncid_w, varid_w));
    }
}

/**
 * copy a piece of a variable and all of it's attributes to another file
 *
 * @param ncid_r netCDF file or group to read
 * @param name name of the variable
 * @param start location to start copying from
 * @param count how many of each dimension to copy
 * @param dimIds dimension IDs from the destination file to attach to the new variable
 * @param ncid_w netCDF file or group to write the variable to
 * @param data write this data to the new variable or copy from read variable if NULL
 */
void copyVariable(int ncid_r, const char* name, size_t* start, size_t* count,
        int* dimIds, int ncid_w, void* data) {

    int status;
    int varid_r;
    int varid_w;
    nc_type xtype;
    int numDims;
    int numAtts;
    size_t chunkSize;
    size_t chunkNelems;
    float chunkPreemption;

    size_t  typeSize;
    int arraySize;
    int i;

    static int localDataSize = 0;
    static char* localData = NULL;

    int shuffle;
    int deflate;
    int deflateLevel;

    status = nc_inq_varid(ncid_r, name, &varid_r);
    if(status == NC_NOERR) {
        NCDIE(nc_inq_var(ncid_r, varid_r, NULL, &xtype, &numDims, NULL, &numAtts));
        NCDIE(nc_get_var_chunk_cache(ncid_r, varid_r,
                &chunkSize, &chunkNelems, &chunkPreemption));
        NCDIE(nc_inq_var_deflate(ncid_r, varid_r,
                &shuffle, &deflate, &deflateLevel));

        NCDIE(nc_def_var(ncid_w, name, xtype, numDims, dimIds, &varid_w));
        NCDIE(nc_set_var_chunk_cache(ncid_w, varid_w,
                chunkSize, chunkNelems, chunkPreemption));
        NCDIE(nc_def_var_deflate(ncid_w, varid_w,
                shuffle, deflate, deflateLevel));

        copyVariableAttributes(ncid_r, varid_r, ncid_w, varid_w);

        NCDIE(nc_inq_type(ncid_r, xtype, NULL, &typeSize));

        if(data==NULL) {
            // calc array size in num of bytes
            arraySize = typeSize;
            for(i=0;i<numDims; i++) {
                arraySize *= count[i];
            }

            // allocate array
            if(arraySize > localDataSize) {
                if(localData)
                    free(localData);
                localDataSize = arraySize;
                localData = (char*) malloc(localDataSize);
                if(localData == NULL) {
                    printf("-E- %s %d: could not allocate data for variable %s\n",
                            __FILE__, __LINE__, name);
                    exit(EXIT_FAILURE);
                }
            }

            NCDIE(nc_get_vara(ncid_r, varid_r, start, count, localData));
            NCDIE(nc_put_var(ncid_w, varid_w, localData));
        } else {
            NCDIE(nc_put_var(ncid_w, varid_w, data));
        }

    } // found the variable
}


/**
 * extract a L2 netCDF file
 *
 * @param infile input file name
 * @param outfile output file name
 * @param spix start pixel (1 based)
 * @param epix ending pixel (1 based)
 * @param sscan start line (1 based)
 * @param escan end line (1 based)
 * @param prodlist product list, comma separated, empty string outputs all products
 * @return 0 = success
 */
int extractNetCDF(const char* infile, const char* outfile, int spix, int epix,
        int sscan, int escan, const char* prodlist) {

    int status;
    idDS ds_id_r, ds_id_w;  // data set file ids for reading and writing
    int rootGroup_r;        // netCDF group for read file root
    int rootGroup_w;        // netCDF group for write file root
    int group_r;            // netCDF group for read file
    int group_w;            // netCDF group for write file
    int subGroup_r;         // netCDF sub group for read file
    int subGroup_w;         // netCDF sub group for write file
    int dim_id;
    int varid;
    size_t start[3] = {0,0,0};
    size_t count[3] = {1,1,1};
    int dimIds[3] = {0,0,0};

    size_t numLines_r; // number of line in read file
    size_t numPixels_r; // number of pixels in read file
    size_t numPixelControlPoints_r; // number of pixel control points in read file
    size_t numLines_w; // number of line in write file
    size_t numPixels_w; // number of pixels in write file
    size_t numPixelControlPoints_w; // number of pixel control points in read file
    size_t numTotalBands; // number of total bands (vis + IR) in file
    size_t numVisibleBands; // number of visible bands in file
    int s_cntl_pt;  // starting pixel control point to write (0 based)
    int e_cntl_pt;  // ending pixel control point to write (0 based)
    int numLinesDimId; // Number_of_Scan_Lines dimension id
    int numPixelsDimId; // Number_of_Scan_Lines dimension id
    int numControlPointsDimId; // Number_of_Scan_Lines dimension id
    int numTotalBandsDimId; // Number_of_Scan_Lines dimension id
    int numVisibleBandsDimId; // Number_of_Scan_Lines dimension id

    int i;
    int numVariables;
    int variableIds[MAX_VARIABLES];
    char name[256];

    // data sets read in
    int* cntl_pt_cols;
    int* cntl_pt_rows;
    float* latArray;
    float* lonArray;
    int* flagArray;

    ds_id_r = openDS(infile);
    if(ds_id_r.fid == FAIL) {
        printf("could not open \"%s\" for reading.\n", infile);
        exit(EXIT_FAILURE);
    }
    if(ds_id_r.fftype != DS_NCDF) {
        printf("could not open \"%s\" is not a netCDF4 file.\n", infile);
        exit(EXIT_FAILURE);
    }
    rootGroup_r = ds_id_r.fid;

    /* Get # of scans and # of pixels */
    /* ------------------------------ */
    status = nc_inq_dimid(rootGroup_r, "number_of_lines", &dim_id);
    if (status) {
        nc_inq_dimid( rootGroup_r, "Number_of_Scan_Lines", &dim_id);
    }
    nc_inq_dimlen(rootGroup_r, dim_id, &numLines_r);
    status = nc_inq_dimid( rootGroup_r, "pixels_per_line", &dim_id);
    if (status) {
        nc_inq_dimid( rootGroup_r, "Pixels_per_Scan_Line", &dim_id);

    }
    nc_inq_dimlen( rootGroup_r, dim_id, &numPixels_r);
    nc_inq_dimid( rootGroup_r, "pixel_control_points", &dim_id);
    nc_inq_dimlen( rootGroup_r, dim_id, &numPixelControlPoints_r);
    status = nc_inq_dimid( rootGroup_r, "number_of_bands", &dim_id);
    if (status) {
        nc_inq_dimid( rootGroup_r, "total_band_number", &dim_id);
    }
    nc_inq_dimlen( rootGroup_r, dim_id, &numTotalBands);
    status = nc_inq_dimid( rootGroup_r, "number_of_reflective_bands", &dim_id);
    if (status){
        nc_inq_dimid( rootGroup_r, "band_number", &dim_id);
    }
    nc_inq_dimlen( rootGroup_r, dim_id, &numVisibleBands);

    char* title = readAttrStr(ds_id_r, "title");
    if (strstr(title, "Level-2") == NULL) {
        printf("\"%s\" is not a L2 file.\n", infile);
        exit(EXIT_FAILURE);
    }

    if (sscan < 1) {
      sscan = 1;
    }
    if(sscan >= numLines_r) {
        printf("sscan needs to be less than number of scans in file.\n");
        exit(EXIT_FAILURE);
    }
    if (escan < 1 || escan > numLines_r) {
      escan = numLines_r;
    }
    if(escan < sscan) {
        printf("escan needs to be greater than sscan.\n");
        exit(EXIT_FAILURE);
    }
    numLines_w = escan-sscan+1;

    if (spix < 1) {
      spix = 1;
    }
    if(spix >= numPixels_r) {
        printf("spix needs to be less than number of pixels in file.\n");
        exit(EXIT_FAILURE);
    }
    if (epix < 1 || epix > numPixels_r) {
      epix = numPixels_r;
    }
    if(epix < spix) {
        printf("epix needs to be greater than spix.\n");
        exit(EXIT_FAILURE);
    }

    // --------------------------------------------------------
    // set to navigation data group
    // --------------------------------------------------------
    nc_inq_ncid(rootGroup_r, "navigation_data", &group_r);

    // read control points
    cntl_pt_cols = malloc(numPixelControlPoints_r * sizeof(int));
    cntl_pt_rows = malloc(numLines_w * sizeof(int));
    if(cntl_pt_cols == NULL || cntl_pt_rows == NULL) {
        printf("could not allocate memory for cntl_pt_cols and cntl_pt_rows arrays\n");
        exit(EXIT_FAILURE);
    }
    start[0] = 0;
    count[0] = numPixelControlPoints_r;
    NCDIE(nc_inq_varid(group_r, "cntl_pt_cols", &varid));
    NCDIE(nc_get_vara_int(group_r, varid, start, count, cntl_pt_cols));

    // find starting pixel control point
    for(i=0; i<numPixelControlPoints_r; i++) {
        if(cntl_pt_cols[i] >= spix) {
            s_cntl_pt=i;
            break;
        }
    }
    if(i>=numPixelControlPoints_r) {
        printf("start pixel control point not found.\n");
        exit(EXIT_FAILURE);
    }

    // adjust spix if not on a control point
    if(cntl_pt_cols[s_cntl_pt] > spix) {
        if(s_cntl_pt > 0)
            s_cntl_pt--;
        spix = cntl_pt_cols[s_cntl_pt];
    }

    // find ending pixel control point
    e_cntl_pt = numPixelControlPoints_r - 1;
    for(i=0; i<numPixelControlPoints_r; i++) {
        if(cntl_pt_cols[i] > epix) {
            e_cntl_pt=i-1;
            break;
        }
    }

    printf("sscan: %d  escan: %d\n", sscan, escan);
    printf("spixl: %d  epixl: %d\n", spix, epix);

    numPixels_w = epix-spix+1;
    numPixelControlPoints_w = e_cntl_pt-s_cntl_pt+1;

    // read lat and lon
    latArray = malloc(numLines_w * numPixelControlPoints_w * sizeof(float));
    lonArray = malloc(numLines_w * numPixelControlPoints_w * sizeof(float));
    if(latArray == NULL || lonArray == NULL) {
        printf("could not allocate memory for lat and lon arrays\n");
        exit(EXIT_FAILURE);
    }
    start[0] = sscan-1;
    start[1] = s_cntl_pt;
    count[0] = numLines_w;
    count[1] = numPixelControlPoints_w;
    NCDIE(nc_inq_varid(group_r, "latitude", &varid));
    NCDIE(nc_get_vara_float(group_r, varid, start, count, latArray));
    NCDIE(nc_inq_varid(group_r, "longitude", &varid));
    NCDIE(nc_get_vara_float(group_r, varid, start, count, lonArray));

    // --------------------------------------------------------
    // set to geophysical data group
    // --------------------------------------------------------
    nc_inq_ncid(rootGroup_r, "geophysical_data", &group_r);

    flagArray = malloc(numLines_w * numPixels_w * sizeof(int));
    if(flagArray == NULL) {
        printf("could not allocate memory for flag array\n");
        exit(EXIT_FAILURE);
    }
    start[0] = sscan-1;
    start[1] = spix-1;
    count[0] = numLines_w;
    count[1] = numPixels_w;
    NCDIE(nc_inq_varid(group_r, "l2_flags", &varid));
    NCDIE(nc_get_vara_int(group_r, varid, start, count, flagArray));


    /* Calculate new navigation metadata */
    /* --------------------------------- */
    #define   GEOBOX_INC 20.0

    int startDone=0;
    int centerDone=0;
    int endDone=0;
    int ccol = numPixelControlPoints_w/2;
    int line, scol, ecol;
    int lineStart;
    float last_lat=0.0;
    int lastGoodLine;
    float geobox[4][100];
    int32 geobox_cnt=0;
    float gring_fval[100];
    int32 gring_ival[100];

    float startCenterLat, startCenterLon;
    float endCenterLat, endCenterLon;
    float geoLatMin, geoLatMax;
    float geoLonMin, geoLonMax;

    // loop to find beginning info
    for(line=0; line<numLines_w; line++) {
        lineStart = numPixelControlPoints_w*line;

        // find good start pixel
        if(!startDone) {
            for(scol=0; scol<ccol; scol++) {
                // check NAVFAIL flag
                if(flagArray[line*numPixels_w + cntl_pt_cols[scol]-1] ^ NAVFAIL) {
                    startDone = 1;
                    geobox[0][geobox_cnt] = lonArray[lineStart+scol];
                    geobox[1][geobox_cnt] = latArray[lineStart+scol];
                    break;
                } // flag good
            } // for col
        } // if not start

        // find good center pixel
        if(!centerDone) {
            // check NAVFAIL flag
            if(flagArray[line*numPixels_w + cntl_pt_cols[ccol]-1] ^ NAVFAIL) {
                centerDone = 1;
                startCenterLon = lonArray[lineStart+ccol];
                startCenterLat = latArray[lineStart+ccol];
                last_lat = startCenterLat;
            } // flag good
        } // if not center

        // find good end pixel
        if(!endDone) {
            for(ecol=numPixelControlPoints_w-1; ecol>=ccol; ecol--) {
                // check NAVFAIL flag
                if(flagArray[line*numPixels_w + cntl_pt_cols[ecol]-1] ^ NAVFAIL) {
                    endDone = 1;
                    geobox[2][geobox_cnt] = lonArray[lineStart+ecol];
                    geobox[3][geobox_cnt] = latArray[lineStart+ecol];
                    break;
                } // flag good
            } // for col
        } // if not start

        if(startDone && centerDone && endDone)
            break;
    }

    // set the min and max lat lon values
    geoLonMin = geoLonMax = geobox[0][geobox_cnt];
    geoLatMin = geoLatMax = geobox[1][geobox_cnt];
    if(geoLonMin > geobox[2][geobox_cnt])
        geoLonMin = geobox[2][geobox_cnt];
    if(geoLatMin > geobox[3][geobox_cnt])
        geoLatMin = geobox[3][geobox_cnt];
    if(geoLonMax < geobox[2][geobox_cnt])
        geoLonMax = geobox[2][geobox_cnt];
    if(geoLatMax < geobox[3][geobox_cnt])
        geoLatMax = geobox[3][geobox_cnt];
    if(geoLonMin > startCenterLon)
        geoLonMin = startCenterLon;
    if(geoLonMax < startCenterLon)
        geoLonMax = startCenterLon;

    geobox_cnt++;

    // loop through the rest of the lines
    for(; line<numLines_w; line++) {
        lineStart = numPixelControlPoints_w*line;

        // find first good pixel on line
        for(scol=0; scol<ccol; scol++) {
            // check NAVFAIL flag
            if(flagArray[line*numPixels_w + cntl_pt_cols[scol]-1] ^ NAVFAIL)
                break;
        }
        if(scol == ccol) // could not find start col, so skip this line
            continue;

        // find last good pixel
        for(ecol=numPixelControlPoints_w-1; ecol>=ccol; ecol--) {
             // check NAVFAIL flag
             if(flagArray[line*numPixels_w + cntl_pt_cols[ecol]-1] ^ NAVFAIL)
                 break;
        }
        if(ecol < ccol) // could not find end col, so skip this line
            continue;

        lastGoodLine = line;

        // set min/max for every line
        if(geoLonMax < lonArray[lineStart+scol])
            geoLonMax = lonArray[lineStart+scol];
        if(geoLonMax < lonArray[lineStart+ecol])
            geoLonMax = lonArray[lineStart+ecol];
        if(geoLonMin > lonArray[lineStart+scol])
            geoLonMin = lonArray[lineStart+scol];
        if(geoLonMin > lonArray[lineStart+ecol])
            geoLonMin = lonArray[lineStart+ecol];

        if(geoLatMax < latArray[lineStart+scol])
            geoLatMax = latArray[lineStart+scol];
        if(geoLatMax < latArray[lineStart+ecol])
            geoLatMax = latArray[lineStart+ecol];
        if(geoLatMin > latArray[lineStart+scol])
            geoLatMin = latArray[lineStart+scol];
        if(geoLatMin > latArray[lineStart+ecol])
            geoLatMin = latArray[lineStart+ecol];

        // load up geobox
        if (fabs(last_lat-latArray[ccol]) > GEOBOX_INC) {
            geobox[0][geobox_cnt] = lonArray[lineStart+scol];
            geobox[1][geobox_cnt] = latArray[lineStart+scol];
            geobox[2][geobox_cnt] = lonArray[lineStart+ecol];
            geobox[3][geobox_cnt] = latArray[lineStart+ecol];
            last_lat = latArray[lineStart + ccol];
            geobox_cnt++;
        }

    } // for lines

    // make sure we add the last line
    lineStart = numPixelControlPoints_w*lastGoodLine;

    // find first good pixel on line
    for(scol=0; scol<ccol; scol++) {
        // check NAVFAIL flag
        if(flagArray[lastGoodLine*numPixels_w + cntl_pt_cols[scol]-1] ^ NAVFAIL)
            break;
    }

    // find last good pixel
    for(ecol=numPixelControlPoints_w-1; ecol>=ccol; ecol--) {
         // check NAVFAIL flag
         if(flagArray[lastGoodLine*numPixels_w + cntl_pt_cols[ecol]-1] ^ NAVFAIL)
             break;
    }

    endCenterLon = lonArray[lineStart+ccol];
    endCenterLat = latArray[lineStart+ccol];

    geobox[0][geobox_cnt] = lonArray[lineStart+scol];
    geobox[1][geobox_cnt] = latArray[lineStart+scol];
    geobox[2][geobox_cnt] = lonArray[lineStart+ecol];
    geobox[3][geobox_cnt] = latArray[lineStart+ecol];
    geobox_cnt++;

    /* Create output netCDF file (delete if it exists) */
    /* ------------------------- */
    if( access( outfile, F_OK ) != -1 ) {
        if(unlink(outfile)) {
            printf("could not delete existing output file %s\n", outfile);
            exit(EXIT_FAILURE);
        }
    }
    ds_id_w = startDS(outfile, DS_NCDF, DS_WRITE, 5);
    rootGroup_w = ds_id_w.fid;

    /* Create dimensions */
    /* ----------------- */
    NCDIE(nc_def_dim(rootGroup_w, "number_of_lines", numLines_w, &numLinesDimId));
    NCDIE(nc_def_dim(rootGroup_w, "pixels_per_line", numPixels_w, &numPixelsDimId));
    NCDIE(nc_def_dim(rootGroup_w, "pixel_control_points", numPixelControlPoints_w, &numControlPointsDimId));
    NCDIE(nc_def_dim(rootGroup_w, "number_of_bands", numTotalBands, &numTotalBandsDimId));
    NCDIE(nc_def_dim(rootGroup_w, "number_of_reflective_bands", numVisibleBands, &numVisibleBandsDimId));

    // --------------------------------------------------------
    // set to sensor_band_parameters group
    // --------------------------------------------------------
    nc_inq_ncid(rootGroup_r, "sensor_band_parameters", &group_r);
    nc_def_grp(rootGroup_w, "sensor_band_parameters", &group_w);

    // copy vcal_gain(visibleBands)
    start[0] = 0;
    count[0] = numVisibleBands;
    dimIds[0] = numVisibleBandsDimId;

    // loop through all of the variables
    status = nc_inq_varids(group_r, &numVariables, variableIds);
    if(status == NC_NOERR) {
        if(numVariables > MAX_VARIABLES) {
            printf("-E- %s %d: too many variables in group \"sensor_band_parameters\"\n",
                    __FILE__, __LINE__);
            exit(0);
        }
        for(i=0; i<numVariables; i++) {
            status = nc_inq_varname(group_r, variableIds[i], name);
            if(strcmp(name, "wavelength") == 0) {
                count[0] = numTotalBands;
                dimIds[0] = numTotalBandsDimId;
                copyVariable(group_r, name, start, count, dimIds, group_w, NULL);
                count[0] = numVisibleBands;
                dimIds[0] = numVisibleBandsDimId;
            } else {
                copyVariable(group_r, name, start, count, dimIds, group_w, NULL);
            }
        }
    }

    // --------------------------------------------------------
    // set to sensor_band_parameters group
    // --------------------------------------------------------
    nc_inq_ncid(rootGroup_r, "scan_line_attributes", &group_r);
    nc_def_grp(rootGroup_w, "scan_line_attributes", &group_w);

    // set up the copy parameters
    start[0] = sscan-1;
    count[0] = numLines_w;
    dimIds[0] = numLinesDimId;

    // loop through all of the variables
    status = nc_inq_varids(group_r, &numVariables, variableIds);
    if(status == NC_NOERR) {
        if(numVariables > MAX_VARIABLES) {
            printf("-E- %s %d: too many variables in group \"scan_line_attributes\"\n",
                    __FILE__, __LINE__);
            exit(0);
        }
        for(i=0; i<numVariables; i++) {
            status = nc_inq_varname(group_r, variableIds[i], name);
            copyVariable(group_r, name, start, count, dimIds, group_w, NULL);
        }
    }

    // --------------------------------------------------------
    // set to geophysical_data group
    // --------------------------------------------------------
    nc_inq_ncid(rootGroup_r, "geophysical_data", &group_r);
    nc_def_grp(rootGroup_w, "geophysical_data", &group_w);

    // set up the copy parameters
    start[0] = sscan - 1;
    start[1] = spix - 1;
    count[0] = numLines_w;
    count[1] = numPixels_w;
    dimIds[0] = numLinesDimId;
    dimIds[1] = numPixelsDimId;

    // loop through all of the variables
    status = nc_inq_varids(group_r, &numVariables, variableIds);
    if(status == NC_NOERR) {
        if(numVariables > MAX_VARIABLES) {
            printf("-E- %s %d: too many variables in group \"scan_line_attributes\"\n",
                    __FILE__, __LINE__);
            exit(0);
        }
        for(i=0; i<numVariables; i++) {
            status = nc_inq_varname(group_r, variableIds[i], name);
            copyVariable(group_r, name, start, count, dimIds, group_w, NULL);
        }
    }

    // --------------------------------------------------------
    // set to navigation data group
    // --------------------------------------------------------
    nc_inq_ncid(rootGroup_r, "navigation_data", &group_r);
    nc_def_grp(rootGroup_w, "navigation_data", &group_w);

    // copy longitude(lines, pix_cntl_points)
    start[0] = sscan - 1;
    start[1] = s_cntl_pt;
    count[0] = numLines_w;
    count[1] = numPixelControlPoints_w;
    dimIds[0] = numLinesDimId;
    dimIds[1] = numControlPointsDimId;
    copyVariable(group_r, "longitude", start, count, dimIds, group_w, lonArray);
    copyVariable(group_r, "latitude", start, count, dimIds, group_w, latArray);

    // copy pixel control point cols
    // modify the cntl_pt_cols to write out
    for(i=0; i<numPixelControlPoints_w; i++) {
        cntl_pt_cols[i] = cntl_pt_cols[s_cntl_pt+i] - cntl_pt_cols[s_cntl_pt] + 1;
    }
    start[0] = 0;
    count[0] = numPixelControlPoints_w;
    dimIds[0] = numControlPointsDimId;
    copyVariable(group_r, "cntl_pt_cols", start, count, dimIds, group_w, cntl_pt_cols);

    // copy pixel control point cols
    // modify the cntl_pt_rows to write out
    cntl_pt_rows = malloc(numLines_w * sizeof(int));
    if(cntl_pt_cols == NULL || cntl_pt_rows == NULL) {
        printf("could not allocate memory for cntl_pt_cols and cntl_pt_rows arrays\n");
        exit(EXIT_FAILURE);
    }
    for(i=0; i<numLines_w; i++) {
        cntl_pt_rows[i] = i + 1;
    }
    start[0] = sscan-1;
    count[0] = numLines_w;
    dimIds[0] = numLinesDimId;
    copyVariable(group_r, "cntl_pt_rows", start, count, dimIds, group_w, cntl_pt_rows);
    copyVariable(group_r, "tilt", start, count, dimIds, group_w, NULL);


    // fill up gring
    int j = 1;
    gring_fval[0] = geobox[0][0];
    for (i=0; i<geobox_cnt; i++) {
        gring_fval[j++] = geobox[2][i];
    }
    for (i=0; i<geobox_cnt-1; i++) {
        gring_fval[j++] = geobox[0][geobox_cnt-1-i];
    }
    NCDIE(nc_put_att_float(group_w, NC_GLOBAL, "gringpointlongitude", NC_FLOAT, j, gring_fval));

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
    NCDIE(nc_put_att_float(group_w, NC_GLOBAL, "gringpointlatitude", NC_FLOAT, j, gring_fval));
    NCDIE(nc_put_att_int(group_w, NC_GLOBAL, "gringpointsequence", NC_INT, j, gring_ival));

    // --------------------------------------------------------
    // set to processing_control group
    // --------------------------------------------------------
    nc_inq_ncid(rootGroup_r, "processing_control", &group_r);
    nc_def_grp(rootGroup_w, "processing_control", &group_w);

    copyGlobalAttributes(group_r, group_w);

    // sub group input_parameters
    nc_inq_ncid(group_r, "input_parameters", &subGroup_r);
    nc_def_grp(group_w, "input_parameters", &subGroup_w);

    copyGlobalAttributes(subGroup_r, subGroup_w);

    // sub group flag_percentages
    nc_inq_ncid(group_r, "flag_percentages", &subGroup_r);
    nc_def_grp(group_w, "flag_percentages", &subGroup_w);

    copyGlobalAttributes(subGroup_r, subGroup_w);

    // --------------------------------------------------------
    // copy global attributes
    // --------------------------------------------------------
    copyGlobalAttributes(rootGroup_r, rootGroup_w);
    
    // write modified global attrbutes
    NCDIE(nc_put_att_float(rootGroup_w, NC_GLOBAL, "start_center_longitude", NC_FLOAT, 1, &startCenterLon));
    NCDIE(nc_put_att_float(rootGroup_w, NC_GLOBAL, "start_center_latitude", NC_FLOAT, 1, &startCenterLat));
    NCDIE(nc_put_att_float(rootGroup_w, NC_GLOBAL, "end_center_longitude", NC_FLOAT, 1, &endCenterLon));
    NCDIE(nc_put_att_float(rootGroup_w, NC_GLOBAL, "end_center_latitude", NC_FLOAT, 1, &endCenterLat));

    NCDIE(nc_put_att_float(rootGroup_w, NC_GLOBAL, "northernmost_latitude", NC_FLOAT, 1, &geoLatMax));
    NCDIE(nc_put_att_float(rootGroup_w, NC_GLOBAL, "southernmost_latitude", NC_FLOAT, 1, &geoLatMin));
    NCDIE(nc_put_att_float(rootGroup_w, NC_GLOBAL, "easternmost_longitude", NC_FLOAT, 1, &geoLonMax));
    NCDIE(nc_put_att_float(rootGroup_w, NC_GLOBAL, "westernmost_longitude", NC_FLOAT, 1, &geoLonMin));

    NCDIE(nc_put_att_float(rootGroup_w, NC_GLOBAL, "geospatial_lat_max", NC_FLOAT, 1, &geoLatMax));
    NCDIE(nc_put_att_float(rootGroup_w, NC_GLOBAL, "geospatial_lat_min", NC_FLOAT, 1, &geoLatMin));
    NCDIE(nc_put_att_float(rootGroup_w, NC_GLOBAL, "geospatial_lon_max", NC_FLOAT, 1, &geoLonMax));
    NCDIE(nc_put_att_float(rootGroup_w, NC_GLOBAL, "geospatial_lon_min", NC_FLOAT, 1, &geoLonMin));
    
    nc_close(rootGroup_r);
    nc_close(rootGroup_w);

    return EXIT_SUCCESS;
}

