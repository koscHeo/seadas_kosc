#include "l12_proto.h"
#define MAX_LAT 89.95

/* DEM interpolation function prototypes */
int interp_nc_height (char  *demfile,
                      float *xlon, float *xlat, float *height);
int interp_dem_height(char  *demfile,
                      float *xlon, float *xlat, float *height);
static int (*interp_height)(char  *demfile,
                            float *xlon, float *xlat, float *height);

/* Terrain correction function prototypes */
int get_nc_height    (char  *demfile,
                      float *lon, float *lat,
                      float *senz, float *sena,
                      float *height);
int get_dem_height   (char  *demfile,
                      float *lon, float *lat,
                      float *senz, float *sena,
                      float *height);
static int (*correct_terrain)(char  *demfile,
                              float *lon, float *lat,
                              float *senz, float *sena,
                              float *height);

/**
 * Load DEM height for one pixel
 * @param[in] demfile
 * @param[in,out] l1str
 * @param[in] ip
 * @param[in] terrain_corrected
 * @return
 */
int get_height(char *demfile, l1str *l1rec, int32_t ip, int terrain_corrected)
{
    static int firstCall = 1;
    int status = 1;
    float *xlon   = &l1rec->lon [ip];
    float *xlat   = &l1rec->lat [ip];
    float *senz   = &l1rec->senz[ip];
    float *sena   = &l1rec->sena[ip];
    float *height = &l1rec->height[ip];

    /* Initial file tests */
    if (firstCall) {
        int ncid;
        int netcdf_dem;
        firstCall = 0;

        /* input file defined? */
        if (demfile == NULL || demfile[0] == 0) {
            fprintf(stderr, "-E- %s line %d: "
                    "Elevation file is NULL.\n",
                    __FILE__, __LINE__);
            return 1;
        }
        printf("Loading DEM info from %s\n",demfile);

        /* test for NetCDF input file */
        status = nc_open(demfile, NC_NOWRITE, &ncid);
        netcdf_dem = (status == NC_NOERR);
        if (netcdf_dem) nc_close(ncid);

        /* set function pointers according to input file type */
        if (netcdf_dem) {
            interp_height = interp_nc_height;
            correct_terrain = get_nc_height;
        } else {
            interp_height = interp_dem_height;
            correct_terrain = get_dem_height;
        }
    }

    /* Interpolate DEM height if terrain correction already done,
       or target too close to poles */
    if ( terrain_corrected 
         || (fabs((double)*xlat) > MAX_LAT) ) {
        status = interp_height(demfile, xlon, xlat, height);
        if (status) {
            fprintf(stderr, "-E- %s line %d: interp_height():\n",
                    __FILE__, __LINE__);
            fprintf(stderr,
                    "xlon=%f xlat=%f height=%f\n",
                    (double)*xlon, (double)*xlat, (double)*height);
        }
    }

    /* Otherwise, do terrain correction */
    else {
        status = correct_terrain(demfile, xlon, xlat, senz, sena, height);
        if (status) {
            fprintf(stderr, "-E- %s line %d: correct_terrain():\n",
                    __FILE__, __LINE__);
            fprintf(stderr,
                    "xlon=%f xlat=%f senz=%f sena=%f height=%f\n",
                    (double)*xlon, (double)*xlat,
                    (double)*senz, (double)*sena, (double)*height);
        }
    }

    return status;
}
