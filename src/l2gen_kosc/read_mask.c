#include "l12_proto.h"
#include "nc_gridutils.h"

/* global variables */
static int landindex = 0;
static const char* landnames[] = {"watermask", "landmask", "z", NULL};
static grid_info_t* landmask_grid = {0};

static int bathindex = 1;
//static const char* bathnames[] = {"height", "z", "depth",
//                                  "elevation", "altitude", NULL};


/* Land Mask tools */
int land_mask_init(char *file) {
    int status;

    /* load structure from NetCDF file*/
    landmask_grid = allocate_gridinfo();
    status = init_gridinfo(file, landnames, landmask_grid);

    /* if not NetCDF, free structure and read binary file */
    if (status != NC_NOERR) {
        free(landmask_grid);
        landmask_grid = NULL;
        status = b128_msk_init(file,landindex);
    }

    return status;
}

int land_mask( float lat, float lon ) {
    double value;
    int status;

    /* get value from NetCDF file*/
    if (landmask_grid != NULL) {
        status = get_bylatlon(landmask_grid, lat, lon, &value);
        if(status) {
            fprintf(stderr, "-W- %s:%d: lat=%f, lon=%f not in landmask file %s\n...assuming it's water.\n",
                    __FILE__, __LINE__, lat, lon, landmask_grid->file);
            return 0;
        }
        return ( (short) value != 1);   // convert from water=1 to land=1
    }

    /* otherwise get from binary file */
    else return b128_msk_get(lat,lon,landindex);
}


/* Bathymetry placeholders */
int bath_mask_init(char *file) {
    int status;
    status = b128_msk_init(file,bathindex);
    return status;
}
int bath_mask( float lat, float lon ) {
    return b128_msk_get(lat,lon,bathindex);
}
