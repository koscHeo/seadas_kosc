#ifndef _FILEHANDLE_H
#define _FILEHANDLE_H

#include <stdio.h>
#include <string.h>
#include "input_struc.h"
#include "hdf.h"
#include <hdf4utils.h>
#include <sensorInfo.h>
#include "l2prod_struc.h"

#define FORWARD      0
#define INVERSE_ZERO 1
#define INVERSE_NLW  2
#define INVERSE_LW   3

#define READ  0
#define WRITE 1

/* Use enum to give each input file format a unique number */
enum fmt_codes {
    FMT_AVIRIS,
    FMT_CLASSAVHRR,
    FMT_CZCSL1A,
    FMT_GOCIL1B,
    FMT_HICOL1B,
    FMT_L1BNCDF,
    FMT_L1HDF,
    FMT_L1XCAL,
    FMT_L2HDF,
    FMT_L2NCDF,
    FMT_L3BIN,
    FMT_L3MAP,
    FMT_MERISCC,
    FMT_MERISL1B,
    FMT_MERISL2,
    FMT_MODISGEO,	// MODIS Geolocation (hdf4)
    FMT_MODISL1B,	// MODIS L1B (hdf4, ocean-color band subset)
    FMT_HMODISL1B,	// MODIS L1B (hdf4, all bands)
    FMT_MOSL1B,
    FMT_OCM2L1B,
    FMT_OCML1B,
    FMT_OCML1BDB,
    FMT_OCTSL1A,
    FMT_OCTSL1B,
    FMT_OLCI,
    FMT_OLIL1B,
    FMT_ORCA,
    FMT_OSMIL1A,
    FMT_PRISM,
    FMT_SEAWIFSL1A,
    FMT_VIIRSGEO,	// VIIRS Geolocation (hdf5)
    FMT_VIIRSGEONC,	// VIIRS Geolocation (NetCDF4)
    FMT_VIIRSL1A,	// VIIRS Level-1A (NetCDF4)
    FMT_VIIRSL1B,	// VIIRS M-band (hdf5)
    FMT_VIIRSL1BNC,	// VIIRS M-band (NetCDF4)
    /* obsolete:
       FMT_L2VAL,
       FMT_MSL12L2,
    */
};

typedef struct filehandle_struct {
    char      name[FILENAME_MAX];
    int32_t   format;
    int32_t   sensorID;
    int32_t   subsensorID;
    char     spatialResolution[10];
    int       modis_subset_compat; /* force modis file to be read as subsetted */
    int32_t   length;
    int32_t   spix;                /* start pixel (0-based)                  */
    int32_t   epix;                /* end pixel (0-based)                    */
    int32_t   dpix;                /* pixel sub-sampling increment           */
    int32_t   npix;
    int32_t   ctl_pt_incr;
    int32_t   nscan;
    int32_t   nbands;
    int32_t   nbandsir;
    int32_t   *bindx;              /* index to closest seawifs band          */
    int32_t   ndets;
    int32_t   mode;
    char      l2prod[PRODSTRLEN];  /* list of L2 products to be included     */
    char      def_l2prod[PRODSTRLEN];  /* list of default L2 products      */
    int32_t   sd_id;               /* hdf file id for the opened output file */
    int32_t   tot_prod;            /* total # of L2 products to be created   */
    char      l2_prod_names[MAXPROD][32];
    l2prodstr* prodptr;        /* array of product structures */
    char      *pro_control;
    char      *input_parms;
    char      *input_files;
    char      *mask_names;
    instr     *input;
    int       percent_cloud;
    int       percent_land;
    int       percent_water;
    char      *calfile;
    char      *geofile;
    float     orbit_node_lon;
    int32_t   orbit_number;
    char      node_crossing_time[32];
    int32_t   flag_cnt[NFLAGS];
    int32     n_inprods;
    int32_t   terrain_corrected;
    int32_t   sv_with_moon;
    int32_t   ocrvc_opt;
    int32_t   grp_id[8];

    int32_t   deflate;

    void      *private_data;
    char      *viirscalparfile;
} filehandle;

#endif
