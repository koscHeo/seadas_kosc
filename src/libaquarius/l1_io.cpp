#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <timeutils.h>
#include <libgen.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>

#include <hdf.h>
#include "hdf5utils.h"
#include "hdf5_Aquarius.h"
#include "tlm_convert.h"
#include "passthebuck.h"

//#define ORBIT_AVG_PERIOD 5872.7983  // Fred Patt 06/17/11
#define ORBIT_AVG_PERIOD 5872.57  // JMG 08/10/11 after orbit lowering

#define DATACENTER "NASA/GSFC Aquarius Data Processing Center"
#define MISSIONCHARACTERISTICS "Nominal orbit: inclination=98.0 (Sun-synchronous); node=6PM (ascending); eccentricity=<0.002; altitude=657 km; ground speed=6.825 km/sec"

#define SENSORCHARACTERISTICS "Number of beams=3; channels per receiver=4; frequency 1.413 GHz; bits per sample=16; instatntaneous field of view=6.5 degrees; science data block period=1.44 sec."

double gpstai2unix( double gpstai);

int isleap(int year)
{
if ( ((year % 400) == 0) || (((year % 4) == 0) && ((year % 100) != 0)) )
        return TRUE;
else
        return FALSE;
}

using namespace std;

namespace Hdf {  

  /*----------------------------------------------------------------- */
  /* Create an HDF5 level1 file                                       */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::createl1( char* l1_filename, int32_t numBlocks,
			       int32_t num_SACD_HKT, double granStart, 
			       char* inputFiles, char* processControl,
			       char *softwareId)
  {
    first_write = true;

    // Set maxscan
    maxscan = numBlocks;

    // granuleStart is in seconds from UTC 6 January 1980

    /* Create the HDF file */
    h5fid = H5Fcreate(l1_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if(h5fid == FAIL) {
       fprintf(stderr,
       "-E- %s line %d: Could not create HDF file, %s .\n",
       __FILE__,__LINE__,l1_filename);
       return(HDF_FUNCTION_ERROR);
    }   
   
    openFlags = H5F_ACC_RDWR;

    grp0 = H5Gopen1(h5fid,"/");
    grp1 = H5Gcreate1(h5fid, "/Block Attributes", 0);
    grp2 = H5Gcreate1(h5fid, "/Raw Aquarius Data", 0);
    grp3 = H5Gcreate1(h5fid, "/SAC-D Telemetry", 0);
    grp4 = H5Gcreate1(h5fid, "/Navigation", 0);
    grp5 = H5Gcreate1(h5fid, "/Converted Telemetry", 0);

    // Initialize private variables
    blknum = 0;
    rad_frmnum = 0;
    atc_frmnum = 0;

    /*                                                                  */
    /* Create the block H5Ds                                            */
    /* ---------------------------------------------------------------- */
    /*                                                                  */

    // Setup dataset create property list
    hid_t proplist = H5Pcreate( H5P_DATASET_CREATE);

    // Block Time

    // Set fill to -1
    double time_fill_value = -1.0;
    H5Pset_fill_value( proplist, H5T_NATIVE_DOUBLE, (void *) &time_fill_value);

    PTB( CreateH5D(
    grp1,                                            /* file id         */
    "blk_sec",                                       /* short name      */
    "Block time, seconds of day",                    /* long name       */
    "seconds",                                       /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_DOUBLE,                               /* HDF number type */
    1,                                               /* rank            */
    numBlocks, 1, 1, 1, 1, 1,                        /* dimension sizes */
    "Number of Blocks", NULL, NULL, NULL, NULL, NULL,/* dimension names */
    proplist) );

    // ATC Frame Number
    int32_t fill_value = -1;
    H5Pset_fill_value( proplist, H5T_STD_I32LE, (void *) &fill_value);

    PTB( CreateH5D(
    grp1,                                            /* file id         */
    "atc_frmnum",                                    /* short name      */
    "ATC Frame Number",                              /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_STD_I32LE,                                   /* HDF number type */
    1,                                               /* rank            */
    numBlocks, 1, 1, 1, 1, 1,                        /* dimension sizes */
    "Number of Blocks", NULL, NULL, NULL, NULL, NULL,/* dimension names */
    proplist) );


    // ATC Sub-Frame Number
    int8_t char_fill_value = -1;
    H5Pset_fill_value( proplist, H5T_NATIVE_UCHAR, (void *) &char_fill_value);

    PTB( CreateH5D(
    grp1,                                            /* file id         */
    "atc_subframe",                                  /* short name      */
    "ATC Sub-Frame Number",                          /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_UCHAR,                                 /* HDF number type */
    1,                                               /* rank            */
    numBlocks, 1, 1, 1, 1, 1,                        /* dimension sizes */
    "Number of Blocks", NULL, NULL, NULL, NULL, NULL,/* dimension names */
    proplist) );

    // Radiometer Frame Number
    H5Pset_fill_value( proplist, H5T_STD_I32LE, (void *) &fill_value);

    PTB( CreateH5D(
    grp1,                                            /* file id         */
    "rad_frmnum",                                    /* short name      */
    "Radiometer Frame Number",                       /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_STD_I32LE,                                 /* HDF number type */
    1,                                               /* rank            */
    numBlocks, 1, 1, 1, 1, 1,                        /* dimension sizes */
    "Number of Blocks", NULL, NULL, NULL, NULL, NULL,/* dimension names */
    proplist) );


    // Radiometer Sub-Frame Number
    H5Pset_fill_value( proplist, H5T_NATIVE_UCHAR, (void *) &char_fill_value);

    PTB( CreateH5D(
    grp1,                                            /* file id         */
    "rad_subframe",                                  /* short name      */
    "Radiometer Sub-Frame Number",                   /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_UCHAR,                                 /* HDF number type */
    1,                                               /* rank            */
    numBlocks, 1, 1, 1, 1, 1,                        /* dimension sizes */
    "Number of Blocks", NULL, NULL, NULL, NULL, NULL,/* dimension names */
    proplist) );

    // Synch and Time Tag
    PTB( CreateH5D(
    grp2,                                            /* file id         */
    "gps_time_tag",                                  /* short name      */
    "Block GPS time tag",                            /* long name       */
    "seconds",                                       /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_STD_I32LE,                                   /* HDF number type */
    1,                                               /* rank            */
    numBlocks, 1, 1, 1, 1, 1,                        /* dimension sizes */
    "Number of Blocks", NULL, NULL, NULL, NULL, NULL,/* dimension names */
    H5P_DEFAULT) );

    PTB( CreateH5D(
    grp2,                                            /* file id         */
    "time_tag_offset",                               /* short name      */
    "Block time offset from GPS",                    /* long name       */
    "62.5 nanosec units",                            /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_STD_I32LE,                                   /* HDF number type */
    1,                                               /* rank            */
    numBlocks, 1, 1, 1, 1, 1,                        /* dimension sizes */
    "Number of Blocks", NULL, NULL, NULL, NULL, NULL,/* dimension names */
    H5P_DEFAULT) );

    PTB( CreateH5D(
    grp2,                                            /* file id         */
    "start_synch",                                   /* short name      */
    "Start-synch word",                              /* long name       */
    "dummy",                                         /* units           */
    0, 0,                                            /* valid range     */
    //    0xefbeadde, 0xefbeadde,                    /* valid range     */
    0, 0,                                            /* slope           */
    H5T_STD_I32LE,                                 /* HDF number type */
    1,                                               /* rank            */
    numBlocks, 1, 1, 1, 1, 1,                        /* dimension sizes */
    "Number of Blocks", NULL, NULL, NULL, NULL, NULL,/* dimension names */
    H5P_DEFAULT) );

    /// Check if it should be LE !!!!!!!!!!!!
    PTB( CreateH5D(
    grp2,                                            /* file id         */
    "checksum",                                      /* short name      */
    "Checksum",                                      /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_USHORT,                               /* HDF number type */
    1,                                               /* rank            */
    numBlocks, 1, 1, 1, 1, 1,                        /* dimension sizes */
    "Number of Blocks", NULL, NULL, NULL, NULL, NULL,/* dimension names */
    H5P_DEFAULT) );


    // Aquarius Raw Telemetry
    PTB( CreateH5D(
    grp2,                                            /* file id         */
    "deploy_tlm",                                    /* short name      */
    "Antenna deployment telemetry",                  /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_UCHAR,                                /* HDF number type */
    2,                                               /* rank            */
    numBlocks, 5, 1, 1, 1, 1,                        /* dimension sizes */
    "Number of Blocks", NULL, NULL, NULL, NULL, NULL,/* dimension names */
    H5P_DEFAULT) );

    PTB( CreateH5D(
    grp2,                                            /* file id         */
    "icds_status",                                   /* short name      */
    "ICDS processing status",                        /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_UCHAR,                                /* HDF number type */
    2,                                               /* rank            */
    numBlocks, 12, 1, 1, 1, 1,                       /* dimension sizes */
    "Number of Blocks", NULL, NULL, NULL, NULL, NULL,/* dimension names */
    H5P_DEFAULT) );

    PTB( CreateH5D(
    grp2,                                            /* file id         */
    "icds_tlm",                                      /* short name      */
    "ICDS engineering telemetry",                    /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_UCHAR,                                /* HDF number type */
    2,                                               /* rank            */
    numBlocks, 24, 1, 1, 1, 1,                       /* dimension sizes */
    "Number of Blocks", NULL, NULL, NULL, NULL, NULL,/* dimension names */
    H5P_DEFAULT) );

    PTB( CreateH5D(
    grp2,                                            /* file id         */
    "temp_tlm",                                      /* short name      */
    "External Temperature Sensors telemetry",        /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_UCHAR,                                /* HDF number type */
    2,                                               /* rank            */
    numBlocks, 70, 1, 1, 1, 1,                       /* dimension sizes */
    "Number of Blocks", NULL, NULL, NULL, NULL, NULL,/* dimension names */
    H5P_DEFAULT) );

    PTB( CreateH5D(
    grp2,                                            /* file id         */
    "apdu_tlm",                                      /* short name      */
    "Aquarius Power Distribution Unit telemetry",    /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_UCHAR,                                /* HDF number type */
    2,                                               /* rank            */
    numBlocks, 5, 1, 1, 1, 1,                        /* dimension sizes */
    "Number of Blocks", NULL, NULL, NULL, NULL, NULL,/* dimension names */
    H5P_DEFAULT) );

    PTB( CreateH5D(
    grp2,                                            /* file id         */
    "atc_tlm",                                       /* short name      */
    "Active Thermal Control Unit telemetry",         /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_UCHAR,                                /* HDF number type */
    2,                                               /* rank            */
    numBlocks, 36, 1, 1, 1, 1,                       /* dimension sizes */
    "Number of Blocks", NULL, NULL, NULL, NULL, NULL,/* dimension names */
    H5P_DEFAULT) );

    PTB( CreateH5D(
    grp2,                                            /* file id         */
    "radiom_rt_tlm",                                 /* short name      */
    "Radiometer real-time telemetry",                /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_UCHAR,                                /* HDF number type */
    3,                                               /* rank            */
    numBlocks/BLOCKS_PER_FRAME,                      /* dimension sizes */
    BLOCKS_PER_FRAME, 
    50, 1, 1, 1,
    "Number of Frames", "Blocks per Frame",          /* dimension names */
    NULL, NULL, NULL, NULL,
    H5P_DEFAULT) );

    PTB( CreateH5D(
    grp2,                                            /* file id         */
    "radiom_nrt_tlm",                                /* short name      */
    "Radiometer non-real-time telemetry",            /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_UCHAR,                                /* HDF number type */
    3,                                               /* rank            */
    numBlocks/BLOCKS_PER_FRAME,                      /* dimension sizes */
    BLOCKS_PER_FRAME, 
    10, 1, 1, 1,
    "Number of Frames", "Blocks per Frame",          /* dimension names */
    NULL, NULL, NULL, NULL,
    H5P_DEFAULT) );

    PTB( CreateH5D(
    grp2,                                            /* file id         */
    "scatter_tlm",                                   /* short name      */
    "Scatterometer telemetry",                       /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_UCHAR,                                /* HDF number type */
    2,                                               /* rank            */
    numBlocks, 37, 1, 1, 1, 1,                       /* dimension sizes */
    "Number of Blocks", NULL, NULL, NULL, NULL, NULL,/* dimension names */
    H5P_DEFAULT) );


    // Raw Radiometer Science Data
    PTB( CreateH5D(
    grp2,                                            /* file id         */
    "radiom_header",                                 /* short name      */
    "Radiometer block header",                       /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_USHORT,                               /* HDF number type */
    2,                                               /* rank            */
    numBlocks/BLOCKS_PER_FRAME,                      /* dimension sizes */
    BLOCKS_PER_FRAME, 1, 1, 1, 1,
    "Number of Frames", "Blocks per Frame",          /* dimension names */
    NULL, NULL, NULL, NULL,
    H5P_DEFAULT) );

    PTB( CreateH5D(
    grp2,                                            /* file id         */
    "radiom_signals",                                /* short name      */
    "Radiometer Antenna Looks",                      /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_USHORT,                               /* HDF number type */
    6,                                               /* rank            */
    numBlocks/BLOCKS_PER_FRAME,                      /* dimension sizes */
    BLOCKS_PER_FRAME, 
    RADIOMETER_SUBCYCLES,
    RADIOMETER_SIGNALS_PER_SUBCYCLE,
    NUMBER_OF_BEAMS, 
    RADIOMETER_POLARIZATIONS,
    "Number of Frames", "Blocks per Frame",          /* dimension names */
    "Radiometer Subcycles",
    "Radiometer Signals per Subcycle",
    "Number of Beams", "Radiometer Polarizations",
    H5P_DEFAULT) );


    PTB( CreateH5D(
    grp2,                                            /* file id         */
    "radiom_cnd",                                    /* short name      */
    "Radiometer CND Looks",                          /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_USHORT,                               /* HDF number type */
    5,                                               /* rank            */
    numBlocks/BLOCKS_PER_FRAME,                      /* dimension sizes */
    BLOCKS_PER_FRAME, 
    RADIOMETER_SUBCYCLES,
    NUMBER_OF_BEAMS, 
    RADIOMETER_POLARIZATIONS, 1,
    "Number of Frames", "Blocks per Frame",          /* dimension names */
    "Radiometer Subcycles",
    "Number of Beams", "Radiometer Polarizations",
    NULL,
    H5P_DEFAULT) );


    PTB( CreateH5D(
    grp2,                                            /* file id         */
    "radiom_lavg",                                   /* short name      */
    "Radiometer Long Accumulations",                 /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_USHORT,                               /* HDF number type */
    5,                                               /* rank            */
    numBlocks/BLOCKS_PER_FRAME,                      /* dimension sizes */
    BLOCKS_PER_FRAME,
    RADIOMETER_LONG_ACCUM, 
    NUMBER_OF_BEAMS, 
    RADIOMETER_POLARIZATIONS, 1,
    "Number of Frames", "Blocks per Frame",          /* dimension names */
    "Radiometer Long Accumulations",
    "Number of Beams", "Radiometer Polarizations",
    NULL,
    H5P_DEFAULT) );


    // Raw Scatterometer Science Data
    PTB( CreateH5D(
    grp2,                                            /* file id         */
    "scatter_headers",                               /* short name      */
    "Scatterometer subcycle headers",                /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_UCHAR,                                /* HDF number type */
    2,                                               /* rank            */
    numBlocks,                                       /* dimension sizes */
    SCATTEROMETER_SUBCYCLES, 1, 1, 1, 1,
    "Number of Blocks", "Scatterometer Subcycles",   /* dimension names */
    NULL, NULL, NULL, NULL,
    H5P_DEFAULT) );

    PTB( CreateH5D(
    grp2,                                            /* file id         */
    "scatter_pwr",                                   /* short name      */
    "Scatterometer Power",                           /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_USHORT,                               /* HDF number type */
    4,                                               /* rank            */
    numBlocks,                                       /* dimension sizes */
    SCATTEROMETER_SUBCYCLES, NUMBER_OF_BEAMS,
    SCATTEROMETER_POLARIZATIONS, 1, 1, 
    "Number of Blocks", "Scatterometer Subcycles",   /* dimension names */
    "Number of Beams", "Scatterometer Polarizations",
    NULL, NULL, 
    H5P_DEFAULT) );

    PTB( CreateH5D(
    grp2,                                            /* file id         */
    "scatter_loop",                                  /* short name      */
    "Scatterometer Loopback Measurements",           /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_USHORT,                               /* HDF number type */
    3,                                               /* rank            */
    numBlocks,                                       /* dimension sizes */
    NUMBER_OF_BEAMS,
    SCATTEROMETER_POLARIZATIONS, 1, 1, 1, 
    "Number of Blocks",                              /* dimension names */
    "Number of Beams", "Scatterometer Polarizations",
    NULL, NULL, NULL,
    H5P_DEFAULT) );

    PTB( CreateH5D(
    grp2,                                            /* file id         */
    "scatter_dc",                                    /* short name      */
    "Scatterometer DC data",                         /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_USHORT,                               /* HDF number type */
    2,                                               /* rank            */
    numBlocks,                                       /* dimension sizes */
    2, 1, 1, 1, 1, 
    "Number of Blocks",                              /* dimension names */
    NULL, NULL, NULL, NULL, NULL,
    H5P_DEFAULT) );

    PTB( CreateH5D(
    grp2,                                            /* file id         */
    "scatter_rfi",                                   /* short name      */
    "Scatterometer RFI flags for H-pol",             /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_UCHAR,                                /* HDF number type */
    2,                                               /* rank            */
    numBlocks,                                       /* dimension sizes */
    NUMBER_OF_BEAMS+1, 1, 1, 1, 1, 
    "Number of Blocks",                              /* dimension names */
    NULL, NULL, NULL, NULL, NULL,
    H5P_DEFAULT) );

    PTB( CreateH5D(
    grp2,                                            /* file id         */
    "pad",                                           /* short name      */
    "Pad",                                           /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_UCHAR,                                /* HDF number type */
    2,                                               /* rank            */
    numBlocks,                                       /* dimension sizes */
    35, 1, 1, 1, 1, 
    "Number of Blocks",                              /* dimension names */
    NULL, NULL, NULL, NULL, NULL,
    H5P_DEFAULT) );


    // SAC-D Telemetry
    PTB( CreateH5D(
    grp3,                                            /* file id         */
    "sacd_hkt",                                       /* short name      */
    "SAC-D raw housekeeping telemetry blocks",       /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_UCHAR,                                /* HDF number type */
    2,                                               /* rank            */
    num_SACD_HKT,                                    /* dimension sizes */
    4000, 1, 1, 1, 1, 
    "Number of SAC-D HKT",                           /* dimension names */
    NULL, NULL, NULL, NULL, NULL,
    H5P_DEFAULT) );


    PTB( CreateH5D(
    grp5,                                            /* file id         */
    "rfe1_analog_tlm",                               /* short name      */
    "RFE1 analog telemetry",                         /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_FLOAT,                                /* HDF number type */
    2,                                               /* rank            */
    numBlocks/BLOCKS_PER_FRAME,                      /* dimension sizes */
    23,
    1, 1, 1, 1,
    "Number of Frames",                              /* dimension names */
    NULL, NULL, NULL, NULL, NULL,
    H5P_DEFAULT) );

    PTB( CreateH5D(
    grp5,                                            /* file id         */
    "rfe2_analog_tlm",                               /* short name      */
    "RFE2 analog telemetry",                         /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_FLOAT,                                /* HDF number type */
    2,                                               /* rank            */
    numBlocks/BLOCKS_PER_FRAME,                      /* dimension sizes */
    23,
    1, 1, 1, 1,
    "Number of Frames",                              /* dimension names */
    NULL, NULL, NULL, NULL, NULL,
    H5P_DEFAULT) );

    PTB( CreateH5D(
    grp5,                                            /* file id         */
    "rfe3_analog_tlm",                               /* short name      */
    "RFE3 analog telemetry",                         /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_FLOAT,                                /* HDF number type */
    2,                                               /* rank            */
    numBlocks/BLOCKS_PER_FRAME,                      /* dimension sizes */
    23,
    1, 1, 1, 1,
    "Number of Frames",                              /* dimension names */
    NULL, NULL, NULL, NULL, NULL,
    H5P_DEFAULT) );



    PTB( CreateH5D(
    grp5,                                            /* file id         */
    "rbe1_analog_tlm",                               /* short name      */
    "RBE1 analog telemetry",                         /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_FLOAT,                                /* HDF number type */
    2,                                               /* rank            */
    numBlocks/BLOCKS_PER_FRAME,                      /* dimension sizes */
    21,
    1, 1, 1, 1,
    "Number of Frames",                              /* dimension names */
    NULL, NULL, NULL, NULL, NULL,
    H5P_DEFAULT) );

    PTB( CreateH5D(
    grp5,                                            /* file id         */
    "rbe2_analog_tlm",                               /* short name      */
    "RBE2 analog telemetry",                         /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_FLOAT,                                /* HDF number type */
    2,                                               /* rank            */
    numBlocks/BLOCKS_PER_FRAME,                      /* dimension sizes */
    21,
    1, 1, 1, 1,
    "Number of Frames",                              /* dimension names */
    NULL, NULL, NULL, NULL, NULL,
    H5P_DEFAULT) );

    PTB( CreateH5D(
    grp5,                                            /* file id         */
    "rbe3_analog_tlm",                               /* short name      */
    "RBE3 analog telemetry",                         /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_FLOAT,                                /* HDF number type */
    2,                                               /* rank            */
    numBlocks/BLOCKS_PER_FRAME,                      /* dimension sizes */
    21,
    1, 1, 1, 1,
    "Number of Frames",                              /* dimension names */
    NULL, NULL, NULL, NULL, NULL,
    H5P_DEFAULT) );


    PTB( CreateH5D(
    grp5,                                            /* file id         */
    "dpu_analog_tlm",                                /* short name      */
    "DPU analog telemetry",                          /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_FLOAT,                                /* HDF number type */
    2,                                               /* rank            */
    numBlocks/BLOCKS_PER_FRAME,                      /* dimension sizes */
    12,
    1, 1, 1, 1,
    "Number of Frames",                              /* dimension names */
    NULL, NULL, NULL, NULL, NULL,
    H5P_DEFAULT) );

    PTB( CreateH5D(
    grp5,                                            /* file id         */
    "dpu_status_tlm",                                /* short name      */
    "DPU discrete telemetry",                        /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_UCHAR,                                /* HDF number type */
    3,                                               /* rank            */
    numBlocks/BLOCKS_PER_FRAME,                      /* dimension sizes */
    BLOCKS_PER_FRAME, 
    13, 1, 1, 1,
    "Number of Frames", "Blocks per Frame",          /* dimension names */
    NULL, NULL, NULL, NULL,
    H5P_DEFAULT) );


    PTB( CreateH5D(
    grp5,                                            /* file id         */
    "radiom_nrt_tlm",                                /* short name      */
    "Radiometer discrete non-real-time telemetry",   /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_USHORT,                               /* HDF number type */
    3,                                               /* rank            */
    numBlocks/BLOCKS_PER_FRAME,                      /* dimension sizes */
    BLOCKS_PER_FRAME, 
    17, 1, 1, 1,
    "Number of Frames", "Blocks per Frame",          /* dimension names */
    NULL, NULL, NULL, NULL,
    H5P_DEFAULT) );


    PTB( CreateH5D(
    grp5,                                            /* file id         */
    "scatter_analog_tlm",                            /* short name      */
    "Scatterometer analog telemetry",                /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_FLOAT,                                /* HDF number type */
    2,                                               /* rank            */
    numBlocks,                                       /* dimension sizes */
    30,
    1, 1, 1, 1,
    "Number of Blocks",                              /* dimension names */
    NULL, NULL, NULL, NULL, NULL,
    H5P_DEFAULT) );

    PTB( CreateH5D(
    grp5,                                            /* file id         */
    "scatter_discrete_tlm",                          /* short name      */
    "Scatterometer discrete telemetry",              /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_UCHAR,                                /* HDF number type */
    2,                                               /* rank            */
    numBlocks,                                       /* dimension sizes */
    8,
    1, 1, 1, 1,
    "Number of Blocks",                              /* dimension names */
    NULL, NULL, NULL, NULL, NULL,
    H5P_DEFAULT) );


    PTB( CreateH5D(
    grp5,                                            /* file id         */
    "ext_temp_analog_tlm",                           /* short name      */
    "External Temperature Sensor analog telemetry",  /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_FLOAT,                                /* HDF number type */
    2,                                               /* rank            */
    numBlocks,                                       /* dimension sizes */
    38,
    1, 1, 1, 1,
    "Number of Blocks",                              /* dimension names */
    NULL, NULL, NULL, NULL, NULL,
    H5P_DEFAULT) );


    PTB( CreateH5D(
    grp5,                                            /* file id         */
    "deploy_analog_tlm",                             /* short name      */
    "Antenna deployment analog telemetry",           /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_FLOAT,                                /* HDF number type */
    2,                                               /* rank            */
    numBlocks,                                       /* dimension sizes */
    4,
    1, 1, 1, 1,
    "Number of Blocks",                              /* dimension names */
    NULL, NULL, NULL, NULL, NULL,
    H5P_DEFAULT) );

    PTB( CreateH5D(
    grp5,                                            /* file id         */
    "deploy_discrete_tlm",                           /* short name      */
    "Antenna deployment discrete telemetry",         /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_STD_I32LE,                                 /* HDF number type */
    2,                                               /* rank            */
    numBlocks,                                       /* dimension sizes */
    7,
    1, 1, 1, 1,
    "Number of Blocks",                              /* dimension names */
    NULL, NULL, NULL, NULL, NULL,
    H5P_DEFAULT) );


    PTB( CreateH5D(
    grp5,                                            /* file id         */
    "icds_analog_tlm",                               /* short name      */
    "ICDS analog telemetry",                         /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_FLOAT,                                /* HDF number type */
    2,                                               /* rank            */
    numBlocks,                                       /* dimension sizes */
    14,
    1, 1, 1, 1,
    "Number of Blocks",                              /* dimension names */
    NULL, NULL, NULL, NULL, NULL,
    H5P_DEFAULT) );

    PTB( CreateH5D(
    grp5,                                            /* file id         */
    "icds_discrete_tlm",                             /* short name      */
    "ICDS discrete telemetry",                       /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_STD_I32LE,                                 /* HDF number type */
    2,                                               /* rank            */
    numBlocks,                                       /* dimension sizes */
    8,
    1, 1, 1, 1,
    "Number of Blocks",                              /* dimension names */
    NULL, NULL, NULL, NULL, NULL,
    H5P_DEFAULT) );


    PTB( CreateH5D(
    grp5,                                            /* file id         */
    "apdu_analog_tlm",                               /* short name      */
    "Aquarius Power Distribution Unit analog telemetry", /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_FLOAT,                                /* HDF number type */
    2,                                               /* rank            */
    numBlocks,                                       /* dimension sizes */
    5,
    1, 1, 1, 1,
    "Number of Blocks",                              /* dimension names */
    NULL, NULL, NULL, NULL, NULL,
    H5P_DEFAULT) );


    PTB( CreateH5D(
    grp5,                                            /* file id         */
    "atc_omt1_analog_tlm",                           /* short name      */
    "ATC OMT1 analog telemetry",                     /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_FLOAT,                                /* HDF number type */
    2,                                               /* rank            */
    numBlocks/BLOCKS_PER_FRAME,                      /* dimension sizes */
    17,
    1, 1, 1, 1,
    "Number of Frames",                              /* dimension names */
    NULL, NULL, NULL, NULL, NULL,
    H5P_DEFAULT) );

    PTB( CreateH5D(
    grp5,                                            /* file id         */
    "atc_omt2_analog_tlm",                           /* short name      */
    "ATC OMT2 analog telemetry",                     /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_FLOAT,                                /* HDF number type */
    2,                                               /* rank            */
    numBlocks/BLOCKS_PER_FRAME,                      /* dimension sizes */
    17,
    1, 1, 1, 1,
    "Number of Frames",                              /* dimension names */
    NULL, NULL, NULL, NULL, NULL,
    H5P_DEFAULT) );

    PTB( CreateH5D(
    grp5,                                            /* file id         */
    "atc_omt3_analog_tlm",                           /* short name      */
    "ATC OMT3 analog telemetry",                     /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_FLOAT,                                /* HDF number type */
    2,                                               /* rank            */
    numBlocks/BLOCKS_PER_FRAME,                      /* dimension sizes */
    17,
    1, 1, 1, 1,
    "Number of Frames",                              /* dimension names */
    NULL, NULL, NULL, NULL, NULL,
    H5P_DEFAULT) );

    PTB( CreateH5D(
    grp5,                                            /* file id         */
    "atc_rbe_analog_tlm",                            /* short name      */
    "ATC RBE analog telemetry",                      /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_FLOAT,                                /* HDF number type */
    2,                                               /* rank            */
    numBlocks/BLOCKS_PER_FRAME,                      /* dimension sizes */
    17,
    1, 1, 1, 1,
    "Number of Frames",                              /* dimension names */
    NULL, NULL, NULL, NULL, NULL,
    H5P_DEFAULT) );


    PTB( CreateH5D(
    grp5,                                            /* file id         */
    "atc_omt1_discrete_tlm",                         /* short name      */
    "ATC OMT1 discrete telemetry",                   /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_UCHAR,                                /* HDF number type */
    2,                                               /* rank            */
    numBlocks/BLOCKS_PER_FRAME,                      /* dimension sizes */
    17,
    1, 1, 1, 1,
    "Number of Frames",                              /* dimension names */
    NULL, NULL, NULL, NULL, NULL,
    H5P_DEFAULT) );

    PTB( CreateH5D(
    grp5,                                            /* file id         */
    "atc_omt2_discrete_tlm",                         /* short name      */
    "ATC OMT2 discrete telemetry",                   /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_UCHAR,                                /* HDF number type */
    2,                                               /* rank            */
    numBlocks/BLOCKS_PER_FRAME,                      /* dimension sizes */
    17,
    1, 1, 1, 1,
    "Number of Frames",                              /* dimension names */
    NULL, NULL, NULL, NULL, NULL,
    H5P_DEFAULT) );

    PTB( CreateH5D(
    grp5,                                            /* file id         */
    "atc_omt3_discrete_tlm",                         /* short name      */
    "ATC OMT3 discrete telemetry",                   /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_UCHAR,                                /* HDF number type */
    2,                                               /* rank            */
    numBlocks/BLOCKS_PER_FRAME,                      /* dimension sizes */
    17,
    1, 1, 1, 1,
    "Number of Frames",                              /* dimension names */
    NULL, NULL, NULL, NULL, NULL,
    H5P_DEFAULT) );

    PTB( CreateH5D(
    grp5,                                            /* file id         */
    "atc_rbe_discrete_tlm",                          /* short name      */
    "ATC RBE discrete telemetry",                    /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_UCHAR,                                /* HDF number type */
    2,                                               /* rank            */
    numBlocks/BLOCKS_PER_FRAME,                      /* dimension sizes */
    17,
    1, 1, 1, 1,
    "Number of Frames",                              /* dimension names */
    NULL, NULL, NULL, NULL, NULL,
    H5P_DEFAULT) );

    H5Pclose( proplist);

    /*                                                                  */
    /* Write out some global attributes                                 */
    /* ---------------------------------------------------------------- */
    /*                                                                  */
    PTB( SetScalarH5A (grp0, "Product Name", H5T_STRING, 
		       basename(l1_filename)));
    PTB( SetScalarH5A (grp0, "Title", H5T_STRING, 
		       (VOIDP) "Aquarius Level 1A Data"));
    PTB( SetScalarH5A (grp0, "Data Center", H5T_STRING, (VOIDP) DATACENTER));
    PTB( SetScalarH5A (grp0, "Mission", H5T_STRING, (VOIDP) "SAC-D Aquarius"));
    PTB( SetScalarH5A (grp0, "Mission Characteristics", H5T_STRING, 
		       (VOIDP) MISSIONCHARACTERISTICS));
    PTB( SetScalarH5A (grp0, "Sensor" , H5T_STRING, (VOIDP) "Aquarius"));
    PTB( SetScalarH5A (grp0, "Data Type", H5T_STRING, (VOIDP) "SCI"));
    PTB( SetScalarH5A (grp0, "Software ID", H5T_STRING, (VOIDP) softwareId));

    // This function converts its input arguments to "YYYYDDDHHMMSSFFF"
    // where YYYY is the year, DDD is the day of the year, HH is the hour,
    // MM is the minute, SS is the second, and FFF is the fraction of a second.
    // The first argument represents the number of seconds elapsed since
    // 1-Jan-1970 00:00:00.000 GMT.  The second argument determines whether
    // the output string represents local time (L) or GMT (G).
    PTB( SetScalarH5A (grp0, "Processing Time", H5T_STRING, 
		       ydhmsf(time(NULL),'G')));

    PTB( SetScalarH5A (grp0, "Input Files", H5T_STRING, (VOIDP) inputFiles));
    PTB( SetScalarH5A (grp0, "Processing Control", H5T_STRING, 
		       (VOIDP) processControl));

    int32_t itemp;

    itemp = NUMBER_OF_BEAMS;
    PTB( SetScalarH5A (grp0, "Number of Beams" , H5T_STD_I32LE, 
		       (VOIDP) &itemp));
    itemp = RADIOMETER_POLARIZATIONS;
    PTB( SetScalarH5A (grp0, "Radiometer Polarizations" , H5T_STD_I32LE, 
		       (VOIDP) &itemp));
    itemp = RADIOMETER_SUBCYCLES;
    PTB( SetScalarH5A (grp0, "Radiometer Subcycles" , H5T_STD_I32LE, 
		       (VOIDP) &itemp));
    itemp = RADIOMETER_SIGNALS_PER_SUBCYCLE;
    PTB( SetScalarH5A (grp0, "Radiometer Signals per Subcycle" , 
		       H5T_STD_I32LE, (VOIDP) &itemp));
    itemp = RADIOMETER_LONG_ACCUM;
    PTB( SetScalarH5A (grp0, "Radiometer Long Accumulations" , 
		       H5T_STD_I32LE, (VOIDP) &itemp));
    itemp = SCATTEROMETER_POLARIZATIONS;
    PTB( SetScalarH5A (grp0, "Scatterometer Polarizations" , 
		       H5T_STD_I32LE, (VOIDP) &itemp));
    itemp = SCATTEROMETER_SUBCYCLES;
    PTB( SetScalarH5A (grp0, "Scatterometer Subcycles" , 
		       H5T_STD_I32LE, (VOIDP) &itemp));

    // granuleStart is in seconds from UTC 6 January 1980
    granuleStart = granStart; // Set granule start in private class variable
    string ydhmsf_str( ydhmsf( gpstai2unix(granuleStart), 'G'));
    PTB( SetScalarH5A (grp0, "Start Time", H5T_STRING, 
		       (void *) ydhmsf_str.c_str()));

    /*
      Can't redefine string attributes.  Might be HDF5 bug
    // Attribute Stubs
    PTB( SetScalarH5A (grp0, "End Time", H5T_STRING, 
		       (void *) ydhmsf_str.c_str()));
    PTB( SetScalarH5A (grp0, "Product Center Time", H5T_STRING, 
		       (void *) ydhmsf_str.c_str()));
    PTB( SetScalarH5A (grp0, "Node Crossing Time", H5T_STRING,
		       (void *) ydhmsf_str.c_str())); 
    */

    istringstream istr;

    istr.str( ydhmsf_str.substr(0,4)); istr >> itemp;
    PTB( SetScalarH5A (grp0, "Start Year", H5T_STD_I32LE, (VOIDP) &itemp));

    istr.clear(); istr.str( ydhmsf_str.substr(4,3)); istr >> itemp;
    PTB( SetScalarH5A (grp0, "Start Day", H5T_STD_I32LE, (VOIDP) &itemp));

    int32_t millisec = get_millisec( &ydhmsf_str);
    PTB( SetScalarH5A (grp0, "Start Millisec", H5T_STD_I32LE, 
		       (VOIDP) &millisec));

    istr.clear(); istr.str( ydhmsf_str.substr(0,4)); istr >> itemp;
    PTB( SetScalarH5A (grp0, "End Year", H5T_STD_I32LE, (VOIDP) &itemp));
    istr.clear(); istr.str( ydhmsf_str.substr(4,3)); istr >> itemp;
    PTB( SetScalarH5A (grp0, "End Day", H5T_STD_I32LE, (VOIDP) &itemp));
    PTB( SetScalarH5A (grp0, "End Millisec", H5T_STD_I32LE, 
		       (VOIDP) &millisec));

    //    PTB( SetScalarH5A (grp0, "Granule Start Time" , 
    //	       H5T_NATIVE_DOUBLE, (VOIDP) &granuleStart));

    return 0;
  }


  /*----------------------------------------------------------------- */
  /* Write orbit times to L1 file                                     */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::writeOrbitTimes(hid_t grp0, 
				     double orbStart, double orbStop)
  {
    string ydhmsf_str;
    ydhmsf_str = ydhmsf( gpstai2unix(orbStart), 'G');

    PTB( SetScalarH5A (grp0, "Orbit Start Time", H5T_STRING, 
		       (void *) ydhmsf_str.c_str()));

    ydhmsf_str = ydhmsf( gpstai2unix(orbStop), 'G');

    PTB( SetScalarH5A (grp0, "Orbit Stop Time", H5T_STRING, 
		       (void *) ydhmsf_str.c_str()));

    //    PTB( SetScalarH5A (grp0, "Orbit Start Time" , 
    //		       H5T_NATIVE_DOUBLE, (VOIDP) &orbStart));

    //PTB( SetScalarH5A (grp0, "Orbit Stop Time" , 
    //		       H5T_NATIVE_DOUBLE, (VOIDP) &orbStop));

    return 0;
  }


  /*----------------------------------------------------------------- */
  /* Write record to L1 file                                          */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::writel1( char* ptr)
  {
    hsize_t start[6], count[6];

    // Block Attributes
    uint32_t ulng;
    unsigned short usht;
    uint32_t tme, nanotme;
    int32_t nanotme0;
    double tmedbl;

    // tme - seconds from 01/06/1980
    // nanotme - 62.5 nanoseconds time units
    // tai = tme - ((nanotme-2) / 16000000)
    memcpy(&tme, &ptr[4], 4);
    memcpy(&nanotme, &ptr[8], 4);
    tme = SWAP_4(tme);

    // 18-Jun-08 Corrected D*62.5 ns to (D-2)*62.5 ns
    // Adam per email from Craig Cheetham, per CTU HLD
    // nanotme0 must be signed to support nanotme = {0,1}  JMG  02/21/12
    nanotme0 = (int32_t) SWAP_4(nanotme) - 2;
    tmedbl = tme - ( nanotme0 / ((double) 16000000.)) + 0.0000001;
    string ydhmsf_str( ydhmsf( gpstai2unix( tmedbl), 'G'));

    uint32_t icds_blknum;
    memcpy(&icds_blknum, &ptr[ICDS_STS_OFFSET+8], 4);
    icds_blknum = SWAP_4(icds_blknum);

    uint16_t rad_header;
    memcpy(&rad_header, &ptr[RADIOMETER_SCI_OFFSET], 2);
    rad_header = SWAP_2(rad_header);

    uint8_t rad_subframe = (ptr[RADIOMETER_SCI_OFFSET+1] & 0xf) % 4;
    uint8_t atc_subframe = ptr[ATC_TLM_OFFSET+1] & 0x3;

    static uint32_t prev_icds_blknum;
    static uint16_t prev_rad_header;
    // Moved to hdf5_Aquarius class JMG 06/05/12
    //    static uint8_t offset_corr;
    if ( first_write) {
      first_icds_blknum = icds_blknum;
      first_atc_subframe = atc_subframe;

      prev_icds_blknum = icds_blknum;
      prev_rad_header = rad_header;

      offset_corr = rad_subframe;
    }

    // Don't write blks outside of L1A filesize
    if ( (icds_blknum - first_icds_blknum) < 0) {
      cout << "Dropping out of bounds blk: " << 
	icds_blknum - first_icds_blknum << endl;
      prev_icds_blknum = icds_blknum;
      return 0;
    }
    if ( (icds_blknum - first_icds_blknum) >= maxscan) {
      cout << "Dropping out of bounds blk: " << 
	icds_blknum - first_icds_blknum << endl;
      prev_icds_blknum = icds_blknum;
      return 0;
    }

    if ( (icds_blknum - prev_icds_blknum) > 1)
      cout << "Jump in icds_blknum from " << prev_icds_blknum << " to "
	   << icds_blknum << endl;

    blknum = icds_blknum - first_icds_blknum;

    // Compute radiometer frame number
    if (prev_rad_header == 0 && rad_header != 0) offset_corr = blknum % 4;
    rad_frmnum = (blknum + offset_corr - rad_subframe) / 4;

    // Compute ATC frame number
    atc_frmnum = (blknum + first_atc_subframe) / 4;

    start[0] = blknum;
    count[0] = 1;

    // blk_sec is seconds from start of day
    // Secs since start of day (mod 86400)
    double blk_sec = get_millisec( &ydhmsf_str) / 1000.;
    PTB( h5d_write(grp1, "blk_sec", &blk_sec, 1, start, count));
    PTB( h5d_write(grp1, "atc_frmnum", &atc_frmnum, 1, start, count));
    PTB( h5d_write(grp1, "atc_subframe", &atc_subframe, 1, start, count));
    PTB( h5d_write(grp1, "rad_frmnum", &rad_frmnum, 1, start, count));
    PTB( h5d_write(grp1, "rad_subframe", &rad_subframe, 1, start, count));

    // Synch and Time Tag
    PTB( h5d_write(grp2, "gps_time_tag", &tme, 1, start, count));
    PTB( h5d_write(grp2, "time_tag_offset", &nanotme0, 1, start, count));

    memcpy(&ulng, &ptr[0], 4); ulng = SWAP_4(ulng);
    PTB( h5d_write(grp2, "start_synch", &ulng, 1, start, count));

    memcpy(&usht, &ptr[2558], 2); usht = SWAP_2(usht);
    PTB( h5d_write(grp2, "checksum", &usht, 1, start, count));

    // Aquarius Raw Telemetry
    start[1] = 0;
    count[1] = 5;
    PTB( h5d_write(grp2, "deploy_tlm", &ptr[12], 2, start, count));
    count[1] = 12;
    PTB( h5d_write(grp2, "icds_status", &ptr[17], 2, start, count));
    count[1] = 24;
    PTB( h5d_write(grp2, "icds_tlm", &ptr[ICDS_TLM_OFFSET], 2, start, count));
    count[1] = 70;
    PTB( h5d_write(grp2, "temp_tlm", &ptr[EXTT_TLM_OFFSET], 2, start, count));
    count[1] = 5;
    PTB( h5d_write(grp2, "apdu_tlm", &ptr[APDU_TLM_OFFSET], 2, start, count));
    count[1] = 36;
    PTB( h5d_write(grp2, "atc_tlm", &ptr[ATC_TLM_OFFSET], 2, start, count));


    ////////////////////////////////////////////////
    //////// Begin Radiometer Raw Telemetry //////// 
    ////////////////////////////////////////////////
    if ( rad_header != 0) {
      // Radiometer Raw Telemetry
      start[0] = rad_frmnum;
      count[0]  = 1;
      start[1] = rad_subframe;
      count[1]  = 1;

      start[2] = 0;
      count[2] = 50;
      PTB( h5d_write(grp2, "radiom_rt_tlm", &ptr[RADIOMETER_SCI_OFFSET+2], 
		     3, start, count));
      count[2] = 10;
      PTB( h5d_write(grp2, "radiom_nrt_tlm", &ptr[RADIOMETER_SCI_OFFSET+2+50], 
		     3, start, count));


      // Radiometer header
      memcpy(&usht, &ptr[RADIOMETER_SCI_OFFSET], 2); usht = SWAP_2(usht);
      PTB( h5d_write(grp2, "radiom_header", &usht, 2, start, count));


      // Swap accum short ints  
      uint32_t naccumbytes = NUMBER_OF_BEAMS*RADIOMETER_POLARIZATIONS*
	((RADIOMETER_SIGNALS_PER_SUBCYCLE+1)*RADIOMETER_SUBCYCLES +
	 RADIOMETER_LONG_ACCUM);

      for (size_t i=0; i<naccumbytes; i++) {
	memcpy(&usht, &ptr[RADIOMETER_ACCUM_OFFSET+2*i], 2);
	usht = SWAP_2(usht);
	//      if (usht > 32767) {
	//cout << "Negative radio count" << endl;
	//cout << i << " " << usht << endl;
	//exit(1);
	//}
	memcpy(&ptr[RADIOMETER_ACCUM_OFFSET+2*i], &usht, 2);
      }

      // Radiometer signals
      start[3] = 0;
      count[3] = RADIOMETER_SIGNALS_PER_SUBCYCLE;
      start[4] = 0;
      count[4] = NUMBER_OF_BEAMS;
      start[5] = 0;
      count[5] = RADIOMETER_POLARIZATIONS;

      for (size_t isubcycl=0; isubcycl<RADIOMETER_SUBCYCLES; isubcycl++) {
	start[2] = isubcycl;
	count[2] = 1;
	PTB( h5d_write(grp2, "radiom_signals", 
		       &ptr[RADIOMETER_ACCUM_OFFSET+24*6*isubcycl], 6, 
		       start, count));
      }


      // Radiometer short accumulations
      start[3] = 0;
      count[3] = NUMBER_OF_BEAMS;
      start[4] = 0;
      count[4] = RADIOMETER_POLARIZATIONS;

      for (size_t isubcycl=0; isubcycl<RADIOMETER_SUBCYCLES; isubcycl++) {
	start[2] = isubcycl;
	count[2] = 1;
	PTB( h5d_write(grp2, "radiom_cnd", 
		       &ptr[RADIOMETER_ACCUM_OFFSET+24*6*isubcycl+24*5], 5, 
		       start, count));
      }


      // Radiometer long accumulations
      start[2] = 0;
      count[2] = RADIOMETER_LONG_ACCUM;
      start[3] = 0;
      count[3] = NUMBER_OF_BEAMS;
      start[4] = 0;
      count[4] = RADIOMETER_POLARIZATIONS;
      PTB( h5d_write(grp2, "radiom_lavg", 
		     &ptr[RADIOMETER_ACCUM_OFFSET+24*6*12], 5, start, count));


      // Un-Swap to recover original big endian buffer
      for (size_t i=0; i<naccumbytes; i++) {
	memcpy(&usht, &ptr[RADIOMETER_ACCUM_OFFSET+2*i], 2); 
	usht = SWAP_2(usht);
	memcpy(&ptr[RADIOMETER_ACCUM_OFFSET+2*i], &usht, 2);
      }
    } // rad_header != 0


    ///////////////////////////////////////////////////
    //////// Begin Scatterometer Raw Telemetry //////// 
    ///////////////////////////////////////////////////

    // scatter_tlm
    start[0] = blknum;
    count[0]  = 1;

    start[1] = 0;
    count[1] = 37;
    PTB( h5d_write(grp2, "scatter_tlm", &ptr[SCATTEROMETER_TLM_OFFSET], 
		   2, start, count));


    // Scatterometer header
    start[1] = 0;
    count[1] = SCATTEROMETER_SUBCYCLES;
    PTB( h5d_write(grp2, "scatter_headers", &ptr[SCATTEROMETER_SCI_OFFSET], 
		   2, start, count));


    uint32_t n_scatter_power = SCATTEROMETER_SUBCYCLES * NUMBER_OF_BEAMS *
      SCATTEROMETER_POLARIZATIONS;
    uint32_t n_scatter_loopback = 
      NUMBER_OF_BEAMS * SCATTEROMETER_POLARIZATIONS;
    uint32_t n_scatter_dc = 2;
    uint32_t n_scatter = n_scatter_power + n_scatter_loopback + n_scatter_dc;

    // Swap scatter bytes
    for (size_t i=0; i<n_scatter; i++) {
      memcpy(&usht, &ptr[SCATTEROMETER_SCI_OFFSET+8+2*i], 2);
      usht = SWAP_2(usht);
      memcpy(&ptr[SCATTEROMETER_SCI_OFFSET+8+2*i], &usht, 2);
    }

    // Scatterometer power
    start[1] = 0;
    count[1] = SCATTEROMETER_SUBCYCLES;
    start[2] = 0;
    count[2] = NUMBER_OF_BEAMS;
    start[3] = 0;
    count[3] = SCATTEROMETER_POLARIZATIONS;
    
    PTB( h5d_write(grp2, "scatter_pwr", &ptr[SCATTEROMETER_SCI_OFFSET+8], 
		   4, start, count));


    // Scatterometer loopback
    start[1] = 0;
    count[1] = NUMBER_OF_BEAMS;
    start[2] = 0;
    count[2] = SCATTEROMETER_POLARIZATIONS;

    PTB( h5d_write(grp2, "scatter_loop", &ptr[SCATTEROMETER_SCI_OFFSET+296], 
		   3, start, count));


    // Scatterometer DC data
    start[1] = 0;
    count[1] = 2;

    PTB( h5d_write(grp2, "scatter_dc", &ptr[SCATTEROMETER_SCI_OFFSET+332], 
		   2, start, count));

    // Un-Swap scatter bytes
    for (size_t i=0; i<n_scatter; i++) {
      memcpy(&usht, &ptr[SCATTEROMETER_SCI_OFFSET+8+2*i], 2);
      usht = SWAP_2(usht);
      memcpy(&ptr[SCATTEROMETER_SCI_OFFSET+8+2*i], &usht, 2);
    }

    // Scatter RFI
    start[1] = 0;
    count[1] = NUMBER_OF_BEAMS+1;
    PTB( h5d_write(grp2, "scatter_rfi", &ptr[2519], 2, start, count));

    // Pad
    start[1] = 0;
    count[1] = 35;
    PTB( h5d_write(grp2, "pad", &ptr[2523], 2, start, count));


    // Converted telemetry
    uint8_t tlm_8;
    unsigned short tlm_16;
    uint32_t tlm_32;
    int32_t s_tlm_32;
    float b2f_08;
    float b2f_10;
    float b2f_12;

    ///////////////////////////////////////////////////
    //////// Begin Converted RAD R/T Telemetry //////// 
    ///////////////////////////////////////////////////

    static float rfe1_analog_tlm[23];
    static float rfe2_analog_tlm[23];
    static float rfe3_analog_tlm[23];
    static float rbe1_analog_tlm[21];
    static float rbe2_analog_tlm[21];
    static float rbe3_analog_tlm[21];
    static float dpu_analog_tlm[12];
    uint8_t dpu_status_tlm[13];
    int16_t nrt_discrete_tlm[17];
    static float extram[6];

    if ( (rad_header & 0xfff0) == 0xfa50) {

      // Swap radiometer R/T HKT
      for (size_t i=0; i<25; i++) {
	memcpy(&usht, &ptr[RADIOMETER_SCI_OFFSET+2+2*i], 2);
	usht = SWAP_2(usht);
	memcpy(&ptr[RADIOMETER_SCI_OFFSET+2+2*i], &usht, 2);
      }

      if (rad_subframe == 0) {

	// DPU_R1P //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+40],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xff);
	dpu_analog_tlm[ 0] = pls_dspl_ohms_12( tlm_16);
	extram[0] = tlm_16 * 250000. / 4096;

	// DPU_R1N //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+42],2);
	tlm_16 = (tlm_16 & 0x0fff) / 1;
	dpu_analog_tlm[ 1] = mns_dspl_ohms_12( tlm_16);
	extram[1] = tlm_16 * 250000. / 4096;

	// DPU_R2P //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+44],2);
	tlm_16 = (tlm_16 & 0xfff0) / 16;
	dpu_analog_tlm[ 2] = pls_dspl_ohms_12( tlm_16);
	extram[2] = tlm_16 * 250000. / 4096;

	// DPU_R2N //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+44],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xf);
	dpu_analog_tlm[ 3] = mns_dspl_ohms_12( tlm_16);
	extram[3] = tlm_16 * 250000. / 4096;

	// DPU_R3 //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+46],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xff);
	dpu_analog_tlm[ 4] = crs_dspl_ohms_12( tlm_16);
	extram[4] = tlm_16 * 250000. / 4096;

	// DPU_R1 //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+48],2);
	tlm_16 = (tlm_16 & 0x0fff) / 1;
	dpu_analog_tlm[ 5] = crs_dspl_ohms_12( tlm_16);
	extram[5] = tlm_16 * 250000. / 4096;

	// DPU_I //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+38],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xf);
	b2f_12 = tlm_16 * 250000 / (1 << 12);
	dpu_analog_tlm[ 11] = b2f_12 * 3.59e-5 - 5.00e-5;


	// RFE1_CBPT //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+14],2);
	tlm_16 = (tlm_16 & 0xffc0) / 64;
	//	cout << crs_dspl_temp_10( tlm_16) << endl;
	rfe1_analog_tlm[ 0] = coarse_temp_10( tlm_16, extram);

	// RFE1_CBT //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+14],4);
	tlm_16 = in32_out16( tlm_32, 10, 0x3f);
	rfe1_analog_tlm[ 1] = coarse_temp_10( tlm_16, extram);

	// RFE1_CNDNDT //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+2],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xf);
	rfe1_analog_tlm[ 2] = 
	  finep_temp_12( tlm_16, extram, 
			 CNST_RFE1CND_A, CNST_RFE_B, CNST_RFE1CND_R0);

	// RFE1_HDL1TN //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+26],2);
	tlm_16 = (tlm_16 & 0xfff0) / 16;
	rfe1_analog_tlm[ 3] = 
	  finen_temp_12( tlm_16, extram,
			 CNST_RFE1HDL_A, CNST_RFE_B, CNST_RFE1HDL_R0);


	// RFE1_HDL1TP //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+24],2);
	tlm_16 = (tlm_16 & 0x0fff) / 1;
	rfe1_analog_tlm[ 4] = 
	  finep_temp_12( tlm_16, extram,
			 CNST_RFE1HDL_A, CNST_RFE_B, CNST_RFE1HDL_R0);

	// RFE1_HDL2T //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+4],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xff);
	rfe1_analog_tlm[ 5] = coarse_temp_12( tlm_16, extram);

	// RFE1_HLNAT //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+10],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xff);
	rfe1_analog_tlm[ 6] = coarse_temp_12( tlm_16, extram);

	// RFE1_VDL1TN //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+22],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xff);
	rfe1_analog_tlm[ 9] = 
	  finen_temp_12( tlm_16, extram,
			 CNST_RFE1VDL_A, CNST_RFE_B, CNST_RFE1VDL_R0);

	// RFE1_VDL1TP //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+20],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xf);
	rfe1_analog_tlm[ 10] = 
	  finep_temp_12( tlm_16, extram,
			 CNST_RFE1VDL_A, CNST_RFE_B, CNST_RFE1VDL_R0);

	// RFE1_VDL2T //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+8],2);
	tlm_16 = (tlm_16 & 0xfff0) / 16;
	rfe1_analog_tlm[ 11] = coarse_temp_12( tlm_16, extram);

	// RFE1_VLNAT //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+12],2);
	tlm_16 = (tlm_16 & 0x0fff) / 1;
	rfe1_analog_tlm[ 12] = coarse_temp_12( tlm_16, extram);

	// RFE1_CNDI //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+2],2);
	tlm_16 = (tlm_16 & 0xfff0) / 16;
	b2f_12 = tlm_16 * 250000 / (1 << 12);
	rfe1_analog_tlm[ 15] = b2f_12 * 2.99902e-5 - 8.79937e-6;

	// RFE1_HNDI //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+6],2);
	tlm_16 = (tlm_16 & 0x0fff) / 1;
	b2f_12 = tlm_16 * 250000 / (1 << 12);
	rfe1_analog_tlm[ 16] = b2f_12 * 3.03703e-5 - 8.91089e-6;

	// RFE1_VNDI //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+8],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xf);
	b2f_12 = tlm_16 * 250000 / (1 << 12);
	rfe1_analog_tlm[ 17] = b2f_12 * 3.03403e-5 - 8.90208e-6;

	// RFE1_12N //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+18],2);
	tlm_16 = (tlm_16 & 0x0ff0) / 16;
	b2f_08 = tlm_16 * 250000 / (1 << 8);
	rfe1_analog_tlm[ 18] = b2f_08 * -7.64366e-5 + 2.24271e-5;

	// RFE1_12P //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+16],4);
	tlm_16 = in32_out16( tlm_32, 8, 0xf);
	b2f_08 = tlm_16 * 250000 / (1 << 8);
	rfe1_analog_tlm[ 19] = b2f_08 * 7.62845e-5 - 2.23825e-5;

	// RFE1_15P //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+16],2);
	tlm_16 = (tlm_16 & 0x0ff0) / 16;
	b2f_08 = tlm_16 * 250000 / (1 << 8);
	rfe1_analog_tlm[ 20] = b2f_08 * 9.58862e-5 - 2.81338e-5;

	// RFE1_5D //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+20],2);
	tlm_16 = (tlm_16 & 0x0ff0) / 16;
	b2f_08 = tlm_16 * 250000 / (1 << 8);
	rfe1_analog_tlm[ 21] = b2f_08 * 3.0674e-5 - 9.00e-6;

	// RFE1_5RF //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+18],4);
	tlm_16 = in32_out16( tlm_32, 8, 0xf);
	b2f_08 = tlm_16 * 250000 / (1 << 8);
	rfe1_analog_tlm[ 22] = b2f_08 * 3.06559e-5 - 8.99469e-6;


	// RFE2_HDL1TN //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+32],2);
	tlm_16 = (tlm_16 & 0xfff0) / 16;
	rfe2_analog_tlm[ 3] = 
	  finen_temp_12( tlm_16, extram, 
			 CNST_RFE2HDL_A, CNST_RFE_B, CNST_RFE2HDL_R0);

	// RFE2_HDL1TP //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+30],2);
	tlm_16 = (tlm_16 & 0x0fff) / 1;
	rfe2_analog_tlm[ 4] = 
	  finep_temp_12( tlm_16, extram,
			 CNST_RFE2HDL_A, CNST_RFE_B, CNST_RFE2HDL_R0);

	// RFE2_VDL1TN //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+28],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xff);
	rfe2_analog_tlm[ 9] = 
	  finen_temp_12( tlm_16, extram,
			 CNST_RFE2VDL_A, CNST_RFE_B, CNST_RFE2VDL_R0);

	// RFE2_VDL1TP //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+26],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xf);
	rfe2_analog_tlm[ 10] = 
	  finep_temp_12( tlm_16, extram,
			 CNST_RFE2VDL_A, CNST_RFE_B, CNST_RFE2VDL_R0);


	// RFE3_HDL1TN //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+38],2);
	tlm_16 = (tlm_16 & 0xfff0) / 16;
	rfe3_analog_tlm[ 3] = 
	  finen_temp_12( tlm_16, extram,
			 CNST_RFE3HDL_A, CNST_RFE_B, CNST_RFE3HDL_R0);

	// RFE3_HDL1TP //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+36],2);
	tlm_16 = (tlm_16 & 0x0fff) / 1;
	rfe3_analog_tlm[ 4] = 
	  finep_temp_12( tlm_16, extram,
			 CNST_RFE3HDL_A, CNST_RFE_B, CNST_RFE3HDL_R0);

	// RFE3_VDL1TN //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+34],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xff);
	rfe3_analog_tlm[ 9] = 
	  finen_temp_12( tlm_16, extram,
			 CNST_RFE3VDL_A, CNST_RFE_B, CNST_RFE3VDL_R0);

	// RFE3_VDL1TP //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+32],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xf);
	rfe3_analog_tlm[ 10] = 
	  finep_temp_12( tlm_16, extram,
			 CNST_RFE3VDL_A, CNST_RFE_B, CNST_RFE3VDL_R0);

      } else if (rad_subframe == 1) {

	// RFE2_CBPT //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+16],2);
	tlm_16 = (tlm_16 & 0x0ffc) / 4;
	rfe2_analog_tlm[ 0] = coarse_temp_10( tlm_16, extram);

	// RFE2_CBT //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+16],4);
	tlm_16 = in32_out16( tlm_32, 10, 0x3);
	rfe2_analog_tlm[ 1] = coarse_temp_10( tlm_16, extram);

	// RFE2_CNDNDT //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+6],2);
	tlm_16 = (tlm_16 & 0xfff0) / 16;
	rfe2_analog_tlm[ 2] = 
	  finep_temp_12( tlm_16, extram,
			 CNST_RFE2CND_A, CNST_RFE_B, CNST_RFE2CND_R0);

	// RFE2_HDL2T //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+6],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xf);
	rfe2_analog_tlm[ 5] = coarse_temp_12( tlm_16, extram);

	// RFE2_HLNAT //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+12],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xf);
	rfe2_analog_tlm[ 6] = coarse_temp_12( tlm_16, extram);

	// RFE2_VDL2T //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+10],2);
	tlm_16 = (tlm_16 & 0x0fff) / 1;
	rfe2_analog_tlm[ 11] = coarse_temp_12( tlm_16, extram);

	// RFE2_VLNAT //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+14],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xff);
	rfe2_analog_tlm[ 12] = coarse_temp_12( tlm_16, extram);

	// RFE2_CNDI //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+4],2);
	tlm_16 = (tlm_16 & 0x0fff) / 1;
	b2f_12 = tlm_16 * 250000 / (1 << 12);
	rfe2_analog_tlm[ 15] = b2f_12 * 2.99259e-5 - 8.78049e-6;

	// RFE2_HNDI //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+8],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xff);
	b2f_12 = tlm_16 * 250000 / (1 << 12);
	rfe2_analog_tlm[ 16] = b2f_12 * 3.03703e-5 - 8.91089e-6;

	// RFE2_VNDI //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+12],2);
	tlm_16 = (tlm_16 & 0xfff0) / 16;
	b2f_12 = tlm_16 * 250000 / (1 << 12);
	rfe2_analog_tlm[ 17] = b2f_12 * 3.03403e-5 - 8.90208e-6;

	// RFE2_12N //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+20],2);
	tlm_16 = (tlm_16 & 0x00ff) / 1;
	b2f_08 = tlm_16 * 250000 / (1 << 8);
	rfe2_analog_tlm[ 18] = b2f_08 * -7.65128e-5 + 2.24495e-5;

	// RFE2_12P //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+20],2);
	tlm_16 = (tlm_16 & 0xff00) / 256;
	b2f_08 = tlm_16 * 250000 / (1 << 8);
	rfe2_analog_tlm[ 19] = b2f_08 * 7.62845e-5 - 2.23825e-5;

	// RFE2_15P //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+18],2);
	tlm_16 = (tlm_16 & 0x00ff) / 1;
	b2f_08 = tlm_16 * 250000 / (1 << 8);
	rfe2_analog_tlm[ 20] = b2f_08 * 9.58263e-5 - 2.81162e-5;

	// RFE2_5D //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+22],2);
	tlm_16 = (tlm_16 & 0x00ff) / 1;
	b2f_08 = tlm_16 * 250000 / (1 << 8);
	rfe2_analog_tlm[ 21] = b2f_08 * 3.0674e-5 - 9.00e-6;

	// RFE2_5RF //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+22],2);
	tlm_16 = (tlm_16 & 0xff00) / 256;
	b2f_08 = tlm_16 * 250000 / (1 << 8);
	rfe2_analog_tlm[ 22] = b2f_08 * 3.06434e-5 - 8.99101e-6;


	// RFE3_CBPT //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+36],2);
	tlm_16 = (tlm_16 & 0xfff0) / 16;
	rfe3_analog_tlm[ 0] = coarse_temp_12( tlm_16, extram);

	// RFE3_CBT //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+36],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xf);
	rfe3_analog_tlm[ 1] = coarse_temp_12( tlm_16, extram);

	// RFE3_CNDNDT //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+24],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xf);
	rfe3_analog_tlm[ 2] = 
	  finep_temp_12( tlm_16, extram,
			 CNST_RFE3CND_A, CNST_RFE_B, CNST_RFE3CND_R0);

	// RFE3_HDL2T //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+26],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xff);
	rfe3_analog_tlm[ 5] = coarse_temp_12( tlm_16, extram);

	// RFE3_HLNAT //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+32],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xff);
	rfe3_analog_tlm[ 6] = coarse_temp_12( tlm_16, extram);

	// RFE3_VDL2T //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+30],2);
	tlm_16 = (tlm_16 & 0xfff0) / 16;
	rfe3_analog_tlm[ 11] = coarse_temp_12( tlm_16, extram);

	// RFE3_VLNAT //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+34],2);
	tlm_16 = (tlm_16 & 0x0fff) / 1;
	rfe3_analog_tlm[ 12] = coarse_temp_12( tlm_16, extram);

	// RFE3_CNDI //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+24],2);
	tlm_16 = (tlm_16 & 0xfff0) / 16;
	b2f_12 = tlm_16 * 250000 / (1 << 12);
	rfe3_analog_tlm[ 15] = b2f_12 * 2.99668e-5 - 8.7925e-6;

	// RFE3_HNDI //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+28],2);
	tlm_16 = (tlm_16 & 0x0fff) / 1;
	b2f_12 = tlm_16 * 250000 / (1 << 12);
	rfe3_analog_tlm[ 16] = b2f_12 * 3.03703E-5 - 8.91089E-6;

	// RFE3_VNDI //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+30],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xf);
	b2f_12 = tlm_16 * 250000 / (1 << 12);
	rfe3_analog_tlm[ 17] = b2f_12 * 3.03403E-5 - 8.90208E-6;

	// RFE3_12N //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+40],2);
	tlm_16 = (tlm_16 & 0x00ff) / 1;
	b2f_08 = tlm_16 * 250000 / (1 << 8);
	rfe3_analog_tlm[ 18] = b2f_08 * -7.64175e-5 + 2.24215e-5;

	// RFE3_12P //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+40],2);
	tlm_16 = (tlm_16 & 0xff00) / 256;
	b2f_08 = tlm_16 * 250000 / (1 << 8);
	rfe3_analog_tlm[ 19] = b2f_08 * 7.62655e-5 - 2.23769e-5;

	// RFE3_15P //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+38],2);
	tlm_16 = (tlm_16 & 0x00ff) / 1;
	b2f_08 = tlm_16 * 250000 / (1 << 8);
	rfe3_analog_tlm[ 20] = b2f_08 * 9.59162e-5 - 2.81426e-5;

	// RFE3_5D //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+42],2);
	tlm_16 = (tlm_16 & 0x00ff) / 1;
	b2f_08 = tlm_16 * 250000 / (1 << 8);
	rfe3_analog_tlm[ 21] = b2f_08 * 3.06801e-5 - 9.0018e-6;

	// RFE3_5RF //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+42],2);
	tlm_16 = (tlm_16 & 0xff00) / 256;
	b2f_08 = tlm_16 * 250000 / (1 << 8);
	rfe3_analog_tlm[ 22] = b2f_08 * 3.06128e-5 - 8.98204e-6;


	// RBE1_D_HT //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+48],2);
	tlm_16 = (tlm_16 & 0x0fff) / 1;
	rbe1_analog_tlm[ 1] = coarse_temp_12( tlm_16, extram);

	// RBE1_D_MT //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+46],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xff);
	rbe1_analog_tlm[ 2] = coarse_temp_12( tlm_16, extram);

	// RBE1_D_PT //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+44],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xf);
	rbe1_analog_tlm[ 3] = coarse_temp_12( tlm_16, extram);

	// RBE1_D_VT //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+44],2);
	tlm_16 = (tlm_16 & 0xfff0) / 16;
	rbe1_analog_tlm[ 4] = coarse_temp_12( tlm_16, extram);


	// RBE1_5RF //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+2],4);
	tlm_16 = in32_out16( tlm_32, 8, 0xf);
	b2f_08 = tlm_16 * 250000 / (1 << 8);
	rbe1_analog_tlm[ 19] = b2f_08 * 3.06495e-5 - 8.99281e-6;


	// DPU_5P //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+2],2);
	tlm_16 = (tlm_16 & 0xfff0) / 16;
	b2f_12 = tlm_16 * 250000 / (1 << 12);
	dpu_analog_tlm[ 10] = b2f_12 * 3.59e-5 - 5.00e-5;

      } else if (rad_subframe == 2) {

	// RFE1_HNDTN //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+6],2);
	tlm_16 = (tlm_16 & 0x0fff) / 1;
	rfe1_analog_tlm[ 7] = 
	  finen_temp_12( tlm_16, extram,
			 CNST_RFE1HND_A, CNST_RFE_B, CNST_RFE1HND_R0);

	// RFE1_HNDTP //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+4],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xff);
	rfe1_analog_tlm[ 8] = 
	  finep_temp_12( tlm_16, extram,
			 CNST_RFE1HND_A, CNST_RFE_B, CNST_RFE1HND_R0);

	// RFE1_VNDTN //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+2],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xf);
	rfe1_analog_tlm[ 13] = 
	  finen_temp_12( tlm_16, extram,
			 CNST_RFE1VND_A, CNST_RFE_B, CNST_RFE1VND_R0);

	// RFE1_VNDTP //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+2],2);
	tlm_16 = (tlm_16 & 0xfff0) / 16;
	rfe1_analog_tlm[ 14] = 
	  finep_temp_12( tlm_16, extram,
			 CNST_RFE1VND_A, CNST_RFE_B, CNST_RFE1VND_R0);

	// RFE2_HNDTN //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+12],2);
	tlm_16 = (tlm_16 & 0x0fff) / 1;
	rfe2_analog_tlm[ 7] = 
	  finen_temp_12( tlm_16, extram,
			 CNST_RFE2HND_A, CNST_RFE_B, CNST_RFE2HND_R0);

	// RFE2_HNDTP //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+10],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xff);
	rfe2_analog_tlm[ 8] = 
	  finep_temp_12( tlm_16, extram,
			 CNST_RFE2HND_A, CNST_RFE_B, CNST_RFE2HND_R0);

	// RFE2_VNDTN //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+8],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xf);
	rfe2_analog_tlm[ 13] = 
	  finen_temp_12( tlm_16, extram,
			 CNST_RFE2VND_A, CNST_RFE_B, CNST_RFE2VND_R0);

	// RFE2_VNDTP //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+8],2);
	tlm_16 = (tlm_16 & 0xfff0) / 16;
	rfe2_analog_tlm[ 14] = 
	  finep_temp_12( tlm_16, extram,
			 CNST_RFE2VND_A, CNST_RFE_B, CNST_RFE2VND_R0);


	// RFE3_HNDTN //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+18],2);
	tlm_16 = (tlm_16 & 0x0fff) / 1;
	rfe3_analog_tlm[ 7] = 
	  finen_temp_12( tlm_16, extram,
			 CNST_RFE3HND_A, CNST_RFE_B, CNST_RFE3HND_R0);

	// RFE3_HNDTP //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+16],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xff);
	rfe3_analog_tlm[ 8] = 
	  finep_temp_12( tlm_16, extram,
			 CNST_RFE3HND_A, CNST_RFE_B, CNST_RFE3HND_R0);

	// RFE3_VNDTN //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+14],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xf);
	rfe3_analog_tlm[ 13] = 
	  finen_temp_12( tlm_16, extram,
			 CNST_RFE3VND_A, CNST_RFE_B, CNST_RFE3VND_R0);

	// RFE3_VNDTP //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+14],2);
	tlm_16 = (tlm_16 & 0xfff0) / 16;
	rfe3_analog_tlm[ 14] = 
	  finep_temp_12( tlm_16, extram,
			 CNST_RFE3VND_A, CNST_RFE_B, CNST_RFE3VND_R0);


	// RBE1_CBT //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+36],4);
	tlm_16 = in32_out16( tlm_32, 10, 0xff);
	rbe1_analog_tlm[ 0] = coarse_temp_10( tlm_16, extram);

	// RBE1_LNA_H1T //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+34],2);
	tlm_16 = (tlm_16 & 0xfff0) / 16;
	rbe1_analog_tlm[ 5] = coarse_temp_12( tlm_16, extram);

	// RBE1_LNA_H2T //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+34],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xf);
	rbe1_analog_tlm[ 6] = coarse_temp_12( tlm_16, extram);

	// RBE1_LNA_M1T //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+30],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xff);
	rbe1_analog_tlm[ 7] = coarse_temp_12( tlm_16, extram);

	// RBE1_LNA_M2T //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+32],2);
	tlm_16 = (tlm_16 & 0x0fff) / 1;
	rbe1_analog_tlm[ 8] = coarse_temp_12( tlm_16, extram);

	// RBE1_LNA_P1T //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+28],2);
	tlm_16 = (tlm_16 & 0xfff0) / 16;
	rbe1_analog_tlm[ 9] = coarse_temp_12( tlm_16, extram);

	// RBE1_LNA_P2T //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+28],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xf);
	rbe1_analog_tlm[ 10] = coarse_temp_12( tlm_16, extram);

	// RBE1_LNA_V1T //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+24],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xff);
	rbe1_analog_tlm[ 11] = coarse_temp_12( tlm_16, extram);

	// RBE1_LNA_V2T //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+26],2);
	tlm_16 = (tlm_16 & 0x0fff) / 1;
	rbe1_analog_tlm[ 12] = coarse_temp_12( tlm_16, extram);

	// RBE1_VCBT //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+38],2);
	tlm_16 = (tlm_16 & 0x3ff0) / 16;
	rbe1_analog_tlm[ 13] = coarse_temp_10( tlm_16, extram);

	// RBE1_12N //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+22],2);
	tlm_16 = (tlm_16 & 0x00ff) / 1;
	b2f_08 = tlm_16 * 250000 / (1 << 8);
	rbe1_analog_tlm[ 14] = b2f_08 * -7.65893e-5 + 2.24719e-5;

	// RBE1_12P //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+22],2);
	tlm_16 = (tlm_16 & 0xff00) / 256;
	b2f_08 = tlm_16 * 250000 / (1 << 8);
	rbe1_analog_tlm[ 15] = b2f_08 * 7.63035e-5 - 2.23881e-5;

	// RBE1_15N //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+20],2);
	tlm_16 = (tlm_16 & 0x00ff) / 1;
	b2f_08 = tlm_16 * 250000 / (1 << 8);
	rbe1_analog_tlm[ 16] = b2f_08 * -9.50542e-5 + 2.78897e-5;

	// RBE1_15P //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+20],2);
	tlm_16 = (tlm_16 & 0xff00) / 256;
	b2f_08 = tlm_16 * 250000 / (1 << 8);
	rbe1_analog_tlm[ 17] = b2f_08 * 9.60664e-5 - 2.81867e-5;

	// RBE1_5 //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+24],2);
	tlm_16 = (tlm_16 & 0xff00) / 256;
	b2f_08 = tlm_16 * 250000 / (1 << 8);
	rbe1_analog_tlm[ 18] = b2f_08 * 3.06372e-5 - 8.98921e-6;


	// RBE2_12N //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+42],2);
	tlm_16 = (tlm_16 & 0x0ff0) / 16;
	b2f_08 = tlm_16 * 250000 / (1 << 8);
	rbe2_analog_tlm[ 14] = b2f_08 * -7.62466e-5 + 2.23714e-5;

	// RBE2_12P //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+40],4);
	tlm_16 = in32_out16( tlm_32, 8, 0xf);
	b2f_08 = tlm_16 * 250000 / (1 << 8);
	rbe2_analog_tlm[ 15] = b2f_08 * 7.63225e-5 - 2.23936e-5;

	// RBE2_15N //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+40],2);
	tlm_16 = (tlm_16 & 0x0ff0) / 16;
	b2f_08 = tlm_16 * 250000 / (1 << 8);
	rbe2_analog_tlm[ 16] = b2f_08 * -9.46144e-5 + 2.77606e-5;

	// RBE2_15P //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+38],4);
	tlm_16 = in32_out16( tlm_32, 8, 0xf);
	b2f_08 = tlm_16 * 250000 / (1 << 8);
	rbe2_analog_tlm[ 17] = b2f_08 * 9.57665e-5 - 2.80987e-5;

	// RBE2_5RF //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+48],2);
	tlm_16 = (tlm_16 & 0x00ff) / 1;
	b2f_08 = tlm_16 * 250000 / (1 << 8);
	rbe2_analog_tlm[ 19] = b2f_08 * 3.06495e-5 - 8.99281e-6;


	// DPU_DT //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+46],2);
	tlm_16 = (tlm_16 & 0xfff0) / 16;
	b2f_12 = tlm_16 * 250000 / (1 << 12);
	dpu_analog_tlm[ 6] = 
	  (sqrt(powf(0.00393858, 2) - 4*(-5.8755e-7) * 
		(1-(b2f_12 * -4.3121e-3 + 2.5859e+3) / 1990.88326) ) 
	   - 0.00393858) / 2 / (-5.8755e-7);

	// DPU_AT2 //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+46],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xf);
	b2f_12 = tlm_16 * 250000 / (1 << 12);
	dpu_analog_tlm[ 7] =  coarse_temp_12( tlm_16, extram);

	// DPU_12N //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+42],4);
	tlm_16 = in32_out16( tlm_32, 10, 0xf);
	b2f_10 = tlm_16 * 250000 / (1 << 10);
	dpu_analog_tlm[ 8] = b2f_10 * -7.24398e-5 + 1.01e-4;

	// DPU_12P //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+44],2);
	tlm_16 = (tlm_16 & 0x03ff) / 1;
	b2f_10 = tlm_16 * 250000 / (1 << 10);
	dpu_analog_tlm[ 9] = b2f_10 * 7.19067e-5 - 1.00e-4;

      } else if (rad_subframe == 3) {

	// RBE2_CBT //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+14],4);
	tlm_16 = in32_out16( tlm_32, 10, 0xff);
	rbe2_analog_tlm[ 0] = coarse_temp_10( tlm_16, extram);

	// RBE2_D_HT //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+40],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xf);
	rbe2_analog_tlm[ 1] = coarse_temp_12( tlm_16, extram);

	// RBE2_D_MT //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+40],2);
	tlm_16 = (tlm_16 & 0xfff0) / 16;
	rbe2_analog_tlm[ 2] = coarse_temp_12( tlm_16, extram);

	// RBE2_D_PT //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+38],2);
	tlm_16 = (tlm_16 & 0x0fff) / 1;
	rbe2_analog_tlm[ 3] = coarse_temp_12( tlm_16, extram);

	// RBE2_D_VT //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+36],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xff);
	rbe2_analog_tlm[ 4] = coarse_temp_12( tlm_16, extram);

	// RBE2_LNA_H1T //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+12],2);
	tlm_16 = (tlm_16 & 0xfff0) / 16;
	rbe2_analog_tlm[ 5] = coarse_temp_12( tlm_16, extram);

	// RBE2_LNA_H2T //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+12],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xf);
	rbe2_analog_tlm[ 6] = coarse_temp_12( tlm_16, extram);

	// RBE2_LNA_M1T //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+8],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xff);
	rbe2_analog_tlm[ 7] = coarse_temp_12( tlm_16, extram);

	// RBE2_LNA_M2T //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+10],2);
	tlm_16 = (tlm_16 & 0x0fff) / 1;
	rbe2_analog_tlm[ 8] = coarse_temp_12( tlm_16, extram);

	// RBE2_LNA_P1T //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+6],2);
	tlm_16 = (tlm_16 & 0xfff0) / 16;
	rbe2_analog_tlm[ 9] = coarse_temp_12( tlm_16, extram);

	// RBE2_LNA_P2T //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+6],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xf);
	rbe2_analog_tlm[ 10] = coarse_temp_12( tlm_16, extram);

	// RBE2_LNA_V1T //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+2],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xff);
	rbe2_analog_tlm[ 11] = coarse_temp_12( tlm_16, extram);

	// RBE2_LNA_V2T //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+4],2);
	tlm_16 = (tlm_16 & 0x0fff) / 1;
	rbe2_analog_tlm[ 12] = coarse_temp_12( tlm_16, extram);

	// RBE2_VCBT //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+16],2);
	tlm_16 = (tlm_16 & 0x3ff0) / 16;
	rbe2_analog_tlm[ 13] = coarse_temp_10( tlm_16, extram);

	// RBE2_5 //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+2],2);
	tlm_16 = (tlm_16 & 0xff00) / 256;
	b2f_08 = tlm_16 * 250000 / (1 << 8);
	rbe2_analog_tlm[ 18] = b2f_08 * 3.0625e-5 - 8.98562e-6;


	// RBE3_CBT //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+34],2);
	tlm_16 = (tlm_16 & 0x0ffc) / 4;
	rbe3_analog_tlm[ 0] = coarse_temp_10( tlm_16, extram);

	// RBE3_D_HT //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+46],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xf);
	rbe3_analog_tlm[ 1] = coarse_temp_12( tlm_16, extram);

	// RBE3_D_MT //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+46],2);
	tlm_16 = (tlm_16 & 0xfff0) / 16;
	rbe3_analog_tlm[ 2] = coarse_temp_12( tlm_16, extram);

	// RBE3_D_PT //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+44],2);
	tlm_16 = (tlm_16 & 0x0fff) / 1;
	rbe3_analog_tlm[ 3] = coarse_temp_12( tlm_16, extram);

	// RBE3_D_VT //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+42],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xff);
	rbe3_analog_tlm[ 4] = coarse_temp_12( tlm_16, extram);

	// RBE3_LNA_H1T //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+30],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xf);
	rbe3_analog_tlm[ 5] = coarse_temp_12( tlm_16, extram);

	// RBE3_LNA_H2T //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+32],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xff);
	rbe3_analog_tlm[ 6] = coarse_temp_12( tlm_16, extram);

	// RBE3_LNA_M1T //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+28],2);
	tlm_16 = (tlm_16 & 0x0fff) / 1;
	rbe3_analog_tlm[ 7] = coarse_temp_12( tlm_16, extram);

	// RBE3_LNA_M2T //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+30],2);
	tlm_16 = (tlm_16 & 0xfff0) / 16;
	rbe3_analog_tlm[ 8] = coarse_temp_12( tlm_16, extram);

	// RBE3_LNA_P1T //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+24],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xf);
	rbe3_analog_tlm[ 9] = coarse_temp_12( tlm_16, extram);

	// RBE3_LNA_P2T //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+26],4);
	tlm_16 = in32_out16( tlm_32, 12, 0xff);
	rbe3_analog_tlm[ 10] = coarse_temp_12( tlm_16, extram);

	// RBE3_LNA_V1T //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+22],2);
	tlm_16 = (tlm_16 & 0x0fff) / 1;
	rbe3_analog_tlm[ 11] = coarse_temp_12( tlm_16, extram);

	// RBE3_LNA_V2T //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+24],2);
	tlm_16 = (tlm_16 & 0xfff0) / 16;
	rbe3_analog_tlm[ 12] = coarse_temp_12( tlm_16, extram);

	// RBE3_VCBT //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+34],4);
	tlm_16 = in32_out16( tlm_32, 10, 0x3);
	rbe3_analog_tlm[ 13] = coarse_temp_10( tlm_16, extram);

	// RBE3_12N //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+20],2);
	tlm_16 = (tlm_16 & 0x0ff0) / 16;
	b2f_08 = tlm_16 * 250000 / (1 << 8);
	rbe3_analog_tlm[ 14] = b2f_08 * -7.63985e-5 + 2.24159e-5;

	// RBE3_12P //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+18],4);
	tlm_16 = in32_out16( tlm_32, 8, 0xf);
	b2f_08 = tlm_16 * 250000 / (1 << 8);
	rbe3_analog_tlm[ 15] = b2f_08 * 7.62466e-5 - 2.23714e-5;

	// RBE3_15N //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+18],2);
	tlm_16 = (tlm_16 & 0x0ff0) / 16;
	b2f_08 = tlm_16 * 250000 / (1 << 8);
	rbe3_analog_tlm[ 16] = b2f_08 * -9.46144e-5 + 2.77606e-5;

	// RBE3_15P //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+16],4);
	tlm_16 = in32_out16( tlm_32, 8, 0xf);
	b2f_08 = tlm_16 * 250000 / (1 << 8);
	rbe3_analog_tlm[ 17] = b2f_08 * 9.57366e-5 - 2.80899e-5;

	// RBE3_5 //
	memcpy(&tlm_32, &ptr[RADIOMETER_SCI_OFFSET+2+20],4);
	tlm_16 = in32_out16( tlm_32, 8, 0xf);
	b2f_08 = tlm_16 * 250000 / (1 << 8);
	rbe3_analog_tlm[ 18] = b2f_08 * 3.06189e-5 - 8.98383e-6;


	// RBE3_5RF //
	memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+48],2);
	tlm_16 = (tlm_16 & 0x00ff) / 1;
	b2f_08 = tlm_16 * 250000 / (1 << 8);
	rbe3_analog_tlm[ 19] = b2f_08 * 3.06311e-5 - 8.98742e-6;

      }

      // Un-Swap radiometer R/T HKT
      for (size_t i=0; i<25; i++) {
	memcpy(&usht, &ptr[RADIOMETER_SCI_OFFSET+2+2*i], 2);
	usht = SWAP_2(usht);
	memcpy(&ptr[RADIOMETER_SCI_OFFSET+2+2*i], &usht, 2);
      }


      // DPU Status Telemetry
      memcpy(&tlm_8, &ptr[RADIOMETER_SCI_OFFSET+2], 1);
 
      // RD_STS_FIFO_ERR
      dpu_status_tlm[0] = (tlm_8 >> 7) & 0x1;

      // RD_STS_FIFO_ERR_VAL
      dpu_status_tlm[1] = (tlm_8 >> 6) & 0x1;

      // RD_STS_SPARE
      dpu_status_tlm[12] = (tlm_8 >> 5) & 0x1;

      // RD_STS_MEM_LDG
      dpu_status_tlm[2] = (tlm_8 >> 4) & 0x1;

      // RD_STS_INTEG_ERR
      dpu_status_tlm[3] = (tlm_8 >> 3) & 0x1;

      // RD_STS_RP_LOST_ERR
      dpu_status_tlm[4] = (tlm_8 >> 2) & 0x1;

      // RD_STS_FS_LOST_ERR
      dpu_status_tlm[5] = (tlm_8 >> 1) & 0x1;

      // RD_STS_CMD_ERR
      dpu_status_tlm[6] = tlm_8 & 0x1;


      memcpy(&tlm_8, &ptr[RADIOMETER_SCI_OFFSET+2+1], 1);

      // RD_STS_GOT_CMD
      dpu_status_tlm[7] = (tlm_8 >> 7) & 0x1;

      // RD_STS_UP_CHK_ERR
      dpu_status_tlm[8] = (tlm_8 >> 6) & 0x1;

      // RD_STS_OP_CHK_ERR
      dpu_status_tlm[9] = (tlm_8 >> 5) & 0x1;

      // RD_STS_OPLUT
      dpu_status_tlm[10] = (tlm_8 >> 2) & 0x7;

      // RD_STS_BLK_CNT
      dpu_status_tlm[11] = tlm_8 & 0x3;

    


      // NRT Discrete Converted Telemetry
      memcpy(&tlm_8, &ptr[RADIOMETER_SCI_OFFSET+2+50+0], 1);

      // RD_SPARE1
      nrt_discrete_tlm[15] = (tlm_8 >> 7) & 0x1;

      // RD_VFA_OVRFLO
      nrt_discrete_tlm[0] = (tlm_8 >> 6) & 0x1;

      // RD_EEPROM_SEL
      nrt_discrete_tlm[1] = (tlm_8 >> 5) & 0x1;

      // RD_UPLUTNO
      nrt_discrete_tlm[2] = (tlm_8 >> 2) & 0x7;

      // RD_LAST_RESET
      nrt_discrete_tlm[3] = tlm_8 & 0x3;


      memcpy(&tlm_8, &ptr[RADIOMETER_SCI_OFFSET+2+50+1], 1);

      // RD_HRCNT
      nrt_discrete_tlm[4] = (tlm_8 >> 4) & 0xf;

      // RD_CMDERR
      nrt_discrete_tlm[5] = tlm_8 & 0xf;


      // RD_CMDCTR
      memcpy(&tlm_8, &ptr[RADIOMETER_SCI_OFFSET+2+50+2], 1);
      nrt_discrete_tlm[6] = tlm_8;

      // RD_BLKCNT
      memcpy(&tlm_8, &ptr[RADIOMETER_SCI_OFFSET+2+50+3], 1);
      nrt_discrete_tlm[7] = tlm_8;


      memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+50+4], 2);
      tlm_16 = SWAP_2(tlm_16);

      // RD_VF1_ERR
      nrt_discrete_tlm[8] = (tlm_16 >> 13) & 0x1;

      // RD_VF2_ERR
      nrt_discrete_tlm[9] = (tlm_16 >> 14) & 0x1;

      // RD_VF3_ERR
      nrt_discrete_tlm[10] = (tlm_16 >> 15) & 0x1;

      // RD_INTEG_STRT
      nrt_discrete_tlm[11] = (tlm_16 >> 7) & 0x3f;

      // RD_INTEG_END
      nrt_discrete_tlm[12] = (tlm_16 & 0x7f);


      memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+50+6], 2);
      tlm_16 = SWAP_2(tlm_16);

      // RD_SPARE2
      nrt_discrete_tlm[16] = tlm_16 >> 15;

      // RD_DUMPADDR
      nrt_discrete_tlm[13] = (tlm_16 & 0x7fff);


      // RD_DUMPDATA
      memcpy(&tlm_16, &ptr[RADIOMETER_SCI_OFFSET+2+50+8], 2);
      tlm_16 = SWAP_2(tlm_16);
      nrt_discrete_tlm[14] = tlm_16;

    } // (rad_header & 0xfff0) == 0xfa50


    ////////////////////////////////////////
    //////// Begin Deploy Telemetry //////// 
    ////////////////////////////////////////
    float deploy_analog_tlm[4];

    // DP_UP_MEC_TMP_PRI_CAL
    memcpy(&tlm_8, &ptr[DPLY_TLM_OFFSET+1],1);
    deploy_analog_tlm[0] = temp_8bit_ext_abr0( tlm_8,
					       CNST_DP_A,
					       CNST_DP_B,
					       CNST_UP_MEC_PRI_R0);

    // DP_UP_MEC_TMP_RED_CAL
    memcpy(&tlm_8, &ptr[DPLY_TLM_OFFSET+2],1);
    deploy_analog_tlm[1] = temp_8bit_ext_abr0( tlm_8,
					       CNST_DP_A,
					       CNST_DP_B,
					       CNST_UP_MEC_RED_R0);

    // DP_LO_MEC_TMP_PRI_CAL
    memcpy(&tlm_8, &ptr[DPLY_TLM_OFFSET+3],1);
    deploy_analog_tlm[2] = temp_8bit_ext_abr0( tlm_8,
					       CNST_DP_A,
					       CNST_DP_B,
					       CNST_LO_MEC_PRI_R0);

    // DP_LO_MEC_TMP_RED_CAL
    memcpy(&tlm_8, &ptr[DPLY_TLM_OFFSET+4],1);
    deploy_analog_tlm[3] = temp_8bit_ext_abr0( tlm_8,
					       CNST_DP_A,
					       CNST_DP_B,
					       CNST_LO_MEC_RED_R0);


    int32_t deploy_discrete_tlm[7];
    // DP_SEP1_STW
    // DP_SEP2_STW
    // DP_LDM_STW
    // DP_UDM_PRI_LATCH
    // DP_UDM_RED_DEPL
    // DP_LDM_PRI_LATCH
    // DP_LDM_RED_DEPL
    memcpy(&tlm_8, &ptr[DPLY_TLM_OFFSET+0],1);
    for (size_t i=0; i<7; i++)
      deploy_discrete_tlm[6-i] = (tlm_8 >> i) & 0x1;


    //////////////////////////////////////
    //////// Begin ICDS Telemetry //////// 
    //////////////////////////////////////

    // ICDS analog telemetry
    float icds_analog_tlm[14];

    // Swap ICDS HKT
    for (size_t i=0; i<12; i++) {
      memcpy(&usht, &ptr[ICDS_TLM_OFFSET+2*i], 2);
      usht = SWAP_2(usht);
      memcpy(&ptr[ICDS_TLM_OFFSET+2*i], &usht, 2);
    }

    // IE_ICDS_P5V
    memcpy(&tlm_8, &ptr[ICDS_TLM_OFFSET+1], 1);
    icds_analog_tlm[0] = ie_volt( tlm_8) * 2;

    // IE_ICDS_P15V
    memcpy(&tlm_8, &ptr[ICDS_TLM_OFFSET+0], 1);
    icds_analog_tlm[1] = ie_volt( tlm_8) * 6;

    // IE_ICDS_N15V
    memcpy(&tlm_8, &ptr[ICDS_TLM_OFFSET+3], 1);
    icds_analog_tlm[2] = ie_volt( tlm_8) * -5;

    // IE_SCI_ADC_5VA
    memcpy(&tlm_8, &ptr[ICDS_TLM_OFFSET+2], 1);
    icds_analog_tlm[3] = ie_volt( tlm_8) * 2;

    // IE_SCI_ADC_VREF
    memcpy(&tlm_8, &ptr[ICDS_TLM_OFFSET+5], 1);
    icds_analog_tlm[4] = ie_volt( tlm_8) * 2;

    // IE_SCI_ADC_TMP
    memcpy(&tlm_8, &ptr[ICDS_TLM_OFFSET+4], 1);
    icds_analog_tlm[5] = temp_8bit_ext_fine( tlm_8,
					     CNST_ICDS_A,
					     CNST_ICDS_B,
					     CNST_SCI_ADC_R0);

    // IE_ICDS_TMP_CHASS
    memcpy(&tlm_8, &ptr[ICDS_TLM_OFFSET+7], 1);
    icds_analog_tlm[6] = temp_8bit_ext_fine( tlm_8,
					     CNST_ICDS_A,
					     CNST_ICDS_B,
					     CNST_ICDS_R0_CHASS);
    
    // IE_ICDS_TMP_RAD6K
    memcpy(&tlm_8, &ptr[ICDS_TLM_OFFSET+6], 1);
    icds_analog_tlm[7] = temp_8bit_rad6k( tlm_8);

    // IE_ICDS_TMP_ATC1
    memcpy(&tlm_8, &ptr[ICDS_TLM_OFFSET+9], 1);
    icds_analog_tlm[8] = temp_8bit_norm_fine( tlm_8,
					     CNST_ICDS_A,
					     CNST_ICDS_B,
					     CNST_ICDS_R0_ATC1);

    //IE_ICDS_TMP_ATC2
    memcpy(&tlm_8, &ptr[ICDS_TLM_OFFSET+8], 1);
    icds_analog_tlm[9] = temp_8bit_norm_fine( tlm_8,
					     CNST_ICDS_A,
					     CNST_ICDS_B,
					     CNST_ICDS_R0_ATC2);

    // Define extram for ET conversions
    uint16_t extram_et[4];

    // IE_ICDS_TLM_CAL_RES1
    memcpy(&tlm_16, &ptr[ICDS_TLM_OFFSET+10], 2);
    icds_analog_tlm[10] = resist_16bit_ext_disp( tlm_16);
    extram_et[0] = tlm_16;

    // IE_ICDS_TLM_CAL_RES2
    memcpy(&tlm_16, &ptr[ICDS_TLM_OFFSET+12], 2);
    icds_analog_tlm[11] = resist_16bit_norm_disp( tlm_16);
    extram_et[1] = tlm_16;

    // IE_ICDS_TLM_CAL_RES3
    memcpy(&tlm_16, &ptr[ICDS_TLM_OFFSET+14], 2);
    icds_analog_tlm[12] = resist_16bit_norm_disp( tlm_16);
    extram_et[2] = tlm_16;

    // IE_ICDS_TLM_CAL_RES4
    memcpy(&tlm_16, &ptr[ICDS_TLM_OFFSET+16], 2);
    icds_analog_tlm[13] = resist_16bit_ext_disp( tlm_16);
    extram_et[3] = tlm_16;


    // ICDS discrete telemetry
    uint32_t icds_discrete_tlm[8];

    // IP_GPS_TIME_LST_BOOT
    memcpy(&tlm_32, &ptr[ICDS_STS_OFFSET], 4);
    icds_discrete_tlm[0] = SWAP_4(tlm_32);

    // IP_CMD_RCVD_CNT
    memcpy(&tlm_16, &ptr[ICDS_STS_OFFSET+4], 2);
    icds_discrete_tlm[1] = SWAP_2(tlm_16);

    // IP_LAST_CMD_RCVD_ID
    memcpy(&tlm_16, &ptr[ICDS_STS_OFFSET+6], 2);
    icds_discrete_tlm[2] = SWAP_2(tlm_16);

    // IP_LST_BLK_NUM
    memcpy(&tlm_32, &ptr[ICDS_STS_OFFSET+8], 4);
    icds_discrete_tlm[3] = SWAP_4(tlm_32);

    // IE_ICDS_ROT_REG_ID
    memcpy(&tlm_8, &ptr[ICDS_TLM_OFFSET+19], 1);
    icds_discrete_tlm[4] = tlm_8;

    // IE_ICDS_ROT_REG_VAL
    memcpy(&tlm_32, &ptr[ICDS_TLM_OFFSET+18], 4);
    tlm_16 = in32_out16( tlm_32, 16, 0xff);
    icds_discrete_tlm[5] = tlm_16;

    // IE_ICDS_SEL_REG_ID
    memcpy(&tlm_8, &ptr[ICDS_TLM_OFFSET+20], 1);
    icds_discrete_tlm[6] = tlm_8;

    // IE_ICDS_SEL_REG_VAL
    memcpy(&tlm_16, &ptr[ICDS_TLM_OFFSET+22], 2);
    icds_discrete_tlm[7] = tlm_16;

    // Un-Swap ICDS HKT
    for (size_t i=0; i<12; i++) {
      memcpy(&usht, &ptr[ICDS_TLM_OFFSET+2*i], 2);
      usht = SWAP_2(usht);
      memcpy(&ptr[ICDS_TLM_OFFSET+2*i], &usht, 2);
    }


    //////////////////////////////////////
    //////// Begin EXTT Telemetry //////// 
    //////////////////////////////////////

    // External Temperature analog telemetry
    float extt_analog_tlm[38];

    // Swap EXTT HKT
    for (size_t i=0; i<35; i++) {
      memcpy(&usht, &ptr[EXTT_TLM_OFFSET+2*i], 2);
      usht = SWAP_2(usht);
      memcpy(&ptr[EXTT_TLM_OFFSET+2*i], &usht, 2);
    }


    memcpy(&tlm_16, &ptr[ICDS_TLM_OFFSET+10], 2);
    //    uint16_t IE_ICDS_TLM_CAL_RES1_DN = SWAP_2( tlm_16);

    memcpy(&tlm_16, &ptr[ICDS_TLM_OFFSET+12], 2);
    uint16_t IE_ICDS_TLM_CAL_RES2_DN = SWAP_2( tlm_16);

    memcpy(&tlm_16, &ptr[ICDS_TLM_OFFSET+14], 2);
    uint16_t IE_ICDS_TLM_CAL_RES3_DN = SWAP_2( tlm_16);

    memcpy(&tlm_16, &ptr[ICDS_TLM_OFFSET+16], 2);
    //uint16_t IE_ICDS_TLM_CAL_RES4_DN = SWAP_2( tlm_16);

    float gain_norm_cmp =
      1.0 * (CNST_CALRES2_NOM-CNST_CALRES3_NOM) /
      (IE_ICDS_TLM_CAL_RES2_DN-IE_ICDS_TLM_CAL_RES3_DN);

    float offset_norm_cmp = 
      CNST_CALRES2_NOM - IE_ICDS_TLM_CAL_RES2_DN * gain_norm_cmp;  

    //    float gain_ext_cmp =
    //1.0 * (CNST_CALRES1_NOM-CNST_CALRES4_NOM) /
    //(IE_ICDS_TLM_CAL_RES1_DN-IE_ICDS_TLM_CAL_RES4_DN);

    //    float offset_ext_cmp = 
    //CNST_CALRES1_NOM - IE_ICDS_TLM_CAL_RES1_DN * gain_ext_cmp;


    // ET_OMT1_H_PROBE_TMP
    memcpy(&tlm_16, &ptr[EXTT_TLM_OFFSET+0], 2);
    /*
    extt_analog_tlm[0] = temp_16bit_norm_disp( tlm_16, 
                                               gain_norm_cmp,
                                               offset_norm_cmp,
                                               CNST_OMT1_H_PROBE_A, 
                                               CNST_OMT1_H_PROBE_B,
                                               CNST_OMT1_H_PROBE_R0);
    */
    extt_analog_tlm[0] = temp_16bit_norm_fine( tlm_16, extram_et, 
					       CNST_OMT1_H_PROBE_A, 
					       CNST_OMT1_H_PROBE_B,
					       CNST_OMT1_H_PROBE_R0);


    // ET_OMT1_V_PROBE_TMP
    memcpy(&tlm_16, &ptr[EXTT_TLM_OFFSET+2], 2);
    extt_analog_tlm[1] = temp_16bit_norm_fine( tlm_16, extram_et,
					       CNST_OMT1_V_PROBE_A, 
					       CNST_OMT1_V_PROBE_B,
					       CNST_OMT1_V_PROBE_R0);
    // ET_OMT2_H_PROBE_TMP
    memcpy(&tlm_16, &ptr[EXTT_TLM_OFFSET+4], 2);
    extt_analog_tlm[2] = temp_16bit_norm_fine( tlm_16, extram_et,
					       CNST_OMT2_H_PROBE_A, 
					       CNST_OMT2_H_PROBE_B,
					       CNST_OMT2_H_PROBE_R0);

    // ET_OMT2_V_PROBE_TMP
    memcpy(&tlm_16, &ptr[EXTT_TLM_OFFSET+6], 2);
    extt_analog_tlm[3] = temp_16bit_norm_fine( tlm_16, extram_et,
					       CNST_OMT2_V_PROBE_A, 
					       CNST_OMT2_V_PROBE_B,
					       CNST_OMT2_V_PROBE_R0);
    // ET_OMT3_H_PROBE_TMP
    memcpy(&tlm_16, &ptr[EXTT_TLM_OFFSET+8], 2);
    extt_analog_tlm[4] = temp_16bit_norm_fine( tlm_16, extram_et,
					       CNST_OMT3_H_PROBE_A, 
					       CNST_OMT3_H_PROBE_B,
					       CNST_OMT3_H_PROBE_R0);
    // ET_OMT3_V_PROBE_TMP
    memcpy(&tlm_16, &ptr[EXTT_TLM_OFFSET+10], 2);
    extt_analog_tlm[5] = temp_16bit_norm_fine( tlm_16, extram_et, 
					       CNST_OMT3_V_PROBE_A, 
					       CNST_OMT3_V_PROBE_B,
					       CNST_OMT3_V_PROBE_R0);

    // ET_SFE_H1_COAX_TMP
    memcpy(&tlm_8, &ptr[EXTT_TLM_OFFSET+13], 1);
    extt_analog_tlm[6] = temp_8bit_norm_abr0( tlm_8, 
					      CNST_SFE_H1_COAX_A, 
					      CNST_SFE_H1_COAX_B,
					      CNST_SFE_H1_COAX_R0);
    // ET_SFE_H2_COAX_TMP
    memcpy(&tlm_8, &ptr[EXTT_TLM_OFFSET+12], 1);
    extt_analog_tlm[7] = temp_8bit_norm_abr0( tlm_8, 
					      CNST_SFE_H2_COAX_A, 
					      CNST_SFE_H2_COAX_B,
					      CNST_SFE_H2_COAX_R0);
    // ET_SFE_H3_COAX_TMP
    memcpy(&tlm_8, &ptr[EXTT_TLM_OFFSET+15], 1);
    extt_analog_tlm[8] = temp_8bit_norm_abr0( tlm_8, 
					      CNST_SFE_H3_COAX_A, 
					      CNST_SFE_H3_COAX_B,
					      CNST_SFE_H3_COAX_R0);
    // ET_SFE_V1_COAX_TMP
    memcpy(&tlm_8, &ptr[EXTT_TLM_OFFSET+14], 1);
    extt_analog_tlm[9] = temp_8bit_norm_abr0( tlm_8, 
					      CNST_SFE_V1_COAX_A, 
					      CNST_SFE_V1_COAX_B,
					      CNST_SFE_V1_COAX_R0);
    // ET_SFE_V2_COAX_TMP
    memcpy(&tlm_8, &ptr[EXTT_TLM_OFFSET+17], 1);
    extt_analog_tlm[10] = temp_8bit_norm_abr0( tlm_8, 
					       CNST_SFE_V2_COAX_A, 
					       CNST_SFE_V2_COAX_B,
					       CNST_SFE_V2_COAX_R0);
    // ET_SFE_V3_COAX_TMP
    memcpy(&tlm_8, &ptr[EXTT_TLM_OFFSET+16], 1);    
    extt_analog_tlm[11] = temp_8bit_norm_abr0( tlm_8, 
					       CNST_SFE_V3_COAX_A, 
					       CNST_SFE_V3_COAX_B,
					       CNST_SFE_V3_COAX_R0);

    // ET_DPLX_1H_TMP
    memcpy(&tlm_16, &ptr[EXTT_TLM_OFFSET+18], 2);
    extt_analog_tlm[12] = temp_16bit_norm_fine( tlm_16, extram_et,
						CNST_DPLX_1H_A, 
						CNST_DPLX_1H_B,
						CNST_DPLX_1H_R0);
    // ET_DPLX_1V_TMP
    memcpy(&tlm_16, &ptr[EXTT_TLM_OFFSET+20], 2);
    extt_analog_tlm[13] = temp_16bit_norm_fine( tlm_16, extram_et,
						CNST_DPLX_1V_A, 
						CNST_DPLX_1V_B,
						CNST_DPLX_1V_R0);
    // ET_DPLX_2H_TMP
    memcpy(&tlm_16, &ptr[EXTT_TLM_OFFSET+22], 2);
    extt_analog_tlm[14] = temp_16bit_norm_fine( tlm_16, extram_et,
						CNST_DPLX_2H_A, 
						CNST_DPLX_2H_B,
						CNST_DPLX_2H_R0);
    // ET_DPLX_2V_TMP
    memcpy(&tlm_16, &ptr[EXTT_TLM_OFFSET+24], 2);
    extt_analog_tlm[15] = temp_16bit_norm_fine( tlm_16, extram_et,
						CNST_DPLX_2V_A, 
						CNST_DPLX_2V_B,
						CNST_DPLX_2V_R0);
    // ET_DPLX_3H_TMP
    memcpy(&tlm_16, &ptr[EXTT_TLM_OFFSET+26], 2);
    extt_analog_tlm[16] = temp_16bit_norm_fine( tlm_16, extram_et,
						CNST_DPLX_3H_A, 
						CNST_DPLX_3H_B,
						CNST_DPLX_3H_R0);
    // ET_DPLX_3V_TMP
    memcpy(&tlm_16, &ptr[EXTT_TLM_OFFSET+28], 2);
    extt_analog_tlm[17] = temp_16bit_norm_fine( tlm_16, extram_et,
						CNST_DPLX_3V_A, 
						CNST_DPLX_3V_B,
						CNST_DPLX_3V_R0);
    // ET_CPLR_1H_TMP
    memcpy(&tlm_16, &ptr[EXTT_TLM_OFFSET+30], 2);
    extt_analog_tlm[18] = temp_16bit_norm_fine( tlm_16, extram_et,
						CNST_CPLR_1H_A, 
						CNST_CPLR_1H_B,
						CNST_CPLR_1H_R0);

    // ET_CPLR_1V_TMP
    memcpy(&tlm_16, &ptr[EXTT_TLM_OFFSET+32], 2);
    extt_analog_tlm[19] = temp_16bit_norm_fine( tlm_16, extram_et,
						CNST_CPLR_1V_A, 
						CNST_CPLR_1V_B,
						CNST_CPLR_1V_R0);
    // ET_CPLR_2H_TMP
    memcpy(&tlm_16, &ptr[EXTT_TLM_OFFSET+34], 2);
    extt_analog_tlm[20] = temp_16bit_norm_fine( tlm_16, extram_et,
						CNST_CPLR_2H_A, 
						CNST_CPLR_2H_B,
						CNST_CPLR_2H_R0);
    // ET_CPLR_2V_TMP
    memcpy(&tlm_16, &ptr[EXTT_TLM_OFFSET+36], 2);
    extt_analog_tlm[21] = temp_16bit_norm_fine( tlm_16, extram_et,
						CNST_CPLR_2V_A, 
						CNST_CPLR_2V_B,
						CNST_CPLR_2V_R0);
    // ET_CPLR_3H_TMP
    memcpy(&tlm_16, &ptr[EXTT_TLM_OFFSET+38], 2);
    extt_analog_tlm[22] = temp_16bit_norm_fine( tlm_16, extram_et, 
						CNST_CPLR_3H_A, 
						CNST_CPLR_3H_B,
						CNST_CPLR_3H_R0);
    // ET_CPLR_3V_TMP
    memcpy(&tlm_16, &ptr[EXTT_TLM_OFFSET+40], 2);
    extt_analog_tlm[23] = temp_16bit_norm_fine( tlm_16, extram_et,
						CNST_CPLR_3V_A, 
						CNST_CPLR_3V_B,
						CNST_CPLR_3V_R0);
    // ET_CND1_TMP
    memcpy(&tlm_16, &ptr[EXTT_TLM_OFFSET+42], 2);
    extt_analog_tlm[24] = temp_16bit_norm_fine( tlm_16, extram_et,
						CNST_CND1_A, 
						CNST_CND1_B,
						CNST_CND1_R0);
    // ET_CND2_TMP
    memcpy(&tlm_16, &ptr[EXTT_TLM_OFFSET+44], 2);
    extt_analog_tlm[25] = temp_16bit_norm_fine( tlm_16, extram_et,
						CNST_CND2_A, 
						CNST_CND2_B,
						CNST_CND2_R0);
    // ET_CND3_TMP
    memcpy(&tlm_16, &ptr[EXTT_TLM_OFFSET+46], 2);
    extt_analog_tlm[26] = temp_16bit_norm_fine( tlm_16, extram_et,
						CNST_CND3_A, 
						CNST_CND3_B,
						CNST_CND3_R0);
    // ET_RFL1_TMP
    memcpy(&tlm_16, &ptr[EXTT_TLM_OFFSET+48], 2);
    extt_analog_tlm[27] = temp_16bit_ext_fine( tlm_16, extram_et,
					       CNST_RFL1_A, 
					       CNST_RFL1_B,
					       CNST_RFL1_R0,
					       CNST_RFL1_H_R);
    // ET_RFL2_TMP
    memcpy(&tlm_16, &ptr[EXTT_TLM_OFFSET+50], 2);
    extt_analog_tlm[28] = temp_16bit_ext_fine( tlm_16, extram_et,
					       CNST_RFL2_A, 
					       CNST_RFL2_B,
					       CNST_RFL2_R0,
					       CNST_RFL2_H_R);
    // ET_RFL3_TMP
    memcpy(&tlm_16, &ptr[EXTT_TLM_OFFSET+52], 2);
    extt_analog_tlm[29] = temp_16bit_ext_fine( tlm_16, extram_et,
					       CNST_RFL3_A, 
					       CNST_RFL3_B,
					       CNST_RFL3_R0,
					       CNST_RFL3_H_R);
    // ET_RFL4_TMP
    memcpy(&tlm_16, &ptr[EXTT_TLM_OFFSET+54], 2);
    extt_analog_tlm[30] = temp_16bit_ext_fine( tlm_16, extram_et,
					       CNST_RFL4_A, 
					       CNST_RFL4_B,
					       CNST_RFL4_R0,
					       CNST_RFL4_H_R);
    // ET_RFL5_TMP
    memcpy(&tlm_16, &ptr[EXTT_TLM_OFFSET+56], 2);
    extt_analog_tlm[31] = temp_16bit_ext_fine( tlm_16, extram_et,
					       CNST_RFL5_A, 
					       CNST_RFL5_B,
					       CNST_RFL5_R0,
					       CNST_RFL5_H_R);
    // ET_RFL6_TMP
    memcpy(&tlm_16, &ptr[EXTT_TLM_OFFSET+58], 2);
    extt_analog_tlm[32] = temp_16bit_ext_fine( tlm_16, extram_et,
					       CNST_RFL6_A, 
					       CNST_RFL6_B,
					       CNST_RFL6_R0,
					       CNST_RFL6_H_R);
    // ET_RFL7_TMP
    memcpy(&tlm_16, &ptr[EXTT_TLM_OFFSET+60], 2);
    extt_analog_tlm[33] = temp_16bit_ext_fine( tlm_16, extram_et,
					       CNST_RFL7_A, 
					       CNST_RFL7_B,
					       CNST_RFL7_R0,
					       CNST_RFL7_H_R);
    // ET_RFL8_TMP
    memcpy(&tlm_16, &ptr[EXTT_TLM_OFFSET+62], 2);
    extt_analog_tlm[34] = temp_16bit_ext_fine( tlm_16, extram_et,
					       CNST_RFL8_A, 
					       CNST_RFL8_B,
					       CNST_RFL8_R0,
					       CNST_RFL8_H_R);
    // ET_FD1_TMP
    memcpy(&tlm_16, &ptr[EXTT_TLM_OFFSET+64], 2);
    extt_analog_tlm[35] = temp_16bit_norm_fine( tlm_16, extram_et,
						CNST_FD1_A, 
						CNST_FD1_B,
						CNST_FD1_R0);
    // ET_FD2_TMP
    memcpy(&tlm_16, &ptr[EXTT_TLM_OFFSET+66], 2);
    extt_analog_tlm[36] = temp_16bit_norm_fine( tlm_16, extram_et,
						CNST_FD2_A, 
						CNST_FD2_B,
						CNST_FD2_R0);
    // ET_FD3_TMP
    memcpy(&tlm_16, &ptr[EXTT_TLM_OFFSET+68], 2);
    extt_analog_tlm[37] = temp_16bit_norm_fine( tlm_16, extram_et,
						CNST_FD3_A, 
						CNST_FD3_B,
						CNST_FD3_R0);

    // Un-Swap EXTT HKT
    for (size_t i=0; i<35; i++) {
      memcpy(&usht, &ptr[EXTT_TLM_OFFSET+2*i], 2);
      usht = SWAP_2(usht);
      memcpy(&ptr[EXTT_TLM_OFFSET+2*i], &usht, 2);
    }


    //////////////////////////////////////
    //////// Begin APDU Telemetry //////// 
    //////////////////////////////////////

    // APDU analog telemetry
    float apdu_analog_tlm[5];

    // AP_ATC_TMP
    memcpy(&tlm_8, &ptr[APDU_TLM_OFFSET+0], 1);
    apdu_analog_tlm[0] = temp_apdu_8bit( tlm_8, CNST_AP_ATC_R0);

    // AP_RFE_TMP
    memcpy(&tlm_8, &ptr[APDU_TLM_OFFSET+1], 1);
    apdu_analog_tlm[1] = temp_apdu_8bit( tlm_8, CNST_AP_RFE_R0);

    // AP_RBE_DPU_TMP
    memcpy(&tlm_8, &ptr[APDU_TLM_OFFSET+2], 1);
    apdu_analog_tlm[2] = temp_apdu_8bit( tlm_8, CNST_AP_RBE_DPU_R0);

    // AP_SCAT_TMP
    memcpy(&tlm_8, &ptr[APDU_TLM_OFFSET+3], 1);
    apdu_analog_tlm[3] = temp_apdu_8bit( tlm_8, CNST_AP_SCAT_R0);

    // AP_ICDS_TMP
    memcpy(&tlm_8, &ptr[APDU_TLM_OFFSET+4], 1);
    apdu_analog_tlm[4] = temp_apdu_8bit( tlm_8, CNST_AP_ICDS_R0);


    /////////////////////////////////////////////////////////
    //////// Begin Scatterometer Converted Telemetry //////// 
    /////////////////////////////////////////////////////////
    float scatter_analog_tlm[30];

    // SC_LBK_B1HV_DC            0       0       16
    memcpy(&tlm_16, &ptr[SCATTEROMETER_TLM_OFFSET+0], 2);
    scatter_analog_tlm[0] = dwnlk_fltg_pt( SWAP_2(tlm_16)) * 1000 / (1 << 9);

    // SC_LBK_B1HV_PWR           2       0       16
    memcpy(&tlm_16, &ptr[SCATTEROMETER_TLM_OFFSET+2], 2);
    scatter_analog_tlm[1] = 
      10 * log10(dwnlk_fltg_pt( SWAP_2(tlm_16)) * 20 / (1 << 18));

    // SC_LVPS_BOX_TEMP          5       0       8
    memcpy(&tlm_8, &ptr[SCATTEROMETER_TLM_OFFSET+5], 1);
    scatter_analog_tlm[2] = 
      temp_8bit_norm_scat( tlm_8, CNST_SCAT_A, CNST_SCAT_B, CNST_LVPS_BOX_R0);

    // SC_SSPA_RF_DCK_TMP        6       0       8
    memcpy(&tlm_8, &ptr[SCATTEROMETER_TLM_OFFSET+6], 1);
    scatter_analog_tlm[3] = 
      temp_8bit_norm_scat( tlm_8, CNST_SCAT_A, CNST_SCAT_B, CNST_SSPA_RF_DCK_R0);

    // SC_SCG_TMP                7       0       8
    memcpy(&tlm_8, &ptr[SCATTEROMETER_TLM_OFFSET+7], 1);
    scatter_analog_tlm[4] = 
      temp_8bit_norm_scat( tlm_8, CNST_SCAT_A, CNST_SCAT_B, CNST_SCG_R0);

    // SC_SBE_LNA_TMP            8       0       8
    memcpy(&tlm_8, &ptr[SCATTEROMETER_TLM_OFFSET+8], 1);
    scatter_analog_tlm[5] = 
      temp_8bit_norm_scat( tlm_8, CNST_SCAT_A, CNST_SCAT_B, CNST_SBE_LNA_R0);

    // SC_SBE_TX_CHN_TMP         9       0       8
    memcpy(&tlm_8, &ptr[SCATTEROMETER_TLM_OFFSET+9], 1);
    scatter_analog_tlm[6] = 
      temp_8bit_norm_scat( tlm_8, CNST_SCAT_A, CNST_SCAT_B, CNST_SBE_TX_CHN_R0);

    // SC_SBE_RX_CHN_TMP         10      0       8
    memcpy(&tlm_8, &ptr[SCATTEROMETER_TLM_OFFSET+10], 1);
    scatter_analog_tlm[7] = 
      temp_8bit_norm_scat( tlm_8, CNST_SCAT_A, CNST_SCAT_B, CNST_SBE_RX_CHN_R0);

    // SC_SFE_TX_LOAD_TMP        11      0       8
    memcpy(&tlm_8, &ptr[SCATTEROMETER_TLM_OFFSET+11], 1);
    scatter_analog_tlm[8] = 
      temp_8bit_norm_scat( tlm_8, CNST_SCAT_A, CNST_SCAT_B, CNST_SFE_TX_LOAD_R0);

    // SC_SBE_STP_ATTEN_TMP      12      0       16
    memcpy(&tlm_16, &ptr[SCATTEROMETER_TLM_OFFSET+12], 2);
    scatter_analog_tlm[9] = temp_16bit_norm_disp( SWAP_2(tlm_16), 
						  gain_norm_cmp,
						  offset_norm_cmp,
						  CNST_SBE_STP_ATTEN_A, 
						  CNST_SBE_STP_ATTEN_B, 
						  CNST_SBE_STP_ATTEN_R0);


    // SC_SFE_LBK_ATTEN_TMP      14      0       16
    memcpy(&tlm_16, &ptr[SCATTEROMETER_TLM_OFFSET+14], 2);
    scatter_analog_tlm[10] = temp_16bit_norm_disp( SWAP_2(tlm_16), 
						   gain_norm_cmp,
						   offset_norm_cmp,
						   CNST_SFE_LBK_ATTEN_A, 
						   CNST_SFE_LBK_ATTEN_B, 
						   CNST_SFE_LBK_ATTEN_R0);

    // SC_SFE_LBK_SW_TMP         16      0       16
    memcpy(&tlm_16, &ptr[SCATTEROMETER_TLM_OFFSET+16], 2);
    scatter_analog_tlm[11] = temp_16bit_norm_disp( SWAP_2(tlm_16), 
						   gain_norm_cmp,
						   offset_norm_cmp,
						   CNST_SFE_LBK_SW_A, 
						   CNST_SFE_LBK_SW_B, 
						   CNST_SFE_LBK_SW_R0);

    // SC_SFE_BM_SW_TMP          18      0       16
    memcpy(&tlm_16, &ptr[SCATTEROMETER_TLM_OFFSET+18], 2);
    scatter_analog_tlm[12] = temp_16bit_norm_disp( SWAP_2(tlm_16), 
						   gain_norm_cmp,
						   offset_norm_cmp,
						   CNST_SFE_BM_SW_A, 
						   CNST_SFE_BM_SW_B, 
						   CNST_SFE_BM_SW_R0);

    // SC_SFE_GND_MON            20      0       8
    memcpy(&tlm_8, &ptr[SCATTEROMETER_TLM_OFFSET+20], 1);
    scatter_analog_tlm[13] = ie_volt( tlm_8);

    // SC_SFE_N15_VLT_MON        21      0       8
    memcpy(&tlm_8, &ptr[SCATTEROMETER_TLM_OFFSET+21], 1);
    scatter_analog_tlm[14] = ie_volt( tlm_8)* -6;

    // SC_SBE_PLM1_TUNG_VLT      22      0       8
    memcpy(&tlm_8, &ptr[SCATTEROMETER_TLM_OFFSET+22], 1);
    scatter_analog_tlm[15] = ie_volt( tlm_8);

    // SC_SBE_PLM2_TUNG_VLT      23      0       8
    memcpy(&tlm_8, &ptr[SCATTEROMETER_TLM_OFFSET+23], 1);
    scatter_analog_tlm[16] = ie_volt( tlm_8);

    // SC_8MHZ_STLO_OUT_PWR_LVL  24      0       8
    memcpy(&tlm_8, &ptr[SCATTEROMETER_TLM_OFFSET+24], 1);
    scatter_analog_tlm[17] = ie_volt( tlm_8);

    // SC_16MHZ_PWR_LVL          25      0       8
    memcpy(&tlm_8, &ptr[SCATTEROMETER_TLM_OFFSET+25], 1);
    scatter_analog_tlm[18] = ie_volt( tlm_8);

    // SC_SBE_PLM1_PWR_LVL       26      0       8
    memcpy(&tlm_8, &ptr[SCATTEROMETER_TLM_OFFSET+26], 1);
    scatter_analog_tlm[19] = ie_volt( tlm_8);

    // SC_SBE_PLM2_PWR_LVL       27      0       8
    memcpy(&tlm_8, &ptr[SCATTEROMETER_TLM_OFFSET+27], 1);
    scatter_analog_tlm[20] = ie_volt( tlm_8);

    // SC_SBE_TX_EXCTR_PWR_MON   28      0       8
    memcpy(&tlm_8, &ptr[SCATTEROMETER_TLM_OFFSET+28], 1);
    scatter_analog_tlm[21] = ie_volt( tlm_8);

    // SC_LVPS_CNVRTR_CURR       29      0       8
    memcpy(&tlm_8, &ptr[SCATTEROMETER_TLM_OFFSET+29], 1);
    scatter_analog_tlm[22] = ie_volt( tlm_8);

    // SC_SSPA_OUTPT_STG_VLT     30      0       8
    memcpy(&tlm_8, &ptr[SCATTEROMETER_TLM_OFFSET+30], 1);
    scatter_analog_tlm[23] = ie_volt( tlm_8) * 10;

    // SC_SSPA_INTRM_STG_VLT     31      0       8
    memcpy(&tlm_8, &ptr[SCATTEROMETER_TLM_OFFSET+31], 1);
    scatter_analog_tlm[24] = ie_volt( tlm_8) * 10;

    // SC_SSPA_INPT_STG_VLT      32      0       8
    memcpy(&tlm_8, &ptr[SCATTEROMETER_TLM_OFFSET+32], 1);
    scatter_analog_tlm[25] = ie_volt( tlm_8) * 5;

    // SC_SFE_P5V_MON            33      0       8
    memcpy(&tlm_8, &ptr[SCATTEROMETER_TLM_OFFSET+33], 1);
    scatter_analog_tlm[26] = ie_volt( tlm_8) * 2;

    // SC_SCG_P5V_MON            34      0       8
    memcpy(&tlm_8, &ptr[SCATTEROMETER_TLM_OFFSET+34], 1);
    scatter_analog_tlm[27] = ie_volt( tlm_8) * 2;

    // SC_SBE_P5V_MON            35      0       8
    memcpy(&tlm_8, &ptr[SCATTEROMETER_TLM_OFFSET+35], 1);
    scatter_analog_tlm[28] = ie_volt( tlm_8) * 2;

    // SC_SBE_P12V_NSW_MON       36      0       8
    memcpy(&tlm_8, &ptr[SCATTEROMETER_TLM_OFFSET+36], 1);
    scatter_analog_tlm[29] = ie_volt( tlm_8) * 4.8;


   uint8_t scatter_discrete_tlm[8];
   memcpy(&tlm_8, &ptr[SCATTEROMETER_TLM_OFFSET+4], 1);

   // SC_STS_FIFO_ERR           4       0       1
   scatter_discrete_tlm[0] = (tlm_8 >> 7) & 0x1;

   // SC_STS_FIFO_ERR_VAL       4       1       1
   scatter_discrete_tlm[1] = (tlm_8 >> 6) & 0x1;

   // SC_STS_SCAT_LVPS_OVR_CUR  4       2       1
   scatter_discrete_tlm[2] = (tlm_8 >> 5) & 0x1;

   // SC_STS_SCAT_LVPS_BUS_VLT  4       3       1
   scatter_discrete_tlm[3] = (tlm_8 >> 4) & 0x1;

   // SC_STS_CG_LTCH_PRM        4       6       1
   scatter_discrete_tlm[4] = (tlm_8 >> 1) & 0x1;

   // SC_STS_CG_LTCH_RED        4       7       1
   scatter_discrete_tlm[5] = tlm_8 & 0x1;

   // SC_SPARE1                 4       4       1
   // SC_SPARE2                 4       5       1
   scatter_discrete_tlm[6] = (tlm_8 >> 3) & 0x1;
   scatter_discrete_tlm[7] = (tlm_8 >> 2) & 0x1;




    /////////////////////////////////////
    //////// Begin ATC Telemetry //////// 
    /////////////////////////////////////

    static float atc_omt1_analog_tlm[17];
    static float atc_omt2_analog_tlm[17];
    static float atc_omt3_analog_tlm[17];
    static float atc_rbe_analog_tlm[17];
    float *atc_analog_tlm[4];
    atc_analog_tlm[0] = atc_omt1_analog_tlm;
    atc_analog_tlm[1] = atc_omt2_analog_tlm;
    atc_analog_tlm[2] = atc_omt3_analog_tlm;
    atc_analog_tlm[3] = atc_rbe_analog_tlm;
    float f32;

    // Swap ATC HKT
    for (size_t i=0; i<18; i++) {
      memcpy(&usht, &ptr[ATC_TLM_OFFSET+2*i], 2);
      usht = SWAP_2(usht);
      memcpy(&ptr[ATC_TLM_OFFSET+2*i], &usht, 2);
    }

    if ( atc_subframe == 0) {
      // AT_PRT1_PRM_CTL_MSBS   2       0       16
      memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+2], 2);
      tlm_32 = tlm_16 * 128;

      // AT_PRT1_PRM_CTL_LSB    32      0       4
      memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+32], 2);
      tlm_32 += (tlm_16 >> 12) * 8;

      // AT_PRT1_PRM_CTL
      f32 = (float) tlm_32 / 128;

      // AT_PRT1_PRM_CTL_RES
      f32 = f32 * CNST_ADC_OHMS - CNST_OMT1PRT1_OFFSET;

      // AT_PRT1_PRM_CTL_TMP
      atc_omt1_analog_tlm[0] = 
	(sqrt( pow(CNST_OMT1PRT1_A, 2) - 4 * CNST_OMT1PRT1_B * 
	       (1 - f32 / CNST_OMT1PRT1_R0)) - CNST_OMT1PRT1_A) / 
	(2 * CNST_OMT1PRT1_B);


      // AT_PRT2_H_PRB_CX_RBE   4       0       16
      memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+4], 2);

      // AT_PRT2_H_PRBCX_RES
      f32 = tlm_16 * CNST_ADC_OHMS - CNST_OMT1PRT2_OFFSET;

      // AT_PRT2_H_PRBCX_TMP
      atc_omt1_analog_tlm[1] = 
	(sqrt( pow(CNST_OMT1PRT2_A, 2) - 4 * CNST_OMT1PRT2_B * 
	       (1 - f32 / CNST_OMT1PRT2_R0)) - CNST_OMT1PRT2_A) / 
	(2 * CNST_OMT1PRT2_B);


      // AT_PRT3_OMT_NCK_SFE    6       0       16
      memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+6], 2);

      // AT_PRT3_OMTNECK_RES
      f32 = tlm_16 * CNST_ADC_OHMS - CNST_OMT1PRT3_OFFSET;

      // AT_PRT3_OMTNECK_TMP
      atc_omt1_analog_tlm[2] = 
	(sqrt( pow(CNST_OMT1PRT3_A, 2) - 4 * CNST_OMT1PRT3_B * 
	       (1 - f32 / CNST_OMT1PRT3_R0) ) - CNST_OMT1PRT3_A) / 
	(2 * CNST_OMT1PRT3_B);


      // AT_PRT4_RED_CTL_MSBS   8       0       16
      memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+8], 2);
      tlm_32 = tlm_16 * 128;

      // AT_PRT4_RED_CTL_LSB    34      0       4
      memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+34], 2);
      tlm_32 += (tlm_16 >> 12) * 8;

      // AT_PRT4_RED_CTL
      f32 = (float) tlm_32 / 128;

      // AT_PRT4_REDUN_CTL_RES
      f32 = f32 * CNST_ADC_OHMS - CNST_OMT1PRT4_OFFSET;

      // AT_PRT4_REDUN_CTL_TMP
      atc_omt1_analog_tlm[3] = 
	(sqrt( pow(CNST_OMT1PRT4_A, 2) - 4 * CNST_OMT1PRT4_B * 
	       (1 - f32 / CNST_OMT1PRT4_R0) ) - CNST_OMT1PRT4_A) / 
	(2 * CNST_OMT1PRT4_B);


      // AT_PRT5_V_PRB_CX_SBE   10      0       16
      memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+10], 2);

      // AT_PRT5_V_PRBCX_RES 
      f32 = tlm_16 * CNST_ADC_OHMS - CNST_OMT1PRT5_OFFSET;

      // AT_PRT5_V_PRBCX_TMP
      atc_omt1_analog_tlm[4] = 
	(sqrt( pow(CNST_OMT1PRT5_A, 2) - 4 * CNST_OMT1PRT5_B * 
	       (1 - f32 / CNST_OMT1PRT5_R0) ) - CNST_OMT1PRT5_A) / 
	(2 * CNST_OMT1PRT5_B);


      // AT_PRT6_OMT_STR_SCG    12      0       16
      memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+12], 2);

      // AT_PRT6_OMT_RES
      f32 = tlm_16 * CNST_ADC_OHMS - CNST_OMT1PRT6_OFFSET;

      // AT_PRT6_OMT_TMP
      atc_omt1_analog_tlm[5] = 
	(sqrt( pow(CNST_OMT1PRT6_A, 2) - 4 * CNST_OMT1PRT6_B * 
	       (1 - f32 / CNST_OMT1PRT6_R0) ) - CNST_OMT1PRT6_A) / 
	(2 * CNST_OMT1PRT6_B);

    } else if ( atc_subframe == 1) {
      // AT_PRT1_PRM_CTL_MSBS   2       0       16
      memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+2], 2);
      tlm_32 = tlm_16 * 128;

      // AT_PRT1_PRM_CTL_LSB    32      0       4
      memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+32], 2);
      tlm_32 += (tlm_16 >> 12) * 8;

      // AT_PRT1_PRM_CTL
      f32 = (float) tlm_32 / 128;

      // AT_PRT1_PRM_CTL_RES
      f32 = f32 * CNST_ADC_OHMS - CNST_OMT2PRT1_OFFSET;

      // AT_PRT1_PRM_CTL_TMP
      atc_omt2_analog_tlm[0] = 
	(sqrt( pow(CNST_OMT2PRT1_A, 2) - 4 * CNST_OMT2PRT1_B * 
	       (1 - f32 / CNST_OMT2PRT1_R0)) - CNST_OMT2PRT1_A) / 
	(2 * CNST_OMT2PRT1_B);


      // AT_PRT2_H_PRB_CX_RBE   4       0       16
      memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+4], 2);

      // AT_PRT2_H_PRBCX_RES
      f32 = tlm_16 * CNST_ADC_OHMS - CNST_OMT2PRT2_OFFSET;

      // AT_PRT2_H_PRBCX_TMP
      atc_omt2_analog_tlm[1] = 
	(sqrt( pow(CNST_OMT2PRT2_A, 2) - 4 * CNST_OMT2PRT2_B * 
	       (1 - f32 / CNST_OMT2PRT2_R0)) - CNST_OMT2PRT2_A) / 
	(2 * CNST_OMT2PRT2_B);


      // AT_PRT3_OMT_NCK_SFE    6       0       16
      memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+6], 2);

      // AT_PRT3_OMTNECK_RES
      f32 = tlm_16 * CNST_ADC_OHMS - CNST_OMT2PRT3_OFFSET;

      // AT_PRT3_OMTNECK_TMP
      atc_omt2_analog_tlm[2] = 
	(sqrt( pow(CNST_OMT2PRT3_A, 2) - 4 * CNST_OMT2PRT3_B * 
	       (1 - f32 / CNST_OMT2PRT3_R0) ) - CNST_OMT2PRT3_A) / 
	(2 * CNST_OMT2PRT3_B);


      // AT_PRT4_RED_CTL_MSBS   8       0       16
      memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+8], 2);
      tlm_32 = tlm_16 * 128;

      // AT_PRT4_RED_CTL_LSB    34      0       4
      memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+34], 2);
      tlm_32 += (tlm_16 >> 12) * 8;

      // AT_PRT4_RED_CTL
      f32 = (float) tlm_32 / 128;

      // AT_PRT4_REDUN_CTL_RES
      f32 = f32 * CNST_ADC_OHMS - CNST_OMT2PRT4_OFFSET;

      // AT_PRT4_REDUN_CTL_TMP
      atc_omt2_analog_tlm[3] = 
	(sqrt( pow(CNST_OMT2PRT4_A, 2) - 4 * CNST_OMT2PRT4_B * 
	       (1 - f32 / CNST_OMT2PRT4_R0) ) - CNST_OMT2PRT4_A) / 
	(2 * CNST_OMT2PRT4_B);


      // AT_PRT5_V_PRB_CX_SBE   10      0       16
      memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+10], 2);

      // AT_PRT5_V_PRBCX_RES 
      f32 = tlm_16 * CNST_ADC_OHMS - CNST_OMT2PRT5_OFFSET;

      // AT_PRT5_V_PRBCX_TMP
      atc_omt2_analog_tlm[4] = 
	(sqrt( pow(CNST_OMT2PRT5_A, 2) - 4 * CNST_OMT2PRT5_B * 
	       (1 - f32 / CNST_OMT2PRT5_R0) ) - CNST_OMT2PRT5_A) / 
	(2 * CNST_OMT2PRT5_B);


      // AT_PRT6_OMT_STR_SCG    12      0       16
      memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+12], 2);

      // AT_PRT6_OMT_RES
      f32 = tlm_16 * CNST_ADC_OHMS - CNST_OMT2PRT6_OFFSET;

      // AT_PRT6_OMT_TMP
      atc_omt2_analog_tlm[5] = 
	(sqrt( pow(CNST_OMT2PRT6_A, 2) - 4 * CNST_OMT2PRT6_B * 
	       (1 - f32 / CNST_OMT2PRT6_R0) ) - CNST_OMT2PRT6_A) / 
	(2 * CNST_OMT2PRT6_B);

    } else if ( atc_subframe == 2) {
      // AT_PRT1_PRM_CTL_MSBS   2       0       16
      memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+2], 2);
      tlm_32 = tlm_16 * 128;

      // AT_PRT1_PRM_CTL_LSB    32      0       4
      memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+32], 2);
      tlm_32 += (tlm_16 >> 12) * 8;

      // AT_PRT1_PRM_CTL
      f32 = (float) tlm_32 / 128;

      // AT_PRT1_PRM_CTL_RES
      f32 = f32 * CNST_ADC_OHMS - CNST_OMT3PRT1_OFFSET;

      // AT_PRT1_PRM_CTL_TMP
      atc_omt3_analog_tlm[0] = 
	(sqrt( pow(CNST_OMT3PRT1_A, 2) - 4 * CNST_OMT3PRT1_B * 
	       (1 - f32 / CNST_OMT3PRT1_R0)) - CNST_OMT3PRT1_A) / 
	(2 * CNST_OMT3PRT1_B);


      // AT_PRT2_H_PRB_CX_RBE   4       0       16
      memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+4], 2);

      // AT_PRT2_H_PRBCX_RES
      f32 = tlm_16 * CNST_ADC_OHMS - CNST_OMT3PRT2_OFFSET;

      // AT_PRT2_H_PRBCX_TMP
      atc_omt3_analog_tlm[1] = 
	(sqrt( pow(CNST_OMT3PRT2_A, 2) - 4 * CNST_OMT3PRT2_B * 
	       (1 - f32 / CNST_OMT3PRT2_R0)) - CNST_OMT3PRT2_A) / 
	(2 * CNST_OMT3PRT2_B);


      // AT_PRT3_OMT_NCK_SFE    6       0       16
      memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+6], 2);

      // AT_PRT3_OMTNECK_RES
      f32 = tlm_16 * CNST_ADC_OHMS - CNST_OMT3PRT3_OFFSET;

      // AT_PRT3_OMTNECK_TMP
      atc_omt3_analog_tlm[2] = 
	(sqrt( pow(CNST_OMT3PRT3_A, 2) - 4 * CNST_OMT3PRT3_B * 
	       (1 - f32 / CNST_OMT3PRT3_R0) ) - CNST_OMT3PRT3_A) / 
	(2 * CNST_OMT3PRT3_B);


      // AT_PRT4_RED_CTL_MSBS   8       0       16
      memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+8], 2);
      tlm_32 = tlm_16 * 128;

      // AT_PRT4_RED_CTL_LSB    34      0       4
      memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+34], 2);
      tlm_32 += (tlm_16 >> 12) * 8;

      // AT_PRT4_RED_CTL
      f32 = (float) tlm_32 / 128;

      // AT_PRT4_REDUN_CTL_RES
      f32 = f32 * CNST_ADC_OHMS - CNST_OMT3PRT4_OFFSET;

      // AT_PRT4_REDUN_CTL_TMP
      atc_omt3_analog_tlm[3] = 
	(sqrt( pow(CNST_OMT3PRT4_A, 2) - 4 * CNST_OMT3PRT4_B * 
	       (1 - f32 / CNST_OMT3PRT4_R0) ) - CNST_OMT3PRT4_A) / 
	(2 * CNST_OMT3PRT4_B);



      // AT_PRT5_V_PRB_CX_SBE   10      0       16
      memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+10], 2);

      // AT_PRT5_V_PRBCX_RES 
      f32 = tlm_16 * CNST_ADC_OHMS - CNST_OMT3PRT5_OFFSET;

      // AT_PRT5_V_PRBCX_TMP
      atc_omt3_analog_tlm[4] = 
	(sqrt( pow(CNST_OMT3PRT5_A, 2) - 4 * CNST_OMT3PRT5_B * 
	       (1 - f32 / CNST_OMT3PRT5_R0) ) - CNST_OMT3PRT5_A) / 
	(2 * CNST_OMT3PRT5_B);


      // AT_PRT6_OMT_STR_SCG    12      0       16
      memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+12], 2);

      // AT_PRT6_OMT_RES
      f32 = tlm_16 * CNST_ADC_OHMS - CNST_OMT3PRT6_OFFSET;

      // AT_PRT6_OMT_TMP
      atc_omt3_analog_tlm[5] = 
	(sqrt( pow(CNST_OMT3PRT6_A, 2) - 4 * CNST_OMT3PRT6_B * 
	       (1 - f32 / CNST_OMT3PRT6_R0) ) - CNST_OMT3PRT6_A) / 
	(2 * CNST_OMT3PRT6_B);

    } else if ( atc_subframe == 3) {
      // AT_PRT1_PRM_CTL_MSBS   2       0       16
      memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+2], 2);
      tlm_32 = tlm_16 * 128;

      // AT_PRT1_PRM_CTL_LSB    32      0       4
      memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+32], 2);
      tlm_32 += (tlm_16 >> 12) * 8;

      // AT_PRT1_PRM_CTL
      f32 = (float) tlm_32 / 128;

      // AT_PRT1_PRM_CTL_RES
      f32 = f32 * CNST_ADC_OHMS - CNST_RBEPRT1_OFFSET;

      // AT_PRT1_PRM_CTL_TMP
      atc_rbe_analog_tlm[0] = 
	(sqrt( pow(CNST_RBEPRT1_A, 2) - 4 * CNST_RBEPRT1_B * 
	       (1 - f32 / CNST_RBEPRT1_R0)) - CNST_RBEPRT1_A) / 
	(2 * CNST_RBEPRT1_B);


      // AT_PRT2_H_PRB_CX_RBE   4       0       16
      memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+4], 2);

      // AT_PRT2_H_PRBCX_RES
      f32 = tlm_16 * CNST_ADC_OHMS - CNST_RBEPRT2_OFFSET;

      // AT_PRT2_RBE_TMP
      atc_rbe_analog_tlm[1] = 
	(sqrt( pow(CNST_RBEPRT2_A, 2) - 4 * CNST_RBEPRT2_B * 
	       (1 - f32 / CNST_RBEPRT2_R0)) - CNST_RBEPRT2_A) / 
	(2 * CNST_RBEPRT2_B);


      // AT_PRT3_OMT_NCK_SFE    6       0       16
      memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+6], 2);

      // AT_PRT3_OMTNECK_RES
      f32 = tlm_16 * CNST_ADC_OHMS - CNST_RBEPRT3_OFFSET;

      // AT_PRT3_SFE_TMP
      atc_rbe_analog_tlm[2] = 
	(sqrt( pow(CNST_RBEPRT3_A, 2) - 4 * CNST_RBEPRT3_B * 
	       (1 - f32 / CNST_RBEPRT3_R0) ) - CNST_RBEPRT3_A) / 
	(2 * CNST_RBEPRT3_B);


      // AT_PRT4_RED_CTL_MSBS   8       0       16
      memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+8], 2);
      tlm_32 = tlm_16 * 128;

      // AT_PRT4_RED_CTL_LSB    34      0       4
      memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+34], 2);
      tlm_32 += (tlm_16 >> 12) * 8;

      // AT_PRT4_RED_CTL
      f32 = (float) tlm_32 / 128;

      // AT_PRT4_REDUN_CTL_RES
      f32 = f32 * CNST_ADC_OHMS - CNST_RBEPRT4_OFFSET;

      // AT_PRT4_REDUN_CTL_TMP
      atc_rbe_analog_tlm[3] = 
	(sqrt( pow(CNST_RBEPRT4_A, 2) - 4 * CNST_RBEPRT4_B * 
	       (1 - f32 / CNST_RBEPRT4_R0) ) - CNST_RBEPRT4_A) / 
	(2 * CNST_RBEPRT4_B);



      // AT_PRT5_V_PRB_CX_SBE   10      0       16
      memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+10], 2);

      // AT_PRT5_V_PRBCX_RES 
      f32 = tlm_16 * CNST_ADC_OHMS - CNST_RBEPRT5_OFFSET;

      // AT_PRT5_SBE_TMP
      atc_rbe_analog_tlm[4] = 
	(sqrt( pow(CNST_RBEPRT5_A, 2) - 4 * CNST_RBEPRT5_B * 
	       (1 - f32 / CNST_RBEPRT5_R0) ) - CNST_RBEPRT5_A) / 
	(2 * CNST_RBEPRT5_B);


      // AT_PRT6_OMT_STR_SCG    12      0       16
      memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+12], 2);

      // AT_PRT6_OMT_RES
      f32 = tlm_16 * CNST_ADC_OHMS - CNST_RBEPRT6_OFFSET;

      // AT_PRT6_SCG_TMP
      atc_rbe_analog_tlm[5] = 
	(sqrt( pow(CNST_RBEPRT6_A, 2) - 4 * CNST_RBEPRT6_B * 
	       (1 - f32 / CNST_RBEPRT6_R0) ) - CNST_RBEPRT6_A) / 
	(2 * CNST_RBEPRT6_B);
    }


    // AT_TMP_CTL_SP_RDBK     14      0       16
    memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+14], 2);
    *(atc_analog_tlm[atc_subframe] + 6) = 
      tlm_16 * AT_CNST_TEMP_GAIN + AT_CNST_TEMP_OFFSET;

    // AT_PID_PRM_RDBK_PRP    16      0       16
    memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+16], 2);
    if ( atc_subframe == 3)
      *(atc_analog_tlm[atc_subframe] + 7) = tlm_16 * AT_CNST_KP_RBE_SCALE;
    else
      *(atc_analog_tlm[atc_subframe] + 7) = tlm_16 * AT_CNST_KP_OMT_SCALE;

    // AT_PID_PRM_RDBK_IGL    18      0       16
    memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+18], 2);
    if ( atc_subframe == 3)
      *(atc_analog_tlm[atc_subframe] + 8) = tlm_16 * AT_CNST_KI_RBE_SCALE;
    else
      *(atc_analog_tlm[atc_subframe] + 8) = tlm_16 * AT_CNST_KI_OMT_SCALE;

    // AT_PID_PRM_RDBK_DRV    20      0       16
    memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+20], 2);
    if ( atc_subframe == 3)
      *(atc_analog_tlm[atc_subframe] + 9) = tlm_16 * AT_CNST_KD_RBE_SCALE;
    else
      *(atc_analog_tlm[atc_subframe] + 9) = tlm_16 * AT_CNST_KD_OMT_SCALE;

    // AT_CMP_PRP             22      0       24
    memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+22], 2);
    tlm_32 = tlm_16 << 8;
    memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+24], 2);
    tlm_32 += tlm_16 >> 8;

    // 24-bit two's complement
    if ( tlm_32 > 0x7fffff)
      s_tlm_32 = tlm_32 - 0x1000000;
    else 
      s_tlm_32 = tlm_32;

    if ( atc_subframe == 3)
      *(atc_analog_tlm[atc_subframe] + 10) = 
	s_tlm_32 * AT_CNST_CMPID_RBE_SCALE;
    else
      *(atc_analog_tlm[atc_subframe] + 10) = 
	s_tlm_32 * AT_CNST_CMPID_OMT_SCALE;

    // AT_CMP_INT             25      0       24
    memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+24], 2);
    tlm_32 = (tlm_16 & 0xff) << 16;
    memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+26], 2);
    tlm_32 += tlm_16;

    // 24-bit two's complement
    if ( tlm_32 > 0x7fffff)
      s_tlm_32 = tlm_32 - 0x1000000;
    else 
      s_tlm_32 = tlm_32;

    if ( atc_subframe == 3)
      *(atc_analog_tlm[atc_subframe] + 11) = 
	s_tlm_32 * AT_CNST_CMPID_RBE_SCALE;
    else
      *(atc_analog_tlm[atc_subframe] + 11) = 
	s_tlm_32 * AT_CNST_CMPID_OMT_SCALE;

    // AT_CMP_DRV             28      0       24
    memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+28], 2);
    tlm_32 = tlm_16 << 8;
    memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET+30], 2);
    tlm_32 += tlm_16 >> 8;

    // 24-bit two's complement
    if ( tlm_32 > 0x7fffff)
      s_tlm_32 = tlm_32 - 0x1000000;
    else 
      s_tlm_32 = tlm_32;

    if ( atc_subframe == 3)
      *(atc_analog_tlm[atc_subframe] + 12) = 
	s_tlm_32 * AT_CNST_CMPID_RBE_SCALE;
    else
      *(atc_analog_tlm[atc_subframe] + 12) = 
	s_tlm_32 * AT_CNST_CMPID_OMT_SCALE;

    // AT_HTR_PWR_OFST_RDBK   32      4       12
    memcpy(&tlm_32, &ptr[ATC_TLM_OFFSET+32], 4);
    tlm_16 = tlm_32 & 0xfff;
    if ( atc_subframe == 3)
      *(atc_analog_tlm[atc_subframe] + 13) = tlm_16 * AT_CNST_PWR_RBE_SCALE;
    else
      *(atc_analog_tlm[atc_subframe] + 13) = tlm_16 * AT_CNST_PWR_OMT_SCALE;
    *(atc_analog_tlm[atc_subframe] + 15) = tlm_16 * AT_CNST_PWR_PRCNT;

    // AT_HTR_TOT_PWR_OUT     34      4       12
    memcpy(&tlm_32, &ptr[ATC_TLM_OFFSET+34], 4);
    tlm_16 = tlm_32 & 0xfff;
    if ( atc_subframe == 3)
      *(atc_analog_tlm[atc_subframe] + 14) = tlm_16 * AT_CNST_PWR_RBE_SCALE;
    else
      *(atc_analog_tlm[atc_subframe] + 14) = tlm_16 * AT_CNST_PWR_OMT_SCALE;

    *(atc_analog_tlm[atc_subframe] + 16) = tlm_16 * AT_CNST_PWR_PRCNT;


    // ATC Discrete Fields
    static unsigned char atc_omt1_discrete_tlm[17];
    static unsigned char atc_omt2_discrete_tlm[17];
    static unsigned char atc_omt3_discrete_tlm[17];
    static unsigned char atc_rbe_discrete_tlm[17];
    unsigned char *atc_discrete_tlm[4];
    atc_discrete_tlm[0] = atc_omt1_discrete_tlm;
    atc_discrete_tlm[1] = atc_omt2_discrete_tlm;
    atc_discrete_tlm[2] = atc_omt3_discrete_tlm;
    atc_discrete_tlm[3] = atc_rbe_discrete_tlm;

    memcpy(&tlm_16, &ptr[ATC_TLM_OFFSET], 2);

    // AT_ST_3LSB_BLK_ID      0       0       3
    *(atc_discrete_tlm[atc_subframe]) = (tlm_16 >> 13) & 0x7;

    // AT_ST_ACTIVE_SNSR      0       3       1
    // AT_ST_CMD_ERR          0       4       1
    // AT_ST_ADC_LTCH         0       5       1
    // AT_ST_ADC1_OVR_MX      0       6       1
    // AT_ST_ADC2_OVR_MX      0       7       1
    // AT_ST_FEW_TMP_FLT_FL   0       8       1
    // AT_ST_MANY_TMP_FLT_FL  0       9       1
    // AT_ST_ADC_RST          0       10      1
    // AT_ST_WD_RST           0       11      1
    // AT_ST_HRD_RST          0       12      1
    // AT_ST_SFT_RST          0       13      1
    for (size_t i=0; i<11; i++)
      *(atc_discrete_tlm[atc_subframe]+11-i) = (tlm_16 >> (i+2)) & 0x1;


    // AT_ST_ADC_TOUT         31      3       1
    // AT_ST_WD_RST_ENA       31      4       1
    // AT_ST_ADC1_PWRD        31      5       1
    // AT_ST_ADC2_PWRD        31      6       1
    // AT_ST_AUTOPRT_ENA      31      7       1
    memcpy(&tlm_8, &ptr[ATC_TLM_OFFSET+30], 1);
    for (size_t i=0; i<5; i++)
      *(atc_discrete_tlm[atc_subframe]+16-i) = (tlm_8 >> i) & 0x1;

    // Un-Swap ATC HKT
    for (size_t i=0; i<18; i++) {
      memcpy(&usht, &ptr[ATC_TLM_OFFSET+2*i], 2);
      usht = SWAP_2(usht);
      memcpy(&ptr[ATC_TLM_OFFSET+2*i], &usht, 2);
    }

    if ( (rad_header & 0xfff0) == 0xfa50) {
      // Write converted telemetry fields
      start[0] = rad_frmnum;
      count[0]  = 1;

      start[1] = 0;
      count[1] = 23;
      //      if ( tmedbl <  851649650  && tmedbl > 851648400) {
      //cout << (long) tmedbl << " " << start[0] << " " << (long) rad_subframe << " " << 
      //  rfe1_analog_tlm[7] << endl;
      //}
      PTB( h5d_write(grp5, "rfe1_analog_tlm", rfe1_analog_tlm, 2, 
		     start, count));
      PTB( h5d_write(grp5, "rfe2_analog_tlm", rfe2_analog_tlm, 2, 
		     start, count));
      PTB( h5d_write(grp5, "rfe3_analog_tlm", rfe3_analog_tlm, 2, 
		     start, count));

      start[1] = 0;
      count[1] = 21;
      PTB( h5d_write(grp5, "rbe1_analog_tlm", rbe1_analog_tlm, 2, 
		     start, count));
      PTB( h5d_write(grp5, "rbe2_analog_tlm", rbe2_analog_tlm, 2, 
		     start, count));
      PTB( h5d_write(grp5, "rbe3_analog_tlm", rbe3_analog_tlm, 2, 
		   start, count));

      start[1] = 0;
      count[1] = 12;
      PTB( h5d_write(grp5, "dpu_analog_tlm", dpu_analog_tlm, 2, 
		     start, count));


      start[1] = rad_subframe;
      count[1]  = 1;

      start[2] = 0;
      count[2] = 13;
      PTB( h5d_write(grp5, "dpu_status_tlm", dpu_status_tlm, 3, 
		     start, count));

      start[2] = 0;
      count[2] = 17;
      PTB( h5d_write(grp5, "radiom_nrt_tlm", nrt_discrete_tlm, 3, 
		     start, count));

    } // (rad_header & 0xfff0) == 0xfa50


    // Zero out arrays if radiometer subframe index = 3
    //    if ( rad_subframe == 3) {
    //memset( rfe1_analog_tlm, 0, 23*sizeof(float));
    //memset( rfe2_analog_tlm, 0, 23*sizeof(float));
    //memset( rfe3_analog_tlm, 0, 23*sizeof(float));
    //memset( rbe1_analog_tlm, 0, 21*sizeof(float));
    //memset( rbe2_analog_tlm, 0, 21*sizeof(float));
    //memset( rbe3_analog_tlm, 0, 21*sizeof(float));
    //memset( dpu_analog_tlm,  0, 12*sizeof(float));
    //}


    // ATC analog
    start[0] = atc_frmnum;
    count[0]  = 1;

    start[1] = 0;
    count[1] = 17;
    PTB( h5d_write(grp5, "atc_omt1_analog_tlm", atc_omt1_analog_tlm, 2, 
		   start, count));
    PTB( h5d_write(grp5, "atc_omt2_analog_tlm", atc_omt2_analog_tlm, 2, 
		   start, count));
    PTB( h5d_write(grp5, "atc_omt3_analog_tlm", atc_omt3_analog_tlm, 2, 
		   start, count));
    PTB( h5d_write(grp5, "atc_rbe_analog_tlm", atc_rbe_analog_tlm, 2, 
		   start, count));
    PTB( h5d_write(grp5, "atc_omt1_discrete_tlm", atc_omt1_discrete_tlm, 2, 
		   start, count));
    PTB( h5d_write(grp5, "atc_omt2_discrete_tlm", atc_omt2_discrete_tlm, 2, 
		   start, count));
    PTB( h5d_write(grp5, "atc_omt3_discrete_tlm", atc_omt3_discrete_tlm, 2, 
		   start, count));
    PTB( h5d_write(grp5, "atc_rbe_discrete_tlm", atc_rbe_discrete_tlm, 2, 
		   start, count));



    start[0] = blknum;
    count[0]  = 1;

    start[1] = 0;
    count[1] = 4;
    PTB( h5d_write(grp5, "deploy_analog_tlm", deploy_analog_tlm, 2, 
		   start, count));

    start[1] = 0;
    count[1] = 7;
    PTB( h5d_write(grp5, "deploy_discrete_tlm", deploy_discrete_tlm, 2, 
		   start, count));

    start[1] = 0;
    count[1] = 14;
    PTB( h5d_write(grp5, "icds_analog_tlm", icds_analog_tlm, 2, 
		   start, count));

    start[1] = 0;
    count[1] = 8;
    PTB( h5d_write(grp5, "icds_discrete_tlm", icds_discrete_tlm, 2, 
		   start, count));

    start[1] = 0;
    count[1] = 38;
    PTB( h5d_write(grp5, "ext_temp_analog_tlm", extt_analog_tlm, 2, 
		   start, count));

    start[1] = 0;
    count[1] = 5;
    PTB( h5d_write(grp5, "apdu_analog_tlm", apdu_analog_tlm, 2, 
		   start, count));

    start[1] = 0;
    count[1] = 30;
    PTB( h5d_write(grp5, "scatter_analog_tlm", scatter_analog_tlm, 2, 
		   start, count));

    start[1] = 0;
    count[1] = 8;
    PTB( h5d_write(grp5, "scatter_discrete_tlm", scatter_discrete_tlm, 2, 
		   start, count));

    prev_icds_blknum = icds_blknum;
    prev_rad_header = rad_header;
    first_write = false;

    granuleStop = tmedbl;

    return 0;
  }

  /*----------------------------------------------------------------- */
  /* Open an existing HDF file                                        */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::openl1( char* l1_filename, unsigned flags)   
  {
    h5fid =  H5Fopen(l1_filename, flags, H5P_DEFAULT);
    if ( h5fid < 0) {
      cout << l1_filename << " not found.";
      exit(1);
    }

    grp0 = H5Gopen1(h5fid,"/");
    if ( grp0 < 0) {
      cout << "Root group not found.";
      exit(1);
    }

    grp1 = H5Gopen1(h5fid, "/Block Attributes");
    if ( grp1 < 0) {
      cout << "'Block Attributes' group not found.";
      exit(1);
    }

    grp2 = H5Gopen1(h5fid, "/Raw Aquarius Data");
    if ( grp2 < 0) {
      cout << "'Raw Aquarius Data' group not found.";
      exit(1);
    }

    grp3 = H5Gopen1(h5fid, "/SAC-D Telemetry");
    if ( grp3 < 0) {
      cout << "'SAC-D Telemetry' group not found.";
      exit(1);
    }

    grp4 = H5Gopen1(h5fid, "/Navigation");
    if ( grp4 < 0) {
      cout << "'/Navigation' group not found.";
      exit(1);
    }

    grp5 = H5Gopen1(h5fid, "/Converted Telemetry");
    if ( grp5 < 0) {
      cout << "'Converted Telemetry' group not found.";
      exit(1);
    }

    hid_t attr;
    attr = H5Aopen_name(grp0, "Number of Blocks");
    H5Aread(attr, H5T_STD_I32LE, &blknum);
    H5Aclose(attr);

    int32_t year, doy, millisec;

    attr = H5Aopen_name(grp0, "Start Year");
    H5Aread(attr, H5T_STD_I32LE, &year);
    H5Aclose(attr);

    attr = H5Aopen_name(grp0, "Start Day");
    H5Aread(attr, H5T_STD_I32LE, &doy);
    H5Aclose(attr);

    attr = H5Aopen_name(grp0, "Start Millisec");
    H5Aread(attr, H5T_STD_I32LE, &millisec);
    H5Aclose(attr);

    // JMG  time is seconds from 01/06/1980
    granuleStart = get_tai( year, doy, (double) millisec); // new

    attr = H5Aopen_name(grp0, "End Year");
    H5Aread(attr, H5T_STD_I32LE, &year);
    H5Aclose(attr);

    attr = H5Aopen_name(grp0, "End Day");
    H5Aread(attr, H5T_STD_I32LE, &doy);
    H5Aclose(attr);

    attr = H5Aopen_name(grp0, "End Millisec");
    H5Aread(attr, H5T_STD_I32LE, &millisec);
    H5Aclose(attr);

    // JMG  time is seconds from 01/06/1980
    granuleStop = get_tai( year, doy, (double) millisec); // new


    /* Save old error handler */
    H5E_auto_t old_func;
    void *old_client_data;
    H5Eget_auto(H5E_DEFAULT , &old_func, &old_client_data);

    char timeString[17];
    hid_t atype = H5Tcopy(H5T_C_S1);
    H5Tset_size(atype, 17);

    H5Eset_auto( H5E_DEFAULT, NULL, NULL); // Turn off error handling
    attr = H5Aopen_name(grp0, "Orbit Start Time");
    H5Eset_auto( H5E_DEFAULT, old_func, old_client_data); // Restore previous error handler
    if (attr == -1) {
      orbitStart = -999;
    } else {
      // time is seconds from 01/06/1980
      H5Aread(attr, atype, timeString);
      orbitStart = get_tai( timeString);
      H5Aclose(attr);
    }

    H5Eset_auto( H5E_DEFAULT, NULL, NULL); // Turn off error handling
    attr = H5Aopen_name(grp0, "Orbit Stop Time");
    H5Eset_auto( H5E_DEFAULT, old_func, old_client_data); // Restore previous error handler
    if (attr == -1) {
      orbitStop = -999;
    } else {
      // time is seconds from 01/06/1980
      H5Aread(attr, atype, timeString);
      orbitStop = get_tai( timeString);
      H5Aclose(attr);
    }

    H5Eset_auto( H5E_DEFAULT, NULL, NULL); // Turn off error handling
    attr = H5Aopen_name(grp0, "Orbit Number");
    H5Eset_auto( H5E_DEFAULT, old_func, old_client_data); // Restore previous error handler
    if (attr == -1) {
      orbitNumber = -999;
    } else {
      // time is seconds from 01/06/1980
      H5Aread(attr, H5T_STD_I32LE, &orbitNumber);
      H5Aclose(attr);
    }

    H5Eset_auto( H5E_DEFAULT, NULL, NULL); // Turn off error handling
    attr = H5Aopen_name(grp0, "Cycle Number");
    H5Eset_auto( H5E_DEFAULT, old_func, old_client_data); // Restore previous error handler
    if (attr == -1) {
      orbitNumber = -999;
    } else {
      // time is seconds from 01/06/1980
      H5Aread(attr, H5T_STD_I32LE, &cycleNumber);
      H5Aclose(attr);
    }

    H5Eset_auto( H5E_DEFAULT, NULL, NULL); // Turn off error handling
    attr = H5Aopen_name(grp0, "Pass Number");
    H5Eset_auto( H5E_DEFAULT, old_func, old_client_data); // Restore previous error handler
    if (attr == -1) {
      orbitNumber = -999;
    } else {
      // time is seconds from 01/06/1980
      H5Aread(attr, H5T_STD_I32LE, &passNumber);
      H5Aclose(attr);
    }

    H5Eset_auto( H5E_DEFAULT, NULL, NULL); // Turn off error handling
    attr = H5Aopen_name(grp0, "Node Crossing Time");
    H5Eset_auto( H5E_DEFAULT, old_func, old_client_data); // Restore previous error handler
    if (attr == -1) {
      nodeCrossingTime = -999;
    } else {
      // time is seconds from 01/06/1980
      H5Aread(attr, atype, timeString);
      if( strcmp( timeString, "N/A") != 0)
	nodeCrossingTime = get_tai( timeString);
      else
	nodeCrossingTime = -999;
      H5Aclose(attr);
    }
    H5Tclose(atype);


    H5Eset_auto( H5E_DEFAULT, NULL, NULL); // Turn off error handling
    attr = H5Aopen_name(grp0, "Orbit Node Longitude");
    H5Eset_auto( H5E_DEFAULT, old_func, old_client_data); // Restore previous error handler
    if (attr == -1) {
      nodeLongitude = -999;
    } else {
      H5Aread(attr, H5T_NATIVE_FLOAT, &nodeLongitude);
      H5Aclose(attr);
    }

    openFlags = flags;

    return 0;
  }


  /*----------------------------------------------------------------- */
  /* Read single radiometer record                                    */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::readl1_radiometer( int32_t blknum, 
					int32_t frmnum, uint8_t subframe, 
					double *blkSec,
					uint16_t *radiomHdr, 
					unsigned short *radiomEar, 
					unsigned short *radiomCnd, 
					unsigned short *radiomCal)
  {

    hsize_t start[6], count[6];

    if (blkSec != NULL) {
      // Block Sec (Seconds from start of day)
      start[0] = blknum;
      count[0] = 1;
      PTB( h5d_read(grp1, "blk_sec", blkSec, 1, start, count));

      // Bail if block second fill value found
      if ( *blkSec == -1)
	return 0;

      // Convert from time of day to TAI
      int32_t year, doy;
      hid_t attr;
      attr = H5Aopen_name(grp0, "Start Year");
      H5Aread(attr, H5T_STD_I32LE, &year);
      H5Aclose(attr);

      attr = H5Aopen_name(grp0, "Start Day");
      H5Aread(attr, H5T_STD_I32LE, &doy);
      H5Aclose(attr);

      // JMG
      // convert blkSec to seconds from 01/06/1980
      *blkSec = get_tai( year, doy, (*blkSec)*1000);
      // Kludge to correct when blkSec & granuleStart are virtually equal.
      if ( (*blkSec + 1.0e-6) < granuleStart) (*blkSec) += 86400.;
      //
    }

    // Radiometer header
    if (radiomHdr != NULL) {
      start[0] = frmnum;
      count[0]  = 1;
      start[1] = subframe;
      count[1]  = 1;
      PTB( h5d_read(grp2, "radiom_header", radiomHdr, 1, start, count));
    }

    // Radiometer signals
    if (radiomEar != NULL) {
      start[0] = frmnum;
      count[0]  = 1;
      start[1] = subframe;
      count[1]  = 1;

      start[3] = 0;
      count[3] = RADIOMETER_SIGNALS_PER_SUBCYCLE;
      start[4] = 0;
      count[4] = NUMBER_OF_BEAMS;
      start[5] = 0;
      count[5] = RADIOMETER_POLARIZATIONS;

      for (size_t isubcycl=0; isubcycl<RADIOMETER_SUBCYCLES; isubcycl++) {
	start[2] = isubcycl;
	count[2] = 1;
	PTB( h5d_read(grp2, "radiom_signals", 
		      &radiomEar[RADIOMETER_SIGNALS_PER_SUBCYCLE*
				 NUMBER_OF_BEAMS*
				 RADIOMETER_POLARIZATIONS*isubcycl], 6, 
		      start, count));
      }
    }

    // Radiometer short accumulations
    if (radiomCnd != NULL) {
      start[3] = 0;
      count[3] = NUMBER_OF_BEAMS;
      start[4] = 0;
      count[4] = RADIOMETER_POLARIZATIONS;

      for (size_t isubcycl=0; isubcycl<RADIOMETER_SUBCYCLES; isubcycl++) {
	start[2] = isubcycl;
	count[2] = 1;
	PTB( h5d_read(grp2, "radiom_cnd", 
		      &radiomCnd[NUMBER_OF_BEAMS*
				 RADIOMETER_POLARIZATIONS*isubcycl], 5, 
		      start, count));
      }
    }

    // Radiometer long accumulations
    if (radiomCal != NULL) {
      start[2] = 0;
      count[2] = RADIOMETER_LONG_ACCUM;
      start[3] = 0;
      count[3] = NUMBER_OF_BEAMS;
      start[4] = 0;
      count[4] = RADIOMETER_POLARIZATIONS;
      PTB( h5d_read(grp2, "radiom_lavg", &radiomCal[0], 5, start, count));
    }

    return 0;
  }


  /*------------------------------------------------------ */
  /* Read calibration temperatures                         */
  /* ----------------------------------------------------- */
  int hdf5_Aquarius::readl1_caltemps( int32_t blknum, 
				      int32_t frmnum, uint8_t subframe, 
				      float *caltemps)
  {
    // Add NaN warnings  JMG  09/02/11
    // Change NaN values to missing (-9999) JMG 09/05/11

    hsize_t start[6], count[6];

    float ext_temps[38];
    float rbe_temps[21];
    float rfe_temps[23];

    // idx values are from tlm_idx.txt in /disk01/aquarius/doc

    uint8_t caltemps_et_idx[29] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,
				   17,18,19,20,33,34,35,36,37,38,39,40};

    uint8_t exttemps_idx[29] = {35,36,37,0,1,2,3,4,5,18,19,20,21,22,23,
				12,13,14,15,16,17,27,28,29,30,31,32,33,34};

    start[0] = blknum;
    count[0] = 1;

    start[1] = 0;
    count[1] = 38;
    PTB( h5d_read(grp5, "ext_temp_analog_tlm", ext_temps, 2, start, count));

    for (size_t i=0; i<29; i++) {
      caltemps[caltemps_et_idx[i]] = ext_temps[exttemps_idx[i]];
      if ( isnan(ext_temps[exttemps_idx[i]])) {
	cout << "Value (0-based): " << setw(2) << right << 
	  ((int) exttemps_idx[i]) << 
	  " in \"ext_temp_analog_tlm\" is NaN for block number " << 
	  setw(4) << blknum << endl;
	caltemps[caltemps_et_idx[i]] = -9999.;
      }
    }


    // RBE temps
    uint8_t caltemps_rbe_idx[4] = {41,42,43,44};
    uint8_t rbetemps_idx[4] = {4,3,2,1};

    start[0] = frmnum;
    count[0] = 1;

    start[1] = 0;
    count[1] = 21;

    PTB( h5d_read(grp5, "rbe1_analog_tlm", rbe_temps, 2, start, count));
    for (size_t i=0; i<4; i++) {
      caltemps[caltemps_rbe_idx[i]] = rbe_temps[rbetemps_idx[i]];
      if ( isnan(rbe_temps[rbetemps_idx[i]])) {
	cout << "Value (0-based): " << setw(2) << right << 
	  ((int) rbetemps_idx[i]) << 
	  " in \"rbe1_analog_tlm\" is NaN for block number " << 
	  setw(4) << blknum << endl;
	caltemps[caltemps_rbe_idx[i]] = -9999.;
      }
    }

    PTB( h5d_read(grp5, "rbe2_analog_tlm", rbe_temps, 2, start, count));
    for (size_t i=0; i<4; i++) {
      caltemps[caltemps_rbe_idx[i]+4] = rbe_temps[rbetemps_idx[i]];
      if ( isnan(rbe_temps[rbetemps_idx[i]])) {
	cout << "Value (0-based): " << setw(2) << right << 
	  ((int) rbetemps_idx[i]) << 
	  " in \"rbe2_analog_tlm\" is NaN for block number " << 
	  setw(4) << blknum << endl;
	caltemps[caltemps_rbe_idx[i]+4] = -9999.;
      }
    }

    PTB( h5d_read(grp5, "rbe3_analog_tlm", rbe_temps, 2, start, count));
    for (size_t i=0; i<4; i++) {
      caltemps[caltemps_rbe_idx[i]+8] = rbe_temps[rbetemps_idx[i]];
      if ( isnan(rbe_temps[rbetemps_idx[i]])) {
	cout << "Value (0-based): " << setw(2) << right << 
	  ((int) rbetemps_idx[i]) << 
	  " in \"rbe3_analog_tlm\" is NaN for block number " << 
	  setw(4) << blknum << endl;
	caltemps[caltemps_rbe_idx[i]+8] = -9999.;
      }
    }


    // RFE temps
    uint8_t caltemps_rfe_idx[3][9] = {{21,22,23,24,70,71,72,73,74},
				      {25,26,27,28,75,76,77,78,79},
				      {29,30,31,32,80,81,82,83,84}};

    uint8_t rfetemps_idx[9] = {3,4,9,10,2,14,8,13,7};

    start[1] = 0;
    count[1] = 23;

    PTB( h5d_read(grp5, "rfe1_analog_tlm", rfe_temps, 2, start, count));
    for (size_t i=0; i<9; i++) {
      caltemps[caltemps_rfe_idx[0][i]] = rfe_temps[rfetemps_idx[i]];
      if ( isnan(rfe_temps[rfetemps_idx[i]])) {
	cout << "Value (0-based): " << setw(2) << right << 
	  ((int) rfetemps_idx[i]) << 
	  " in \"rfe1_analog_tlm\" is NaN for block number " << 
	  setw(4) << blknum << endl;
	caltemps[caltemps_rfe_idx[0][i]] = -9999.;
      }
    }

    PTB( h5d_read(grp5, "rfe2_analog_tlm", rfe_temps, 2, start, count));
    for (size_t i=0; i<9; i++) {
      caltemps[caltemps_rfe_idx[1][i]] = rfe_temps[rfetemps_idx[i]];
      if ( isnan(rfe_temps[rfetemps_idx[i]])) {
	cout << "Value (0-based): " << setw(2) << right << 
	  ((int) rfetemps_idx[i]) << 
	  " in \"rfe2_analog_tlm\" is NaN for block number " << 
	  setw(4) << blknum << endl;
	caltemps[caltemps_rfe_idx[1][i]] = -9999.;
      }
    }

    PTB( h5d_read(grp5, "rfe3_analog_tlm", rfe_temps, 2, start, count));
    for (size_t i=0; i<9; i++) {
      caltemps[caltemps_rfe_idx[2][i]] = rfe_temps[rfetemps_idx[i]];
      if ( isnan(rfe_temps[rfetemps_idx[i]])) {
	cout << "Value (0-based): " << setw(2) << right << 
	  ((int) rfetemps_idx[i]) << 
	  " in \"rfe3_analog_tlm\" is NaN for block number " << 
	  setw(4) << blknum << endl;
	caltemps[caltemps_rfe_idx[2][i]] = -9999.;
      }
    }

    return 0;
  }


  /*----------------------------------------------------------------- */
  /* Read block attribute fields                                      */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::readl1_frame_info( int32_t blknum, 
					int32_t *atc_frmnum, 
					uint8_t *atc_subframe,
					int32_t *rad_frmnum, 
					uint8_t *rad_subframe)
  {

    hsize_t start[6], count[6];

    start[0] = blknum;
    count[0] = 1;

    if ( atc_frmnum != NULL) 
      PTB( h5d_read(grp1, "atc_frmnum", atc_frmnum, 1, start, count));
    if ( atc_subframe != NULL) 
      PTB( h5d_read(grp1, "atc_subframe", atc_subframe, 1, start, count));
 
    if ( rad_frmnum != NULL) 
      PTB( h5d_read(grp1, "rad_frmnum", rad_frmnum, 1, start, count));
    if ( rad_subframe != NULL) 
      PTB( h5d_read(grp1, "rad_subframe", rad_subframe, 1, start, count));

    return 0;
  }


  /*----------------------------------------------------------------- */
  /* Write frame block attribute fields                               */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::writel1_frame_info( int32_t blknum, 
					 int32_t *atc_frmnum, 
					 uint8_t *atc_subframe,
					 int32_t *rad_frmnum, 
					 uint8_t *rad_subframe)
  {

    hsize_t start[6], count[6];

    start[0] = blknum;
    count[0] = 1;

    PTB( h5d_write(grp1, "atc_frmnum", atc_frmnum, 1, start, count));
    PTB( h5d_write(grp1, "atc_subframe", atc_subframe, 1, start, count));
 
    PTB( h5d_write(grp1, "rad_frmnum", rad_frmnum, 1, start, count));
    PTB( h5d_write(grp1, "rad_subframe", rad_subframe, 1, start, count));

    return 0;
  }


  /*----------------------------------------------------------------- */
  /* Read icds block number                                           */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::readl1_icds_blknum( int32_t blknum, int32_t *icds_blknum)
  {
    hsize_t start[6], count[6];

    start[0] = blknum;
    count[0] = 1;

    start[1] = 3;
    count[1] = 1;

    PTB( h5d_read(grp5, "icds_discrete_tlm", icds_blknum, 2, start, count));

    return 0;
  }


  /*----------------------------------------------------------------- */
  /* Read dpu status tlm                                              */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::readl1_dpu_status_tlm( int32_t rad_frmnum, 
					    int32_t rad_subframe,
					    int32_t status_byte,
					    uint8_t *dpu_status_tlm)
  {
    hsize_t start[6], count[6];

    start[0] = rad_frmnum;
    count[0] = 1;

    start[1] = rad_subframe;
    count[1] = 1;

    start[2] = status_byte;
    count[2] = 1;

    PTB( h5d_read(grp5, "dpu_status_tlm", dpu_status_tlm, 3, start, count));

    return 0;
  }


  /*----------------------------------------------------------------- */
  /* Read dpu status tlm                                              */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::readl1_radiom_nrt_tlm( int32_t rad_frmnum, 
					    int32_t rad_subframe,
					    int32_t status_byte,
					    uint8_t *radiom_nrt_tlm)
  {
    hsize_t start[6], count[6];

    start[0] = rad_frmnum;
    count[0] = 1;

    start[1] = rad_subframe;
    count[1] = 1;

    start[2] = status_byte;
    count[2] = 1;

    PTB( h5d_read(grp5, "radiom_nrt_tlm", radiom_nrt_tlm, 3, start, count));

    return 0;
  }


  /*----------------------------------------------------------------- */
  /* Create orbital fields in L1 file                                 */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::createl1_eph( int32_t numOrbVec)
  {
    // Orbit Time
    PTB( CreateH5D(
    grp4,                                            /* file id         */
    "orb_time",                                      /* short name      */
    "Time tag of orbit vectors",                     /* long name       */
    "seconds",                                       /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_DOUBLE,                               /* HDF number type */
    1,                                               /* rank            */
    numOrbVec, 1, 1, 1, 1, 1,                        /* dimension sizes */
    "Number of Orbit Vectors", 
    NULL, NULL, NULL, NULL, NULL,                    /* dimension names */
    H5P_DEFAULT) );

    // Orbit Position J2000
    PTB( CreateH5D(
    grp4,                                            /* file id         */
    "orb_pos",                                       /* short name      */
    "Orbital position vector",                       /* long name       */
    "meters",                                        /* units           */
    -7100000., 7100000.,                             /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_DOUBLE,                               /* HDF number type */
    2,                                               /* rank            */
    numOrbVec, NUMBER_OF_BEAMS, 1, 1, 1, 1,          /* dimension sizes */
    "Number of Orbit Vectors", 
    NULL, NULL, NULL, NULL, NULL,                    /* dimension names */
    H5P_DEFAULT) );

    // Orbit Velocity J2000
    PTB( CreateH5D(
    grp4,                                            /* file id         */
    "orb_vel",                                       /* short name      */
    "Orbital velocity vector",                       /* long name       */
    "meters per second",                             /* units           */
    -7600., 7600.,                                   /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_DOUBLE,                               /* HDF number type */
    2,                                               /* rank            */
    numOrbVec, 3, 1, 1, 1, 1,                        /* dimension sizes */
    "Number of Orbit Vectors", 
    NULL, NULL, NULL, NULL, NULL,                    /* dimension names */
    H5P_DEFAULT) );

    return 0;
  }



  /*----------------------------------------------------------------- */
  /* Write (all) orbital records to L1 file                           */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::writel1_eph( int32_t numOrbVec, 
				  double *time, double *pos, double *vel)
  {
    hsize_t start[6], count[6];

    this->createl1_eph( numOrbVec);

    start[0] = 0;
    count[0] = numOrbVec;
    PTB( h5d_write(grp4, "orb_time", time, 1, start, count));

    start[1] = 0;
    count[1] = 3;
    PTB( h5d_write(grp4, "orb_pos", pos, 2, start, count));
    PTB( h5d_write(grp4, "orb_vel", vel, 2, start, count));


    PTB( SetScalarH5A (grp0, "Number of Orbit Vectors" , H5T_STD_I32LE, 
		       (VOIDP) &numOrbVec));

    return 0;
  }



  /*----------------------------------------------------------------- */
  /* Set number of Orbit Vectors                                      */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::setOrbVec( int32_t numOrbVec) 
  {
    PTB( SetScalarH5A (grp0, "Number of Orbit Vectors" , H5T_STD_I32LE, 
		       (VOIDP) &numOrbVec));
 
    return 0;
  }



  /*----------------------------------------------------------------- */
  /* Create attitude/quaternion fields in L1 file                     */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::createl1_att( int32_t numAttSamp)
  {
    // Attitude Time
    PTB( CreateH5D(
    grp4,                                            /* file id         */
    "att_time",                                      /* short name      */
    "Time tag of attitude data",                     /* long name       */
    "seconds",                                       /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_DOUBLE,                               /* HDF number type */
    1,                                               /* rank            */
    numAttSamp, 1, 1, 1, 1, 1,                       /* dimension sizes */
    "Number of Attitude Samples", 
    NULL, NULL, NULL, NULL, NULL,                    /* dimension names */
    H5P_DEFAULT) );

    // Attitude Angles (roll, pitch, yaw)
    PTB( CreateH5D(
    grp4,                                            /* file id         */
    "att_ang",                                       /* short name      */
    "Spacecraft roll, pitch, yaw",                   /* long name       */
    "degrees",                                       /* units           */
    -180., 180.,                                     /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_DOUBLE,                               /* HDF number type */
    2,                                               /* rank            */
    numAttSamp, 3, 1, 1, 1, 1,                       /* dimension sizes */
    "Number of Attitude Samples", 
    NULL, NULL, NULL, NULL, NULL,                    /* dimension names */
    H5P_DEFAULT) );

    // Quaternions
    PTB( CreateH5D(
    grp4,                                            /* file id         */
    "quaternion",                                    /* short name      */
    "ECI-to-spacecraft quaternion",                  /* long name       */
    "",                                              /* units           */
    -1., 1.,                                         /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_DOUBLE,                               /* HDF number type */
    2,                                               /* rank            */
    numAttSamp, 4, 1, 1, 1, 1,                       /* dimension sizes */
    "Number of Attitude Samples", 
    NULL, NULL, NULL, NULL, NULL,                    /* dimension names */
    H5P_DEFAULT) );

    // Attitude Flags
    PTB( CreateH5D(
    grp4,                                            /* file id         */
    "att_flags",                                     /* short name      */
    "Attitude flags",                                /* long name       */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    H5T_STD_I32LE,                                   /* HDF number type */
    1,                                               /* rank            */
    numAttSamp, 1, 1, 1, 1, 1,                       /* dimension sizes */
    "Number of Attitude Samples", 
    NULL, NULL, NULL, NULL, NULL,                    /* dimension names */
    H5P_DEFAULT) );

    return 0;
  }



  /*----------------------------------------------------------------- */
  /* Write (all) attitude/quaternion records to L1 file               */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::writel1_att( int32_t numAttSamp, double *time, 
				  double *ang, double *quat, uint32_t *flags)
  {
    hsize_t start[6], count[6];

    this->createl1_att( numAttSamp);

    start[0] = 0;
    count[0] = numAttSamp;
    PTB( h5d_write(grp4, "att_time", time, 1, start, count));

    start[1] = 0;
    count[1] = 3;
    PTB( h5d_write(grp4, "att_ang", ang, 2, start, count));

    start[1] = 0;
    count[1] = 4;
    PTB( h5d_write(grp4, "quaternion", quat, 2, start, count));

    PTB( h5d_write(grp4, "att_flags", flags, 1, start, count));

    PTB( SetScalarH5A (grp0, "Number of Attitude Samples" , H5T_STD_I32LE, 
		       (VOIDP) &numAttSamp));

    return 0;
  }


  /*----------------------------------------------------------------- */
  /* Set number of Attitude Samples                                   */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::setAttSamp( int32_t numAttSamp) 
  {
    PTB( SetScalarH5A (grp0, "Number of Attitude Samples" , H5T_STD_I32LE, 
		       (VOIDP) &numAttSamp));

    return 0;
  }



  /*----------------------------------------------------------------- */
  /* Get number of orbital records in L1 file                         */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::readl1_eph( uint32_t *numOrbVec)
  {
    H5E_auto_t old_func;
    void *old_client_data;
    H5Eget_auto(H5E_DEFAULT , &old_func, &old_client_data);

    H5Eset_auto( H5E_DEFAULT, NULL, NULL); // Turn off error handling
    hid_t attr = H5Aopen_name(grp0, "Number of Orbit Vectors");
    H5Eset_auto( H5E_DEFAULT, old_func, old_client_data); // Restore previous error handler

    if ( attr < 0) return -1;

    H5Aread(attr, H5T_STD_I32LE, numOrbVec);
    H5Aclose(attr);

    return 0;
  }


  /*----------------------------------------------------------------- */
  /* Read orbital records from L1 file                                */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::readl1_eph( double *time, double *pos, double *vel)
  {
    hsize_t start[6], count[6];
    uint32_t numOrbVec;

    hid_t attr = H5Aopen_name(grp0, "Number of Orbit Vectors");
    H5Aread(attr, H5T_STD_I32LE, &numOrbVec);
    H5Aclose(attr);

    start[0] = 0;
    count[0] = numOrbVec;
    if (time != NULL)
      PTB( h5d_read(grp4, "orb_time", time, 1, start, count));

    // Convert from time of day to TAI
    int32_t year, doy;
    attr = H5Aopen_name(grp0, "Start Year");
    H5Aread(attr, H5T_STD_I32LE, &year);
    H5Aclose(attr);

    attr = H5Aopen_name(grp0, "Start Day");
    H5Aread(attr, H5T_STD_I32LE, &doy);
    H5Aclose(attr);

    // JMG
    // time is seconds from 01/06/1980
    // Day changes with leap seconds now correctly handled
    // Equivalent change made to readl1_att()
    // 07/05/12  JMG
    // Leap year fix 07/14/12  JMG
    double prevTime = -1;
    for (size_t i=0; i<numOrbVec; i++) {
      if ( time[i] < prevTime) {
	if ( (isleap(year) == 1 && doy == 366) ||
	     (isleap(year) == 0 && doy == 365)) {
	  doy = 1;
	  year++;
	} else {
	  doy++;
	}
      }
      prevTime = time[i];
      time[i] = get_tai( year, doy, time[i]*1000);
    }
    //    for (size_t i=0; i<numOrbVec; i++)
    // cout << "EPH: " << i << " " << (int32_t) time[i] << " " << endl;

    start[1] = 0;
    count[1] = 3;
    if (pos != NULL)
      PTB( h5d_read(grp4, "orb_pos", pos, 2, start, count));
    if (vel != NULL)
      PTB( h5d_read(grp4, "orb_vel", vel, 2, start, count));

    return 0;
  }



  /*----------------------------------------------------------------- */
  /* Get number of attitude samples in L1 file                        */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::readl1_att( uint32_t *numAttSamp)
  {
    hid_t attr = H5Aopen_name(grp0, "Number of Attitude Samples");
    H5Aread(attr, H5T_STD_I32LE, numAttSamp);
    H5Aclose(attr);

    return 0;
  }


  /*----------------------------------------------------------------- */
  /* Read attitude records from L1 file                               */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::readl1_att( double *time, double *ang, double *quat)
  {
    hsize_t start[6], count[6];
    uint32_t numAttSamp;

    hid_t attr = H5Aopen_name(grp0, "Number of Attitude Samples");
    H5Aread(attr, H5T_STD_I32LE, &numAttSamp);
    H5Aclose(attr);

    start[0] = 0;
    count[0] = numAttSamp;
    PTB( h5d_read(grp4, "att_time", time, 1, start, count));

    // Convert from time of day to TAI
    int32_t year, doy;
    attr = H5Aopen_name(grp0, "Start Year");
    H5Aread(attr, H5T_STD_I32LE, &year);
    H5Aclose(attr);

    attr = H5Aopen_name(grp0, "Start Day");
    H5Aread(attr, H5T_STD_I32LE, &doy);
    H5Aclose(attr);

    // JMG
    // time is seconds from 01/06/1980
    double prevTime = -1;
    for (size_t i=0; i<numAttSamp; i++) {
      if ( time[i] < prevTime) {
	if ( (isleap( (int) year) == 1 && doy == 366) ||
	     (isleap(year) == 0 && doy == 365)) {
	  doy = 1;
	  year++;
	} else {
	  doy++;
	}
      }
      prevTime = time[i];
      time[i] = get_tai( year, doy, time[i]*1000);
    }

    //    for (size_t i=0; i<numAttSamp; i++)
    //cout << "ATT: " << i << " " << (int32_t) time[i] << " " << endl;

    start[1] = 0;
    count[1] = 3;
    if (ang != NULL)
      PTB( h5d_read(grp4, "att_ang", ang, 2, start, count));

    if (quat != NULL) {
      count[1] = 4;
      PTB( h5d_read(grp4, "quaternion", quat, 2, start, count));
    }

    return 0;
  }


  /*----------------------------------------------------------------- */
  /* Copy L1A dataset                                                 */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::copyl1( int32_t recIdx, int32_t subfrmIdx, int32_t wrtIdx,
			     hid_t grp[2], char *dsName)
  {
#define COPYBUF 60000
    hsize_t start[6], count[6];
    static char buffer[COPYBUF];
    uint32_t rank;

    // Get dimention sizes (stored in count)
    PTB( h5d_read(grp[0], dsName, NULL, 0, NULL, count));

    // Compute start & count
    start[0] = recIdx;
    count[0] = 1;
    rank = 1;
    if (subfrmIdx != -1) {
      start[1] = subfrmIdx;
      count[1] = 1;
      rank = 2;
    }

    for (size_t i=rank; i<6; i++) {
      if (count[i] == 0) break;
      start[i] = 0;
      count[i] = count[i];
      rank++;
    }


    // Check if enough space (kludge)
    int32_t bufsz = count[0];
    for (size_t i=1; i<rank; i++) bufsz *= count[i];
    if (8*bufsz > COPYBUF) {
      printf("Potential buffer overflow in copyl1\n");
      exit (110);
    }

    PTB( h5d_read(grp[0],  dsName, buffer, rank, start, count));
    start[0] = wrtIdx;
    PTB( h5d_write(grp[1], dsName, buffer, rank, start, count));

    return 0;
  }


  /*----------------------------------------------------------------- */
  /* Get number of HKT records samples in L1 file                     */
  /*----------------------------------------------------------------- */
  int hdf5_Aquarius::readl1_hkt( uint32_t *numHktRec)
  {
    hsize_t start[2], count[2];

    count[0] = 1;

    start[1] = 580+816;
    count[1] = 4;

    uint32_t i=0;
    uint32_t gps_tme;

    uint8_t gps_hkt[4];
    uint8_t tmp[2];

    while (1) {
      start[0] = i;

      PTB( h5d_read(grp3, "sacd_hkt", (VOIDP) gps_hkt, 2, start, count));

      memcpy( tmp, &gps_hkt[0], 2);
      memcpy( &gps_hkt[0], &gps_hkt[2], 2);
      memcpy( &gps_hkt[2], tmp, 2);
      memcpy( &gps_tme, gps_hkt, 4);

      if ( gps_tme == 0) break; else i++;

      *numHktRec = i;
    }

    return 0;
  }


  /*----------------------------------------------------------------- */
  /* Get HKT time                                                     */
  /*----------------------------------------------------------------- */
  int hdf5_Aquarius::readl1_hkttme( uint32_t numHktRec, uint32_t *time)
  {
    hsize_t start[2], count[2];

    count[0] = 1;

    start[1] = 580+816;
    count[1] = 4;

    uint32_t gps_tme;

    uint8_t gps_hkt[4];
    uint8_t tmp[2];

    for (size_t i=0; i<numHktRec; i++) {
      start[0] = i;

      PTB( h5d_read(grp3, "sacd_hkt", (VOIDP) gps_hkt, 2, start, count));

      memcpy( tmp, &gps_hkt[0], 2);
      memcpy( &gps_hkt[0], &gps_hkt[2], 2);
      memcpy( &gps_hkt[2], tmp, 2);
      memcpy( &gps_tme, gps_hkt, 4);

      time[i] = gps_tme;
    }

    return 0;
  }


  /*----------------------------------------------------------------- */
  /* Get HKT time                                                     */
  /*----------------------------------------------------------------- */
  int hdf5_Aquarius::readl1_hktacsmode( uint32_t numHktSamp, uint8_t *acsmode)
  {
    hsize_t start[2], count[2];

    start[0] = 0;
    count[0] = numHktSamp;

    start[1] = 580;
    count[1] = 1;

    PTB( h5d_read(grp3, "sacd_hkt", (VOIDP) acsmode, 2, start, count));

    return 0;
  }


  /*----------------------------------------------------------------- */
  /* Write single HKT record to L1 file                               */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::writel1_hkt( uint32_t recnum, uint8_t *hktrec)
  {
    hsize_t start[2], count[2];

    start[0] = recnum;
    count[0] = 1;

    start[1] = 0;
    count[1] = 4000;

    PTB( h5d_write(grp3, "sacd_hkt", (VOIDP) hktrec, 2, start, count));

    return 0;
  }



  /*----------------------------------------------------------------- */
  /* Close L1 file                                                    */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::closel1()
  {

    if (grp0 != -1) H5Gclose(grp0);
    if (grp1 != -1) H5Gclose(grp1);
    if (grp2 != -1) H5Gclose(grp2);
    if (grp3 != -1) H5Gclose(grp3);
    if (grp4 != -1) H5Gclose(grp4);
    if (grp5 != -1) H5Gclose(grp5);

    H5Fclose(h5fid);

    return 0;
  }


  int hdf5_Aquarius::closel1( double granuleStop)
  {
    if (openFlags == H5F_ACC_RDWR) {
      // Add 1 to blknum to get Number of Blocks
      blknum++;

      PTB( SetScalarH5A (grp0, "Number of Blocks" , H5T_STD_I32LE, 
			 (VOIDP) &blknum));

      PTB( SetScalarH5A (grp0, "Number of RAD Frames" , H5T_STD_I32LE, 
			 (VOIDP) &rad_frmnum));

      PTB( SetScalarH5A (grp0, "Number of ATC Frames" , H5T_STD_I32LE, 
			 (VOIDP) &atc_frmnum));

      // Seconds from UTC 6 January 1980
      string ydhmsf_str( ydhmsf( gpstai2unix( granuleStop), 'G'));
      PTB( SetScalarH5A (grp0, "End Time", H5T_STRING, 
			 (void *) ydhmsf_str.c_str()));

      istringstream istr;
      int32_t itemp;
      ydhmsf_str = ydhmsf( gpstai2unix( granuleStop), 'G');
      int32_t millisec = get_millisec( &ydhmsf_str);

      istr.clear(); istr.str( ydhmsf_str.substr(0,4)); istr >> itemp;
      PTB( SetScalarH5A (grp0, "End Year", H5T_STD_I32LE, (VOIDP) &itemp));
      istr.clear(); istr.str( ydhmsf_str.substr(4,3)); istr >> itemp;
      PTB( SetScalarH5A (grp0, "End Day", H5T_STD_I32LE, (VOIDP) &itemp));
      PTB( SetScalarH5A (grp0, "End Millisec", H5T_STD_I32LE, 
			 (VOIDP) &millisec));

      //      PTB( SetScalarH5A (grp0, "Granule Stop Time" , 
      //		 H5T_NATIVE_DOUBLE, (VOIDP) &granuleStop));
    }
    if (grp0 != -1) H5Gclose(grp0);
    if (grp1 != -1) H5Gclose(grp1);
    if (grp2 != -1) H5Gclose(grp2);
    if (grp3 != -1) H5Gclose(grp3);
    if (grp4 != -1) H5Gclose(grp4);
    if (grp5 != -1) H5Gclose(grp5);

    H5Fclose(h5fid);

    return 0;
  }


  int hdf5_Aquarius::closel1( double granuleStart, double granuleStop,
			      double nodeCrossingTime, float nodeLongitude,
			      int32_t *oplut_accum, float percent_non_default,
			      int32_t missingBlks)
  {
    if (openFlags == H5F_ACC_RDWR) {
      PTB( SetScalarH5A (grp0, "Number of Blocks" , H5T_STD_I32LE, 
			 (VOIDP) &blknum));

      PTB( SetScalarH5A (grp0, "Number of RAD Frames" , H5T_STD_I32LE, 
			 (VOIDP) &rad_frmnum));

      PTB( SetScalarH5A (grp0, "Number of ATC Frames" , H5T_STD_I32LE, 
			 (VOIDP) &atc_frmnum));

      // Seconds from UTC 6 January 1980
      string ydhmsf_str( ydhmsf( gpstai2unix( granuleStop), 'G'));
      PTB( SetScalarH5A (grp0, "End Time", H5T_STRING, 
			 (void *) ydhmsf_str.c_str()));

      ydhmsf_str = ydhmsf( gpstai2unix( 0.5*(granuleStart+granuleStop)), 'G');
      PTB( SetScalarH5A (grp0, "Product Center Time", H5T_STRING, 
			 (void *) ydhmsf_str.c_str()));

      if ( nodeCrossingTime > 0) {
	ydhmsf_str = ydhmsf( gpstai2unix( nodeCrossingTime), 'G');
	PTB( SetScalarH5A (grp0, "Node Crossing Time", H5T_STRING,
			   (void *) ydhmsf_str.c_str())); 

	char *ocvarroot_str = getenv("OCVARROOT");
	if (ocvarroot_str == 0x0) {
	  printf("Environment variable OCVARROOT not defined.\n");
	  exit(1);
	}

	double nodeCrossTimes[32];
	uint32_t nodeCrossOrbits[32];
	string orbitCrossing = ocvarroot_str;
	orbitCrossing += "/aquarius/equatorial_crossing_times.txt";
	ifstream orbitCrossingFile;
	orbitCrossingFile.open( orbitCrossing.c_str());
	if ( !orbitCrossingFile.is_open()) {
	  cout << orbitCrossing.c_str() << " not found." << endl;
	  exit(1);
	}
	int32_t orbitNumber;
	double orbitTime;
	for (size_t i=0; i<32; i++) {
	  orbitCrossingFile >> nodeCrossOrbits[i];
	  orbitCrossingFile >> nodeCrossTimes[i];
	  if ( nodeCrossOrbits[i] == 0) break;
	}

	// Change to hard-coded parameters JMG  09/08/11
	// Make more robust  JMG  09/27/11
	int32_t orbit;
	for (size_t i=0; i<32; i++) {
	  if ( nodeCrossTimes[i] > nodeCrossingTime) break;

	  orbit = (int32_t) ((nodeCrossingTime + (ORBIT_AVG_PERIOD/2) - 
			      nodeCrossTimes[i]) / 
			     ORBIT_AVG_PERIOD) + nodeCrossOrbits[i];
	}

	PTB( SetScalarH5A (grp0, "Orbit Number", H5T_STD_I32LE, 
			   (VOIDP) &orbit));

	// cycle & pass 1-based
	int32_t cycle = 1 + (orbit - FIRSTORBIT) / NORBITCYCLE ;
	int32_t pass  = orbit - (cycle - 1) * NORBITCYCLE - FIRSTORBIT + 1;
	PTB( SetScalarH5A (grp0, "Cycle Number", H5T_STD_I32LE, 
			   (VOIDP) &cycle));
	PTB( SetScalarH5A (grp0, "Pass Number", H5T_STD_I32LE, 
			   (VOIDP) &pass));
      } else {
	PTB( SetScalarH5A (grp0, "Node Crossing Time", H5T_STRING,
			   (void *) "N/A")); 
      }
      istringstream istr;

      int32_t itemp;
      ydhmsf_str = ydhmsf( gpstai2unix( granuleStop), 'G');
      int32_t millisec = get_millisec( &ydhmsf_str);

      istr.clear(); istr.str( ydhmsf_str.substr(0,4)); istr >> itemp;
      PTB( SetScalarH5A (grp0, "End Year", H5T_STD_I32LE, (VOIDP) &itemp));
      istr.clear(); istr.str( ydhmsf_str.substr(4,3)); istr >> itemp;
      PTB( SetScalarH5A (grp0, "End Day", H5T_STD_I32LE, (VOIDP) &itemp));
      PTB( SetScalarH5A (grp0, "End Millisec", H5T_STD_I32LE, 
			 (VOIDP) &millisec));

      PTB( SetScalarH5A (grp0, "Orbit Node Longitude", H5T_NATIVE_FLOAT, 
			 (VOIDP) &nodeLongitude));

      PTB(h5a_set(grp0, "Radiometer LUTs Block Count", H5T_STD_I32LE, 8, 
		  (VOIDP) oplut_accum));

      PTB( SetScalarH5A (grp0, "Percent Non-default Radiometer LUTs", 
			 H5T_NATIVE_FLOAT, 
			 (VOIDP) &percent_non_default));

      PTB( SetScalarH5A (grp0, "Missing Blocks", H5T_STD_I32LE,
			 (VOIDP) &missingBlks));
    }

    if (grp0 != -1) H5Gclose(grp0);
    if (grp1 != -1) H5Gclose(grp1);
    if (grp2 != -1) H5Gclose(grp2);
    if (grp3 != -1) H5Gclose(grp3);
    if (grp4 != -1) H5Gclose(grp4);
    if (grp5 != -1) H5Gclose(grp5);

    H5Fclose(h5fid);

    return 0;
  }

} // namespace Hdf
