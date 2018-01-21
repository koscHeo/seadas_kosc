#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <libgen.h>


#include <timeutils.h>
#include <hdf.h>
#include "hdf5utils.h"
#include "hdf5_Aquarius.h"

#include "passthebuck.h"

#define DATACENTER "NASA/GSFC Aquarius Data Processing Center"
#define MISSIONCHARACTERISTICS "Nominal orbit: inclination=98.0 (Sun-synchronous); node=6PM (ascending); eccentricity=<0.002; altitude=657 km; ground speed=6.825 km/sec"

#define SENSORCHARACTERISTICS "Number of beams=3; channels per receiver=4; frequency 1.413 GHz; bits per sample=16; instatntaneous field of view=6.5 degrees; science data block period=1.44 sec."

double gpstai2unix( double gpstai);

namespace Hdf {

  /*----------------------------------------------------------------- */
  /* Create an HDF5 level2 file                                       */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::createl2( char* l2_filename, int32_t numBlocks,
			       int32_t nprod, 
			       int8_t *l2activeprodflag,
			       double granStart,   
			       double nodeCrossingTime, float nodeLongitude,
			       int32_t orbitNumber, int32_t cycleNumber, 
			       int32_t passNumber,
			       char* inputFiles, char* processControl,
			       char* processControlScat,
			       char *softver, char *pversion, bool iopt_l1b)


  {
    // Set object variable blknum
    blknum = numBlocks;

    /* Create the HDF file */
    h5fid = H5Fcreate(l2_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if(h5fid == FAIL) {
      fprintf(stderr,
	      "-E- %s line %d: Could not create HDF file, %s .\n",
	      __FILE__,__LINE__,l2_filename);
      return(HDF_FUNCTION_ERROR);
    }   

    openFlags = H5F_ACC_RDWR;
  
    grp0 = H5Gopen1(h5fid,"/");
    if ( iopt_l1b != 1) {
      grp1 = H5Gcreate1(h5fid, "/Block Attributes", 0);
    } else {
      grp1 = -1;
    }
    if ( iopt_l1b != 1) {
      grp2 = H5Gcreate1(h5fid, "/Aquarius Flags", 0);
    } else {
      grp2 = -1;
    }
    if (! iopt_l1b) {
      grp3 = H5Gcreate1(h5fid, "/Aquarius Data", 0);
    } else {
      grp3 = -1;
    }
    grp4 = H5Gcreate1(h5fid, "/Navigation", 0);
    grp5 = H5Gcreate1(h5fid, "/Converted Telemetry", 0);
    if ( iopt_l1b) {
      grp6 = H5Gcreate1(h5fid, "/Calibration", 0);
    } else {
      grp6 = -1;
    }

    //Write nProd object variable
    nProd = nprod;

    // Setup dataset create property list
    hid_t proplist = H5Pcreate( H5P_DATASET_CREATE);

    // Set fill to -1
    double time_fill_value = -1.0;
    H5Pset_fill_value( proplist, H5T_NATIVE_DOUBLE, (void *) &time_fill_value);

    if ( iopt_l1b != 1) {
      PTB( CreateH5D(
      grp1,                                            /* file id         */
      "sec",                                           /* short name      */
      "Block time, seconds of day",                    /* long name       */
      "",                                              /* standard name   */
      "seconds",                                       /* units           */
      0, 0,                                            /* valid range     */
      0, 0,                                            /* slope           */
      -9999,                                           /* fill value      */
      H5T_NATIVE_DOUBLE,                               /* HDF number type */
      1,                                               /* rank            */
      numBlocks, 1, 1, 1, 1, 1,                        /* dimension sizes */
      "Number of Blocks", NULL, NULL, NULL, NULL, NULL,/* dimension names */
      proplist) );

      PTB( CreateH5D(
      grp1,                                            /* file id         */
      "secGPS",                                        /* short name      */
      "Block time, GPS time",                          /* long name       */
      "",                                              /* standard name   */
      "seconds",                                       /* units           */
      0, 0,                                            /* valid range     */
      0, 0,                                            /* slope           */
      -9999,                                           /* fill value      */
      H5T_NATIVE_DOUBLE,                               /* HDF number type */
      1,                                               /* rank            */
      numBlocks, 1, 1, 1, 1, 1,                        /* dimension sizes */
      "Number of Blocks", NULL, NULL, NULL, NULL, NULL,/* dimension names */
      proplist) );

      PTB( CreateH5D(
      grp1,                                            /* file id         */
      "solar xray flux",                               /* short name      */
      "Solar xray flux (0.1 - 0.8 nanometers)",        /* long name       */
      "",                                              /* standard name   */
      "watts/m2",                                      /* units           */
      0, 0,                                            /* valid range     */
      0, 0,                                            /* slope           */
      -9999,                                           /* fill value      */
      H5T_NATIVE_FLOAT,                                /* HDF number type */
      1,                                               /* rank            */
      numBlocks, 1, 1, 1, 1, 1,                        /* dimension sizes */
      "Number of Blocks", NULL, NULL, NULL, NULL, NULL,/* dimension names */
      proplist) );
    }

    if ( ! iopt_l1b) {
      PTB( CreateH5D(
       grp1,                                            /* file id         */
       "rad_samples",                                   /* short name      */
       "Number of radiometer samples per average",      /* long name       */
       "",                                              /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       H5T_NATIVE_USHORT,                               /* HDF number type */
       3,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, RADIOMETER_POLARIZATIONS, 
       1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       "Number of Radiometer Polarizations", 
       NULL, NULL, NULL,
       H5P_DEFAULT) );

      PTB( CreateH5D(
       grp1,                                            /* file id         */
       "scat_samples",                                  /* short name      */
       "Number of scatterometer samples per average",   /* long name       */
       "",                                              /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       H5T_NATIVE_USHORT,                               /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TAV] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_TaV",                                       /* short name      */
       "Radiometer Ta V polarization",                  /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }


    if ( l2activeprodflag[TAH] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_TaH",                                       /* short name      */
       "Radiometer Ta H polarization",                  /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TA3] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_Ta3",                                       /* short name      */
       "Radiometer Ta 3rd Stokes",                      /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TFV] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_TfV",                                       /* short name      */
       "Radiometer Ta V polarization (rfi filtered)",   /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }


    if ( l2activeprodflag[TFH] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_TfH",                                       /* short name      */
       "Radiometer Ta H polarization (rfi filtered)",   /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TF3] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_Tf3",                                       /* short name      */
       "Radiometer Ta 3rd Stokes (rfi filtered)",       /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }


    if ( l2activeprodflag[TAV0] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_TaV0",                                      /* short name      */
       "Radiometer Ta V polar (no exp drift)",          /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }


    if ( l2activeprodflag[TAH0] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_TaH0",                                      /* short name      */
       "Radiometer Ta H polar (no exp drift)",          /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TA30] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_Ta30",                                      /* short name      */
       "Radiometer Ta 3rd Stokes (no exp drift)",       /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TFV0] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_TfV0",                                      /* short name      */
       "Radiometer Ta V (rfi filtered/no exp drift)",   /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }


    if ( l2activeprodflag[TFH0] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_TfH0",                                      /* short name      */
       "Radiometer Ta H (rfi filtered/no exp drift)",   /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TF30] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_Tf30",                                      /* short name      */
       "Radiometer Ta 3 (rfi filtered/no exp drift)",   /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TOIV] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_toi_V",                                     /* short name      */
       "Radiometer TOI Tb V polarization",              /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TOIH] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_toi_H",                                     /* short name      */
       "Radiometer TOI Tb H polarization",              /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TOI3] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_toi_3",                                     /* short name      */
       "Radiometer TOI Tb 3rd Stokes",                  /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TOAV] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_toa_V_nolc",                                /* short name      */
       "Radiometer TOA Tb V polarization (nolc)",       /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TOAH] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_toa_H_nolc",                                /* short name      */
       "Radiometer TOA Tb H polarization (nolc)",       /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TOAVLC] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_toa_V",                                     /* short name      */
       "Radiometer TOA Tb V polarization",              /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TOAHLC] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_toa_H",                                     /* short name      */
       "Radiometer TOA Tb H polarization",              /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[FARTAH] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_far_TaH",                                   /* short name      */
       "Radiometer Faraday Angle",                      /* long name       */
       "",                                              /* standard name   */
       "Degrees",                                       /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TAGALDV] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_galact_Ta_dir_V",                           /* short name      */
       "Radiometer Galactic Direct Corr V polar",       /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TAGALDH] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_galact_Ta_dir_H",                           /* short name      */
       "Radiometer Galactic Direct Corr H polar",       /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TAGALD3] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_galact_Ta_dir_3",                           /* short name      */
       "Radiometer Galactic Direct Corr 3rd Stokes",    /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TAGALRV] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_galact_Ta_ref_V",                           /* short name      */
       "Radiometer Galactic Reflect Corr V polar",      /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TAGALRH] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_galact_Ta_ref_H",                           /* short name      */
       "Radiometer Galactic Reflect Corr H polar",      /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TAGALR3] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_galact_Ta_ref_3",                           /* short name      */
       "Radiometer Galactic Reflect Corr 3rd Stokes",   /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TAGALRGOV] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_galact_Ta_ref_GO_V",                        /* short name      */
       "Radiometer Gal Reflect Corr V polar (GO)",      /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TAGALRGOH] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_galact_Ta_ref_GO_H",                        /* short name      */
       "Radiometer Gal Reflect Corr H polar (GO)",      /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[DTAGALV] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_galact_dTa_V",                              /* short name      */
       "Empirical Gal Reflect Adj V polar",             /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[DTAGALH] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_galact_dTa_H",                              /* short name      */
       "Empirical Gal Reflect Adj H polar",             /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TASUNDV] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_solar_Ta_dir_V",                            /* short name      */
       "Radiometer Solar Direct Corr V polar",          /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TASUNDH] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_solar_Ta_dir_H",                            /* short name      */
       "Radiometer Solar Direct Corr H polar",          /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TASUND3] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_solar_Ta_dir_3",                            /* short name      */
       "Radiometer Solar Direct Corr 3rd Stokes",       /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TASUNRV] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_solar_Ta_ref_V",                            /* short name      */
       "Radiometer Solar Reflect Corr V polar",         /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TASUNRH] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_solar_Ta_ref_H",                            /* short name      */
       "Radiometer Solar Reflect Corr H polar",         /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TASUNR3] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_solar_Ta_ref_3",                            /* short name      */
       "Radiometer Solar Reflect Corr 3rd Stokes",      /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }


    if ( l2activeprodflag[TASUNBV] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_solar_Ta_bak_V",                            /* short name      */
       "Radiometer Solar Back Scattered V polar",       /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TASUNBH] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_solar_Ta_bak_H",                            /* short name      */
       "Radiometer Solar Back Scattered H polar",       /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TASUNB3] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_solar_Ta_bak_3",                            /* short name      */
       "Radiometer Solar Back Scattered 3rd Stokes",    /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TAMONRV] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_moon_Ta_ref_V",                             /* short name      */
       "Radiometer Lunar Reflect Corr V polar",         /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TAMONRH] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_moon_Ta_ref_H",                             /* short name      */
       "Radiometer Lunar Reflect Corr H polar",         /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TAMONR3] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_moon_Ta_ref_3",                             /* short name      */
       "Radiometer Lunar Reflect Corr 3rd Stokes",      /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TBV] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_TbV_nolc",                                  /* short name      */
       "Earth surface Tb V polarization (no lc)",       /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TBH] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_TbH_nolc",                                  /* short name      */
       "Earth surface Tb H polarization (nolc)",        /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TBVRC] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_TbV_rc_nolc",                               /* short name      */
       "Tb V polarization (roughness corr, nolc)",      /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TBHRC] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_TbH_rc_nolc",                               /* short name      */
       "Tb H polarization (roughness corr, nolc)",      /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TBVRCLC] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_TbV_rc",                                    /* short name      */
       "Tb V polarization (rough corr)",                /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TBHRCLC] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_TbH_rc",                                    /* short name      */
       "Tb H polarization (rough corr)",                /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TBVNORCLC] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_TbV",                                       /* short name      */
       "Tb V polarization (no rough corr)",             /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TBHNORCLC] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_TbH",                                       /* short name      */
       "Tb H polarization (no rough corr)",             /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TBCON] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_Tb_consistency_nolc",                       /* short name      */
       "Tb consistency check (nolc)",                   /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TBCONLC] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_Tb_consistency",                            /* short name      */
       "Tb consistency check",                          /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[SSS] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "SSS_nolc",                                      /* short name      */
       "Sea Surface Salinity (nolc)",                   /* long name       */
       "sea_surface_salinity",                          /* standard name   */
       "PSU",                                           /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[SSSLC] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "SSS",                                           /* short name      */
       "Sea Surface Salinity",                          /* long name       */
       "sea_surface_salinity",                          /* standard name   */
       "PSU",                                           /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[SSSERRRAN] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "SSS_unc_ran",                                   /* short name      */
       "Sea Surface Salinity uncertainty (random)",     /* long name       */
       "",                                              /* standard name   */
       "PSU",                                           /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[SSSERRSYS] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "SSS_unc_sys",                                   /* short name      */
       "Sea Surface Salinity uncertainty (systematic)", /* long name       */
       "",                                              /* standard name   */
       "PSU",                                           /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[SSSADJ] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "SSS_bias_adj",                                  /* short name      */
       "Sea Surface Salinity bias adjusted",            /* long name       */
       "",                                              /* standard name   */
       "PSU",                                           /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[DENSITY] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "density" ,                                      /* short name      */
       "Derived Surface Density based on TEOS-10",      /* long name       */
       "sea_surface_density",                           /* standard name   */
       "kg/m^3",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[SPICE] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "spiciness" ,                                    /* short name      */
       "Derived Surface Spiciness based on TEOS-10",    /* long name       */
       "",                                              /* standard name   */
       "kg/m^3",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[EXPTAV] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_exp_TaV",                                   /* short name      */
       "Radiometer Ta V (expected)" ,                   /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }


    if ( l2activeprodflag[EXPTAH] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_exp_TaH",                                   /* short name      */
       "Radiometer Ta H (expected)",                    /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }


    if ( l2activeprodflag[EXPTA3] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_exp_Ta3",                                   /* short name      */
       "Radiometer Ta 3rd Stokes (expected)",           /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }


    if ( l2activeprodflag[EXPTAV_HHH] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_exp_TaV_hhh",                               /* short name      */
       "Radiometer Ta V (expected/HHH winds)" ,         /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }


    if ( l2activeprodflag[EXPTAH_HHH] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_exp_TaH_hhh",                               /* short name      */
       "Radiometer Ta H (expected/HHH winds)",          /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }


    if ( l2activeprodflag[EXPTA3_HHH] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_exp_Ta3_hhh",                               /* short name      */
       "Radiometer Ta 3rd Stokes (expected/HHH winds)", /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }


    if ( l2activeprodflag[EXPTBV] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_exp_TbV",                                   /* short name      */
       "Radiometer Tb V (expected)" ,                   /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[EXPTBH] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_exp_TbH",                                   /* short name      */
       "Radiometer Tb H (expected)",                    /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[EXPTBV0] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_exp_TbV0",                                  /* short name      */
       "Radiometer Tb V (expected smooth)" ,            /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[EXPTBH0] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_exp_TbH0",                                  /* short name      */
       "Radiometer Tb H (expected smooth)",             /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[HHWINDSPD] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_hh_wind_speed",                             /* short name      */
       "Radiometer HH wind speed",                      /* long name       */
       "",                                              /* standard name   */
       "meters/sec",                                    /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[HHHWINDSPD] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_hhh_wind_speed",                            /* short name      */
       "Radiometer HHH wind speed",                     /* long name       */
       "",                                              /* standard name   */
       "meters/sec",                                    /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[DTBSSTWSPDV] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_dtb_sst_wspd_V",                            /* short name      */
       "Radiometer SST bias emissivity correction V",   /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[DTBSSTWSPDH] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_dtb_sst_wspd_H",                            /* short name      */
       "Radiometer SST bias emissivity correction H",   /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[VVANT] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "scat_VV_ant",                                   /* short name      */
       "ANT Scatterometer NRCS for VV polarization",    /* long name       */
       "",                                              /* standard name   */
       "db",                                            /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[HHANT] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "scat_HH_ant",                                   /* short name      */
       "ANT Scatterometer NRCS for HH polarization",    /* long name       */
       "",                                              /* standard name   */
       "db",                                            /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[HVANT] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "scat_HV_ant",                                   /* short name      */
       "ANT Scatterometer NRCS for HV polarization",    /* long name       */
       "",                                              /* standard name   */
       "db",                                            /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[VHANT] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "scat_VH_ant",                                   /* short name      */
       "ANT Scatterometer NRCS for VH polarization",    /* long name       */
       "",                                              /* standard name   */
       "db",                                            /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[VVTOA] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "scat_VV_toa",                                   /* short name      */
       "TOA Scatterometer NRCS for VV polarization",    /* long name       */
       "",                                              /* standard name   */
       "db",                                            /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[HHTOA] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "scat_HH_toa",                                   /* short name      */
       "TOA Scatterometer NRCS for HH polarization",    /* long name       */
       "",                                              /* standard name   */
       "db",                                            /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[HVTOA] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "scat_HV_toa",                                   /* short name      */
       "TOA Scatterometer NRCS for HV polarization",    /* long name       */
       "",                                              /* standard name   */
       "db",                                            /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[VHTOA] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "scat_VH_toa",                                   /* short name      */
       "TOA Scatterometer NRCS for VH polarization",    /* long name       */
       "",                                              /* standard name   */
       "db",                                            /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[VVEXP] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "scat_VV_exp",                                   /* short name      */
       "Expected Sigma0 for VV polarization",           /* long name       */
       "",                                              /* standard name   */
       "db",                                            /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[HHEXP] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "scat_HH_exp",                                   /* short name      */
       "Expected Sigma0 for HH polarization",           /* long name       */
       "",                                              /* standard name   */
       "db",                                            /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                             /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[HVEXP] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "scat_HV_exp",                                   /* short name      */
       "Expected Sigma0 for HV polarization",           /* long name       */
       "",                                              /* standard name   */
       "db",                                            /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[VHEXP] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "scat_VH_exp",                                   /* short name      */
       "Expected Sigma0 for VH polarization",           /* long name       */
       "",                                              /* standard name   */
       "db",                                            /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TOTTOA] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "scat_tot_toa",                                  /* short name      */
       "TOA Scatterometer (Total)",                     /* long name       */
       "",                                              /* standard name   */
       "db",                                            /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }



    if ( l2activeprodflag[VVKPCANT] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "Kpc_VV_ant",                                    /* short name      */
       "Kpc statistical uncertainty for ANT VV NRCS",   /* long name       */
       "",                                              /* standard name   */
       "",                                              /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[HHKPCANT] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "Kpc_HH_ant",                                    /* short name      */
       "Kpc statistical uncertainty for ANT HH NRCS",   /* long name       */
       "",                                              /* standard name   */
       "",                                              /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[HVKPCANT] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "Kpc_HV_ant",                                    /* short name      */
       "Kpc statistical uncertainty for ANT HV NRCS",   /* long name       */
       "",                                              /* standard name   */
       "",                                              /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[VHKPCANT] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "Kpc_VH_ant",                                    /* short name      */
       "Kpc statistical uncertainty for ANT VH NRCS",   /* long name       */
       "",                                              /* standard name   */
       "",                                              /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }


    if ( l2activeprodflag[VVKPCTOA] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "Kpc_VV_toa",                                    /* short name      */
       "Kpc statistical uncertainty for TOA VV NRCS",   /* long name       */
       "",                                              /* standard name   */
       "",                                              /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[HHKPCTOA] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "Kpc_HH_toa",                                    /* short name      */
       "Kpc statistical uncertainty for TOA HH NRCS",   /* long name       */
       "",                                              /* standard name   */
       "",                                              /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[HVKPCTOA] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "Kpc_HV_toa",                                    /* short name      */
       "Kpc statistical uncertainty for TOA HV NRCS",   /* long name       */
       "",                                              /* standard name   */
       "",                                              /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[VHKPCTOA] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "Kpc_VH_toa",                                    /* short name      */
       "Kpc statistical uncertainty for TOA VH NRCS",   /* long name       */
       "",                                              /* standard name   */
       "",                                              /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[TOTKPCTOA] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "Kpc_total",                                     /* short name      */
       "Statistical uncertainty for total power NRCS",  /* long name       */
       "",                                              /* standard name   */
       "",                                              /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[SWINDUNC] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "wind_uncertainty",                              /* short name      */
       "Estimated wind speed error",                    /* long name       */
       "",                                              /* standard name   */
       "meters/sec",                                    /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }


    if ( l2activeprodflag[SWINDSPD] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "scat_wind_speed",                               /* short name      */
       "Scatterometer Wind Speed",                      /* long name       */
       "",                                              /* standard name   */
       "meters/sec",                                    /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }


    if ( l2activeprodflag[ESURFV] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "scat_esurf_V",                                  /* short name      */
       "excess surface scatterometer emissivity (V pol)",  /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[ESURFH] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "scat_esurf_H",                                  /* short name      */
       "excess surface scatterometer emissivity (H pol)",  /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[ESURFVUNC] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "scat_esurf_V_uncertainty",                      /* short name      */
       "Uncertainty in surface emmissivity (V-pol)",    /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[ESURFHUNC] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "scat_esurf_H_uncertainty",                      /* short name      */
       "Uncertainty in surface emmissivity (H-pol)",    /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }



    if ( l2activeprodflag[AWINDSPD] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "anc_wind_speed",                                /* short name      */
       "Ancillary Wind Speed 10m above surface",        /* long name       */
       "",                                              /* standard name   */
       "meters/sec",                                    /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }


    if ( l2activeprodflag[AWINDDIR] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "anc_wind_dir",                                  /* short name      */
       "Ancillary Wind Direction 10m above surface",    /* long name       */
       "",                                              /* standard name   */
       "degrees",                                       /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[ASURTEMP] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "anc_surface_temp",                              /* short name      */
       "Surface Temperature",                           /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[ASUBTEMP] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "anc_subsurf_temp",                              /* short name      */
       "Sub-Surface Temperature",                       /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[ASURP] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "anc_surface_pressure",                          /* short name      */
       "Surface Pressure",                              /* long name       */
       "",                                              /* standard name   */
       "Pascals",                                       /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[ACWAT] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "anc_cwat",                                      /* short name      */
       "Cloud Water",                                   /* long name       */
       "",                                              /* standard name   */
       "kg m-2",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[ASM] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "anc_sm",                                        /* short name      */
       "Soil Moisture",                                 /* long name       */
       "",                                              /* standard name   */
       "fraction",                                      /* units           */
       0, 1,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[ASWE] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "anc_swe",                                       /* short name      */
       "Snow Water Equivalent",                         /* long name       */
       "",                                              /* standard name   */
       "kg/m^2",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[ASSS] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "anc_SSS",                                       /* short name      */
       "Sea surface salinity (HYCOM)",                  /* long name       */
       "",                                              /* standard name   */
       "PSU",                                           /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[ATRANS] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "anc_trans",                                     /* short name      */
       "Atmospheric Transmittance",                     /* long name       */
       "",                                              /* standard name   */
       "",                                              /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[ATBUP] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "anc_Tb_up",                                     /* short name      */
       "Upwelling atmospheric brightness temperature",  /* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[ATBDW] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "anc_Tb_dw",                                     /* short name      */
       "Downwelling atmospheric brightness temperature",/* long name       */
       "",                                              /* standard name   */
       "Kelvin",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[ASWH] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "anc_swh",                                       /* short name      */
       "Significant wave height",                       /* long name       */
       "",                                              /* standard name   */
       "meters",                                        /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[RLANDFRC] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_land_frac",                                 /* short name      */
       "Fraction of land contamination (radiometer)",   /* long name       */
       "",                                              /* standard name   */
       "",                                              /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[RICEFRC] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_ice_frac",                                  /* short name      */
       "Fraction of ice contamination (radiometer)",    /* long name       */
       "",                                              /* standard name   */
       "",                                              /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[RGICE] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "rad_gice",                                      /* short name      */
       "Gain-smoothed ice fraction (radiometer)",       /* long name       */
       "",                                              /* standard name   */
       "",                                              /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

   if ( l2activeprodflag[SLANDFRC] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "scat_land_frac",                                /* short name      */
       "Fraction of land contamination (scatterometer)",/* long name       */
       "",                                              /* standard name   */
       "",                                              /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    if ( l2activeprodflag[SICEFRC] == 1) {
      PTB( CreateH5D(
       grp3,                                            /* file id         */
       "scat_ice_frac",                                 /* short name      */
       "Fraction of ice contamination (scatterometer)", /* long name       */
       "",                                              /* standard name   */
       "",                                              /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -999,                                            /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
    }
 
     if ( ! iopt_l1b) {
      PTB( CreateH5D(
       grp2,                                            /* file id         */
       "radiometer_flags",                              /* short name      */
       "Radiometer data quality flags",                 /* long name       */
       "",                                              /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       H5T_STD_U32LE,                                   /* HDF number type */
       3,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 4, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       "Max Radiometer Flags", NULL, NULL, NULL,
       H5P_DEFAULT) );

      PTB( CreateH5D(
       grp2,                                            /* file id         */
       "rad_rfi_flags",                                 /* short name      */
       "Radiometer RFI flags",                          /* long name       */
       "",                                              /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       H5T_NATIVE_UCHAR,                                /* HDF number type */
       4,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 4, RADIOMETER_SUBCYCLES, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       "Number of Polarizations", 
       "Number of Radiometer Subcycles", NULL, NULL,
       H5P_DEFAULT) );

      PTB( CreateH5D(
       grp2,                                            /* file id         */
       "scatterometer_flags",                           /* short name      */
       "Scatterometer data quality flags",              /* long name       */
       "",                                              /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       H5T_STD_U32LE,                                   /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );

      PTB( CreateH5D(
       grp2,                                            /* file id         */
       "scat_rfi_flags",                                /* short name      */
       "Scatterometer RFI flags",                       /* long name       */
       "",                                              /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       H5T_NATIVE_UCHAR,                                /* HDF number type */
       3,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 2, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       "Number of Polarizations", NULL, NULL, NULL,
       H5P_DEFAULT) );
    }

    PTB( CreateH5D(
     grp5,                                            /* file id         */
     "rad_caltemps",                                  /* short name      */
     "Radiometer calibration temperatures",           /* long name       */
     "",                                              /* standard name   */
     "",                                              /* units           */
     0, 0,                                            /* valid range     */
     0, 0,                                            /* slope           */
     -9999,                                           /* fill value      */
     H5T_NATIVE_FLOAT,                                /* HDF number type */
     2,                                               /* rank            */
     numBlocks,                                       /* dimension sizes */
     NUMCALTEMPS , 1, 1, 1, 1, 
     "Number of Blocks", "Number of CalTemps",        /* dimension names */
     NULL, NULL, NULL, NULL,
     H5P_DEFAULT) );


   if ( ! iopt_l1b) {
      PTB( CreateH5D(
       grp5,                                            /* file id         */
       "rad_gvv",                                       /* short name      */
       "Radiometer VV gain",                            /* long name       */
       "",                                              /* standard name   */
       "",                                              /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );

      PTB( CreateH5D(
       grp5,                                            /* file id         */
       "rad_ghh"     ,                                  /* short name      */
       "Radiometer HH gain",                            /* long name       */
       "",                                              /* standard name   */
       "",                                              /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );

      PTB( CreateH5D(
       grp5,                                            /* file id         */
       "rad_gpp",                                       /* short name      */
       "Radiometer PP gain",                            /* long name       */
       "",                                              /* standard name   */
       "",                                              /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );

      PTB( CreateH5D(
       grp5,                                            /* file id         */
       "rad_gmm"     ,                                  /* short name      */
       "Radiometer MM gain",                            /* long name       */
       "",                                              /* standard name   */
       "",                                              /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
 
      PTB( CreateH5D(
       grp5,                                            /* file id         */
       "rad_ov",                                        /* short name      */
       "Radiometer V offset",                           /* long name       */
       "",                                              /* standard name   */
       "",                                              /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );

      PTB( CreateH5D(
       grp5,                                            /* file id         */
       "rad_oh",                                        /* short name      */
       "Radiometer H offset",                           /* long name       */
       "",                                              /* standard name   */
       "",                                              /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );

      PTB( CreateH5D(
       grp5,                                            /* file id         */
       "rad_op",                                        /* short name      */
       "Radiometer P offset",                           /* long name       */
       "",                                              /* standard name   */
       "",                                              /* units           */
       0, 0,                                            /* valid range     */
       0, 0,                                            /* slope           */
       -9999,                                           /* fill value      */
       H5T_NATIVE_FLOAT,                                /* HDF number type */
       2,                                               /* rank            */
       numBlocks,                                       /* dimension sizes */
       NUMBER_OF_BEAMS, 1, 1, 1, 1, 
       "Number of Blocks", "Number of Beams",           /* dimension names */
       NULL, NULL, NULL, NULL,
       H5P_DEFAULT) );
  
       PTB( CreateH5D(
        grp5,                                            /* file id         */
        "rad_om",                                        /* short name      */
        "Radiometer M offset",                           /* long name       */
       "",                                              /* standard name   */
        "",                                              /* units           */
        0, 0,                                            /* valid range     */
        0, 0,                                            /* slope           */
        -9999,                                           /* fill value      */
        H5T_NATIVE_FLOAT,                                /* HDF number type */
        2,                                               /* rank            */
        numBlocks,                                       /* dimension sizes */
        NUMBER_OF_BEAMS, 1, 1, 1, 1, 
        "Number of Blocks", "Number of Beams",           /* dimension names */
        NULL, NULL, NULL, NULL,
        H5P_DEFAULT) );
    }

    H5Pclose( proplist);


    // Write flag attributes
    hid_t dataset;
    if ( ! iopt_l1b) {
      dataset = H5Dopen1(grp2, "radiometer_flags");
      PTB( SetScalarH5A (dataset, "RFI contamination", (hid_t) H5T_STRING, (VOIDP) "RFI"));
      PTB( SetScalarH5A (dataset, "Rain in main beam", (hid_t) H5T_STRING, (VOIDP) "RAIN"));
      PTB( SetScalarH5A (dataset, "Land contamination", (hid_t) H5T_STRING, (VOIDP) "LAND"));
      PTB( SetScalarH5A (dataset, "Sea ice contamination", (hid_t) H5T_STRING, (VOIDP) "ICE"));
      PTB( SetScalarH5A (dataset, "Wind/foam contamination", (hid_t) H5T_STRING, (VOIDP) "WIND"));
      PTB( SetScalarH5A (dataset, "Unusual brighness temperature", (hid_t) H5T_STRING, (VOIDP) "TEMP"));
      PTB( SetScalarH5A (dataset, "Direct solar flux contamination", (hid_t) H5T_STRING, (VOIDP) "FLUXD"));
      PTB( SetScalarH5A (dataset, "Reflected solar flux contamination", (hid_t) H5T_STRING, (VOIDP) "FLUXR"));
      PTB( SetScalarH5A (dataset, "Sun glint", (hid_t) H5T_STRING, (VOIDP) "SUNGLINT"));
      PTB( SetScalarH5A (dataset, "Moon contamination", (hid_t) H5T_STRING, (VOIDP) "MOON"));
      PTB( SetScalarH5A (dataset, "Galactic contamination", (hid_t) H5T_STRING, (VOIDP) "GALACTIC"));
      PTB( SetScalarH5A (dataset, "Non-nominal navigation", (hid_t) H5T_STRING, (VOIDP) "NAV"));
      PTB( SetScalarH5A (dataset, "SA overflow", (hid_t) H5T_STRING, (VOIDP) "SAOVERFLOW"));
      PTB( SetScalarH5A (dataset, "Roughness correction failure", (hid_t) H5T_STRING, (VOIDP) "ROUGH"));
      PTB( SetScalarH5A (dataset, "Solar flare contamination", (hid_t) H5T_STRING, (VOIDP) "FLARE"));
      PTB( SetScalarH5A (dataset, "Pointing anomaly", (hid_t) H5T_STRING, (VOIDP) "POINTING"));
      PTB( SetScalarH5A (dataset, "Tb consistency", (hid_t) H5T_STRING, (VOIDP) "TBCONS"));
      PTB( SetScalarH5A (dataset, "Cold water", (hid_t) H5T_STRING, (VOIDP) "COLDWATER"));
      PTB( SetScalarH5A (dataset, "RFI level", (hid_t) H5T_STRING, (VOIDP) "TFTADIFF"));
      PTB( SetScalarH5A (dataset, "Moon/Galaxy contamination", (hid_t) H5T_STRING, (VOIDP) "REFL_1STOKES"));
      PTB( SetScalarH5A (dataset, "RFI regional contamination", (hid_t) H5T_STRING, (VOIDP) "RFI_REGION"));

      dataset = H5Dopen1(grp2, "scatterometer_flags");
      PTB( SetScalarH5A (dataset, "RFI corruption of signal", (hid_t) H5T_STRING, (VOIDP) "RFI"));
      PTB( SetScalarH5A (dataset, "Rain in main beam", (hid_t) H5T_STRING, (VOIDP) "RAIN"));
      PTB( SetScalarH5A (dataset, "Negative TOA sigma0", (hid_t) H5T_STRING, (VOIDP) "NEGSIG"));
      PTB( SetScalarH5A (dataset, "Non-nominal attitude", (hid_t) H5T_STRING, (VOIDP) "BADATT"));
      PTB( SetScalarH5A (dataset, "Faraday rotation removal", (hid_t) H5T_STRING, (VOIDP) "FARADAY"));
    }

    /*                                                                  */
    /* Write out some global attributes                                 */
    /* ---------------------------------------------------------------- */
    /*                                                                  */
    PTB( SetScalarH5A (grp0, "Product Name", H5T_STRING, 
		       basename(l2_filename)));
    PTB( SetScalarH5A (grp0, "Title", H5T_STRING, 
		       (VOIDP) "Aquarius Level 2 Data"));
    PTB( SetScalarH5A (grp0, "Data Center", H5T_STRING, (VOIDP) DATACENTER));
    PTB( SetScalarH5A (grp0, "Mission", H5T_STRING, (VOIDP) "SAC-D Aquarius"));
    PTB( SetScalarH5A (grp0, "Mission Characteristics", H5T_STRING, 
		       (VOIDP) MISSIONCHARACTERISTICS));
    PTB( SetScalarH5A (grp0, "Sensor" , H5T_STRING, (VOIDP) "Aquarius"));
    PTB( SetScalarH5A (grp0, "Sensor Characteristics" , H5T_STRING, 
		       (VOIDP) SENSORCHARACTERISTICS));
    PTB( SetScalarH5A (grp0, "Data Type", H5T_STRING, (VOIDP) "SCI"));
    PTB( SetScalarH5A (grp0, "Software ID", H5T_STRING, (VOIDP) softver));
    PTB( SetScalarH5A (grp0, "Processing Version", H5T_STRING, 
		       (VOIDP) pversion));
    PTB( SetScalarH5A (grp0, "Processing Time", H5T_STRING, 
		       ydhmsf(time(NULL),'G')));
    PTB( SetScalarH5A (grp0, "Input Files", H5T_STRING, (VOIDP) inputFiles));
    PTB( SetScalarH5A (grp0, "Processing Control", H5T_STRING, 
		       (VOIDP) processControl));

    PTB( SetScalarH5A (grp0, "Scatterometer Processing Control", H5T_STRING, 
		       (VOIDP) processControlScat));

    PTB( SetScalarH5A (grp0, "_lastModified", H5T_STRING, 
		       ydhmsf(time(NULL),'G')));
    PTB( SetScalarH5A (grp0, "Conventions", H5T_STRING, (VOIDP) "CF-1.6")); 
    PTB( SetScalarH5A (grp0, "institution", H5T_STRING, (VOIDP) "NASA/GSFC OBPG")); 

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
    //    itemp = RADIOMETER_LONG_ACCUM;
    //PTB( SetScalarH5A (grp0, "Radiometer Long Accumulations" , 
    //		       H5T_STD_I32LE, (VOIDP) &itemp));
    itemp = SCATTEROMETER_POLARIZATIONS;
    PTB( SetScalarH5A (grp0, "Scatterometer Polarizations" , 
		       H5T_STD_I32LE, (VOIDP) &itemp));
    itemp = SCATTEROMETER_SUBCYCLES;
    PTB( SetScalarH5A (grp0, "Scatterometer Subcycles" , 
		       H5T_STD_I32LE, (VOIDP) &itemp));
    //    itemp = 2;
    //PTB( SetScalarH5A (grp0, "Flag Polarizations" , 
    //		       H5T_STD_I32LE, (VOIDP) &itemp));

    // granuleStart is in seconds from UTC 6 January 1980
    granuleStart = granStart; // Set granule start in private class variable
    string ydhmsf_str( ydhmsf( gpstai2unix( granuleStart), 'G'));
    PTB( SetScalarH5A (grp0, "Start Time", H5T_STRING, 
		       (void *) ydhmsf_str.c_str()));

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

    ydhmsf_str = ydhmsf( gpstai2unix(  0.5*(granuleStart+granuleStop)), 'G');
    // Removed  JMG  09/09/11
    //    PTB( SetScalarH5A (grp0, "Product Center Time", H5T_STRING, 
    //		       (void *) ydhmsf_str.c_str()));

    ydhmsf_str = ydhmsf( gpstai2unix(  nodeCrossingTime), 'G');
    PTB( SetScalarH5A (grp0, "Node Crossing Time", H5T_STRING,
		       (void *) ydhmsf_str.c_str())); 

    PTB( SetScalarH5A (grp0, "Orbit Node Longitude", H5T_NATIVE_FLOAT, 
		       (VOIDP) &nodeLongitude));


    PTB( SetScalarH5A (grp0, "Latitude Units", H5T_STRING,
		       (void *) "degrees North")); 

    PTB( SetScalarH5A (grp0, "Longitude Units", H5T_STRING,
		       (void *) "degrees East")); 

    PTB( SetScalarH5A (grp0, "Orbit Number", H5T_STD_I32LE,
		       (VOIDP) &orbitNumber));

    PTB( SetScalarH5A (grp0, "Cycle Number", H5T_STD_I32LE,
		       (VOIDP) &cycleNumber));
    
    PTB( SetScalarH5A (grp0, "Pass Number", H5T_STD_I32LE,
		       (VOIDP) &passNumber));

    return 0;
  }


  /*----------------------------------------------------------------- */
  /* Write navigation records to L2 file                              */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::writel2_nav( int32_t nBlks, bool *inout, 
				  double *pos, double *vel,
				  double *ang,
				  float *clat, float *clon,
				  float *celtht, float *celphi, 
				  float *suntht, float *sunphi, 
				  float *sunglt, float *moonglt, 
				  float *glxlat, float *glxlon,
				  float *cellatfoot, float *cellonfoot,
				  double *sund, double *sunr, double *moond,
				  double *zang,
				  double *sclon, double *sclat, double *scalt,
				  float *scat_clon, float *scat_clat,
				  float *scat_elon, float *scat_elat,
				  float *scat_polarization_roll,
				  uint8_t *acsmode,
				  bool browse, bool iopt_l1b)
  {
    hsize_t start[6], count[6];

    int32_t nActiveBlks = 0;
    for (int32_t iblk=0; iblk<nBlks; iblk++) {
      if ( inout[iblk] == true) nActiveBlks++;
    } 

    // Orbit Position J2000
    PTB( CreateH5D(
    grp4,                                            /* file id         */
    "orb_pos",                                       /* short name      */
    "Orbital position vector",                       /* long name       */
    "",                                              /* standard name   */
    "meters",                                        /* units           */
    -7100000., 7100000.,                             /* valid range     */
    0, 0,                                            /* slope           */
    -9999,                                           /* fill value      */
    H5T_NATIVE_DOUBLE,                               /* HDF number type */
    2,                                               /* rank            */
    nActiveBlks, NUMBER_OF_BEAMS, 1, 1, 1, 1 ,       /* dimension sizes */
    "Number of Blocks", "Number of Beams", 
    NULL, NULL, NULL, NULL,                          /* dimension names */
    H5P_DEFAULT) );

    // Orbit Velocity J2000
    PTB( CreateH5D(
    grp4,                                            /* file id         */
    "orb_vel",                                       /* short name      */
    "Orbital velocity vector",                       /* long name       */
    "",                                              /* standard name   */
    "meters per second",                             /* units           */
    -7600., 7600.,                                   /* valid range     */
    0, 0,                                            /* slope           */
    -9999,                                           /* fill value      */
    H5T_NATIVE_DOUBLE,                               /* HDF number type */
    2,                                               /* rank            */
    nActiveBlks, NUMBER_OF_BEAMS, 1, 1, 1, 1 ,       /* dimension sizes */
    "Number of Blocks", "Number of Beams", 
    NULL, NULL, NULL, NULL,                          /* dimension names */
    H5P_DEFAULT) );

    // Attitude Angles (roll, pitch, yaw)
    PTB( CreateH5D(
    grp4,                                            /* file id         */
    "att_ang",                                       /* short name      */
    "Spacecraft roll, pitch, yaw",                   /* long name       */
    "",                                              /* standard name   */
    "degrees",                                       /* units           */
    -180., 180.,                                     /* valid range     */
    0, 0,                                            /* slope           */
    -9999,                                           /* fill value      */
    H5T_NATIVE_DOUBLE,                               /* HDF number type */
    2,                                               /* rank            */
    nActiveBlks, NUMBER_OF_BEAMS, 1, 1, 1, 1 ,       /* dimension sizes */
    "Number of Blocks", "Number of Beams", 
    NULL, NULL, NULL, NULL,                          /* dimension names */
    H5P_DEFAULT) );

    // Beam center latitude
    PTB( CreateH5D(
    grp4,                                            /* file id         */
    "beam_clat",                                     /* short name      */
    "Beam Center Latitude",                          /* long name       */
    "",                                              /* standard name   */
    "degrees",                                       /* units           */
    -90., 90.,                                       /* valid range     */
    0, 0,                                            /* slope           */
    -9999,                                           /* fill value      */
    H5T_NATIVE_FLOAT,                                /* HDF number type */
    2,                                               /* rank            */
    nActiveBlks, NUMBER_OF_BEAMS, 1, 1, 1, 1 ,       /* dimension sizes */
    "Number of Blocks", "Number of Beams", 
    NULL, NULL, NULL, NULL,                          /* dimension names */
    H5P_DEFAULT) );


    // Beam center longitude
    PTB( CreateH5D(
    grp4,                                            /* file id         */
    "beam_clon",                                     /* short name      */
    "Beam Center Longitude",                         /* long name       */
    "",                                              /* standard name   */
    "degrees",                                       /* units           */
    -180., 180.,                                     /* valid range     */
    0, 0,                                            /* slope           */
    -9999,                                           /* fill value      */
    H5T_NATIVE_FLOAT,                                /* HDF number type */
    2,                                               /* rank            */
    nActiveBlks, NUMBER_OF_BEAMS, 1, 1, 1, 1 ,       /* dimension sizes */
    "Number of Blocks", "Number of Beams", 
    NULL, NULL, NULL, NULL,                          /* dimension names */
    H5P_DEFAULT) );

    // ACS control mode
    PTB( CreateH5D(
    grp4,                                            /* file id         */
    "acs_mode",                                      /* short name      */
    "ACS Control Mode",                              /* long name       */
    "",                                              /* units           */
    3,6,                                             /* valid range     */
    0, 0,                                            /* slope           */
    H5T_NATIVE_UCHAR,                                /* HDF number type */
    1,                                               /* rank            */
    nActiveBlks, 1, 1, 1, 1, 1 ,                     /* dimension sizes */
    "Number of Blocks", NULL, 
    NULL, NULL, NULL, NULL,                          /* dimension names */
    H5P_DEFAULT) );

    if ( ! browse) {
      // Boresight earth incidence angle
      PTB( CreateH5D(
      grp4,                                            /* file id         */
      "celtht",                                        /* short name      */
      "Boresight Earth Incidence Angle",               /* long name       */
      "",                                              /* standard name   */
      "degrees",                                       /* units           */
      -180., 180.,                                     /* valid range     */
      0, 0,                                            /* slope           */
      -9999,                                           /* fill value      */
      H5T_NATIVE_FLOAT,                                /* HDF number type */
      2,                                               /* rank            */
      nActiveBlks, NUMBER_OF_BEAMS, 1, 1, 1, 1 ,       /* dimension sizes */
      "Number of Blocks", "Number of Beams", 
      NULL, NULL, NULL, NULL,                          /* dimension names */
      H5P_DEFAULT) );

      // Boresight earth projection angle
      PTB( CreateH5D(
      grp4,                                            /* file id         */
      "celphi",                                        /* short name      */
      "Boresight Earth Azimuth Angle",                 /* long name       */
      "",                                              /* standard name   */
      "degrees",                                       /* units           */
      -180., 180.,                                     /* valid range     */
      0, 0,                                            /* slope           */
      -9999,                                           /* fill value      */
      H5T_NATIVE_FLOAT,                                /* HDF number type */
      2,                                               /* rank            */
      nActiveBlks, NUMBER_OF_BEAMS, 1, 1, 1, 1 ,       /* dimension sizes */
      "Number of Blocks", "Number of Beams", 
      NULL, NULL, NULL, NULL,                          /* dimension names */
      H5P_DEFAULT) );

      // Sun vector earth incidence angle
      PTB( CreateH5D(
      grp4,                                            /* file id         */
      "suntht",                                        /* short name      */
      "Sun Vector Earth Incidence Angle",              /* long name       */
      "",                                              /* standard name   */
      "degrees",                                       /* units           */
      -180., 180.,                                     /* valid range     */
      0, 0,                                            /* slope           */
      -9999,                                           /* fill value      */
      H5T_NATIVE_FLOAT,                                /* HDF number type */
      2,                                               /* rank            */
      nActiveBlks, NUMBER_OF_BEAMS, 1, 1, 1, 1 ,       /* dimension sizes */
      "Number of Blocks", "Number of Beams", 
      NULL, NULL, NULL, NULL,                          /* dimension names */
      H5P_DEFAULT) );

      // Sun vector earth projection angle
      PTB( CreateH5D(
      grp4,                                            /* file id         */
      "sunphi",                                        /* short name      */
      "Sun Vector Earth Azimuth Angle",                /* long name       */
      "",                                              /* standard name   */
      "degrees",                                       /* units           */
      -180., 180.,                                     /* valid range     */
      0, 0,                                            /* slope           */
      -9999,                                           /* fill value      */
      H5T_NATIVE_FLOAT,                                /* HDF number type */
      2,                                               /* rank            */
      nActiveBlks, NUMBER_OF_BEAMS, 1, 1, 1, 1 ,       /* dimension sizes */
      "Number of Blocks", "Number of Beams", 
      NULL, NULL, NULL, NULL,                          /* dimension names */
      H5P_DEFAULT) );

      // Sun glint angle
      PTB( CreateH5D(
      grp4,                                            /* file id         */
      "sunglt",                                        /* short name      */
      "Sun Glint Angle",                               /* long name       */
      "",                                              /* standard name   */
      "degrees",                                       /* units           */
      -180., 180.,                                     /* valid range     */
      0, 0,                                            /* slope           */
      -9999,                                           /* fill value      */
      H5T_NATIVE_FLOAT,                                /* HDF number type */
      2,                                               /* rank            */
      nActiveBlks, NUMBER_OF_BEAMS, 1, 1, 1, 1 ,       /* dimension sizes */
      "Number of Blocks", "Number of Beams", 
      NULL, NULL, NULL, NULL,                          /* dimension names */
      H5P_DEFAULT) );

      // Moon glint angle
      PTB( CreateH5D(
      grp4,                                            /* file id         */
      "moonglt",                                       /* short name      */
      "Moon Glint Angle",                              /* long name       */
      "",                                              /* standard name   */
      "degrees",                                       /* units           */
      -180., 180.,                                     /* valid range     */
      0, 0,                                            /* slope           */
      -9999,                                           /* fill value      */
      H5T_NATIVE_FLOAT,                                /* HDF number type */
      2,                                               /* rank            */
      nActiveBlks, NUMBER_OF_BEAMS, 1, 1, 1, 1 ,       /* dimension sizes */
      "Number of Blocks", "Number of Beams", 
      NULL, NULL, NULL, NULL,                          /* dimension names */
      H5P_DEFAULT) );

      // Galaxy declination (J2000)
      PTB( CreateH5D(
      grp4,                                            /* file id         */
      "glxlat",                                        /* short name      */
      "Galaxy Declination (J2000)",                    /* long name       */
      "",                                              /* standard name   */
      "degrees",                                       /* units           */
      -90., 90.,                                       /* valid range     */
      0, 0,                                            /* slope           */
      -9999,                                           /* fill value      */
      H5T_NATIVE_FLOAT,                                /* HDF number type */
      2,                                               /* rank            */
      nActiveBlks, NUMBER_OF_BEAMS, 1, 1, 1, 1 ,       /* dimension sizes */
      "Number of Blocks", "Number of Beams", 
      NULL, NULL, NULL, NULL,                          /* dimension names */
      H5P_DEFAULT) );

      // Galaxy right ascention (J2000)
      PTB( CreateH5D(
      grp4,                                            /* file id         */
      "glxlon",                                        /* short name      */
      "Galaxy Right Ascention (J2000)",                /* long name       */
      "",                                              /* standard name   */
      "degrees",                                       /* units           */
      -180., 180.,                                     /* valid range     */
      0, 0,                                            /* slope           */
      -9999,                                           /* fill value      */
      H5T_NATIVE_FLOAT,                                /* HDF number type */
      2,                                               /* rank            */
      nActiveBlks, NUMBER_OF_BEAMS, 1, 1, 1, 1 ,       /* dimension sizes */
      "Number of Blocks", "Number of Beams", 
      NULL, NULL, NULL, NULL,                          /* dimension names */
      H5P_DEFAULT) );

      // Geodetic latitudes 3-db footprint
      PTB( CreateH5D(
      grp4,                                            /* file id         */
      "cellatfoot",                                    /* short name      */
      "Geodectic Latitudes (3 dB)",                    /* long name       */
      "",                                              /* standard name   */
      "degrees",                                       /* units           */
      -90., 90.,                                       /* valid range     */
      0, 0,                                            /* slope           */
      -9999,                                           /* fill value      */
      H5T_NATIVE_FLOAT,                                /* HDF number type */
      3,                                               /* rank            */
      nActiveBlks, NUMBER_OF_BEAMS, 4, 1, 1, 1 ,       /* dimension sizes */
      "Number of Blocks", "Number of Beams", 
      NULL, NULL, NULL, NULL,                          /* dimension names */
      H5P_DEFAULT) );

      // East longitudes 3-db footprint
      PTB( CreateH5D(
      grp4,                                            /* file id         */
      "cellonfoot",                                    /* short name      */
      "East Longitudes (3 dB)",                        /* long name       */
      "",                                              /* standard name   */
      "degrees",                                       /* units           */
      -180., 180.,                                     /* valid range     */
      0, 0,                                            /* slope           */ 
      -9999,                                           /* fill value      */
      H5T_NATIVE_FLOAT,                                /* HDF number type */
      3,                                               /* rank            */
      nActiveBlks, NUMBER_OF_BEAMS, 4, 1, 1, 1 ,       /* dimension sizes */
      "Number of Blocks", "Number of Beams", 
      NULL, NULL, NULL, NULL,                          /* dimension names */
      H5P_DEFAULT) );

      // Earth-to-sun unit vector in eci coordinates
      PTB( CreateH5D(
      grp4,                                            /* file id         */
      "sund",                                          /* short name      */
      "Earth-to-Sun unit vector (eci)",                /* long name       */
      "",                                              /* standard name   */
      "",                                              /* units           */
      0, 0,                                            /* valid range     */
      0, 0,                                            /* slope           */
      -9999,                                           /* fill value      */
      H5T_NATIVE_DOUBLE,                               /* HDF number type */
      2,                                               /* rank            */
      nActiveBlks, 3, 1, 1, 1, 1 ,                     /* dimension sizes */
      "Number of Blocks", 
      NULL, NULL, NULL, NULL, NULL,                    /* dimension names */
      H5P_DEFAULT) );

      // Spacecraft-sun reflection point on earth in eci coordinates
      PTB( CreateH5D(
      grp4,                                            /* file id         */
      "sunr",                                          /* short name      */
      "Sun reflection unit vector (eci)",              /* long name       */
      "",                                              /* standard name   */
      "",                                              /* units           */
      0, 0,                                            /* valid range     */
      0, 0,                                            /* slope           */
      -9999,                                           /* fill value      */
      H5T_NATIVE_DOUBLE,                               /* HDF number type */
      2,                                               /* rank            */
      nActiveBlks, 3, 1, 1, 1, 1 ,                     /* dimension sizes */
      "Number of Blocks", 
      NULL, NULL, NULL, NULL, NULL,                    /* dimension names */
      H5P_DEFAULT) );

      // Earth-to-moon unit vector in eci coordinates
      PTB( CreateH5D(
      grp4,                                            /* file id         */
      "moond",                                         /* short name      */
      "Earth-to-Moon unit vector (eci)",               /* long name       */
      "",                                              /* standard name   */
      "",                                              /* units           */
      0, 0,                                            /* valid range     */
      0, 0,                                            /* slope           */
      -9999,                                           /* fill value      */
      H5T_NATIVE_DOUBLE,                               /* HDF number type */
      2,                                               /* rank            */
      nActiveBlks, 3, 1, 1, 1, 1 ,                     /* dimension sizes */
      "Number of Blocks", 
      NULL, NULL, NULL, NULL, NULL,                    /* dimension names */
      H5P_DEFAULT) );

      // Intra-orbit angle
      PTB( CreateH5D(
      grp4,                                            /* file id         */
      "zang",                                          /* short name      */
      "Intra-Orbit Angle",                             /* long name       */
      "",                                              /* standard name   */
      "degrees",                                       /* units           */
      -180., 180.,                                     /* valid range     */
      0, 0,                                            /* slope           */
      -9999,                                           /* fill value      */
      H5T_NATIVE_DOUBLE,                               /* HDF number type */
      1,                                               /* rank            */
      nActiveBlks, 1, 1, 1, 1, 1 ,                     /* dimension sizes */
      "Number of Blocks",
      NULL, NULL, NULL, NULL, NULL,                    /* dimension names */
      H5P_DEFAULT) );

      // S/C longitude
      PTB( CreateH5D(
      grp4,                                            /* file id         */
      "sclon",                                         /* short name      */
      "Spacecraft nadir point longitude",              /* long name       */
      "",                                              /* standard name   */
      "degrees",                                       /* units           */
      -180., 180.,                                     /* valid range     */
      0, 0,                                            /* slope           */
      -9999,                                           /* fill value      */
      H5T_NATIVE_DOUBLE,                               /* HDF number type */
      1,                                               /* rank            */
      nActiveBlks, 1, 1, 1, 1, 1 ,                     /* dimension sizes */
      "Number of Blocks",
      NULL, NULL, NULL, NULL, NULL,                    /* dimension names */
      H5P_DEFAULT) );

      // S/C latitude
      PTB( CreateH5D(
      grp4,                                            /* file id         */
      "sclat",                                         /* short name      */
      "Spacecraft nadir point latitude",               /* long name       */
      "",                                              /* standard name   */
      "degrees",                                       /* units           */
      -90., 90.,                                       /* valid range     */
      0, 0,                                            /* slope           */
      -9999,                                           /* fill value      */
      H5T_NATIVE_DOUBLE,                               /* HDF number type */
      1,                                               /* rank            */
      nActiveBlks, 1, 1, 1, 1, 1 ,                     /* dimension sizes */
      "Number of Blocks",
      NULL, NULL, NULL, NULL, NULL,                    /* dimension names */
      H5P_DEFAULT) );

      // S/C altitude
      PTB( CreateH5D(
      grp4,                                            /* file id         */
      "scalt",                                         /* short name      */
      "Spacecraft altitude",                           /* long name       */
      "",                                              /* standard name   */
      "meters" ,                                       /* units           */
      0, 0,                                            /* valid range     */
      0, 0,                                            /* slope           */
      -9999,                                           /* fill value      */
      H5T_NATIVE_DOUBLE,                               /* HDF number type */
      1,                                               /* rank            */
      nActiveBlks, 1, 1, 1, 1, 1 ,                     /* dimension sizes */
      "Number of Blocks",
      NULL, NULL, NULL, NULL, NULL,                    /* dimension names */
      H5P_DEFAULT) );
    }

    if ( ! iopt_l1b) {
      // Scatterometer beam center latitude
      PTB( CreateH5D(
      grp4,                                            /* file id         */
      "scat_beam_clat",                                /* short name      */
      "Scatterometer Beam Center Latitude",            /* long name       */
      "",                                              /* standard name   */
      "degrees",                                       /* units           */
      -90., 90.,                                       /* valid range     */
      0, 0,                                            /* slope           */
      -9999,                                           /* fill value      */
      H5T_NATIVE_FLOAT,                                /* HDF number type */
      2,                                               /* rank            */
      nActiveBlks, NUMBER_OF_BEAMS, 1, 1, 1, 1 ,       /* dimension sizes */
      "Number of Blocks", "Number of Beams", 
      NULL, NULL, NULL, NULL,                          /* dimension names */
      H5P_DEFAULT) );

      // Scatterometer beam center longitude
      PTB( CreateH5D(
      grp4,                                            /* file id         */
      "scat_beam_clon",                                /* short name      */
      "Scatterometer Beam Center Longitude",           /* long name       */
      "",                                              /* standard name   */
      "degrees",                                       /* units           */
      -180., 180.,                                     /* valid range     */
      0, 0,                                            /* slope           */
      -9999,                                           /* fill value      */
      H5T_NATIVE_FLOAT,                                /* HDF number type */
      2,                                               /* rank            */
      nActiveBlks, NUMBER_OF_BEAMS, 1, 1, 1, 1 ,       /* dimension sizes */
      "Number of Blocks", "Number of Beams", 
      NULL, NULL, NULL, NULL,                          /* dimension names */
      H5P_DEFAULT) );

      // Scatterometer polarization roll angle
      PTB( CreateH5D(
      grp4,                                            /* file id         */
      "scat_polarization_roll",                        /* short name      */
      "Scatterometer Polarization Roll Angle",         /* long name       */
      "",                                              /* standard name   */
      "degrees",                                       /* units           */
      -180., 180.,                                     /* valid range     */
      0, 0,                                            /* slope           */
      -9999,                                           /* fill value      */
      H5T_NATIVE_FLOAT,                                /* HDF number type */
      2,                                               /* rank            */
      nActiveBlks, NUMBER_OF_BEAMS, 1, 1, 1, 1 ,       /* dimension sizes */
      "Number of Blocks", "Number of Beams", 
      NULL, NULL, NULL, NULL,                          /* dimension names */
      H5P_DEFAULT) );
    }

    if ( ! browse) {
      // Scatterometer latitude footprint
      PTB( CreateH5D(
      grp4,                                            /* file id         */
      "scat_latfoot",                                  /* short name      */
      "Scatterometer Latitude Footprint",              /* long name       */
      "",                                              /* standard name   */
      "degrees",                                       /* units           */
      -90., 90.,                                       /* valid range     */
      0, 0,                                            /* slope           */
      -9999,                                           /* fill value      */
      H5T_NATIVE_FLOAT,                                /* HDF number type */
      3,                                               /* rank            */
      nActiveBlks, NUMBER_OF_BEAMS, 4, 1, 1, 1 ,       /* dimension sizes */
      "Number of Blocks", "Number of Beams", 
      NULL, NULL, NULL, NULL,                          /* dimension names */
      H5P_DEFAULT) );

      // Scatterometer longitude footprint
      PTB( CreateH5D(
      grp4,                                            /* file id         */
      "scat_lonfoot",                                  /* short name      */
      "Scatterometer Longitude Footprint",             /* long name       */
      "",                                              /* standard name   */
      "degrees",                                       /* units           */
      -180., 180.,                                     /* valid range     */
      0, 0,                                            /* slope           */
      -9999,                                           /* fill value      */
      H5T_NATIVE_FLOAT,                                /* HDF number type */
      3,                                               /* rank            */
      nActiveBlks, NUMBER_OF_BEAMS, 4, 1, 1, 1 ,       /* dimension sizes */
      "Number of Blocks", "Number of Beams", 
      NULL, NULL, NULL, NULL,                          /* dimension names */
      H5P_DEFAULT) );
    } // if (! browse)

    if ( nBlks == nActiveBlks) {
      start[0] = 0;
      count[0] = blknum;

      start[1] = 0;
      count[1] = 3;
      PTB( h5d_write(grp4, "orb_pos", pos, 2, start, count));
      PTB( h5d_write(grp4, "orb_vel", vel, 2, start, count));
      PTB( h5d_write(grp4, "att_ang", ang, 2, start, count));
      
      PTB( h5d_write(grp4, "beam_clat", clat, 2, start, count));
      PTB( h5d_write(grp4, "beam_clon", clon, 2, start, count));

      PTB( h5d_write(grp4, "acs_mode", acsmode, 1, start, count));

      if ( ! browse) {
	PTB( h5d_write(grp4, "celtht", celtht, 2, start, count));
	PTB( h5d_write(grp4, "celphi", celphi, 2, start, count));

	PTB( h5d_write(grp4, "suntht", suntht, 2, start, count));
	PTB( h5d_write(grp4, "sunphi", sunphi, 2, start, count));

	PTB( h5d_write(grp4, "sunglt", sunglt, 2, start, count));
	PTB( h5d_write(grp4, "moonglt", moonglt, 2, start, count));

	PTB( h5d_write(grp4, "glxlat", glxlat, 2, start, count));
	PTB( h5d_write(grp4, "glxlon", glxlon, 2, start, count));
      }

      if ( ! iopt_l1b) {
	PTB( h5d_write(grp4, "scat_beam_clat", scat_clat, 2, start, count));
	PTB( h5d_write(grp4, "scat_beam_clon", scat_clon, 2, start, count));
	PTB( h5d_write(grp4, "scat_polarization_roll", scat_polarization_roll, 2, start, count));
      }

      start[2] = 0;
      count[2] = 4;

      if ( ! browse) {
	PTB( h5d_write(grp4, "cellatfoot", cellatfoot, 3, start, count));
	PTB( h5d_write(grp4, "cellonfoot", cellonfoot, 3, start, count));

	PTB( h5d_write(grp4, "scat_latfoot", scat_elat, 3, start, count));
	PTB( h5d_write(grp4, "scat_lonfoot", scat_elon, 3, start, count));

	PTB( h5d_write(grp4, "zang", zang, 1, start, count));

	PTB( h5d_write(grp4, "sclon", sclon, 1, start, count));
	PTB( h5d_write(grp4, "sclat", sclat, 1, start, count));
	PTB( h5d_write(grp4, "scalt", scalt, 1, start, count));
      }

    } else {
      int i=0;
      for (int32_t iblk=0; iblk<nBlks; iblk++) {
	if ( inout[iblk] == true) {
	  start[0] = i;
	  count[0] = 1;

	  start[1] = 0;
	  count[1] = 3;

	  PTB( h5d_write(grp4, "orb_pos", &pos[3*iblk], 2, start, count));
	  PTB( h5d_write(grp4, "orb_vel", &vel[3*iblk], 2, start, count));
	  PTB( h5d_write(grp4, "att_ang", &ang[3*iblk], 2, start, count));
      
	  PTB( h5d_write(grp4, "beam_clat", &clat[3*iblk], 2, 
			 start, count));
	  PTB( h5d_write(grp4, "beam_clon", &clon[3*iblk], 2, 
			 start, count));

	  PTB( h5d_write(grp4, "acs_mode", &acsmode[iblk], 1, start, count));

	  if ( ! browse) {
	    PTB( h5d_write(grp4, "celtht", &celtht[3*iblk], 2, start, count));
	    PTB( h5d_write(grp4, "celphi", &celphi[3*iblk], 2, start, count));

	    PTB( h5d_write(grp4, "suntht", &suntht[3*iblk], 2, start, count));
	    PTB( h5d_write(grp4, "sunphi", &sunphi[3*iblk], 2, start, count));

	    PTB( h5d_write(grp4, "sunglt", &sunglt[3*iblk], 2, start, count));
	    PTB( h5d_write(grp4, "moonglt", &moonglt[3*iblk], 2, start, count));

	    PTB( h5d_write(grp4, "glxlat", &glxlat[3*iblk], 2, start, count));
	    PTB( h5d_write(grp4, "glxlon", &glxlon[3*iblk], 2, start, count));
	  }

	  if ( ! iopt_l1b) {
	    PTB( h5d_write(grp4, "scat_beam_clat", &scat_clat[3*iblk], 2, 
			   start, count));
	    PTB( h5d_write(grp4, "scat_beam_clon", &scat_clon[3*iblk], 2, 
			   start, count));
	    PTB( h5d_write(grp4, "scat_polarization_roll", 
			   &scat_polarization_roll[3*iblk], 2, start, count));
	  }

	  start[2] = 0;
	  count[2] = 4;

	  if ( ! browse) {
	    PTB( h5d_write(grp4, "cellatfoot", &cellatfoot[4*3*iblk], 3, 
			   start, count));
	    PTB( h5d_write(grp4, "cellonfoot", &cellonfoot[4*3*iblk], 3, 
			   start, count));

	    PTB( h5d_write(grp4, "scat_latfoot", &scat_elat[4*3*iblk], 3, 
			   start, count));
	    PTB( h5d_write(grp4, "scat_lonfoot", &scat_elon[4*3*iblk], 3, 
			   start, count));


	    PTB( h5d_write(grp4, "sund", &sund[3*iblk], 2, start, count));
	    PTB( h5d_write(grp4, "sunr", &sunr[3*iblk], 2, start, count));
	    PTB( h5d_write(grp4, "moond", &moond[3*iblk], 2, start, count));

	    PTB( h5d_write(grp4, "zang", &zang[iblk], 1, start, count));

	    PTB( h5d_write(grp4, "sclon", &sclon[iblk], 1, start, count));
	    PTB( h5d_write(grp4, "sclat", &sclat[iblk], 1, start, count));
	    PTB( h5d_write(grp4, "scalt", &scalt[iblk], 1, start, count));
	  }
	  i++;
	}
      } 
    }

    //  PTB( SetScalarH5A (grp0, "Number of Orbit Vectors" , H5T_STD_I32LE, 
    //	       (VOIDP) &blknum));


    return 0;
  }



  /*----------------------------------------------------------------- */
  /* Read navigation records from L2 file                             */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::readl2_nav( float *beam_clat, float *beam_clon,
				 double *zang)
  {
    hsize_t start[6], count[6];

    start[0] = 0;
    count[0] = blknum;

    start[1] = 0;
    count[1] = 3;

    PTB( h5d_read(grp4, "beam_clat", beam_clat, 2, start, count));
    PTB( h5d_read(grp4, "beam_clon", beam_clon, 2, start, count));
    PTB( h5d_read(grp4, "zang", zang, 1, start, count));

    return 0;
  }


  /*----------------------------------------------------------------- */
  /* Read navigation records (with footprints) from L2 file           */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::readl2_nav( float *beam_clat, float *beam_clon,
				 float *cellatfoot, float *cellonfoot)
  {
    hsize_t start[6], count[6];

    start[0] = 0;
    count[0] = blknum;

    start[1] = 0;
    count[1] = 3;

    PTB( h5d_read(grp4, "beam_clat", beam_clat, 2, start, count));
    PTB( h5d_read(grp4, "beam_clon", beam_clon, 2, start, count));

    start[2] = 0;
    count[2] = 4;

    PTB( h5d_read(grp4, "cellatfoot", cellatfoot, 3, start, count));
    PTB( h5d_read(grp4, "cellonfoot", cellonfoot, 3, start, count));

    return 0;
  }


  /*----------------------------------------------------------------- */
  /* Open an existing HDF file                                        */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::openl2( char* l2_filename, unsigned flags,
			     int32_t *numBlks, 
			     double *granStart, double *granStop)
{
    h5fid =  H5Fopen(l2_filename, flags, H5P_DEFAULT);
    if ( h5fid < 0) {
      cout << l2_filename << " not found.";
      exit(1);
    }

    grp0 = H5Gopen1(h5fid,"/");
    grp1 = H5Gopen1(h5fid, "/Block Attributes");
    grp2 = H5Gopen1(h5fid, "/Aquarius Flags");
    grp3 = H5Gopen1(h5fid, "/Aquarius Data");;
    grp4 = H5Gopen1(h5fid, "/Navigation");
    //    grp5 = H5Gopen1(h5fid, "/Converted Telemetry");
    grp5 = -1;
    grp6 = -1;

    hid_t attr = H5Aopen_name(grp0, "Number of Blocks");
    H5Aread(attr, H5T_STD_I32LE, &blknum);
    H5Aclose(attr);

    /* Save old error handler */
    H5E_auto_t old_func;
    void *old_client_data;
    H5Eget_auto(H5E_DEFAULT , &old_func, &old_client_data);

    char timeString[17];
    hid_t atype = H5Tcopy(H5T_C_S1);
    H5Tset_size(atype, 17);

    H5Eset_auto( H5E_DEFAULT, NULL, NULL); // Turn off error handling
    attr = H5Aopen_name(grp0, "Start Time");
    H5Eset_auto( H5E_DEFAULT, old_func, old_client_data); // Restore previous error handler
    if (attr == -1) {
      granuleStart = -999;
    } else {
      H5Aread(attr, atype, timeString);
      granuleStart = get_tai( timeString);
      H5Aclose(attr);
    }

    H5Eset_auto( H5E_DEFAULT, NULL, NULL); // Turn off error handling
    attr = H5Aopen_name(grp0, "End Time");
    H5Eset_auto( H5E_DEFAULT, old_func, old_client_data); // Restore previous error handler
    if (attr == -1) {
      granuleStop = -999;
    } else {
      H5Aread(attr, atype, timeString);
      granuleStop = get_tai( timeString);
      H5Aclose(attr);
    }

    //    H5Eset_auto( H5E_DEFAULT, NULL, NULL); // Turn off error handling
    //attr = H5Aopen_name(grp0, "Orbit Start Time");
    //H5Eset_auto( H5E_DEFAULT, old_func, old_client_data); // Restore previous error handler
    //if (attr == -1) {
    //  orbitStart = -999;
    //} else {
    //  H5Aread(attr, atype, timeString);
    //  orbitStart = get_tai( timeString);
    //  H5Aclose(attr);
    //}

    //H5Eset_auto( H5E_DEFAULT, NULL, NULL); // Turn off error handling
    //attr = H5Aopen_name(grp0, "Orbit Stop Time");
    //H5Eset_auto( H5E_DEFAULT, old_func, old_client_data); // Restore previous error handler
    //if (attr == -1) {
    //  orbitStop = -999;
    //} else {
    //  H5Aread(attr, atype, timeString);
    //  orbitStop = get_tai( timeString);
    //  H5Aclose(attr);
    //}
    H5Tclose(atype);

    openFlags = flags;
    *numBlks = blknum;

    *granStart = granuleStart;
    *granStop  = granuleStop;

    /*
      Gran Start/Stop now corresponds with Orbit Start/Stop 05/13/09
    *orbStart = orbitStart;
    *orbStop  = orbitStop;
    */
    return 0;
  }


  /*----------------------------------------------------------------- */
  /* Write data/flag record to L2 file                                */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::writel2( int32_t iblk, hid_t grp,
			      char *prodName, VOIDP data)
  {
    hsize_t start[6], count[6];
   
    start[0] = iblk;
    count[0] = 1;
    start[1] = 0;
    count[1] = NUMBER_OF_BEAMS;
    PTB( h5d_write(grp, prodName, data, 2, start, count));

    return 0;
  }


  /*----------------------------------------------------------------- */
  /* Write rad flag record to L2 file                                 */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::writel2radflag( int32_t iblk, hid_t grp, VOIDP data)
  {
    hsize_t start[3], count[3];

    start[0] = iblk;
    count[0] = 1;

    start[1] = 0;
    count[1] = NUMBER_OF_BEAMS;

    start[2] = 0;
    count[2] = 4;

    PTB( h5d_write(grp,  "radiometer_flags", data, 3, start, count));

    return 0;
  }


  /*----------------------------------------------------------------- */
  /* Write rad flag record to L2 SM file                                 */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::writel2radflagsm( int32_t iblk, hid_t grp, VOIDP data)
  {
    hsize_t start[2], count[2];

    start[0] = iblk;
    count[0] = 1;

    start[1] = 0;
    count[1] = NUMBER_OF_BEAMS;

    PTB( h5d_write(grp,  "radiometer_flags", data, 2, start, count));

    return 0;
  }


  /*----------------------------------------------------------------- */
  /* Write rad RFI flag record to L2 file                             */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::writel2radrfi( int32_t iblk, hid_t grp, VOIDP data)
  {
    hsize_t start[4], count[4];

    start[0] = iblk;
    count[0] = 1;

    start[1] = 0;
    count[1] = NUMBER_OF_BEAMS;

    start[2] = 0;
    count[2] = RADIOMETER_POLARIZATIONS;

    start[3] = 0;
    count[3] = RADIOMETER_SUBCYCLES;

    PTB( h5d_write(grp, "rad_rfi_flags", data, 4, start, count));

    return 0;
  }


  /*----------------------------------------------------------------- */
  /* Write scat RFI flag record to L2 file                            */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::writel2scatrfi( int32_t iblk, hid_t grp, VOIDP data)
  {
    hsize_t start[3], count[3];
   
    start[0] = iblk;
    count[0] = 1;
    start[1] = 0;
    count[1] = NUMBER_OF_BEAMS;
    start[2] = 0;
    count[2] = 2;
    PTB( h5d_write(grp, "scat_rfi_flags", data, 3, start, count));

    return 0;
  }


  /*----------------------------------------------------------------- */
  /* Write rad/scat samples to L2 file                                */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::writel2samples( int32_t iblk, hid_t grp, 
				     char *prodName, VOIDP data)
  {
    hsize_t start[3], count[3];
   
    start[0] = iblk;
    count[0] = 1;
    start[1] = 0;
    count[1] = NUMBER_OF_BEAMS;

    if ( strcmp( prodName, "rad_samples") == 0) {
      start[2] = 0;
      count[2] = RADIOMETER_POLARIZATIONS;

      PTB( h5d_write(grp, prodName, data, 3, start, count));
    } else {
      PTB( h5d_write(grp, prodName, data, 2, start, count));
    }

    return 0;
  }


  /*----------------------------------------------------------------- */
  /* Write sec record to L2 file                                      */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::writel2sec( int32_t iblk, hid_t grp, VOIDP data)
  {
    hsize_t start[6], count[6];
   
    double blkSec;
    memcpy( &blkSec, data, sizeof(double));

    string ydhmsf_str( ydhmsf( gpstai2unix(  blkSec), 'G'));

    // Secs since start of day (mod 86400)
    double secStartOfDay = get_millisec( &ydhmsf_str) / 1000.;

    // Convert to midblock time
    secStartOfDay += (double) 0.72000e0;

    start[0] = iblk;
    count[0] = 1;
    PTB( h5d_write(grp, "sec", &secStartOfDay, 1, start, count));

    PTB( h5d_write(grp, "secGPS", &blkSec, 1, start, count));

    return 0;
  }


  /*----------------------------------------------------------------- */
  /* Write sec record to L2 file                                      */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::writeSolarXrayFlux( int32_t iblk, hid_t grp, VOIDP data)
  {
    hsize_t start[6], count[6];
   
    start[0] = iblk;
    count[0] = 1;
    PTB( h5d_write(grp, "solar xray flux", data, 1, start, count));

    return 0;
  }


  /*----------------------------------------------------------------- */
  /* Write caltemps to L2 file                                        */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::writel2_caltemps( int32_t iblk,  VOIDP data)
  {
    hsize_t start[6], count[6];
   
    start[0] = iblk;
    count[0] = 1;
    start[1] = 0;
    count[1] = NUMCALTEMPS;
    PTB( h5d_write(grp5, "rad_caltemps", data, 2, start, count));

    return 0;
  }


  /*----------------------------------------------------------------- */
  /* Write gainOff to L2 file                                         */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::writel2_gainoff( int32_t iblk, gainoff *gainOff,
				      bool iopt_l1b)
  {
    hsize_t start[6], count[6];
    start[0] = iblk;
    count[0] = 1;
    start[1] = 0;
    count[1] = NUMBER_OF_BEAMS;

    if ( ! iopt_l1b) {
      PTB( h5d_write(grp5, "rad_gvv", gainOff->gvv[iblk], 2, start, count));
      PTB( h5d_write(grp5, "rad_ghh", gainOff->ghh[iblk], 2, start, count));
      PTB( h5d_write(grp5, "rad_gpp", gainOff->gpp[iblk], 2, start, count));
      PTB( h5d_write(grp5, "rad_gmm", gainOff->gmm[iblk], 2, start, count));

      PTB( h5d_write(grp5, "rad_ov", gainOff->ov[iblk], 2, start, count));
      PTB( h5d_write(grp5, "rad_oh", gainOff->oh[iblk], 2, start, count));
      PTB( h5d_write(grp5, "rad_op", gainOff->op[iblk], 2, start, count));
      PTB( h5d_write(grp5, "rad_om", gainOff->om[iblk], 2, start, count));
      //PTB( h5d_write(grp5, "rad_gpv", gainOff->gpv[iblk], 2, start, count));
      //PTB( h5d_write(grp5, "rad_gph", gainOff->gph[iblk], 2, start, count));
      //PTB( h5d_write(grp5, "rad_gmv", gainOff->gmv[iblk], 2, start, count));
      //PTB( h5d_write(grp5, "rad_gmh", gainOff->gmh[iblk], 2, start, count));
      //PTB( h5d_write(grp5, "rad_gpU", gainOff->gpU[iblk], 2, start, count));
      //PTB( h5d_write(grp5, "rad_gmU", gainOff->gmU[iblk], 2, start, count));
      //PTB( h5d_write(grp5, "rad_hpU", gainOff->hpU[iblk], 2, start, count));
      //PTB( h5d_write(grp5, "rad_hmU", gainOff->hmU[iblk], 2, start, count));
    }

    return 0;
  }


  /*----------------------------------------------------------------- */
  /* Read sec record from L2 file                                     */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::readl2sec( int32_t iblk, hid_t grp, VOIDP data)
  {
    hsize_t start[6], count[6];
   
    start[0] = iblk;
    count[0] = 1;
    PTB( h5d_read(grp, "sec", data, 1, start, count));

    double blkSec;
    memcpy( &blkSec, data, sizeof(double));

    // Convert from time of day to TAI
    int32_t year, doy;
    hid_t attr;

    char timeString[17];
    hid_t atype = H5Tcopy(H5T_C_S1);
    H5Tset_size(atype, 17);
    attr = H5Aopen_name(grp0, "Start Time");
    H5Aread(attr, atype, timeString);
    H5Aclose(attr);
    H5Tclose(atype);

    istringstream istr;
    istr.str( string(timeString).substr(0,4));
    istr >> year;
    istr.clear();
    istr.str( string(timeString).substr(4,3));
    istr >> doy;

    /*
    attr = H5Aopen_name(grp0, "Start Year");
    H5Aread(attr, H5T_STD_I32LE, &year);
    H5Aclose(attr);

    attr = H5Aopen_name(grp0, "Start Day");
    H5Aread(attr, H5T_STD_I32LE, &doy);
    H5Aclose(attr);
    */

    // blkSec is TAI80
    if ( blkSec != -1) {
      blkSec = get_tai( year, doy, blkSec*1000);
      if ( blkSec < granuleStart) blkSec += 86400.;
    }

    memcpy( data, &blkSec, sizeof(double));

    return 0;
  }


  /*----------------------------------------------------------------- */
  /* Read record from L2 file                                         */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::readl2( int32_t iblk, hid_t grp, 
			     char *prodName, VOIDP data)
  {
    hsize_t start[6], count[6];
   
    start[0] = iblk;
    count[0] = 1;
    start[1] = 0;
    count[1] = NUMBER_OF_BEAMS;
    PTB( h5d_read(grp, prodName, data, 2, start, count));

    return 0;
  }


  /*----------------------------------------------------------------- */
  /* Read radiometer flags from L2 file                               */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::readl2radflag( int32_t iblk, hid_t grp, VOIDP data)
  {
    hsize_t start[6], count[6];
   
    start[0] = iblk;
    count[0] = 1;
    start[1] = 0;
    count[1] = NUMBER_OF_BEAMS;
    start[2] = 0;
    count[2] = 4;
    PTB( h5d_read(grp, "radiometer_flags", data, 3, start, count));

    return 0;
  }

  /*----------------------------------------------------------------- */
  /* Read radiometer flags from L2 SM file                            */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::readl2radflagsm( int32_t iblk, hid_t grp, VOIDP data)
  {
    hsize_t start[6], count[6];
   
    start[0] = iblk;
    count[0] = 1;
    start[1] = 0;
    count[1] = NUMBER_OF_BEAMS;
    PTB( h5d_read(grp, "radiometer_flags", data, 2, start, count));

    return 0;
  }


  /*----------------------------------------------------------------- */
  /* Close L2 file                                                    */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::closel2()
  {

    if (openFlags == H5F_ACC_RDWR) {
      PTB( SetScalarH5A (grp0, "Number of Blocks" , H5T_STD_I32LE, 
			 (VOIDP) &blknum));
    }

    if (grp0 != -1) H5Gclose(grp0);
    if (grp1 != -1) H5Gclose(grp1);
    if (grp2 != -1) H5Gclose(grp2);
    if (grp3 != -1) H5Gclose(grp3);
    if (grp4 != -1) H5Gclose(grp4);
    if (grp5 != -1) H5Gclose(grp5);
    if (grp6 != -1) H5Gclose(grp6);

    H5Fclose(h5fid);

    return 0;
  }


  int hdf5_Aquarius::closel2( double granuleStop, float percentWater,
			      float percentRFI, float *solar_attr,
			      string *rad_anc_files, string *scat_anc_files,
			      string *rad_calfiles, string *rad_tables,
			      char *anomaly_status, float *c_deltaTND, 
			      const char *radLimitsStr, const char *nominalNav,
			      float *offset_corr)
  {
    if (openFlags == H5F_ACC_RDWR) {
      PTB( SetScalarH5A (grp0, "Number of Blocks" , H5T_STD_I32LE, 
			 (VOIDP) &blknum));

      // Seconds from UTC 6 January 1980
      string ydhmsf_str( ydhmsf( gpstai2unix( granuleStop), 'G'));
      PTB( SetScalarH5A (grp0, "End Time", H5T_STRING, 
			 (void *) ydhmsf_str.c_str()));

      istringstream istr;
      int32_t itemp;
      ydhmsf_str = ydhmsf( gpstai2unix( granuleStop), 'G');
      int32_t millisec = get_millisec( &ydhmsf_str);

      istr.clear(); istr.str( ydhmsf_str.substr(0,4)); istr >> itemp;
      PTB( SetScalarH5A (grp0, "End Year", H5T_STD_I32LE, 
			 (VOIDP) &itemp));
      istr.clear(); istr.str( ydhmsf_str.substr(4,3)); istr >> itemp;
      PTB( SetScalarH5A (grp0, "End Day", H5T_STD_I32LE, (VOIDP) &itemp));
      PTB( SetScalarH5A (grp0, "End Millisec", H5T_STD_I32LE, 
			 (VOIDP) &millisec));


      PTB( SetScalarH5A (grp0, "Percent Water", H5T_NATIVE_FLOAT, 
			 (VOIDP) &percentWater));

      PTB( SetScalarH5A (grp0, "Percent RFI", H5T_NATIVE_FLOAT, 
			 (VOIDP) &percentRFI));

      PTB( SetScalarH5A (grp0, "Mean Solar 1415 MHz Flux", H5T_NATIVE_FLOAT, 
			 (VOIDP) &solar_attr[1]));

      // PTB( SetScalarH5A (grp0, "Peak Solar 1415 MHz Flux", H5T_NATIVE_FLOAT, 
      //		 (VOIDP) &solar_attr[0]));

      // PTB( SetScalarH5A (grp0, "Mean Solar Xray Flux", H5T_NATIVE_FLOAT, 
      //		 (VOIDP) &solar_attr[3]));

      // PTB( SetScalarH5A (grp0, "Peak Solar Xray Flux", H5T_NATIVE_FLOAT, 
      //		 (VOIDP) &solar_attr[2]));

      //      PTB( SetScalarH5A (grp0, "Cold Sky Calibration", H5T_STRING,
      //		 (void *) coldSky)); 


      // Generate radiometer ancillary file attributes
      string str;
      char rad_anc_name[8][80];
      for (int32_t i=0; i<3; i++) {
	if ( rad_anc_files[i].compare("") != 0) {

	  // Get ancillary filenames from ancillary "y" file attributes
	  hid_t h5fid =  
	    H5Fopen(rad_anc_files[i].c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	  hid_t grp = H5Gopen1(h5fid, "/");
	  hid_t attr = H5Aopen_name(grp, "Input Ancillary Files");
	  hid_t type = H5Tcopy(H5T_C_S1);
	  H5Tset_size(type, 80);
	  H5Aread(attr, type, rad_anc_name);
          H5Aclose(attr);
	  H5Gclose(grp);
	  H5Fclose(h5fid);
          H5Tclose(type);

	  // Generate attribute
	  str.assign(basename((char*) rad_anc_files[i].c_str()));
	  str.append(":");
	  for (int32_t j=0; j<8; j++) {
	    str.append(rad_anc_name[j]);
	    if ( j < 7) str.append(",");
	  }

	  // Write attribute
	  stringstream ss;
	  ss << (i + 1);
	  string rad_anc_fldname = "RAD_ANCILLARY_FILE" + ss.str();
	  PTB( SetScalarH5A (grp0, rad_anc_fldname.c_str(), H5T_STRING, 
			     (void *) str.c_str()));
	}
      }

      // Write scatterometer ancillary files attribute
      PTB( SetScalarH5A (grp0, "Scatterometer Ancillary Files", 
			 H5T_STRING, (void *) scat_anc_files->c_str()));

      // Write radiometer cal files attribute
      PTB( SetScalarH5A (grp0, "Radiometer Calibration Files", 
			 H5T_STRING, (void *) rad_calfiles->c_str()));

      // Write radiometer data tables attribute
      PTB( SetScalarH5A (grp0, "Radiometer Data Tables", 
			 H5T_STRING, (void *) rad_tables->c_str()));

      // Write scatterometer prt/dict coefficient files
      str.assign(basename(getenv("AQUARIUS_ATC_PRT")));
      str.append(",");
      str.append(basename(getenv("AQUARIUS_EXT_PRT")));
      str.append(",");
      str.append(basename(getenv("AQUARIUS_SCAT_PRT")));
      str.append(",");
      str.append(basename(getenv("AQUARIUS_RAD_PRT")));
      str.append(",");
      str.append(basename(getenv("AQUARIUS_CMD_DICT")));
      PTB( SetScalarH5A (grp0, "Scatterometer Coefficient Files", 
			 H5T_STRING, (void *) str.c_str()));
    }

    if ( c_deltaTND[0] != 0.0) {
      float delta[3];
      delta[0] = c_deltaTND[0];
      delta[1] = c_deltaTND[2];
      delta[2] = c_deltaTND[4];
      PTB(h5a_set(grp0, "Delta TND V coefficient", H5T_NATIVE_FLOAT, 3, 
		  (VOIDP) delta));

      delta[0] = c_deltaTND[1];
      delta[1] = c_deltaTND[3];
      delta[2] = c_deltaTND[5];
      PTB(h5a_set(grp0, "Delta TND H coefficient", H5T_NATIVE_FLOAT, 3, 
		  (VOIDP) delta));
    } // if ( c_deltaTND[0] != 0.0)

    PTB( SetScalarH5A (grp0, "Radiometer Flag Limits", 
		       H5T_STRING, (void *) radLimitsStr));

    PTB( SetScalarH5A (grp0, "Nominal Navigation", 
		       H5T_STRING, (void *) nominalNav));

    PTB(h5a_set(grp0, "Radiometer Offset Correction", H5T_NATIVE_FLOAT, 6, 
		(VOIDP) offset_corr));

    if ( strlen(anomaly_status) != 0) {
      PTB( SetScalarH5A (grp0, "Anomaly Status String", 
			 H5T_STRING, (void *) anomaly_status));
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



  /*----------------------------------------------------------------- */
  /* Create an HDF5 level1 file                                       */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::createl2sm( char* l2_filename, int32_t numBlocks,
				 char* starttime, char* endtime,   
				 double nodeCrossingTime, float nodeLongitude,
				 int32_t orbitNumber, int32_t cycleNumber, 
				 int32_t passNumber,
				 char* inputFiles, char* processControl,
				 char *softver, char *pversion, 
                                 char *ancillary_files,
				 char *nominalNav, char *anomaly_status,
                                 char *input_L2_file)
  {
    // Set object variable blknum
    blknum = numBlocks;

    /* Create the HDF file */
    h5fid = H5Fcreate(l2_filename, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    if(h5fid == FAIL) {
      fprintf(stderr,
	      "-E- %s line %d: Could not create HDF file, %s .\n",
	      __FILE__,__LINE__,l2_filename);
      return(HDF_FUNCTION_ERROR);
    }   

    openFlags = H5F_ACC_RDWR;
  
    grp0 = H5Gopen1(h5fid,"/");
    grp1 = H5Gcreate1(h5fid, "/Block Attributes", 0);
    grp2 = H5Gcreate1(h5fid, "/Aquarius Flags", 0);
    grp3 = H5Gcreate1(h5fid, "/Aquarius Data", 0);
    grp4 = H5Gcreate1(h5fid, "/Navigation", 0);
    grp5 = -1;
    grp6 = -1;

    // Setup dataset create property list
    hid_t proplist = H5Pcreate( H5P_DATASET_CREATE);

    // Set fill to -1
    double time_fill_value = -1.0;
    H5Pset_fill_value( proplist, H5T_NATIVE_DOUBLE, (void *) &time_fill_value);

    PTB( CreateH5D(
     grp1,                                            /* file id         */
     "sec",                                           /* short name      */
     "Block time, seconds of day",                    /* long name       */
     "",                                              /* standard name   */
     "seconds",                                       /* units           */
     0, 0,                                            /* valid range     */
     0, 0,                                            /* slope           */
     -9999,                                           /* fill value      */
     H5T_NATIVE_DOUBLE,                               /* HDF number type */
     1,                                               /* rank            */
     numBlocks, 1, 1, 1, 1, 1,                        /* dimension sizes */
     "Number of Blocks", NULL, NULL, NULL, NULL, NULL,/* dimension names */
     proplist) );

    PTB( CreateH5D(
     grp1,                                            /* file id         */
     "secGPS",                                        /* short name      */
     "Block time, GPS time",                          /* long name       */
     "",                                              /* standard name   */
     "seconds",                                       /* units           */
     0, 0,                                            /* valid range     */
     0, 0,                                            /* slope           */
     -9999,                                           /* fill value      */
     H5T_NATIVE_DOUBLE,                               /* HDF number type */
     1,                                               /* rank            */
     numBlocks, 1, 1, 1, 1, 1,                        /* dimension sizes */
     "Number of Blocks", NULL, NULL, NULL, NULL, NULL,/* dimension names */
     proplist) );


    PTB( CreateH5D(
     grp3,                                            /* file id         */
     "rad_sm",                                        /* short name      */
     " Volumetric Soil Moisture",                     /* long name       */
     "",                                              /* standard name   */
     "m3/m3",                                         /* units           */
     0, 1,                                            /* valid range     */
     0, 0,                                            /* slope           */
     -9999,                                           /* fill value      */
     H5T_NATIVE_FLOAT,                                /* HDF number type */
     2,                                               /* rank            */
     numBlocks,                                       /* dimension sizes */
     NUMBER_OF_BEAMS, 1, 1, 1, 1, 
     "Number of Blocks", "Number of Beams",           /* dimension names */
     NULL, NULL, NULL, NULL,
     H5P_DEFAULT) );


    PTB( CreateH5D(
     grp3,                                            /* file id         */
     "rad_land_frac",                                 /* short name      */
     "Fraction of land contamination (radiometer)",   /* long name       */
     "",                                              /* standard name   */
     "",                                              /* units           */
     0, 0,                                            /* valid range     */
     0, 0,                                            /* slope           */
     -9999,                                           /* fill value      */
     H5T_NATIVE_FLOAT,                                /* HDF number type */
     2,                                               /* rank            */
     numBlocks,                                       /* dimension sizes */
     NUMBER_OF_BEAMS, 1, 1, 1, 1, 
     "Number of Blocks", "Number of Beams",           /* dimension names */
     NULL, NULL, NULL, NULL,
     H5P_DEFAULT) );

    PTB( CreateH5D(
     grp3,                                            /* file id         */
     "rad_ice_frac",                                  /* short name      */
     "Fraction of ice contamination (radiometer)",    /* long name       */
     "",                                              /* standard name   */
     "",                                              /* units           */
     0, 0,                                            /* valid range     */
     0, 0,                                            /* slope           */
     -9999,                                           /* fill value      */
     H5T_NATIVE_FLOAT,                                /* HDF number type */
     2,                                               /* rank            */
     numBlocks,                                       /* dimension sizes */
     NUMBER_OF_BEAMS, 1, 1, 1, 1, 
     "Number of Blocks", "Number of Beams",           /* dimension names */
     NULL, NULL, NULL, NULL,
     H5P_DEFAULT) );

    PTB( CreateH5D(
     grp3,                                            /* file id         */
     "rad_TbV",                                       /* short name      */
     "Earth surface Tb V polarization",               /* long name       */
     "",                                              /* standard name   */
     "Kelvin",                                        /* units           */
     0, 0,                                            /* valid range     */
     0, 0,                                            /* slope           */
     -9999,                                           /* fill value      */
     H5T_NATIVE_FLOAT,                                /* HDF number type */
     2,                                               /* rank            */
     numBlocks,                                       /* dimension sizes */
     NUMBER_OF_BEAMS, 1, 1, 1, 1, 
     "Number of Blocks", "Number of Beams",           /* dimension names */
     NULL, NULL, NULL, NULL,
     H5P_DEFAULT) );

    PTB( CreateH5D(
     grp3,                                            /* file id         */
     "rad_TbH",                                       /* short name      */
     "Earth surface Tb H polarization",               /* long name       */
     "",                                              /* standard name   */
     "Kelvin",                                        /* units           */
     0, 0,                                            /* valid range     */
     0, 0,                                            /* slope           */
     -9999,                                           /* fill value      */
     H5T_NATIVE_FLOAT,                                /* HDF number type */
     2,                                               /* rank            */
     numBlocks,                                       /* dimension sizes */
     NUMBER_OF_BEAMS, 1, 1, 1, 1, 
     "Number of Blocks", "Number of Beams",           /* dimension names */
     NULL, NULL, NULL, NULL,
     H5P_DEFAULT) );

    PTB( CreateH5D(
     grp3,                                            /* file id         */
     "anc_surface_temp",                              /* short name      */
     "Surface Temperature",                           /* long name       */
     "",                                              /* standard name   */
     "Kelvin",                                        /* units           */
     0, 0,                                            /* valid range     */
     0, 0,                                            /* slope           */
     -999,                                            /* fill value      */
     H5T_NATIVE_FLOAT,                                /* HDF number type */
     2,                                               /* rank            */
     numBlocks,                                       /* dimension sizes */
     NUMBER_OF_BEAMS, 1, 1, 1, 1, 
     "Number of Blocks", "Number of Beams",           /* dimension names */
     NULL, NULL, NULL, NULL,
     H5P_DEFAULT) );

    PTB( CreateH5D(
     grp3,                                            /* file id         */
     "anc_subsurf_temp",                              /* short name      */
     "Sub-Surface Temperature",                       /* long name       */
     "",                                              /* standard name   */
     "Kelvin",                                        /* units           */
     0, 0,                                            /* valid range     */
     0, 0,                                            /* slope           */
     -999,                                            /* fill value      */
     H5T_NATIVE_FLOAT,                                /* HDF number type */
     2,                                               /* rank            */
     numBlocks,                                       /* dimension sizes */
     NUMBER_OF_BEAMS, 1, 1, 1, 1, 
     "Number of Blocks", "Number of Beams",           /* dimension names */
     NULL, NULL, NULL, NULL,
     H5P_DEFAULT) );

    PTB( CreateH5D(
     grp3,                                            /* file id         */
     "anc_swe",                                       /* short name      */
     "Snow Water Equivalent",                         /* long name       */
     "",                                              /* standard name   */
     "kg/m^2",                                        /* units           */
     0, 0,                                            /* valid range     */
     0, 0,                                            /* slope           */
     -999,                                            /* fill value      */
     H5T_NATIVE_FLOAT,                                /* HDF number type */
     2,                                               /* rank            */
     numBlocks,                                       /* dimension sizes */
     NUMBER_OF_BEAMS, 1, 1, 1, 1, 
     "Number of Blocks", "Number of Beams",           /* dimension names */
     NULL, NULL, NULL, NULL,
     H5P_DEFAULT) );

    PTB( CreateH5D(
     grp3,                                            /* file id         */
     "anc_sm",                                        /* short name      */
     "Soil Moisture",                                 /* long name       */
     "",                                              /* standard name   */
     "",                                              /* units           */
     0, 1,                                            /* valid range     */
     0, 0,                                            /* slope           */
     -999,                                            /* fill value      */
     H5T_NATIVE_FLOAT,                                /* HDF number type */
     2,                                               /* rank            */
     numBlocks,                                       /* dimension sizes */
     NUMBER_OF_BEAMS, 1, 1, 1, 1, 
     "Number of Blocks", "Number of Beams",           /* dimension names */
     NULL, NULL, NULL, NULL,
     H5P_DEFAULT) );

    PTB( CreateH5D(
     grp3,                                            /* file id         */
     "scat_VV_toa",                                   /* short name      */
     "TOA Scatterometer NRCS for VV polarization",    /* long name       */
     "",                                              /* standard name   */
     "db",                                            /* units           */
     0, 0,                                            /* valid range     */
     0, 0,                                            /* slope           */
     -999,                                            /* fill value      */
     H5T_NATIVE_FLOAT,                                /* HDF number type */
     2,                                               /* rank            */
     numBlocks,                                       /* dimension sizes */
     NUMBER_OF_BEAMS, 1, 1, 1, 1, 
     "Number of Blocks", "Number of Beams",           /* dimension names */
     NULL, NULL, NULL, NULL,
     H5P_DEFAULT) );

    PTB( CreateH5D(
      grp3,                                            /* file id         */
      "scat_HH_toa",                                   /* short name      */
      "TOA Scatterometer NRCS for HH polarization",    /* long name       */
      "",                                              /* standard name   */
      "db",                                            /* units           */
      0, 0,                                            /* valid range     */
      0, 0,                                            /* slope           */
      -999,                                            /* fill value      */
      H5T_NATIVE_FLOAT,                                /* HDF number type */
      2,                                               /* rank            */
      numBlocks,                                       /* dimension sizes */
      NUMBER_OF_BEAMS, 1, 1, 1, 1, 
      "Number of Blocks", "Number of Beams",           /* dimension names */
      NULL, NULL, NULL, NULL,
      H5P_DEFAULT) );

    PTB( CreateH5D(
     grp3,                                            /* file id         */
     "scat_HV_toa",                                   /* short name      */
     "TOA Scatterometer NRCS for HV polarization",    /* long name       */
     "",                                              /* standard name   */
     "db",                                            /* units           */
     0, 0,                                            /* valid range     */
     0, 0,                                            /* slope           */
     -999,                                            /* fill value      */
     H5T_NATIVE_FLOAT,                                /* HDF number type */
     2,                                               /* rank            */
     numBlocks,                                       /* dimension sizes */
     NUMBER_OF_BEAMS, 1, 1, 1, 1, 
     "Number of Blocks", "Number of Beams",           /* dimension names */
     NULL, NULL, NULL, NULL,
     H5P_DEFAULT) );

    PTB( CreateH5D(
     grp3,                                            /* file id         */
     "scat_VH_toa",                                   /* short name      */
     "TOA Scatterometer NRCS for VH polarization",    /* long name       */
     "",                                              /* standard name   */
     "db",                                            /* units           */
     0, 0,                                            /* valid range     */
     0, 0,                                            /* slope           */
     -999,                                            /* fill value      */
     H5T_NATIVE_FLOAT,                                /* HDF number type */
     2,                                               /* rank            */
     numBlocks,                                       /* dimension sizes */
     NUMBER_OF_BEAMS, 1, 1, 1, 1, 
     "Number of Blocks", "Number of Beams",           /* dimension names */
     NULL, NULL, NULL, NULL,
     H5P_DEFAULT) );


    PTB( CreateH5D(
     grp2,                                            /* file id         */
     "radiometer_flags",                              /* short name      */
     "Radiometer data quality flags",                 /* long name       */
     "",                                              /* units           */
     0, 0,                                            /* valid range     */
     0, 0,                                            /* slope           */
     H5T_STD_U32LE,                                   /* HDF number type */
     2,                                               /* rank            */
     numBlocks,                                       /* dimension sizes */
     NUMBER_OF_BEAMS, 1, 1, 1, 1, 
     "Number of Blocks", "Number of Beams",           /* dimension names */
     NULL, NULL, NULL, NULL,
     H5P_DEFAULT) );


    H5Pclose( proplist);


    // Write flag attributes
    hid_t dataset;
    dataset = H5Dopen1(grp2, "radiometer_flags");
    PTB( SetScalarH5A (dataset, "No SM retrieval", (hid_t) H5T_STRING, (VOIDP) "SMRET"));
    PTB( SetScalarH5A (dataset, "Brighness Temp", (hid_t) H5T_STRING, (VOIDP) "TB"));
    PTB( SetScalarH5A (dataset, "Orbit Manuevers", (hid_t) H5T_STRING, (VOIDP) "ORBIT"));
    PTB( SetScalarH5A (dataset, "RFI", (hid_t) H5T_STRING, (VOIDP) "RFI"));
    PTB( SetScalarH5A (dataset, "Surface Temp", (hid_t) H5T_STRING, (VOIDP) "TSURF"));
    PTB( SetScalarH5A (dataset, "Frozen Ground", (hid_t) H5T_STRING, (VOIDP) "FROZ"));
    PTB( SetScalarH5A (dataset, "Snow", (hid_t) H5T_STRING, (VOIDP) "SNOW"));
    PTB( SetScalarH5A (dataset, "Ice", (hid_t) H5T_STRING, (VOIDP) "ICE"));
    PTB( SetScalarH5A (dataset, "NDVI", (hid_t) H5T_STRING, (VOIDP) "NDVI"));
    PTB( SetScalarH5A (dataset, "Vegetation", (hid_t) H5T_STRING, (VOIDP) "VEG"));
    PTB( SetScalarH5A (dataset, "Urban", (hid_t) H5T_STRING, (VOIDP) "URBAN"));
    PTB( SetScalarH5A (dataset, "Soil", (hid_t) H5T_STRING, (VOIDP) "SOIL"));
    PTB( SetScalarH5A (dataset, "Water", (hid_t) H5T_STRING, (VOIDP) "WATER"));
    
    /*                                                                  */
    /* Write out some global attributes                                 */
    /* ---------------------------------------------------------------- */
    /*                                                                  */
    PTB( SetScalarH5A (grp0, "Product Name", H5T_STRING, 
		       basename(l2_filename)));
    PTB( SetScalarH5A (grp0, "Title", H5T_STRING, 
		       (VOIDP) "Aquarius Level-2 Soil Moisture Data"));
    PTB( SetScalarH5A (grp0, "Data Center", H5T_STRING, (VOIDP) DATACENTER));
    PTB( SetScalarH5A (grp0, "Mission", H5T_STRING, (VOIDP) "SAC-D Aquarius"));
    PTB( SetScalarH5A (grp0, "Mission Characteristics", H5T_STRING, 
		       (VOIDP) MISSIONCHARACTERISTICS));
    PTB( SetScalarH5A (grp0, "Sensor" , H5T_STRING, (VOIDP) "Aquarius"));
    PTB( SetScalarH5A (grp0, "Sensor Characteristics" , H5T_STRING, 
		       (VOIDP) SENSORCHARACTERISTICS));
    PTB( SetScalarH5A (grp0, "Data Type", H5T_STRING, (VOIDP) "SM"));
    PTB( SetScalarH5A (grp0, "Software ID", H5T_STRING, (VOIDP) softver));
    PTB( SetScalarH5A (grp0, "Processing Version", H5T_STRING, 
		       (VOIDP) pversion));
    PTB( SetScalarH5A (grp0, "Processing Time", H5T_STRING, 
		       ydhmsf(time(NULL),'G')));
    PTB( SetScalarH5A (grp0, "Input Files", H5T_STRING, (VOIDP) inputFiles));
    PTB( SetScalarH5A (grp0, "Input L2 File", H5T_STRING, 
                       (VOIDP) input_L2_file));
    PTB( SetScalarH5A (grp0, "Processing Control", H5T_STRING, 
		       (VOIDP) processControl));

    PTB( SetScalarH5A (grp0, "_lastModified", H5T_STRING, 
		       ydhmsf(time(NULL),'G')));
    PTB( SetScalarH5A (grp0, "Conventions", H5T_STRING, (VOIDP) "CF-1.6")); 
    PTB( SetScalarH5A (grp0, "institution", H5T_STRING, 
                       (VOIDP) "NASA/GSFC OBPG")); 

    int32_t itemp;

    itemp = NUMBER_OF_BEAMS;
    PTB( SetScalarH5A (grp0, "Number of Beams" , H5T_STD_I32LE, 
		       (VOIDP) &itemp));

    PTB( SetScalarH5A (grp0, "Start Time", H5T_STRING, 
		       (void *) starttime));

    istringstream istr;
    string ydhmsf_str( starttime);

    istr.str( ydhmsf_str.substr(0,4)); istr >> itemp;
    PTB( SetScalarH5A (grp0, "Start Year", H5T_STD_I32LE, (VOIDP) &itemp));

    istr.clear(); istr.str( ydhmsf_str.substr(4,3)); istr >> itemp;
    PTB( SetScalarH5A (grp0, "Start Day", H5T_STD_I32LE, (VOIDP) &itemp));

    int32_t millisec = get_millisec( &ydhmsf_str);
    PTB( SetScalarH5A (grp0, "Start Millisec", H5T_STD_I32LE, 
		       (VOIDP) &millisec));

    PTB( SetScalarH5A (grp0, "End Time", H5T_STRING, 
		       (void *) endtime));
    ydhmsf_str.assign( endtime);

    istr.clear(); istr.str( ydhmsf_str.substr(0,4)); istr >> itemp;
    PTB( SetScalarH5A (grp0, "End Year", H5T_STD_I32LE, (VOIDP) &itemp));
    istr.clear(); istr.str( ydhmsf_str.substr(4,3)); istr >> itemp;
    PTB( SetScalarH5A (grp0, "End Day", H5T_STD_I32LE, (VOIDP) &itemp));
    PTB( SetScalarH5A (grp0, "End Millisec", H5T_STD_I32LE, 
		       (VOIDP) &millisec));

    ydhmsf_str = ydhmsf( gpstai2unix(  nodeCrossingTime), 'G');
    PTB( SetScalarH5A (grp0, "Node Crossing Time", H5T_STRING,
    		       (void *) ydhmsf_str.c_str())); 

    PTB( SetScalarH5A (grp0, "Orbit Node Longitude", H5T_NATIVE_FLOAT, 
    		       (VOIDP) &nodeLongitude));


    PTB( SetScalarH5A (grp0, "Latitude Units", H5T_STRING,
		       (void *) "degrees North")); 

    PTB( SetScalarH5A (grp0, "Longitude Units", H5T_STRING,
		       (void *) "degrees East")); 

    PTB( SetScalarH5A (grp0, "Orbit Number", H5T_STD_I32LE,
    		       (VOIDP) &orbitNumber));

    PTB( SetScalarH5A (grp0, "Cycle Number", H5T_STD_I32LE,
    		       (VOIDP) &cycleNumber));
    
    PTB( SetScalarH5A (grp0, "Pass Number", H5T_STD_I32LE,
    		       (VOIDP) &passNumber));

    PTB( SetScalarH5A (grp0, "Nominal Navigation", 
		       H5T_STRING, (void *) nominalNav));

    if ( strcmp(anomaly_status, "N/A") != 0) {
      PTB( SetScalarH5A (grp0, "Anomaly Status String", 
			 H5T_STRING, (void *) anomaly_status));
    }

    PTB( SetScalarH5A (grp0, "Ancillary Files",
 		       H5T_STRING, (void *) ancillary_files));

    return 0;
  }


  /*----------------------------------------------------------------- */
  /* Write navigation records to L2 SM file                           */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::writel2_navsm( int32_t nBlks, float *clat, float *clon,
				    double *zang, double *ang)
  {
    hsize_t start[6], count[6];

    // Attitude Angles (roll, pitch, yaw)
    PTB( CreateH5D(
    grp4,                                            /* file id         */
    "att_ang",                                       /* short name      */
    "Spacecraft roll, pitch, yaw",                   /* long name       */
    "degrees",                                       /* standard name   */
    "",                                              /* units           */
    -180., 180.,                                     /* valid range     */
    0, 0,                                            /* slope           */
    -9999,                                           /* fill value      */
    H5T_NATIVE_DOUBLE,                               /* HDF number type */
    2,                                               /* rank            */
    nBlks, NUMBER_OF_BEAMS, 1, 1, 1, 1 ,             /* dimension sizes */
    "Number of Blocks", "Number of Beams", 
    NULL, NULL, NULL, NULL,                          /* dimension names */
    H5P_DEFAULT) );

    // Beam center latitude
    PTB( CreateH5D(
    grp4,                                            /* file id         */
    "beam_clat",                                     /* short name      */
    "Beam Center Latitude",                          /* long name       */
    "",                                              /* standard name   */
    "degrees",                                       /* units           */
    -90., 90.,                                       /* valid range     */
    0, 0,                                            /* slope           */
    -9999,                                           /* fill value      */
    H5T_NATIVE_FLOAT,                                /* HDF number type */
    2,                                               /* rank            */
    nBlks, NUMBER_OF_BEAMS, 1, 1, 1, 1 ,             /* dimension sizes */
    "Number of Blocks", "Number of Beams", 
    NULL, NULL, NULL, NULL,                          /* dimension names */
    H5P_DEFAULT) );


    // Beam center longitude
    PTB( CreateH5D(
    grp4,                                            /* file id         */
    "beam_clon",                                     /* short name      */
    "Beam Center Longitude",                         /* long name       */
    "",                                              /* standard name   */
    "degrees",                                       /* units           */
    -180., 180.,                                     /* valid range     */
    0, 0,                                            /* slope           */
    -9999,                                           /* fill value      */
    H5T_NATIVE_FLOAT,                                /* HDF number type */
    2,                                               /* rank            */
    nBlks, NUMBER_OF_BEAMS, 1, 1, 1, 1 ,             /* dimension sizes */
    "Number of Blocks", "Number of Beams", 
    NULL, NULL, NULL, NULL,                          /* dimension names */
    H5P_DEFAULT) );

    // Intra-orbit angle
    PTB( CreateH5D(
    grp4,                                            /* file id         */
    "zang",                                          /* short name      */
    "Intra-Orbit Angle",                             /* long name       */
    "",                                              /* standard name   */
    "degrees",                                       /* units           */
    -180., 180.,                                     /* valid range     */
    0, 0,                                            /* slope           */
    -9999,                                           /* fill value      */
    H5T_NATIVE_DOUBLE,                               /* HDF number type */
    1,                                               /* rank            */
    nBlks, 1, 1, 1, 1, 1 ,                           /* dimension sizes */
    "Number of Blocks",
    NULL, NULL, NULL, NULL, NULL,                    /* dimension names */
    H5P_DEFAULT) );


    for (int32_t iblk=0; iblk<nBlks; iblk++) {
      start[0] = iblk;
      count[0] = 1;

      start[1] = 0;
      count[1] = 3;

      PTB( h5d_write(grp4, "beam_clat", &clat[3*iblk], 2, 
		     start, count));
      PTB( h5d_write(grp4, "beam_clon", &clon[3*iblk], 2, 
		     start, count));
      PTB( h5d_write(grp4, "zang", &zang[iblk], 1, start, count));

      PTB( h5d_write(grp4, "att_ang", &ang[3*iblk], 2, start, count)); 
    }

    return 0;
  }


  /*----------------------------------------------------------------- */
  /* Write navigation records to L2 SM file                           */
  /* ---------------------------------------------------------------- */
  int hdf5_Aquarius::writel2_s_acc( int32_t nBlks, 
                                    const char *name, const char *description,
                                    float *s_acc)
  {
    hsize_t start[6]={0,0,0,0,0,0};
    hsize_t count[6];

    // real(4), dimension(n_Sacc, n_subcyc, npol, n_rad, max_cyc) :: S_acc_raw
    //#define RADIOMETER_SUBCYCLES 12
    //#define RADIOMETER_SIGNALS_PER_SUBCYCLE 5
    //#define NUMBER_OF_BEAMS 3
    //#define RADIOMETER_POLARIZATIONS 4

    PTB( CreateH5D(
    grp6,                                            /* file id         */
    name,                                            /* short name      */
    description,                                     /* long name       */
    "",                                              /* standard name   */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    -9999,                                           /* fill value      */
    H5T_NATIVE_FLOAT,                                /* HDF number type */
    5,                                               /* rank            */
    nBlks, NUMBER_OF_BEAMS, 
    RADIOMETER_POLARIZATIONS, RADIOMETER_SUBCYCLES,
    RADIOMETER_SIGNALS_PER_SUBCYCLE+1, 1,            /* dimension sizes */
    "Number of Blocks", "Number of Beams", 
    "Number of Radiometer Polarizations", 
    "Number of Radiometer Subcycles", 
    "Radiometer Signals per Subcycle+1" , NULL,      /* dimension names */
    H5P_DEFAULT) );

    count[0] = nBlks;
    count[1] = NUMBER_OF_BEAMS;
    count[2] = RADIOMETER_POLARIZATIONS;
    count[3] = RADIOMETER_SUBCYCLES;
    count[4] = RADIOMETER_SIGNALS_PER_SUBCYCLE+1;

    PTB( h5d_write(grp6, name, s_acc, 5, start, count));

    return 0;
  }


  /*-------------------------------------------- */
  /* Write long accumulation calibration records */
  /* ------------------------------------------- */
  int hdf5_Aquarius::writel2_l_acc( int32_t nBlks,
                                    const char *name, const char *description,
                                    float *l_acc)
  {
    hsize_t start[6]={0,0,0,0,0,0};
    hsize_t count[6];

    PTB( CreateH5D(
    grp6,                                            /* file id         */
    name,                                            /* short name      */
    description,                                     /* long name       */
    "",                                              /* standard name   */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    -9999,                                           /* fill value      */
    H5T_NATIVE_FLOAT,                                /* HDF number type */
    4,                                               /* rank            */
    nBlks, NUMBER_OF_BEAMS, 
    RADIOMETER_POLARIZATIONS, 
    RADIOMETER_LONG_ACCUM, 1, 1,                     /* dimension sizes */
    "Number of Blocks", "Number of Beams", 
    "Number of Radiometer Polarizations", 
    "Radiometer Long Accumulations", 
    NULL, NULL,                                      /* dimension names */
    H5P_DEFAULT) );

    count[0] = nBlks;
    count[1] = NUMBER_OF_BEAMS;
    count[2] = RADIOMETER_POLARIZATIONS;
    count[3] = RADIOMETER_LONG_ACCUM;

    PTB( h5d_write(grp6, name, l_acc, 4, start, count));

    return 0;
  }


  int hdf5_Aquarius::writel2_sec_TA( int32_t nBlks, 
                                     Hdf::hdf5_Aquarius *l1afile, 
                                     float *TA)
  {
    PTB( CreateH5D(
    grp6,                                            /* file id         */
    "sec",                                           /* short name      */
    "Block time, seconds of day",                    /* long name       */
    "",                                              /* standard name   */
    "seconds",                                       /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    -9999,                                           /* fill value      */
    H5T_NATIVE_DOUBLE,                               /* HDF number type */
    1,                                               /* rank            */
    nBlks, 1, 1, 1, 1, 1,                            /* dimension sizes */
    "Number of Blocks", NULL, NULL, NULL, NULL, NULL, /* dimension names */
    H5P_DEFAULT) );

    PTB( CreateH5D(
    grp6,                                            /* file id         */
    "secGPS",                                        /* short name      */
    "Block time, GPS time",                          /* long name       */
    "",                                              /* standard name   */
    "seconds",                                       /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    -9999,                                           /* fill value      */
    H5T_NATIVE_DOUBLE,                               /* HDF number type */
    1,                                               /* rank            */
    nBlks, 1, 1, 1, 1, 1,                            /* dimension sizes */
    "Number of Blocks", NULL, NULL, NULL, NULL, NULL, /* dimension names */
    H5P_DEFAULT) );


    double blkSec;
    for (int32_t iblk=0; iblk<nBlks; iblk++) {
      l1afile->readl1_radiometer(iblk, 0, 0, &blkSec, NULL,
                                 NULL, NULL, NULL);
     this->writel2sec(iblk, grp6, (VOIDP) &blkSec);
    }


    PTB( CreateH5D(
    grp6,                                            /* file id         */
    "TA",                                            /* short name      */
    "Radiometer Antenna Temperatures",               /* long name       */
    "",                                              /* standard name   */
    "Kelvin",                                        /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    -9999,                                           /* fill value      */
    H5T_NATIVE_FLOAT,                                /* HDF number type */
    3,                                               /* rank            */
    nBlks,                                           /* dimension sizes */
    NUMBER_OF_BEAMS, 3, 1, 1, 1, 
    "Number of Blocks", "Number of Beams",           /* dimension names */
    "Number of Stokes polarizations", NULL, NULL, NULL,
    H5P_DEFAULT) );


    hsize_t start[6]={0,0,0,0,0,0};
    hsize_t count[6];

    count[0] = nBlks;
    count[1] = NUMBER_OF_BEAMS;
    count[2] = 3;

    PTB( h5d_write(grp6, "TA", TA, 3, start, count));

    return 0;
  }


  /*-------------------------------------------- */
  /* Write long accumulation calibration records */
  /* ------------------------------------------- */
  int hdf5_Aquarius::writel2_drift_corr( int32_t nBlks, 
                                         float *TA_hat_drift_corr)
  {
    hsize_t start[6]={0,0,0,0,0,0};
    hsize_t count[6];

    PTB( CreateH5D(
    grp6,                                            /* file id         */
    "TA_hat_drift_corr",                             /* short name      */
    "TA_hat instrument drift correction",            /* long name       */
    "",                                              /* standard name   */
    "",                                              /* units           */
    0, 0,                                            /* valid range     */
    0, 0,                                            /* slope           */
    -9999,                                           /* fill value      */
    H5T_NATIVE_FLOAT,                                /* HDF number type */
    3,                                               /* rank            */
    nBlks, NUMBER_OF_BEAMS, 2, 1, 1, 1,              /* dimension sizes */
    "Number of Blocks", "Number of Beams", 
    "Number of VH Polarizations", 
    NULL, NULL, NULL,                                /* dimension names */
    H5P_DEFAULT) );

    count[0] = nBlks;
    count[1] = NUMBER_OF_BEAMS;
    count[2] = 2;

    PTB( h5d_write(grp6, "TA_hat_drift_corr", TA_hat_drift_corr, 
                   3, start, count));

    return 0;
  }

} // namespace Hdf


