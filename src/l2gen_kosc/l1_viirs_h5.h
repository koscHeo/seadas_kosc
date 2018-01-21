#ifndef  _L1_VIIRS_H5_H
#define  _L1_VIIRS_H5_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "passthebuck.h"
#include "l1_struc.h"
#include "filehandle.h"
#include "h5io.h"

/*  define fill condition values  */
#define NA_FLOAT32_FILL   -999.9  /* Algorithm Exclusions  */
#define NA_UINT16_FILL  65535
#define MISS_FLOAT32_FILL  -999.8  /* Missing at Time of Processing  */
#define MISS_UINT16_FILL  65534
#define ONBOARD_PT_FLOAT32_FILL  -999.7  /*  Onboard Pixel Trimi (bow tie) */
#define ONBOARD_PT_UINT16_FILL  65533
#define ONGROUND_PT_FLOAT32_FILL  -999.6  /* On-ground Pixel Trim */
#define ONGROUND_PT_UINT16_FILL  65532
#define ERR_FLOAT32_FILL  -999.5  /* Cannot Calculate */
#define ERR_UINT16_FILL  65531
#define ELINT_FLOAT32_FILL  -999.4  /* Ellipsoid Intersection Failed */
#define ELINT_UINT16_FILL  65530
#define VDNE_FLOAT32_FILL  -999.3  /* Value Does Not Exist */
#define VDNE_UINT16_FILL  65529
#define SOUB_FLOAT32_FILL  -999.2  /* Scaling Out Of Bounds */
#define SOUB_UINT16_FILL  65528

int closel1_viirs_h5 (filehandle *l1file);
int openl1_viirs_h5  (filehandle *l1file);
int readl1_viirs_h5  (filehandle *l1file, int32 recnum, l1str *l1rec, int lonlat);
int gen_sdr_suite( char * );
int set_f_cal_corr( h5io_str *, filehandle *, int64_t );
int rd_vir_f_tbl( char *, int64_t, int );
int viirs_u58_yds( int64_t, short *, short *, double * );

#endif
