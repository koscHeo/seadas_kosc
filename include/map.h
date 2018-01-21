#ifndef MAP_H /* avoid re-inclusion */
#define MAP_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "hdf.h"
#include "mfhdf.h"
#include "palette.h"

#define MAXVAL     255

#define NOMATCH_ERR  -2

/*  Time Conversion Constants  */
#define MSECHOUR 3600000
#define MSECMIN  60000
#define MSECSEC  1000

/* define meta_l3m structure */
typedef struct{
 	char	*pname;
	char	*title;
	char 	*dcenter;
	char 	*station;
	float32	station_lat;
	float32 station_lon;
	char	*mission;
	char 	*mission_char;
	char	*sensor;
	char	*sensor_name;
	char	*sensor_char;
	char	*prodtype;
	char	*replace;
	char	*softid;
	char	*softnm;
	char	*softver;
	char	*ptime;
	char	*infiles;
	char	*proccon;
	char	*proclog;
	char    *flag_use;
	uint8	eng_q_use[4];
	int16   bin_syear;
	int16   bin_sday;
	int16   bin_eyear;
	int16   bin_eday;
 	char    *stime;
	char	*etime;
	int16 	syear;
	int16 	sday;
	int32 	smsec;
	int16	eyear;
	int16	eday;
	int32	emsec;
	int32   orbit;
	int32   start_orb;
	int32   end_orb;
	char	*mapproj;
	char	*lat_units;
	char	*lon_units;
	float32 nlat;
	float32 slat;
	float32 elon;
	float32 wlon;
	float32 lat_step;
	float32 lon_step;
	int32	nbins;
	int32	nrows;
	int32	ncols;
	char    *parameter;
	char 	*measure;
	char    *units;
	char    *scaling_type;
	char 	*scaling_eqn;
	char	*base;
	float32 slope;
	float32 intercept;
	float32 scaled_data_min;
	float32 scaled_data_max;
	char   *scaled_data_type;
	float32 data_min;
	float32 data_max;
  int32   bintype; 
   }meta_struct;

#endif /* MAP_H */

