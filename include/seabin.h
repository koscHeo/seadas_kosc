#ifndef SEABIN_H /* avoid re-inclusion */
#define SEABIN_H

#include <math.h>
#include <ctype.h>
#include "hdf.h"
#include "mfhdf.h"
#include "meta_l3b.h"
#include "seaproto.h"

/*  Grid Parameters  */
#define	SEABIN		"SeaWiFS Binning Algorithm - first cut"
#define NBINS		4320
#define RADIUS		6378.137
#define	MAXNORTH	90.0f
#define	MAXSOUTH	-90.0f
#define	SEAMLON		-180.0f
/*#define NUMROWS		2160*/
/*#define	VSIZE		0.08333333f*/
#define GRPNAME		"Level-3 Binned Data"
#define GRPCLASS	"PlanetaryGrid"
#define GEOMNAME	"SEAGrid"
#define GEOMCLASS	"Geometry"
#define NDXNAME		"BinIndex"
#define NDXCLASS	"Index"
#define MSTRNAME	"BinList"
#define MSTRCLASS	"DataMain"
#define SLVCLASS	"DataSubordinate"
#define SNAME0          "nLw_412"
#define SNAME1          "nLw_443"
#define SNAME2          "nLw_490"
#define SNAME3          "nLw_510"
#define SNAME4          "nLw_555"
#define SNAME5          "La_670"
#define SNAME6          "CZCS_pigment"
#define SNAME7          "chlor_a"
#define SNAME8          "K_490"
#define SNAME9          "chlor_a_K_490"
#define SNAME10         "eps_78"
#define SNAME11         "tau_865"
#define GEOMSIZE        ((3*sizeof(int32))+(4*sizeof(float64)))
#define NDXSIZE         ((5*sizeof(int32))+(2*sizeof(float64)))
#define MSTRSIZE        ((3*sizeof(int16))+(2*sizeof(int32))+(sizeof(float32))+(sizeof(int8)))
/*
#define MSTRSIZE        ((4*sizeof(int16))+(1*sizeof(int32))+(sizeof(float32))+(sizeof(int8)))
#define MSTRSIZE        ((3*sizeof(int16))+(1*sizeof(int32))+(sizeof(float32)))
*/
#define SLVSIZE         (2*sizeof(float32))

/*  Global (file-level) Attributes  */
#define TITLE_VAL	 "SeaWiFS Level-3 Binned Data"
#define SENSORNM_VAL	 "SeaWiFS"
#define DCENTER_VAL	 "NASA/GSFC SeaWiFS Data Processing Center"
#define MISSION_VAL	 "SeaStar SeaWiFS"
#define MSNCHAR_VAL	 "Nominal orbit: inclination = 98.2 (Sun-synchronous); node = 12 noon local (descending); eccentricity = <0.002; altitude = 705 km; ground speed = 6.75 km/sec"
#define SENSOR_VAL	 "Sea-viewing Wide Field-of-view Sensor (SeaWiFS)"
#define SNSCHAR_VAL        "Number of bands = 8; number of active bands = 8; wavelengths per band (nm) = 412, 443, 490, 510, 555, 670, 765, 865; bits per pixel = 10; instantaneous field-of-view = 1.5835 mrad; pixels per scan = 1285; scan rate = 6/sec; sample rate = 7710/sec"

#define LATUNITS_VAL 	"degrees North"
#define LONUNITS_VAL 	"degrees East"
#define STATION_VAL  	"Wallops Flight Facility"
#define ATTR     "Attr0.0"

/*  File Types  */
#define SCENE   1
#define DAILY   2
#define OTHER   3

/*  File Type Strings  */
#define SSPAN  "Scene"
#define DSPAN  "Daily"
#define D8SPAN "8-day"
#define MSPAN  "Monthly"
#define YSPAN  "Yearly"

/*  Registration Constants  */
#define UNKNOWN	0
#define NW	1
#define NORTH	2
#define NE	3
#define WEST	4
#define CENTER	5
#define EAST	6
#define SW	7
#define SOUTH	8
#define SE	9

/*  Misc.  */
#define YES	1
#define NO	0

#define RADCONV 0.017453292

#define MAXVAL     255
#define MAXDIMS 3
#define MAX_OUT 16

#define NPARMS  25

/*#define TOTBINS 5940422*/

/*  Buffer Sizes (bins)  */
#ifndef SMALL
#define SCNINIT		150000
#define DAYINIT		500000
#define BIGINIT		1000000
#define BLKSIZE		2500
#else
#define SCNINIT		1000
#define DAYINIT		5000
#define BIGINIT		10000
#define BLKSIZE		1000
#endif

/*  Time Conversion Constants  */
#define MSECHOUR 3600000
#define MSECMIN  60000
#define MSECSEC  1000

/*  Specific to level 3 binned data read routines */

#define MAX_IN  100
#define BUFSZ   1000 

#define FSIZE           "File size"

#define PRODUCT_NAME    "Product Name"
#define TITLE           "title"
#define SENSOR_NAME     "sensor"
#define DATA_CENTER     "project"
#define MISSION         "mission"
#define UNITS           "Units"
#define PROD_TYPE       "Product Type"
#define PVERSION        "Processing Version"
#define REPLACE         "Replacement Flag"
#define SOFT_NAME       "Software Name"
#define SOFT_VER        "Software Version"
#define PTIME           "Processing Time"
#define PROC_CON        "Processing Control"
#define INPUT_PARMS     "Input Parameters"
#define FLAG_NAMES	"L2 Flag Names"
#define FLAG_USE	"L2 Flag Usage"
#define ENG_Q_USE	"L2 Engineering Quality Usage"
#define INFILES         "Input Files"

#define STIME           "time_coverage_start"
#define END_TIME        "time_coverage_end"
#define ORBIT		"Orbit"
#define S_ORBIT		"Start Orbit"
#define E_ORBIT		"End Orbit"
#define LAT_UNITS       "geospatial_lat_units"
#define LON_UNITS       "geospatial_lon_units"

#define DATABINS        "Data Bins"
#define PCT_DATABINS    "Percent Data Bins"
#endif /* SEABIN_H */



