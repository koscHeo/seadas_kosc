#ifndef READL2SCAN_H
#define READL2SCAN_H

#include <dfutils.h>

#define byte unsigned char

#define MAXNFILES 800 /* Increase to 544 09/25/06 JMG - 800 WDR */

#ifdef __cplusplus
extern "C" {
#endif

typedef struct l2prod_struct {

  int32    fileindex;

  char     filename[128];
  char     oformat[32];

  int32    nrec;
  int32    nsamp;

  int16    syear;
  int16    sday;
  int32    smsec;
  int16    eyear;
  int16    eday;
  int32    emsec;
  int32    orbit;
  char     dtype[8];

  int32    ntilts;
  int16    tilt_flags[20];
  int16    tilt_ranges[2][20];

  char     *flagnames;
  int32    flagmask;

  int32    year;
  int32    day;
  int32    msec;

  int32    *year_cache;
  int32    *day_cache;
  int32    *msec_cache;

  float32  *geoloc;
  float32  *longitude;
  float32  *latitude;
  float32  *geonav[6];

  float32  *lon_cntl;
  float32  *lat_cntl;
  float32  *cntl_pnts;
  float32  *cntl_pnts_cache;
  float32  *spline_arr;

  int32    nprod;
  char     *prodname[1000];
  float32  bv_unscaled[1000];
  int16    bv_scaled[1000];

  float32  *l2_data;
  int32    *l2_flags;

  byte     eng_qual[4];
  byte     s_flags[4];
  int32    nflag[8];

  int32    geointerp;

  byte    *mside;
  byte    *detnum;
  int32    *pixnum;


} l2_prod;


typedef struct meta_l2Struct {
        char    *product_name; /* ATTR Product name(file name)         */
        char    *title;
        char    *data_center;   /* ATTR data_center, processing center  */
        char    *mission;       /* ATTR mission                         */
        char    *mission_char;  /* ATTR Mission Characteristics         */
        char    *sensor_name;   /* ATTR sensor name                     */
        char    *sensor;        /* ATTR sensor                          */
        char    *sensor_char;   /* ATTR instrumentInformation          */
        char    *sw_id;         /* ATTR Software ID                     */
        char    *infiles;       /* ATTR Input files                     */
        char    *stime;         /* ATTR Start time                      */
        char    *etime;         /* ATTR End time                        */
        char    *ctime;         /* ATTR scene center time               */
        char    *ntime;         /* ATTR Node crossing time              */
        char    *snode;         /* ATTR Start Node                      */
        char    *enode;         /* ATTR End Node                        */
        int     orbnum;         /* ATTR orbit number                    */
        char    *norad1;        /* ATTR NORAD elements, first line      */
        char    *norad2;        /* ATTR NORAD elements, second line     */
        int     pix_start;      /* ATTR LAC Pixel Start Number          */
        int     pix_sub;        /* ATTR LAC Pixel Subsampling           */
        int     ncrec;          /* ATTR scene center scan line          */
        int     nfrec;          /* ATTR number of filled scan line      */
        byte    ff_mis;         /* ATTR FF missing frames               */
        byte    sd_mis;         /* ATTR SDPS missing frames             */
        float   flags_pc[32];   /* MFSD % data for each quality flag    */
        char    *lat_units;     /* ATTR Latitude units                  */
        char    *lon_units;     /* ATTR Longitude units                 */
        float   northlat;       /* ATTR Northernmost latitude           */
        float   southlat;       /* ATTR Southernmost latitude           */
        float   westlon;        /* ATTR Westernmost longitude           */
        float   eastlon;        /* ATTR Easternmost longitude           */
        float   startclat;      /* ATTR Start Center Latitude           */
        float   startclon;      /* ATTR Start Center Longitude          */
        float   endclat;        /* ATTR End Center Latitude             */
        float   endclon;        /* ATTR End Center Longitude            */
        float   nodel;          /* ATTR Orbit node longitude            */
        int     ntilts;         /* MFSD Sensor Tilt                     */
        /* Calibration Vgroup */
        short   entry_year;
        short   entry_day;
        short   ref_year;
        short   ref_day;
        short   ref_minute;
} meta_l2Type;


/* Prototypes */
void free_rowgroup_cache();
void init_rowgroup_cache();
int32 get_dtype(int32 dtype, ds_format_t fileformat);
int32 openL2(char *, char *, l2_prod *);
int32 reopenL2(int32, l2_prod *);
int32 readL2(l2_prod *l2_str, int32 ifile, int32 recnum, int32 iprod, 
             unsigned char *scan_in_rowgroup);
int32 readlonlat(l2_prod *l2_str, int32 ifile, int32 *start, int32 *edges, 
                 unsigned char *scan_in_rowgroup);
int32 closeL2(l2_prod *, int32);
int32 freeL2(l2_prod *);
int32 findprod(l2_prod *, char* );
int32 readL2meta(meta_l2Type *, int32);
int32 freeL2meta(meta_l2Type *);
int32 getL3units(l2_prod *, int32, char *, char *);
#ifdef __cplusplus
}
#endif
#endif
