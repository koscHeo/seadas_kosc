#ifndef METAL3B_H /* avoid re-inclusion */
#define METAL3B_H

#include <dfutils.h>
#include <hdf5.h>

#define SM_ATTRSZ	1024
#define MD_ATTRSZ 10000
#define LG_ATTRSZ MAX_ORDER        /* MAX_ORDER is defined in HDF as 65535 */

#ifdef  __cplusplus
extern "C" {
#endif

typedef struct {
    char product_name[SM_ATTRSZ];
    char title[SM_ATTRSZ];
    int sensorID;
    char sensor_name[SM_ATTRSZ];
    char data_center[SM_ATTRSZ];
    char mission[SM_ATTRSZ];
    char mission_char[SM_ATTRSZ];
    char sensor[SM_ATTRSZ];
    char sensor_char[SM_ATTRSZ];
    char station[SM_ATTRSZ];
    float32 station_lat;
    float32 station_lon;
    char units[SM_ATTRSZ];
    char prod_type[SM_ATTRSZ];
    char pversion[SM_ATTRSZ];
    char replace[SM_ATTRSZ];
    char soft_name[SM_ATTRSZ];
    char soft_ver[SM_ATTRSZ];
    char ptime[SM_ATTRSZ];
    char proc_con[MD_ATTRSZ];
    char input_parms[MAX_ORDER];
    char flag_names[SM_ATTRSZ];
    char infiles[LG_ATTRSZ];
    double startTime;
    double endTime;
    int32 orbit;
    int32 start_orb;
    int32 end_orb;
    char lat_units[SM_ATTRSZ];
    char lon_units[SM_ATTRSZ];
    int32 data_bins;
    float32 pct_databins;
    char binning_scheme[SM_ATTRSZ];
    float32 north;
    float32 south;
    float32 east;
    float32 west;
    double geospatial_resolution;
    char bin_resolution[32];
} meta_l3bType;

void calculate_temporal_range(meta_l3bType *meta_l3b);
int write_l3b_meta_hdf4(int32 sd_id, meta_l3bType *meta_l3b);
int write_l3b_meta_netcdf4(idDS ds_id, meta_l3bType *meta_l3b);

int32 rdattr(int32 sdfid, char *attr_name, void *buf);
int32 getattrsz(int32 id, char *attr_name, int32 *nt, int32 *count);
int32 read_l3b_meta_hdf4(int32 sdfid, meta_l3bType *meta_l3b);
int32 read_l3b_meta_hdf5(hid_t grp0, meta_l3bType *meta_l3b);
int32 read_l3b_meta_netcdf4(int ncid, meta_l3bType *meta_l3b);


#ifdef  __cplusplus
}
#endif

#endif  /* METAL3B_H */ 
