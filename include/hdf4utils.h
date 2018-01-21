#ifndef  HDF4UTILS_H
#define  HDF4UTILS_H

#include <netcdf.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <hdf.h>
#include <mfhdf.h>
#include <productInfo.h>
#include <passthebuck.h>

#ifdef __cplusplus
extern "C" {
#endif

/// return the sizeof dtype in bytes
int32 hdf_sizeof(int32 dtype);

char *GetFileDesc(const char *filename);

int CreateSDS(
        int32 sd_id, 
        const char *sname,    /* short name */
        const char *lname,    /* long name */
        const char *standard_name, /* NetCDF standard name (not set if passed NULL or "") */
        const char *units,    /* units (not set if passed NULL or "") */
        double low,     /* low end of valid range */
        double high,    /* high end of range (no range set if low >= high) */
        float slope,    /* scale factor (not set if 0)  */
        float offset,   /* scaling offset (not set if 0)  */
        int32 nt,       /* HDF number type */
        int32 rank,     /* number of dimensions (must be <= 3) */
        int32 d0,       /* size of 1st dimension */
        int32 d1,       /* size of 2nd dimension */
        int32 d2,       /* size of 3rd dimension (1 if rank < 2) */
        const char *dn0,      /* name of 1st dimension */
        const char *dn1,      /* name of 2nd dimension (NULL if rank < 2) */
        const char *dn2       /* name of 3rd dimension (NULL if rank < 3) */
        );
int sd_create(int32 id, const char *nam, int32 typ, int32 rank, int32 d0,
        int32 d1, int32 d2, int32 *sds_id);
int sd_select(int32 sd_id, const char *name, int32 *sds_id);
int sd_endaccess(int32 id);
int sd_setdimnames(int32 id, const char *d0, const char *d1, const char *d2);
int sd_setdimname(int32 sds_id, int32 dim_number, const char *name);
int sd_setattr(int32 id, const char *nam, int32 typ, int32 cnt, const void* data);
int sd_readdata(int32 sd_id, const char *name, VOIDP data, int32 s0, int32 s1,
        int32 s2, int32 e0, int32 e1, int32 e2);
int sd_writedata(int32 sd_id, const char *name, const VOIDP data, int32 s0,
        int32 s1, int32 s2, int32 e0, int32 e1, int32 e2);
int AddSdsToVgroup(int32 sd_id, int32 v_id, const char *name);
int v_attach(int32 h_id, int32 *v_id);
int getDims(int32 fileID, const char sdsname[], int32 dims[]);
int get_type(int32 fileID, const char sdsname[], int32 *dtype);
int rdSDS(int32 fileID, const char sdsname[], int32 start1, int32 start2,
         int32 edges1, int32 edges2, VOIDP array_data);
int getHDFattr(int32 fileID, const char attrname[], const char sdsname[],
         VOIDP data);
int32 read_SDS(int32 sdfid, const char *sds_name, void *buffer);
intn attach_vdata(int32 fid, const char *sname);

intn rdvdata(int32 vskey, const char *fields, int32 start, int32 nelt, 
        unsigned char *databuf);

#define READ_GLBL_ATTR(nam,ptr) {                                       \
  if(SDreadattr(sd_id,SDfindattr(sd_id,(nam)),(VOIDP)(ptr))){           \
    fprintf(stderr,                                                     \
    "-E- %s line %d: Could not get global attribute, %s.\n",            \
    __FILE__,__LINE__,(nam));                                           \
  }                                                                     \
}

#define READ_GLBL_ATTR_E(nam,ptr) {                                     \
    if(SDreadattr(sd_id,SDfindattr(sd_id,(nam)),(VOIDP)(ptr))){         \
    fprintf(stderr,                                                     \
    "-E- %s line %d: Could not get global attribute, %s.\n",            \
    __FILE__,__LINE__,(nam));                                           \
    exit(1);                                                            \
  }                                                                     \
}

#define READ_SDS(nam,ptr,s0,s1,s2,e0,e1,e2) {                           \
  int32 start[3];                                                       \
  int32 edge[3];                                                        \
  edge[0]=(e0); edge[1]=(e1); edge[2]=(e2);                             \
  start[0]=(s0); start[1]=(s1); start[2]=(s2);                          \
  if(SDreaddata(SDselect(sd_id, SDnametoindex(sd_id, (nam))),           \
  start, NULL, edge, (VOIDP)(ptr)) == FAIL){                            \
    fprintf(stderr,"-E- %s line %d: Could not read SDS, %s.\n",         \
    __FILE__,__LINE__,(nam));                                           \
  }                                                                     \
}

#define READ_SDS_E(nam,ptr,s0,s1,s2,e0,e1,e2) {                         \
  int32 start[3];                                                       \
  int32 edge[3];                                                        \
  edge[0]=(e0); edge[1]=(e1); edge[2]=(e2);                             \
  start[0]=(s0); start[1]=(s1); start[2]=(s2);                          \
  if(SDreaddata(SDselect(sd_id, SDnametoindex(sd_id, (nam))),           \
  start, NULL, edge, (VOIDP)(ptr)) == FAIL){                            \
    fprintf(stderr,"-E- %s line %d: Could not read SDS, %s.\n",         \
    __FILE__,__LINE__,(nam));                                           \
    exit(1);                                                            \
  }                                                                     \
}

#define READ_SDS_ID(sd_id,nam,ptr,s0,s1,s2,s3,e0,e1,e2,e3) {            \
  int32 start[4];                                                       \
  int32 edge[4];                                                        \
  edge[0]=(e0); edge[1]=(e1); edge[2]=(e2);  edge[3]=(e3);              \
  start[0]=(s0); start[1]=(s1); start[2]=(s2); start[3]=(s3);           \
  if(SDreaddata(SDselect((sd_id), SDnametoindex((sd_id), (nam))),       \
  start, NULL, edge, (VOIDP)(ptr)) == FAIL){                            \
    fprintf(stderr,"-E- %s line %d: Could not read SDS, %s.\n",         \
    __FILE__,__LINE__,(nam));                                           \
  }                                                                     \
}

#define READ_SDS_ID_E(sd_id,nam,ptr,s0,s1,s2,s3,e0,e1,e2,e3) {          \
  int32 start[4];                                                       \
  int32 edge[4];                                                        \
  edge[0]=(e0); edge[1]=(e1); edge[2]=(e2);  edge[3]=(e3);              \
  start[0]=(s0); start[1]=(s1); start[2]=(s2); start[3]=(s3);           \
  if(SDreaddata(SDselect((sd_id), SDnametoindex((sd_id), (nam))),       \
  start, NULL, edge, (VOIDP)(ptr)) == FAIL){                            \
    fprintf(stderr,"-E- %s line %d: Could not read SDS, %s.\n",         \
    __FILE__,__LINE__,(nam));                                           \
    exit(1);                                                            \
  }                                                                     \
}


/*
 * below here are the functions that used to be in hdf4_utils.h
 */

typedef struct att_info_struct {
    int32_t id;        /* HDF4 identifier for parent object  */
    int32_t index;     /* index of attribute within parent */
    char name[H4_MAX_NC_NAME];  /* name of the attribute */
    int32_t nvals;              /* number of values */
    int32_t ntype;     /* number type for stored data (hntdefs.h) */
    void *data;        /* pointer to data */
} att_struct;

typedef struct sds_info_struct {
    int32_t id;        /* HDF4 Science Data Set identifier */
    char name[H4_MAX_NC_NAME];       /* name of the HDF4 data set */
    int32_t ndims;                   /* number of dimensions */
    int32_t dimlen[H4_MAX_VAR_DIMS]; /* dimension lengths */
    int32_t ntype;     /* number type for stored data (hntdefs.h) */
    int32_t nvals;     /* number of values in data */
    void *data;        /* pointer to data */
    int32_t natts;     /* number of attributes */
    att_struct *atts;  /* pointer to array of attributes */
} sds_struct;

const char *hdf_typename(const int32_t ntype);
char *fmt_hdf_val(const void *array, const int32_t i, const int32_t ntype);
void fopen_warn(const char *filename, const char *file, int32_t const line);
void fopen_err (const char *filename, const char *file, int32_t const line);

int init_sds_byname(int32_t fileid, const char *sdsname, sds_struct *sds);
int readall_sds(sds_struct *sds);
void free_sds(sds_struct *sds);
void print_sds_info(const sds_struct sds);

int load_att_byname(int32_t obj_id, const char *attname, att_struct *att);
void print_att_vals(const att_struct att);
void print_att_info(const att_struct att);

#define EOSMETALEN 32768
int parse_odl(const char *odltext, const char *object, char *value);
int get_hdfeos_meta(int32_t sd_id, char *attribute, char *name, char *data);

#define read_att(obj_id,attname,valptr)		\
    SDreadattr(obj_id,				\
               SDfindattr(obj_id,attname),	\
               (void *) valptr)


#ifdef __cplusplus
}
#endif


#endif
