/* 
 * File:   dfutils.h
 * Author: dshea
 *
 * Created on May 11, 2015, 10:44 AM
 */

#ifndef DFUTILS_H
#define	DFUTILS_H

#include <hdf4utils.h>
#include <netcdfutils.h>
#include <l2prod_struc.h>


#ifdef	__cplusplus
extern "C" {
#endif

typedef enum {
    DS_HDF,     // use a HDF file
    DS_NCDF     // use a netCDF file
} ds_format_t;

typedef enum {
    DS_READ,     // open file for reading
    DS_WRITE     // open/create file for writing
} ds_access_t;

typedef struct {
  int32 fid;
  int32 sid;
  ds_format_t fftype;
  int32 deflate;
} idDS;


int s2u(const char *in, char *out);
int getProdlist( const char *fname, char **prodlist, int32 *l2_flags_type);

idDS startDS(const char *filename, ds_format_t format,
        ds_access_t accessmode, int32 deflate);
idDS openDS(const char *filename);
int endDS( idDS ds_id);

#if 0
int createDS(
         idDS ds_id,
         const char *sname,    /* short name */
         const char *lname,    /* long name */
         const char *standard_name, /* NetCDF standard name (not set if passed NULL or "") */
	 char *reference;
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
#endif
         
int createDS(
	     idDS ds_id, int sensorId,
	     const char *sname,    /* short name */
	     int32 dm[3], /* dimension sizes */
	     const char dm_name[3][80] /* dimension names */
	     );

int32 selectDS( idDS ds_id, const char *l2_prod_names);
int32 checkDS( idDS ds_id, const char *l2_prod_name);
int readDS(idDS ds_id, const char *name, int32 *start, int32 *stride,
        int32 *count, VOIDP data);
int writeDS(idDS ds_id, const char *name, const VOIDP data, int32 s0, int32 s1,
        int32 s2, int32 e0, int32 e1, int32 e2);
int endaccessDS(idDS ds_id);

int fileInfo( idDS ds_id, int32 *n_datasets, int32 *n_globalattr);
int getDimsDS(idDS ds_id, const char sdsname[], int32 dims[]);
int getTypeDS( idDS ds_id, const char sdsname[H4_MAX_NC_NAME], int32 *dtype);

int setAttr( idDS ds_id, const char *nam, int32 typ, int32 cnt, const VOIDP data);
int8_t findAttr( idDS ds_id, const char *nam);
int readAttr( idDS ds_id, const char *nam, VOIDP data);
char* readAttrStr(idDS ds_id, const char *name);
int infoAttr( idDS ds_id, const char *nam, int32 *dtype, int32 *count);

int SetChrGA(idDS ds_id, const char *name, const char *value);
int SetF32GA(idDS ds_id, const char *name, float32 value);
int SetF64GA(idDS ds_id, const char *name, float64 value);
int SetI8GA(idDS ds_id, const char *name, uint8 value);
int SetI16GA(idDS ds_id, const char *name, int16 value);
int SetI32GA(idDS ds_id, const char *name, int32 value);

int16_t getDataTypeInt(productInfo_t *p_info);
int16_t *float2int16(float fbuf[], int32_t spix, int32_t npix, int32_t incr,
        float slope, float offset);
int8_t *float2int8(float fbuf[], int32_t spix, int32_t npix, int32_t incr,
        float slope, float offset);
uint16_t *float2uint16(float fbuf[], int32_t spix, int32_t npix, int32_t incr,
        float slope, float offset);
uint8_t *float2uint8(float fbuf[], int32_t spix, int32_t npix, int32_t incr,
        float slope, float offset);
    
void *scale_sds(float *data, l2prodstr *p, int nbands);
float *unscale_sds(VOIDP *data, l2prodstr *p, int32_t spix, int32_t npix, int incr);

// void *scale_sds(float *data, productInfo_t *p, int32_t nbands, int32_t npix);
// float *unscale_sds(void *data, productInfo_t *p, int32_t spix, int32_t npix,
//        int incr);

// netcdf specific routines
int CreateNCDF(
	       idDS ds_id, 
	       const char *sname,    /* short name */
	       const char *lname,    /* long name */
	       const char *standard_name, /* NetCDF standard name (not set if passed NULL or "") */
	       const char *reference,
	       const char *comment,
	       const char *units,    /* units (not set if passed NULL or "") */
	       double low,     /* low end of valid range */
	       double high,    /* high end of range (no range set if low >= high) */
	       float scale_factor,    /* scale factor (not set if 0)  */
	       float add_offset,   /* scaling offset (not set if 0)  */
	       int32_t fillValue,       /* fill value */
	       int32_t nt,       /* NCDF number type */
	       int32_t rank,     /* number of dimensions (must be <= 3) */
	       int32_t dimids[3]/* dimension ids */
	       );


#ifdef	__cplusplus
}
#endif

#endif	/* DFUTILS_H */

