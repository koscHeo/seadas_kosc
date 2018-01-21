#ifndef bin_util_h
#define bin_util_h

#pragma GCC diagnostic ignored "-Wpadded"
#ifndef PI
#define PI  3.141592653589793
#endif

#include <sstream>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <sstream>

#include "hdf.h"
#include "mfhdf.h"
#include "hdf5.h"

#define MAXNPROD 256
#define MAXNVDATA MAXNPROD+3


#define SWAP_4(x) ( ((x) << 24) | \
         (((x) << 8) & 0x00ff0000) | \
         (((x) >> 8) & 0x0000ff00) | \
         ((x) >> 24) )

using namespace std;

namespace Hdf {

  typedef struct binListStruct {
    int32 bin_num;
    int16 nobs;
    int16 nscenes;
    int16 time_rec;
    float32 weights;
    uchar8 sel_cat;
    int32 flags_set;
    float32 lat;
    float32 lon;
  } binListStructure;

  typedef struct binIndexStruct {
    int32 row_num; // 0-based
    float64 vsize;
    float64 hsize;
    int32 basebin;
    int32 beg;
    int32 ext;
    int32 numbin;
  } binIndexStructure;

  typedef struct binIndexStruct_cdf4 {
    int32 basebin;
    int32 beg;
    int32 ext;
    int32 numbin;
  } binIndexStructure_cdf4;

  typedef struct binListStruct_hdf5 {
    int32_t bin_num;
    short nobs;
    short nscenes;
    float weights;
    int64_t flags_set;
  } binListStructure_hdf5;

  typedef struct binListStruct_cdf4 {
    uint32_t bin_num;
    short nobs;
    short nscenes;
    float weights;
    float time_rec;
  } binListStructure_cdf4;

  int create_vdata( int32 file_id, int32 vg_id, 
		    int32 *vdata_id, const char *vdata_name, 
		    const char *class_name,
		    int32 n_flds, char const * const fldname[], int32 type[],
		    int32 noext, int32 *aid);

  int32 write_vdata( int vdata_id, int32 n_recs_to_write, void *data);

  int read_binList( int n_elem, int32 vdata_id_binlist, 
		    binListStruct* binList);

  int write_binList( int n_elem, int32 vdata_id_binlist, 
		     binListStruct* binList);

  int write_prodData( int n_elem, int32 vdata_id_proddata, float32 *data, 
		      binListStruct* binList);
  int copy_prodData( int n_elem, int32 *binsToCopy, 
		     char const * const fldname3[],
		     int32 in_vdata_id_proddata, int32 out_vdata_id_proddata);
  int create_compound( hid_t group_id, const char *dataset_name,
		       hid_t *dataset_id, hid_t *type_id, size_t typesize, 
		       int32_t n_flds, char const * const fldname[], 
		       size_t offset[],
		       hid_t type[], hid_t *filespace, hid_t dataspace);
}

#endif
