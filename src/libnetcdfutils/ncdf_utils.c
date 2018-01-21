/* =========================================================== */
/* Module ncdf_utils.c                                         */
/*                                                             */
/* NCDF4 I/O utilities.                                        */
/*                                                             */ 
/* Written By:                                                 */
/*     Joel Gales, Futurtech                                   */
/*                                                             */
/* Modification History:                                       */
/*     Joel Gales, Futuretech, OBPG Project, Nov 2013.         */
/*           Add support for CF-compliant metadata             */
/* =========================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h> 
#include <time.h>
#include <math.h>
#include <netcdf.h>
#include <netcdfutils.h>

void
report_err(const int stat, const int line, const char *file) {
  if (stat != NC_NOERR) {
    (void)fprintf(stderr,"line %d of %s: %s\n", line, file, nc_strerror(stat));
    fflush(stderr);
  }
}

void
check_err(const int stat, const int line, const char *file) {
  if (stat != NC_NOERR) {
    (void)fprintf(stderr,"line %d of %s: %s\n", line, file, nc_strerror(stat));
    fflush(stderr);
    exit(1);
  }
}

int writeBinList_nc( int32_t deflate, int32_t grpid, int32_t nbins_to_write, 
		     const void *data) {
  // Write BinList DataSet
  int status;
  int varid;

  typedef struct {
    uint32_t binnum;
    int16_t  nobs;
    int16_t  nscenes;
    float    weights;
    float    time_rec;
  }  binListStruct;
  binListStruct binList;

  nc_type binListType;

  int binListDim;

  static size_t startp;
  size_t countp;

  if ( nc_inq_varid( grpid, "BinList", &varid) != NC_NOERR) {
    status = nc_def_compound( grpid, sizeof(binListStruct), 
			      "binListType", &binListType);
    check_err(status,__LINE__,__FILE__);

    status = nc_insert_compound( grpid, binListType, "bin_num", 
				 NC_COMPOUND_OFFSET( binListStruct, binnum), 
				 NC_UINT);
    check_err(status,__LINE__,__FILE__);

    status = nc_insert_compound( grpid, binListType, "nobs", 
				 NC_COMPOUND_OFFSET( binListStruct, nobs), 
				 NC_SHORT);
    check_err(status,__LINE__,__FILE__);

    status = nc_insert_compound( grpid, binListType, "nscenes", 
				 NC_COMPOUND_OFFSET( binListStruct, nscenes), 
				 NC_SHORT);
    check_err(status,__LINE__,__FILE__);

    status = nc_insert_compound( grpid, binListType, "weights", 
				 NC_COMPOUND_OFFSET( binListStruct, weights), 
				 NC_FLOAT);
    check_err(status,__LINE__,__FILE__);

    status = nc_insert_compound( grpid, binListType, "time_rec", 
				 NC_COMPOUND_OFFSET( binListStruct, time_rec), 
				 NC_FLOAT);
    check_err(status,__LINE__,__FILE__);

    status = nc_def_dim( grpid, "binListDim", NC_UNLIMITED, &binListDim);
    check_err(status,__LINE__,__FILE__);

    status = nc_def_var( grpid, "BinList", binListType, 1, &binListDim, 
			 &varid);
    check_err(status,__LINE__,__FILE__);


    /* First set chunking */
    size_t chunksize = 256;
    status = nc_def_var_chunking( grpid, varid, NC_CHUNKED, &chunksize);
    check_err(status,__LINE__,__FILE__);
    if (status != NC_NOERR) exit(1);

    if ( deflate > 0) {
      /* Now we can set compression */
      status = nc_def_var_deflate( grpid, varid, NC_SHUFFLE, 1, deflate);
      check_err(status,__LINE__,__FILE__);
      if (status != NC_NOERR) exit(1);
    } 

    startp = 0;
  }

  countp = nbins_to_write;
  status = nc_put_vara( grpid, varid, &startp, &countp, data);
  check_err(status,__LINE__,__FILE__);
  startp += countp;

  return 0;
}


int writeBinData_nc( int32_t deflate, int32_t grpid, int32_t nbins_to_write, 
		     int32_t iprod, const char *prodname, const void *data) {
  int status=0;
  int varid;

  static nc_type binDataType=-1;

  static int binDataDim;

  static size_t startp[256];
  size_t countp;
  int dimidsp;
  if ( nc_inq_varid( grpid, prodname, &varid) != NC_NOERR) {
    if ( binDataType == -1) {
      status = nc_def_compound( grpid, 8, 
				"binDataType", &binDataType);
      check_err(status,__LINE__,__FILE__);
      status = nc_insert_compound( grpid, binDataType, "sum", 0, NC_FLOAT);
      check_err(status,__LINE__,__FILE__);

      status = nc_insert_compound( grpid, binDataType, "sum_squared", 4, 
				   NC_FLOAT);
      check_err(status,__LINE__,__FILE__);

      status = nc_def_dim( grpid, "binDataDim", NC_UNLIMITED, &binDataDim);
      check_err(status,__LINE__,__FILE__);
    }

    status = nc_def_var( grpid, prodname, binDataType, 1, &binDataDim, 
			 &varid);
    if ( status != NC_NOERR) {
      (void) fprintf(stderr,"line %d of %s: %s\n", __LINE__, __FILE__, 
		     nc_strerror(status));
      (void) fprintf(stderr,"product %s\n", prodname);
      fflush(stderr);

      //exit(1);
    }

    /* First set chunking */
    size_t chunksize = 256;
    status = nc_def_var_chunking( grpid, varid, NC_CHUNKED, &chunksize);
    check_err(status,__LINE__,__FILE__);
    if (status != NC_NOERR) exit(1);

    if ( deflate > 0) {
      /* Now we can set compression */
      status = nc_def_var_deflate( grpid, varid, NC_SHUFFLE, 1, deflate);
      check_err(status,__LINE__,__FILE__);
      if (status != NC_NOERR) exit(1);
    } 

    startp[iprod] = 0;

    }

  countp = nbins_to_write;
  status = nc_put_vara( grpid, varid, &startp[iprod], &countp, data);
  check_err(status,__LINE__,__FILE__);
  startp[iprod] += countp;

  return 0;
}


int writeBinIndex_nc( int32_t grpid, int32_t n_write, const void *data) {
  // Write BinIndex DataSet
  int status;
  int varid;

  nc_type binIndexType;

  int binIndexDim;

  static size_t startp;
  size_t countp;

  if ( nc_inq_varid( grpid, "BinIndex", &varid) != NC_NOERR) {
    status = nc_def_compound( grpid, 16, "binIndexType", &binIndexType);
    check_err(status,__LINE__,__FILE__);

    status = nc_insert_compound( grpid, binIndexType, "start_num",
				 0, NC_UINT);
    check_err(status,__LINE__,__FILE__);

    status = nc_insert_compound( grpid, binIndexType, "begin",
				 4, NC_UINT);
    check_err(status,__LINE__,__FILE__);

    status = nc_insert_compound( grpid, binIndexType, "extent",
				 8, NC_UINT);
    check_err(status,__LINE__,__FILE__);

    status = nc_insert_compound( grpid, binIndexType, "max",
				 12, NC_UINT);
    check_err(status,__LINE__,__FILE__);

    status = nc_def_dim( grpid, "binIndexDim", NC_UNLIMITED, &binIndexDim);
    check_err(status,__LINE__,__FILE__);

    status = nc_def_var( grpid, "BinIndex", binIndexType, 1, &binIndexDim, 
			 &varid);
    check_err(status,__LINE__,__FILE__);

    /* First set chunking */
    size_t chunksize = 256;
    status = nc_def_var_chunking( grpid, varid, NC_CHUNKED, &chunksize);
    check_err(status,__LINE__,__FILE__);
    if (status != NC_NOERR) exit(1);

//    if ( deflate > 0) {
//      /* Now we can set compression */
//      status = nc_def_var_deflate( grpid, varid, NC_SHUFFLE, 1, deflate);
//      check_err(status,__LINE__,__FILE__);
//      if (status != NC_NOERR) exit(1);
//    }

    startp = 0;
  }

  countp = n_write;
  status = nc_put_vara( grpid, varid, &startp, &countp, data);
  check_err(status,__LINE__,__FILE__);
  startp += countp;

  return 0;
}


int writeQuality_nc( int32_t deflate, int32_t grpid, int32_t nbins_to_write, const void *data) {
  // Write Quality DataSet
  int status;
  int varid;

  nc_type qualityType;

  int qualityDim;

  static size_t startp;
  size_t countp;

  if ( nc_inq_varid( grpid, "qual_l3", &varid) != NC_NOERR) {
    //    status = nc_def_compound( grpid, 1, "qualityType", &qualityType);
    //check_err(status,__LINE__,__FILE__);

    //   status = nc_insert_compound( grpid, qualityType, "qual_l3",
    //				 0, NC_UINT);
    //check_err(status,__LINE__,__FILE__);

    status = nc_def_dim( grpid, "qualityDim", NC_UNLIMITED, &qualityDim);
    check_err(status,__LINE__,__FILE__);

    //    status = nc_def_var( grpid, "qual_l3", qualityType, 1, &qualityDim, 
    //			 &varid);

    status = nc_def_var( grpid, "qual_l3", NC_BYTE, 1, &qualityDim, 
			 &varid);
    check_err(status,__LINE__,__FILE__);
    
    /* First set chunking */
    size_t chunksize = 256;
    status = nc_def_var_chunking( grpid, varid, NC_CHUNKED, &chunksize);
    check_err(status,__LINE__,__FILE__);
    if (status != NC_NOERR) exit(1);

    if ( deflate > 0) {
      /* Now we can set compression */
      status = nc_def_var_deflate( grpid, varid, NC_SHUFFLE, 1, deflate);
      check_err(status,__LINE__,__FILE__);
      if (status != NC_NOERR) exit(1);
    } 

    startp = 0;
  }

  countp = nbins_to_write;
  status = nc_put_vara( grpid, varid, &startp, &countp, data);
  check_err(status,__LINE__,__FILE__);
  startp += countp;

  return 0;
}


