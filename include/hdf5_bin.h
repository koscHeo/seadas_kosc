#ifndef hdf5_bin_h
#define hdf5_bin_h

#include "bin_util.h"
#include "meta_l3b.h"

namespace Hdf {
 
  class hdf5_bin {
    hid_t h5fid;
    hid_t grp0, grp1;
    hid_t access_mode;

    int32 n_datasets;
    hid_t h5table_id[3][MAXNPROD];

    int n_data_prod;

    int32 bindata_idx;
    int32 binlist_idx;
    int32 binindex_idx;
    char proddata_name[MAXNPROD][80];

    bool active_data_table[MAXNPROD];

    char *product_array[MAXNPROD];

    int32  *numbin;
    int32  *basebin;
    float *latbin;

    binListStruct_hdf5 *binList;

  public:
    hdf5_bin();
    ~hdf5_bin();

    int open( const char* l3b_filename);
    int create( const char* l3b_filename, int32 nrows);

    int readBinIndex( int row_num_to_read);
    int readBinList( int nbins_to_read);
    int readQual( unsigned char* qual, int nbins_to_read);
    int readSums( float32* sums, int nbins_to_read, int iprod);

    int writeBinList( hsize_t nbins_to_write);

    int writeSums( float* sums, hsize_t nbins_to_write, const char *prodname);

    int write( const char *product_list, hsize_t nwrite, float *data,
	       binListStruct_hdf5* binList);
    int close();

    int nprod() { return n_data_prod;}
    int get_numbin( int irow) { return numbin[irow];}
    int get_basebin( int irow) { return basebin[irow];}

    hid_t get_index_table() { return h5table_id[0][binindex_idx];}
    hid_t get_list_table() { return h5table_id[0][binlist_idx];}
    hid_t get_data_table( int i) { return h5table_id[0][bindata_idx+i];}
    hid_t get_grp0() { return grp0;}

    int get_bin_num( int kbin) { return binList[kbin].bin_num;}
    int get_nobs( int kbin) { return binList[kbin].nobs;}
    int get_nscenes( int kbin) { return binList[kbin].nscenes;}
    float get_weights( int kbin) { return binList[kbin].weights;}

    char *get_prodname( int i) { return proddata_name[i];}

    int copymeta( int32 nfiles, Hdf::hdf5_bin *input_binfile[]);

    int set_bin_num( int offset, int bin_num)
    { binList[offset].bin_num = bin_num; return 0;}
    int inc_nobs( int offset, int nobs)
    { binList[offset].nobs += nobs; return 0;}
    int inc_nscenes( int offset, int nscenes)
    { binList[offset].nscenes += nscenes; return 0;}
    int inc_weights( int offset, float weights)
    { binList[offset].weights += weights; return 0;}
    int set_weights( int offset, float weights)
    { binList[offset].weights = weights; return 0;}

    int clear_binlist() 
    { memset( binList, 0, 2*nrows*sizeof(binListStruct_hdf5)); return 0;}
    int copy_binlist( int src, int dest)
    { memcpy( &binList[dest], &binList[src], sizeof(binListStruct_hdf5));
      return 0;}
  
    int setDataPtr( int nbins_to_read)
    { binDataPtr += nbins_to_read; return 0;}

    int incNumRec( int n_write)
    { n_data_records += n_write; return 0;}

    int32 totbins;
    int32 nrows;
    int32 nbins;
    hsize_t n_data_records;
    hsize_t n_active_prod;
    bool active_data_prod[MAXNVDATA];

    hsize_t binListPtr;
    hsize_t binDataPtr;

    binIndexStruct binIndex;

    meta_l3bType meta_l3b;
 };

}

#endif

