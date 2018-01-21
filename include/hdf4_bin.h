#ifndef hdf4_bin_h
#define hdf4_bin_h

#include "bin_util.h"
#include "meta_l3b.h"

namespace Hdf {

  class hdf4_bin {
    int32 file_id;
    int32 sd_id;
    int32 vg_id;
    int32 access_mode;

    int32 n_data_prod;
    int32 vdata_id[MAXNVDATA];
    int32 seagrid_idx;
    int32 binindex_idx;
    int32 binlist_idx;
    int32 bin_ptr;
    int32 bindata_idx;
    int32 binqual_idx;

    char proddata_name[MAXNPROD][80];
    char *product_array[MAXNPROD];
    int32 aid[MAXNPROD];
    int32  *numbin;
    int32  *basebin;
    float32 *latbin;
    bool hasQual;

    binListStruct *binList;

  public:
    hdf4_bin();
    ~hdf4_bin();

    int create( char* l3b_filename, int32 nrows);
    int open( char* l3b_filename);

    int readBinIndex( int row_num_to_read);
    int read( float32* data, binListStruct* binList);
    int read( float32* data, float32* var, binListStruct* binList);

    int read( float32* data, binListStruct* binList, int nbins_to_read);
    int read( float32* data, float32* var, binListStruct* binList, 
		int nbins_to_read);

    int readBinList( int nbins_to_read);
    int readBinList();
    int readQual( uint8* qual, int nbins_to_read);
    int readSums( float32* sums, int nbins_to_read, int iprod);

    int nprod() { return n_data_prod;}
    int get_numbin( int irow) { return numbin[irow];}
    int get_basebin( int irow) { return basebin[irow];}

    int get_bin_num( int kbin) { return binList[kbin].bin_num;}
    int get_nobs( int kbin) { return binList[kbin].nobs;}
    int get_nscenes( int kbin) { return binList[kbin].nscenes;}
    float get_weights( int kbin) { return binList[kbin].weights;}

    char *get_prodname( int i) { return proddata_name[i];}

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
    { memset( binList, 0, 2*nrows*sizeof(binListStruct)); return 0;}
    int copy_binlist( int src, int dest)
    { memcpy( &binList[dest], &binList[src], sizeof(binListStruct)); return 0;}

    int write( char *product_list, int32 nwrite, float32 *data, 
	       binListStruct* binList);
    int writeBinList( int32 nbins_to_write);
    int writeQual( uint8* qual, int nbins_to_write);
    int writeSums( float32* sums, int nbins_to_write, char *prodname);

    int copy( char *product_list, int32 nwrite, int32 *binsToCopy,
	      Hdf::binListStruct *inBinList, Hdf::hdf4_bin *input_binfile);
    int copymeta( int32 nfiles, Hdf::hdf4_bin *input_binfile[]);

    int close();

    bool has_qual();

    int32 nrows;
    int32 n_data_records;
    int32 n_active_prod;
    bool active_data_prod[MAXNVDATA];
    int32  totbins;
    int32 noext;

    binIndexStruct binIndex;
    meta_l3bType meta_l3b;
 };

}

#endif

