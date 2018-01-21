#ifndef hdf_bin_h
#define hdf_bin_h

#pragma GCC diagnostic ignored "-Wpadded"
#include "bin_util.h"
#include "meta_l3b.h"


namespace Hdf {
  class hdf_bin;


  hdf_bin* openBinObject(const char* binFileName);


  class hdf_bin {

  protected:
    int32 n_data_prod;
    char proddata_name[MAXNPROD][80];
    char *product_array[MAXNPROD];

    binIndexStruct binIndex;
    int32 *numbin;
    int32 *basebin;
    float *latbin;

    size_t binListPtr;
    size_t lastBinListPtr;

  public:
    hdf_bin();
    virtual ~hdf_bin();

    virtual int query();
    virtual int query( char* product_list);
    virtual int query( char ***prod_array);
    virtual int get_prodname( int iprod, char *prodname);

    virtual const char* getProdName(int prodNum) const;
    virtual int getProdIndex(const char *prodname) const;
    virtual const char* getActiveProdName(int prodNum) const;

    virtual int read( char* product_list);

    virtual int get_beg()=0;
    virtual int get_ext()=0;

    virtual int open( const char* l3b_filename)=0;
    virtual int create( const char* l3b_filename, int32 nrows)=0;

    virtual int readBinIndex( int row_num_to_read)=0;
    virtual int readBinList( int nbins_to_read)=0;
    virtual int readBinList( int nbins_to_read, int32 list_reset_ptr)=0;
    virtual int readBinList()=0;

    virtual int writeBinList( int32 nbins_to_write)=0;
    virtual int readQual( uint8* qual, int32 nbins_to_read)=0;
    virtual int readQual( uint8* qual, int32 nbins_to_read, 
    			  int32 row_num_to_read)=0;

    /**
     * Read Bin File Product Data
     * @param sums array to place sum and sumSquares for each product
     * @param nbins_to_read number of consecutive bins to read
     * @param iprod product index or -1 to read all active products set by read()
     * @return 0 if OK
     */
    virtual int readSums( float32* sums, int32 nbins_to_read, int iprod)=0;
    virtual int readSums( float32* sums, int32* listOfBins, int32 nbins_to_read, 
			  int iprod)=0;
    virtual int writeQual( uint8* qual, int32 nbins_to_write)=0;
    virtual int writeSums( float32* sums, int32 nbins_to_write, 
			   const char *prodname)=0;

    virtual int get_numbin( int irow)=0;
    virtual int get_basebin( int irow)=0;

    virtual int get_bin_num( int kbin)=0;
    virtual int get_nobs( int kbin)=0;
    virtual int get_nscenes( int kbin)=0; 
    virtual float get_weights( int kbin)=0;
    virtual float get_time_rec( int kbin)=0;

    virtual hid_t get_index_table()=0;
    virtual hid_t get_list_table()=0;
    virtual hid_t get_data_table( int i)=0;

    virtual int clear_binlist()=0;
    virtual int copy_binlist( int src, int dest)=0;
    virtual int set_bin_num( int offset, int bin_num)=0;
    virtual int inc_nobs( int offset, int nobs)=0;
    virtual int inc_nscenes( int offset, int nscenes)=0;
    virtual int inc_weights( int offset, float weights)=0;
    virtual int set_weights( int offset, float weights)=0;
    virtual int inc_time_rec( int offset, float time_rec)=0;
    virtual bool has_qual()=0;
    virtual int setDataPtr( int nbins_to_read)=0;
    virtual int setDataPtrAbsolute(int32 recordNum)=0;
    virtual int incNumRec( int n_write)=0;
    virtual int close()=0;

    virtual int32 nprod() { return n_data_prod;}
    virtual int32 get_list_ptr() { return binListPtr;}

    virtual int copymeta( int32 nfiles, Hdf::hdf_bin *input_binfile[]);

    int32 totbins;
    int32 nrows;

    bool active_data_prod[MAXNVDATA];
    int32 n_data_records;
    int32 n_active_prod;
    bool isHDF5;
    bool isCDF4;
    bool hasQual;
    bool hasNoext;

    uint32_t deflate;

    meta_l3bType meta_l3b;
  };


  class hdf4_bin: public hdf_bin {
    int32 file_id;
    int32 sd_id;
    int32 vg_id;
    int32 access_mode;

    int32 vdata_id[MAXNVDATA];
    int32 seagrid_idx;
    int32 binindex_idx;
    int32 binlist_idx;
    int32 bin_ptr;
    int32 list_reset_ptr;
    int32 bindata_idx;
    int32 binqual_idx;

    int32 aid[MAXNPROD];

    binListStruct *binList;
    int lastNumBins;

  public:
    hdf4_bin();
    virtual ~hdf4_bin();

    int get_beg() { return binIndex.beg;}
    int get_ext() { return binIndex.ext;}

    int create( const char* l3b_filename, int32 nrows);
    int open( const char* l3b_filename);

    int read( char* product_list) { return hdf_bin::read( product_list);}
    int copymeta( int32 nfiles, Hdf::hdf4_bin *input_binfile[]) 
    {return hdf_bin::copymeta( nfiles, (Hdf::hdf_bin **) input_binfile);}

    int readBinIndex( int row_num_to_read);
    int read( float32* data, binListStruct* binList);
    int read( float32* data, float32* var, binListStruct* binList);
    int read( float32* data, binListStruct* binList, int nbins_to_read);
    int read( float32* data, float32* var, binListStruct* binList, 
    	      int nbins_to_read);

    int readBinList( int nbins_to_read);
    int readBinList();
    int readBinList( int nbins_to_read, int32 list_reset_ptr);
    int readQual( uint8* qual, int32 nbins_to_read);
    int readQual( uint8* qual, int32 nbins_to_read, int32 row_num_to_read);

    int readSums( float32* sums, int32 nbins_to_read, int iprod);
    int readSums( float32* sums, int32* listOfBins, int32 nbins_to_read, 
		  int iprod);
    int get_numbin( int irow) { return numbin[irow];}
    int get_basebin( int irow) { return basebin[irow];}

    int get_bin_num( int kbin) { return binList[kbin].bin_num;}
    int get_nobs( int kbin) { return binList[kbin].nobs;}
    int get_nscenes( int kbin) { return binList[kbin].nscenes;}
    float get_weights( int kbin) { return binList[kbin].weights;}
    float get_time_rec( int kbin) { return binList[kbin].time_rec;}

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
    int inc_time_rec( int offset, float time_rec)
    { binList[offset].time_rec += time_rec; return 0;}

    int clear_binlist() 
    { memset( binList, 0, 2*nrows*sizeof(binListStruct)); return 0;}

    int copy_binlist( int src, int dest)
    { memcpy( &binList[dest], &binList[src], sizeof(binListStruct)); return 0;}

    int write( char *product_list, int32 nwrite, float32 *data, 
	       binListStruct* binList);
    int writeBinList( int32 nbins_to_write);
    int writeQual( uint8* qual, int32 nbins_to_write);
    int writeSums( float32* sums, int32 nbins_to_write, const char *prodname);

    int copy( char *product_list, int32 nwrite, int32 *binsToCopy,
	      Hdf::binListStruct *inBinList, Hdf::hdf4_bin *input_binfile);

    int close();

    bool has_qual();

    int setDataPtr( int nbins_to_read) {return 0;}
    int setDataPtrAbsolute(int32 recordNum);

    int incNumRec( int n_write) {return 0;}

    hid_t get_index_table() { return 0;}
    hid_t get_list_table() { return 0;}
    hid_t get_data_table( int i) { return 0;}

    int32 noext;

 };


  class hdf5_bin: public hdf_bin { 
    hid_t h5fid;
    hid_t grp0, grp1;
    hid_t access_mode;

    int32 n_datasets;
    hid_t h5table_id[3][MAXNPROD];

    int32 bindata_idx;
    int32 binlist_idx;
    int32 binindex_idx;

    binListStruct_hdf5 *binList;

    int lastNumBins;

  public:
    hdf5_bin();
    ~hdf5_bin();

    int get_beg() { return binIndex.beg;}
    int get_ext() { return binIndex.ext;}

    int open( const char* l3b_filename);
    int create( const char* l3b_filename, int32 nrows);

    int readBinIndex( int row_num_to_read);
    int readBinList( int nbins_to_read);
    int readBinList( int nbins_to_read, int32 list_reset_ptr);
    int readBinList( );

    int readQual( unsigned char* qual, int32 nbins_to_read);
    int readQual( unsigned char* qual, int32 nbins_to_read, 
		  int32 row_num_to_read);
    int readSums( float32* sums, int32 nbins_to_read, int iprod);
    int readSums( float32* sums, int32* listOfBins, int32 nbins_to_read, 
		  int iprod);

    int writeBinList( int32 nbins_to_write);

    int writeSums( float* sums, int32 nbins_to_write, const char *prodname);
    int writeQual( uint8* qual, int32 nbins_to_write) {return 0;}

    int write( const char *product_list, hsize_t nwrite, float *data,
	       binListStruct_hdf5* binList);
    int close();

    bool has_qual() {return false;}

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
    float get_time_rec( int kbin) { return 0.0;}

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
    int inc_time_rec( int offset, float time_rec)
    { return 0;}

    int clear_binlist() 
    { memset( binList, 0, 2*nrows*sizeof(binListStruct_hdf5)); return 0;}
    int copy_binlist( int src, int dest)
    { memcpy( &binList[dest], &binList[src], sizeof(binListStruct_hdf5));
      return 0;}
  
    int setDataPtr( int nbins_to_read)
    { binDataPtr += nbins_to_read; return 0;}
    int setDataPtrAbsolute(int32 recordNum)
    { binDataPtr = recordNum; return 0; }

    int incNumRec( int n_write)
    { n_data_records += n_write; return 0;}
 
    hsize_t binDataPtr;
 };


  class cdf4_bin: public hdf_bin { 
    int ncid;
    int grp0, grp1;
    int access_mode;

    int n_datasets;

    int bindata_idx;
    int binqual_idx;
    int binlist_idx;
    int binindex_idx;

    binListStruct_cdf4 *binList;
    int lastNumBins;

  protected:
    binIndexStruct_cdf4 binIndex;

  public:
    cdf4_bin();
    ~cdf4_bin();

    int get_beg() { return binIndex.beg;}
    int get_ext() { return binIndex.ext;}

    int open( const char* l3b_filename);
    int create( const char* l3b_filename, int32 nrows);

    int readBinIndex( int row_num_to_read);
    int readBinList( int nbins_to_read);
    int readBinList( int nbins_to_read, int32 list_reset_ptr);
    int readBinList();

    int readQual( unsigned char* qual, int32 nbins_to_read);
    int readQual( unsigned char* qual, int32 nbins_to_read, 
		  int32 row_num_to_read);
    int readSums( float32* sums, int32 nbins_to_read, int iprod);
    int readSums( float32* sums, int32* listOfBins, int32 nbins_to_read, 
		  int iprod);

    int writeBinList( int32 nbins_to_write);

    int writeSums( float* sums, int32 nbins_to_write, const char *prodname);
    int writeQual( uint8* qual, int32 nbins_to_write);

    int write( const char *product_list, hsize_t nwrite, float *data,
	       binListStruct_cdf4* binList);
    int close();

    bool has_qual();

    int get_numbin( int irow) { return numbin[irow];}
    int get_basebin( int irow) { return basebin[irow];}

    hid_t get_grp0() { return grp0;}
    hid_t get_index_table() { return 0;}
    hid_t get_list_table() { return 0;}
    hid_t get_data_table( int i) { return 0;}

    int get_bin_num( int kbin) { return binList[kbin].bin_num;}
    int get_nobs( int kbin) { return binList[kbin].nobs;}
    int get_nscenes( int kbin) { return binList[kbin].nscenes;}
    float get_weights( int kbin) { return binList[kbin].weights;}
    float get_time_rec( int kbin) { return binList[kbin].time_rec;}

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
    int inc_time_rec( int offset, float time_rec)
    { binList[offset].time_rec += time_rec; return 0;}

    int clear_binlist() 
    { memset( binList, 0, 2*nrows*sizeof(binListStruct_cdf4)); return 0;}
    int copy_binlist( int src, int dest)
    { memcpy( &binList[dest], &binList[src], sizeof(binListStruct_cdf4));
      return 0;}
  
    int setDataPtr( int nbins_to_read)
        { binDataPtr += nbins_to_read; return 0;}
    int setDataPtrAbsolute(int32 recordNum)
        { binDataPtr = recordNum; return 0; }

    int incNumRec( int n_write)
    { n_data_records += n_write; return 0;}
 
    size_t binDataPtr;
 };

}

#endif



