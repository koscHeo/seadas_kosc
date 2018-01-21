#ifndef GETL3B_H /* avoid re-inclusion */
#define GETL3B_H

#include "seaproto.h"

#ifdef	__cplusplus
extern "C" {
#endif

  intn get_fsize(char *prod_type);

  intn get_bins(int32 fid);

  intn rd_common_data(int32 vskey, int32 prod_ID, int32 start, int32 nelts);

  intn  getindex(int32 fst_call, int32 prod_ID, int32 vskey, int32 sbin,
		int32 nrec, int32 *begin, int32 *extent, int32 *begin_index, 
		int32 *sum_extent, int32 *bins_infile, int32 *indexlist, 
		int32 *binno);

  intn  out_data(int32 fst_call, int32 fid, int32 sdfid, int32 prod_ID, 
	l3b_prod *parm_opt, int32 sbin, int32 nrec, int32 binlist_key,  int32 *begin, 
	int32 *extent, int32 *p_vskeys, int32 *last_bin, int32 *binno, 
	int16 *nobs, int16 *time_rec, int16 *nscenes, float32 *weights, 
	int8 *sel_cat, int32 *flags_set, float32 *l3b_data);

  intn  alloc_parm_bufs(int32 prod_ID, int32 nbins);

  intn  free_parm_bufs(int32 prod_ID);

  intn  rd_data(int32 prod_ID, int32 binlist_key, int32 *p_vskey, int32 start, 
                int32 nelts, l3b_prod *parm_opt);

  intn  get_parm_vskeys(int32 fid, int32 prod_ID, int32 *p_vskey);

  intn  attach_slave(int32 fid, char *sname, int32 parm, int32 *p_vskeys); 

  intn  get_index(int32 sbin, int32 prod_ID, int32 bufsz);


#ifdef	__cplusplus
}
#endif


#endif /* GETL3B_H  */

