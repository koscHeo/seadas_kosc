#ifndef SEAPROTOI_H  /* avoid re-inclusion */
#define SEAPROTOI_H

/*  Prototypes  */
int init_bins( int nrows, float seam_lon ); 
void binxf_bin2ll( int idx, float *lat, float *lon, float *dellat, 
		float *dellon);
int32 get_coords(int32 *begin, int32 *lastdatabin, int32 lastrec, float32 *nlat,
                        float32 *slat, float32 *elon, float32 *wlon);
int32 setupgrid(int32 fileid, int32 fid);
int32 closegrid(int32 fileid, int32 gridid);
int32 writegeom(int32 fid, int32 nbins);
int32 setupindex(int32 fid);
intn writeindex(int32 ndxid, int32 *start_num, int32 *begin, int32 *extent, 
		int32 *maxbin);
intn closeindex(int32 ndxid);
int32 setupmaster(int32 fid, int32 fsize);
int32 *setupslaves(int32 fid, int32 fsize, char *l3b_path, l3b_prod *parm_opt);
intn closedata(int32 fileid, int32 mstrid, int32 *slvid, l3b_prod *parm_opt);
intn buffbins(int32 *fstcall, int32 fileid, int32 mstrid, int32 *slvid, 
	      int32 nrec, int32 *binno, int16 *nobs, int16 *ttag, 
	      int16 *nscenes, float32 *weights, int8 *sel_cat, 
	      int32 *flags_set, float32 *l3b_data, l3b_prod *parm_opt);
intn flushbuff(int32 fileid, int32 mstrid, int32 *slvid, l3b_prod *parm_opt);
int32 writeattr(int32 fid, char *name, int32 nt, intn count, VOIDP data);

#endif /* SEAPROTOI_H */
