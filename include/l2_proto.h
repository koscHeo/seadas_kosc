#ifndef L2_PROTO_H_
#define L2_PROTO_H_

extern int get_l2_open
	PROTO((char *l2_path,
	short int *syear,short int *sday,int *smsec,
	short int *eyear,short int *eday,int *emsec,
	int *orbit,char **dtype,int *nrec,int *nsamp, char **flag_names,
	char **mskflg,
	int *ntilts, short int **tilt_flags, short int **tilt_ranges,
	float **tilt_lats, float **tilt_lons, meta_l2Type *meta_l2));

#ifdef GET_ALL_REC_L2
extern int get_l2_record
	PROTO((int prod_ID, int recno, int *msec, byte **eng_qual,
	byte **s_flags, float **fl2_data, short int **l2_flags,
	float **geoloc, int *nav_flag, rec_l2Type *rec_l2));

#endif /* GET_ALL_REC_L2 */

#ifndef GET_ALL_REC_L2
/* -------------------------------------------------------------------------- */
extern int get_l2_record
	PROTO((int prod_ID, int recno, int *msec, byte **eng_qual,
	byte **s_flags, float **fl2_data, short int **l2_flags,
	float **geoloc, int *nav_flag, browse_rec_l2Type *brec_l2));

/* -------------------------------------------------------------------------- */
#endif /* !GET_ALL_REC_L2 */

extern int get_l2_close
	PROTO((int prod_ID));

extern int set_l2get_buffer_blksize
	PROTO((int prod_ID,int next_blk_size));

extern int toggle_l2get_alt_counts_method
	PROTO((int prod_ID));


#endif /* L2_PROTO_H_ */
