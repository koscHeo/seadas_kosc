#ifndef SEAPROTO_H
#define SEAPROTO_H

#ifndef METAL3B_H /* avoid re-inclusion */
#include "meta_l3b.h"
#endif

#ifdef	__cplusplus
extern "C" {
#endif

typedef struct l3b_prod_struct {

  char     *l3b_prodlist;
  char     *l3b_units;
  char     *prodname[100];
  int32    code[100];

} l3b_prod;



/*  Prototypes  */
int32 put_l3b_open(char *l3b_path, char *replaces, int32 fsize, 
		   char *prod_type, char *ptime, int32 orbit, 
		   int32 start_orb, int32 end_orb, char *proc_con,
		   char *soft_name, char *soft_ver,    
		   char *input_parms, l3b_prod *parm_opt, 
		   meta_l3bType *meta_l3b);
intn put_l3b_record(int32 file_id, int32 nrec, int32 *binno, int16 *nobs,
                    int16 *time_rec, int16 *nscenes, float32 *weights,
                    int8 *sel_cat, int32 *flags_set, float32 *l3b_data,
                    l3b_prod *parm_opt);
intn put_l3b_close(int32 prod_ID, int16 bin_syear, int16 bin_sday, 
			int16 bin_eyear, int16 bin_eday, int16 syear, 
			int16 sday, int32 smsec, int16 eyear, int16 eday, 
			int32 emsec, char *infiles,char *flag_names,
			char *flag_use, uint8 *eng_q_use,
                        l3b_prod *parm_opt);

intn get_l3b_open(char *l3b_path, int32 prod_ID, int16 *bin_syear, int16 *bin_sday,
                int16 *bin_eyear, int16 *bin_eday,  int16 *syear, int16 *sday,
                int32 *smsec, int16 *eyear, int16 *eday, int32 *emsec,
		char *infiles, int32 *start_orb, int32 *end_orb, 
		char *flag_names, char *flag_use, uint8 *eng_q_use, 
		int32 *fsize, char *prod_type, int32 *nbins, int32 *nrows, 
		int32 *max_row, float32 *bin_hgt, float32 *seam_lon, 
		int32 *ibinr, int32 *nbinr, int32 *irecr, int32 *nrecr, 
		meta_l3bType *meta_l3b, l3b_prod *parm_opt);

intn get_l3b_record(int32 prod_ID, int32 sbin, int32 nrec, l3b_prod *parm_opt, 
                int32 *binno, int16 *nobs, int16 *time_rec, int16 *nscenes,
                float32 *weights, int8 *sel_cat, int32 *flags_set, 
		float32 *l3b_data);

intn get_l3b_close(int32 prod_ID);



int32 get_beg_ext(int32 n_bins_write, int32 *binnum_data,
			     int32 *basebin, int32 nrows,
			     int32 *beg, int32 *ext);

int32 wr_vdata(char *outname, int32 fileid_w, int32 vgid, char *name, 
	       char *class1, int32 n_flds, int32 n_recs_to_write,
	       char *fldname[], int32 type[], int32 noext, uint8 *data,
	       int32 verbose);


/*****************************/
/* prototypes from libbin    */
/*****************************/

void bin_init(int32 nrow, int32 **nbin, int32 **bbin, float **lbin,
			 int32 *tbin);
void bin2ll(int32 bin, float *lat, float *lon);
void ll2bin(float lat, float lon, int32 *bin);
void ll2rc(float lat, float lon, int32 *row, int32 *col);
void rc2ll(int32 row, int32 col, float *lat, float *lon);
void rc2bin(int32 row, int32 col, int32 *bin);
void bin2rc(int32 bin, int32 *row, int32 *col);
void old_bin2ll(int32 bin, float *lat, float *lon);

#ifdef	__cplusplus
}
#endif


#endif /* SEAPROTO_H */




