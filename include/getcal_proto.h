#ifndef CALPROTOI_H /* avoid re-inclusion */
#define CALPROTOI_H

#include "l1a.h"
#include "mfhdf.h"

  void read_caltable(char *cal_path);

  int32_t l1b_rad (int syear, int sday, int32_t smsec,int32_t msec,
                    char *dtype, int32_t nsta, int32_t ninc, int32_t npix,
                    float *dark_mean, short *gain, short *tdi,
                    short *scan_temp, float *inst_temp, int mside,
                    short *l1a_data, float *l1b_data, cal_mod_struc *cal_mod);

  int32 get_cal(char *cal_path, int16 syear, int16 sday, int16 eday,
	      int32 msec, int32 npix, int32 nsta, int32 ninc,
	      char *dtype, int16 *tdi, int16 *cal_year, int16 *cal_day,
              int16 *ref_year, int16 *ref_day, int16 *ref_min,
	      float32 fp_temps[256][8], float32 scan_mod[2][1285],
              float64 *tfactor_const, float64 *tfactor_linear_1,
              float64 *tfactor_exponential_1, float64 *tfactor_linear_2,
              float64 *tfactor_exponential_2,float64 *cal_offset,
	      float64 *inst_tcorr, float64 *inst_tref, float64 *fp_tcorr,
	      float64 *fp_terf, float64 *mside1_const, 
	      float64 *mside1_linear_1, float64 *mside1_exponential_1, 
              float64 *mside1_linear_2, float64 *mside1_exponential_2, 
              float64 *mside2_const, float64 *mside2_linear_1, 
              float64 *mside2_exponential_1, float64 *mside2_linear_2, 
              float64 *mside2_exponential_2, float32 counts[8][4][5], 
              float32 rads[8][4][5]); 

  int32 read_parm_data(int32 fid, int32 sdfid, int32 index, int32 idoffs[8][16],
                     float32 gains[8][16], float32 fp_temps[256][8],
                     float32 scan_mod[2][1285], float64 *tfactor_const, 
                     float64 *tfactor_linear_1, float64 *tfactor_exponential_1,
                     float64 *tfactor_linear_2, float64 *tfactor_exponential_2,
		     float64 *cal_offset, float64 *inst_tcorr, 
		     float64 *inst_tref, float64 *fp_tcorr, float64 *fp_tref,
                     float64 *mside1_const, float64 *mside1_linear_1, 
                     float64 *mside1_exponential_1, float64 *mside1_linear_2, 
                     float64 *mside1_exponential_2, float64 *mside2_const, 
                     float64 *mside2_linear_1, float64 *mside2_exponential_1, 
                     float64 *mside2_linear_2, float64 *mside2_exponential_2, 
                     int16 tdi_list[256][4]);

  void calc_knees(int16 *tdi, int16 tdi_list[256][4], int32 idoffs[8][16],
                	float32 gains[8][16], float32 counts[8][4][5],
                	float32 rads[8][4][5]);

int32 read_time_data(int32 fid, int16 *year, int16 *day, int32 *msec, 
			int32 *elts);

int32 get_index(int32 fid, int16 syear, int16 sday, int16 eday, int32 msec,
			int16 *cal_year, int16 *cal_day);

void setup_scanmod(int32 npix, int32 nsta, int32 ninc,
		   float32 scan_mod[2][1285]);

void sort_srads(float32 *srads, int32 *oindex);

int32 get_ref_time(int32 sdfid,
			int16 *ref_year, int16 *ref_day, int16 *ref_min);


int32 get_tindex(int32 fid, int16 syear, int16 sday, int16 eday, int32 msec,
		int16 *cal_year, int16 *cal_day);


#endif  /* CALPROTOI_H */
