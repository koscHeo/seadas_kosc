#ifndef CALPROTOI_H /* avoid re-inclusion */
#define CALPROTOI_H

#include "mfhdf.h"

int32 get_cal_osmi(char *cal_path, int16 syear, int16 sday, int16 eday, int32 msec, 
		   int16 *cal_year, int16 *cal_day, int32 *cal_msec,
		   float32 *eoffset, float32 *egain, float32 *temp, 
		   float32 *mirror, float32 *t_const,
		   float32 *t_linear, float32 *t_quadratic, 
		   float32 *slope, float32 *dc, 
		   float32 *sm);


int32 get_index_osmi(int32 fid, int16 syear, int16 sday, int16 eday, int32 msec,
		     int16 *cal_year, int16 *cal_day, int32 *cal_msec);

int32 read_parm_data_osmi(int32 fid, int32 sdfid, int32 index, float32 *eoffset, 
			  float32 *egain, float32 *temp, float32 *mirror,  
			  float32 *t_const, float32 *t_linear, float32 *t_quadratic, 
			  float32 *slope, float32 *dc, float32 *sm);

#endif  /* CALPROTOI_H */
