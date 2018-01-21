#ifndef ANCPROTO_H     /* avoid re-inclusion */
#define ANCPROTO_H

typedef struct time_change {
	float64 start;
	float64 inc;
} timech;


  intn get_ancillary(float32 *in_lat, float32 *in_lon,  int16 cnt, int16 syear,
                        int16 sday, int16 eday, int32 time, char *filename1, 
			char *filename2, char *filename3, char *anc_cor_file, 
                        int16 parm_flag, float32 *interp, int16  *qcflag);

  int32 set_files (int16 parm_flag, int16 syear, int16 sday, int16 eday, 
	int32 msec, char *filename1, char *filename2, char *filename3, 
	char *file1, char *file2, timech *dtime1, timech *dtime2, intn *toms );

  int32 ck_files_in_buf(int16 PA_flag, int16 parm_flag, char *f1, char *f2,
                        int16 month, int16 *error_flag, int16 *read_flag);

  intn  read_climatology(char *file1, int16 parm_flag, int16 month, 
			int32 *dims, float32 *lat_buf, 
			float32 *lon_buf, void *parm_buf);

  intn check_on_TOMS(int16 parm_flag, char *file1, char *file2,
	  char *filename1, char *filename2, char *filename3, float64 s_jd1, 
	  float64 s_jd2, float64 s_jd3, float64 e_jd1, float64 e_jd2, 
	  float64 e_jd3, float64 d_jd, timech *dt1, timech *dt2 ); 

  intn  read_NRT(char *file1, char *file2, char *anc_cor_file, 
                int16 parm_flag, int16 read_flag, int32 *dims, float32 *lat_buf,
		float32 *lon_buf, void *parm_buf1, void *parm_buf2,
		int8 *qc_buf1, int8 *qc_buf2);

  void resiz_anc( int16 *data, int8 *qc, int32 nlat, int32 nlon1, int32 nlon2,
                float32 * lons1, float32 * lons2 );

  void extract_data_pts(int16 PA_flag, int16 parm_flag, timech DTime1, timech DTime2,
	int16 nsamp, float32 *in_latlist, float32 *in_lonlist, 
	float32 *lat_bufp, float32 *lon_bufp, int32 *dims, int8 *qc_buf1, 
	int8 *qc_buf2, void *parm_buf1, void *parm_buf2, intn toms, 
	float32 *out_lat_list, float32 *out_lon_list, int16 *qcflag, float32 *interp);

  void gregor(int16 jday, int16 year, int16 *month, int16 *day);

  void interpolate(int16 PA_flag, int16 parm_flag, float64 DT1, float64 DT2,
         	   	float32 in_lat, float32 in_lon, float32 *lat_list, 
                   	float32 *lon_list, void *data_p1, void *data_p2,
                    int8 *qc1, int8 *qc2, float32 *intpdata, int32 *int_qc );

  intn get_clim_data (int32 fid, int32 sdfid, int16 parm_flag, 
		     int16 month, int32 *dims, 
		     float32 *lat_buf, float32 *lon_buf, void *parm_buf);

  intn get_NRT_data(int32 fid, int32 sdfid, char *anc_cor_file, 
                    int16 parm_flag, int32 *dims, float32 *lat_buf, 
		    float32 *lon_buf, void *parm_buf, int8 *qc_buf);

  intn openHDF(char *infile, int32 *sdfid, int32 *fid);

  intn rdsds(int32 sdsid, int32 *dimsizes, void *data_buf);

  intn rdlatlon(int32 fid, int32 *dims, float32 *lat_buf, 
                float32 *lon_buf);

  intn get_refs(int32 fid, int32 sdfid, int32 vid, char *parm_label, 
			char *QC_label, int32 *geom_id, int32 *parm_sdsid, 
			int32 *QC_sdsid, int32 *dims); 

  intn get_clim_refs(int32 fid, int32 sdfid, int32 vid, int16 month,
        char *obs_label, int32 *geom_id, int32 *parm_sdsid, int32 *dims);

  int32 anc_get_time(char *filename, float64 *s_jd, float64 *e_jd);

  intn closeHDF(int32 sdfid, int32 fid);

#endif /* ANCPROTO_H */
