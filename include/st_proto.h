#ifndef ST_PROTO_H
#define ST_PROTO_H

  int32 stray_light_gac(int32 *initial, float32 Ltyp_frac, float32 Styp_frac,
		int32 nscans, int32 nsamples, int32 scan_no, int16 gn, 
		float32 *rads, float32 *l1b_data, int32 *sl_scan, 
		int32 *sl_flag);

  int32 stray_light_lac(int32 *initial, float32 Ltyp_frac, float32 Styp_frac,
		int32 nscans, int32 nsamples, int32 scan_no, int16 gn, 
		float32 *rads, float32 *l1b_data, int32 *sl_scan, 
		int32 *sl_flag);

  int32 stray_light_corr(int32 *initial, float32 Ltyp_frac, float32 Styp_frac,
		int32 nscans, int32 nsamples, int32 scan_no, char *dtype, 
		int16 gn, float32 *rads, float32 *l1b_data, int32 *sl_scan, 
		int16 *l2_flags, int32 *AS_pixels);

#endif  /* ST_PROTO_H */

