#ifndef CALL1APROTO_OSMI_H
#define CALL1APROTO_OSMI_H

int32 calibrate_l1a_osmi(char     *cal_path, 
                     int16    syear, 
                     int16    sday, 
                     int32    smsec, 
		     int16    eday,  
                     int32    msec, 
                     char     *dtype, 		/* "","SOL","TDI","IGC","LUN" */
                     int32    st_samp, 		/* 1 */
                     int32    nsamp, 		/* 1044 */
                     int32    fpixel,		/* 0-95 */
                     int16    gain[4], 		/* quadrant gain settings from modified ISD */
                     int16    offset,		/* quadrant offset settings from mod ISD */
                     int16    scan_temp, 	/* sensor temps from mod ISD */
                     int16    *l1a_data, 
		     float32  *l1b_data, 
                     cal_mod_struc *cal_mod);

int jul2jday(int year, int yday);

int cal2jday(int year, int month, int mday);

#endif /* CALL1APROTO_H */
