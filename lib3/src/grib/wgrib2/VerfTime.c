#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "grb2.h"
#include "wgrib2.h"
#include "fnlist.h"

static char *months = "janfebmaraprmayjunjulaugsepoctnovdec";

/*
 * HEADER:400:vt:inv:0:verf time = reference_time + forecast_time, -v2 for alt format
 * 9/2006  w.ebisuzaki
 * 1/2007  check error code on verftime
 */
int f_vt(ARG0) {

    int year, month, day, hour, minute, second;

    if (mode >= 0) {
        if (verftime(sec, &year, &month, &day, &hour, &minute, &second) == 0) {
	    if (mode != 2) {
	        sprintf(inv_out,"%4.4d%2.2d%2.2d%2.2d", year,month,day,hour);
	    }
	    else {
               sprintf(inv_out,"vt=%2.2dZ%2.2d%c%c%c%4.4d", hour,day,months[month*3-3],
		months[month*3-2], months[month*3-1], year);
	    }
        }
        else {
           sprintf(inv_out,"vt=?");
        }
    }
    return 0;
}

/*
 * HEADER:400:VT:inv:0:verf time = reference_time + forecast_time (YYYYMMDDHHMMSS)
 */
int f_VT(ARG0) {

    int year, month, day, hour, minute, second;

    if (mode >= 0) {
        if (verftime(sec, &year, &month, &day, &hour, &minute, &second) == 0) {
            sprintf(inv_out,"vt=%4.4d%2.2d%2.2d%2.2d%2.2d%2.2d", year,month,day,hour,minute,second);
        }
        else {
            sprintf(inv_out,"vt=?");
       }
    }
    return 0;
}


/*
           Returns the verification time: reference_time + forecast_time
 * 9/2006  w. ebisuzaki
 * 1/2007  w. ebisuzaki return error code
 * 11/2007  w. ebisuzaki fixed code for non forecasts
 */

int verftime(unsigned char **sec, int *year, int *month, int *day, int *hour, int *minute, int *second) {

    unsigned char *p;
    int units;
    unsigned int dtime;

    p = sec[1];
    *year = (p[12] << 8) | p[13];
    *month = p[14];
    *day = p[15];
    *hour = p[16];
    *minute = p[17];
    *second = p[18];

    units = code_table_4_4(sec);
    if (units == -1) return 1;

    dtime = pds_fcst_time(sec);

    if (dtime == 0xffffffff) return 0;

    return add_time(year, month, day, hour, minute, second, dtime, units);
}


