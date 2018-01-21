#include <stdio.h>
#include "grb2.h"
#include "wgrib2.h"
#include "fnlist.h"

static char *months = "janfebmaraprmayjunjulaugsepoctnovdec";

/*
 * HEADER:100:t:inv:0:time YYYYMMDDHH, -v2 for alt format
 */

int f_t(ARG0) {

    unsigned char *p;
    char *m;
    if (mode >= 0) {
	p = sec[1];
	if (mode != 2) {
	   sprintf(inv_out,"d=%4.4d%2.2d%2.2d%2.2d", (p[12]<<8)+p[13], p[14],p[15],p[16]);
	}
	else {
	   m = months + (p[14]-1)*3;
	   sprintf(inv_out,"%2.2dZ%2.2d%c%c%c%4.4d", p[16],p[15],m[0],m[1],m[2],(p[12]<<8)+p[13]);
	}
    }
    return 0;
}

/*
 * HEADER:101:T:inv:0:time YYYYMMDDHHMMSS
 */

int f_T(ARG0) {
    unsigned char *p;
    if (mode >= 0) {
	p = sec[1];
	sprintf(inv_out,"D=%4.4d%2.2d%2.2d%2.2d%2.2d%2.2d", (p[12]<<8)+p[13], 
	p[14],p[15],p[16],p[17],p[18]);
    }
    return 0;
}
/*
 * HEADER:103:YY:inv:0:year
 */
 
int f_YY(ARG0) {
    int year;
    if (mode >= 0) {
        year = (sec[1][12]<<8)+sec[1][13];
        sprintf(inv_out, mode == 0 ? "YY=%4.4d" : "%4.4d", year);
    }
    return 0;
}

/*
 * HEADER:105:MM:inv:0:month
 */
 
int f_MM(ARG0) {
    if (mode >= 0) {
        sprintf(inv_out,"MM=%2.2d",sec[1][14]);
    }
    return 0;
} 

/*
 * HEADER:110:RT:inv:0:Reference Time
 */
 
int f_RT(ARG0) {
    if (mode >= 0) {
	switch (sec[1][11]) {
	case 0: sprintf(inv_out,"RT=analysis");
		break;
	case 1: sprintf(inv_out,"RT=Start of fcst");
		break;
	case 2: sprintf(inv_out,"RT=Verf time of fcst");
		break;
	case 3: sprintf(inv_out,"RT=Obs time");
		break;
	default:
		break;
	}
    }
    return 0;
}

/*
 * HEADER:110:center:inv:0:center 
 */

int f_center(ARG0) {
    int ctr;
    char *string, tmp[20];
    
    if (mode >= 0) {
        ctr = GB2_Center(sec);
        switch (ctr) {
#include "code_table0.dat"
	}
	sprintf(inv_out,"center=%s", string);
    }
    return 0;
}

/*
 * HEADER:110:subcenter:inv:0:subcenter 
 */

int f_subcenter(ARG0) {
    int ctr, subctr;
    char *string, tmp[20];
    
    if (mode >= 0) {
        ctr = GB2_Center(sec);
        subctr = GB2_Subcenter(sec);
	string = NULL;
	if (ctr == 7) {
           switch (subctr) {
#include "ncep_tableC.dat"
	   }
	}
	if (string == NULL) sprintf(inv_out,"subcenter=%d", subctr);
	else sprintf(inv_out,"subcenter=%s", string);
    }
    return 0;
}
