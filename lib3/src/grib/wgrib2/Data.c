#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "grb2.h"
#include "wgrib2.h"
#include "fnlist.h"

/*
 * Data.c
 *
 *  Some routines that examine the data
 *
 * 2006: Public Domain: Wesley Ebisuzaki
 * 1/2007 cleanup M. Schwarb
 *
 */



extern int decode, latlon;
extern int file_append;
extern float *lat, *lon;

/*
 * HEADER:100:stats:inv:0:statistical summary of data values
 */

int f_stats(ARG0) {
    double sum, min, max;
    int first, n;
    unsigned int i;

    if (mode == -1) {
        decode = 1;
    }
    else if (mode >= 0) {
	first = 1;
	n = 0;
	sum =UNDEFINED;
	min = max = UNDEFINED;
	for (i = 0; i < ndata; i++) {
	    if (!UNDEFINED_VAL(data[i])) {
		n++;
		if (first) {
	 	    max = min = data[i];
		    first = 0;
		    sum = data[i];
		}
		else {
		    max = max < data[i]? data[i] : max;
		    min = min > data[i]? data[i] : min;
		    sum += data[i];
		}
	    }
	}
	if (n) sum /= n;
        sprintf(inv_out,"ndata=%u:undef=%u:mean=%lg:min=%lg:max=%lg", ndata, ndata-n, sum, min, max);
    }
    return 0;
}
/*
 * HEADER:100:max:inv:0:print maximum value
 */

int f_max(ARG0) {
    double mx;
    int ok;
    unsigned int i;

    if (mode == -1) {
        decode = 1;
    }
    else if (mode >= 0) {
        mx = 0.0;
        ok = 0;
        for (i = 0; i < ndata; i++) {
	    if (!UNDEFINED_VAL(data[i])) {
                if (ok) {
		    mx = mx > data[i] ? mx : data[i];
                }
                else {
                    ok = 1;
                    mx = data[i];
               }
	    }
	}
        if (ok) sprintf(inv_out,"max=%lg",mx);
	else sprintf(inv_out,"max=undefined");

    }
    return 0;
}

/*
 * HEADER:100:min:inv:0:print minimum value
 */

int f_min(ARG0) {
    double mx;
    int ok;
    unsigned int i;

    if (mode == -1) {
        decode = 1;
    }
    else if (mode >= 0) {
        mx = 0.0;
        ok = 0;
        for (i = 0; i < ndata; i++) {
            if (!UNDEFINED_VAL(data[i])) {
                if (ok) {
                    mx = mx < data[i] ? mx : data[i];
                }
                else {
                    ok = 1;
                    mx = data[i];
               }
            }
        }
        if (ok) sprintf(inv_out,"min=%lg",mx);
        else sprintf(inv_out,"min=undefined");

    }
    return 0;
}
