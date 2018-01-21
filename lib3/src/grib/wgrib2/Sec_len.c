#include <stdio.h>
#include <string.h>
#include "grb2.h"
#include "wgrib2.h"
#include "fnlist.h"

/*
 * length of various sections
 * public domain 2007: Wesley Ebisuzaki
 */

/*
 * HEADER:700:Sec_len:inv:0:length of various grib sections
 */


int f_Sec_len(ARG0) {
    int i;
    if (mode >= 0) {
	sprintf(inv_out,"sec len ");
        inv_out += strlen(inv_out);
        for (i = 1; i <= 7; i++) {
	    if (sec[i] == NULL) {
		sprintf(inv_out, " %d=NULL", i);
	    }
	    else {
		sprintf(inv_out, " %d=%d", i, uint4(sec[i]));
	    }
            inv_out += strlen(inv_out);
	}
    }
    return 0;
}
