#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "grb2.h"
#include "wgrib2.h"
#include "fnlist.h"

/*
 * this file contains nice to know values 
 */

/*
 * HEADER:-1:pds_fcst_time:inv:0:fcst_time(1) in units given by pds
 */
int f_pds_fcst_time(ARG0) {
    unsigned int p;
    if (mode >= 0) {
        p = pds_fcst_time(sec);
	if (p != 0xffffffff) sprintf(inv_out,"pds_fcst_time1=%u", p);
    }
    return 0;
}

unsigned int pds_fcst_time(unsigned char **sec) {
    int p;
    p = GB2_ProdDefTemplateNo(sec);
    if (p <= 10 || p == 1000 || p == 1001 || p == 1002 || p == 1100 || p == 1101)  {
	return uint4(sec[4]+18);
    }
    return 0xffffffff;
}
