#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "grb2.h"
#include "wgrib2.h"
#include "fnlist.h"

extern char *item_deliminator;

/*
 * HEADER:100:s:inv:0:simple inventory
 */

/*
 * this is a simple macro .. see how easy it is!
 * would be more complicated if functions used static variables
 * minor complication if need to set decode or latlon flags
 */

int f_s(ARG0) {

    if (mode >= 0) {
	f_t(CALL_ARG0);
	strcat(inv_out,item_deliminator);
	inv_out += strlen(inv_out);

	f_var(CALL_ARG0);
	strcat(inv_out,item_deliminator);
	inv_out += strlen(inv_out);

	f_lev(CALL_ARG0);
	strcat(inv_out,item_deliminator);
	inv_out += strlen(inv_out);

	f_ftime(CALL_ARG0);
	strcat(inv_out,item_deliminator);
	inv_out += strlen(inv_out);

	f_ens(CALL_ARG0);
	inv_out += strlen(inv_out);
    }
    return 0;
}

/*
 * HEADER:100:verf:inv:0:simple inventory using verification time
 */
int f_verf(ARG0) {

    if (mode >= 0) {
        f_vt(CALL_ARG0);
        strcat(inv_out,item_deliminator);
        inv_out += strlen(inv_out);

        f_var(CALL_ARG0);
        strcat(inv_out,item_deliminator);
        inv_out += strlen(inv_out);

        f_lev(CALL_ARG0);
        strcat(inv_out,item_deliminator);
        inv_out += strlen(inv_out);

        f_ftime(CALL_ARG0);
        inv_out += strlen(inv_out);
        strcat(inv_out,item_deliminator);

	f_ens(CALL_ARG0);
	inv_out += strlen(inv_out);
    }
    return 0;
}

