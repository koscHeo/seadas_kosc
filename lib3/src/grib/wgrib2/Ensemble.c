#include <stdio.h>
#include <string.h>
#include "grb2.h"
#include "wgrib2.h"
#include "fnlist.h"

/*
 * 2/2007 Public Domain: Wesley Ebisuzaki
 */

/*
 * HEADER:200:ens:inv:0:ensemble information
 */

int f_ens(ARG0) {
    int pdt, n, type, typefcst;
    if (mode >= 0) {
        pdt = code_table_4_0(sec);
	typefcst = code_table_4_7(sec);
        if (pdt == 1 || pdt == 11) {
	    type = code_table_4_6(sec);
	    n = sec[4][35];
	    switch(type) {
	        case 0: sprintf(inv_out,"ens=hi-res ctl"); break;
	        case 1: sprintf(inv_out,"ens=low-res ctl"); break;
	        case 2: sprintf(inv_out,"ens=+%d",n); break;
	        case 3: sprintf(inv_out,"ens=-%d",n); break;
	        default: sprintf(inv_out,"ens=?"); break;
	    }
	    inv_out += strlen(inv_out);
	    if (typefcst >= 0) {
		*inv_out++=' ';
		*inv_out=0;
	    }
	}
	if (typefcst >= 0) {
	    switch(typefcst) {
	        case 0: sprintf(inv_out,"ens-mean"); break;
	        case 1: sprintf(inv_out,"wt ens-mean"); break;
	        case 2: sprintf(inv_out,"std dev"); break;
	        case 3: sprintf(inv_out,"normalized std dev"); break;
	        case 4: sprintf(inv_out,"spread"); break;
	        case 5: sprintf(inv_out,"large anom index"); break;
	        case 6: sprintf(inv_out,"cluster mean"); break;
		default: sprintf(inv_out,"unknown derived fcst"); break;
	    }
	    inv_out += strlen(inv_out);
	}
    }
    return 0;
}

/*
 * HEADER:200:N_ens:inv:0:number of ensemble members
 */
int f_N_ens(ARG0) {
    int pdt, n;
    if (mode >= 0) {
        pdt = code_table_4_0(sec);
	n = -1;
        if (pdt == 2 || pdt == 12 || pdt == 13 || pdt == 14) n = sec[4][35];
	if (n > 0) sprintf(inv_out,"%d ens members", n); 
    }
    return 0;
}
