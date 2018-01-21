#include <stdio.h>
#include <stdlib.h>
#include "grb2.h"
#include "wgrib2.h"
#include "fnlist.h"

extern enum input_type input;
extern int header, dump_rec, dump_submsg;
extern int mode;

extern int file_append;
extern char *item_deliminator;

/*
 * HEADER:100:i:setup:0:read Inventory from stdin
 */
int f_i(ARG0) {

    if (mode == -1) input = inv_mode;
    return 0;
}

/*
 * HEADER:100:v0:misc:0:not verbose (v=0)
 */

int f_v0(ARG0) {
    set_mode(0);
    return 0;
}

/*
 * HEADER:100:v:misc:0:verbose (v=1)
 */

int f_v(ARG0) {
    set_mode(1);
    return 0;
}

/*
 * HEADER:100:v2:misc:0:really verbose (v=2)
 */

int f_v2(ARG0) {
    set_mode(2);
    return 0;
}


/*
 * HEADER:100:header:misc:0:f77 header or nx-ny header in text output
 */

int f_header(ARG0) {
    header = 1;
    return 0;
}

/*
 * HEADER:100:no_header:misc:0:no f77 header or nx-ny header in text output
 */

int f_no_header(ARG0) {
    header = 0;
    return 0;
}

/*
 * HEADER:100:nl:inv:0:inserts new line into inventory
 */

int f_nl(ARG0) {
    if (mode >= 0) sprintf(inv_out, "\n");
    return 0;
}

/*
 * HEADER:100:print:inv:1:inserts string into inventory
 */

int f_print(ARG1) {
    if (mode >= 0) sprintf(inv_out,"%s", arg1);
    return 0;
}

/*
 * HEADER:100:colon:misc:1:replace item deliminator (:) with X
 */

int f_colon(ARG1) {
    item_deliminator = arg1;
    return 0;
}

/*
 * HEADER:100:one_line:misc:0:puts all on one line (makes into inventory format)
 */

int f_one_line(ARG0) {
    nl = " ";
    return 0;
}


/*
 * HEADER:100:append:setup:0:append mode, write to existing output files
 */
int f_append(ARG0) {
    if (mode == -1) file_append = 1;
    return 0;
}
/*
 * HEADER:100:no_append:setup:0:not append mode, write to new output files (default)
 */
int f_no_append(ARG0) {
    if (mode == -1) file_append = 0;
    return 0;
}



