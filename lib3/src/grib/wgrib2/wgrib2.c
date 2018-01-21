/* wgrib2 main module:  w. ebisuzaki
 *
 * 1/2007 mods M. Schwarb: unsigned int ndata
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <ctype.h>

#include "grb2.h"
#include "wgrib2.h"
#include "fnlist.h"
#include "grib2.h"

/* #define DEBUG */

/* global variables .. can be modified by funtions */

int mode=0;             /* -2=finalize, -1=initialize,  0 .. N is verbosity mode */
int header=1;           /* file header flag */
int flush_mode = 0;	/* flush of output 1 = yes */
#ifdef USE_REGEX
int match = 0;
#endif
int last_message = 0;	/* last message to process if set */



enum input_type input = all_mode;
enum output_order_type output_order = wesn;

char inv_out[INV_BUFFER]; 		/* inv functions write to this buffer */

int only_submsg = 0;    /* if only_submsg > 0 .. only process submsg number (only_submsg) */



/* output file variables */

int file_append = 0;
int dump_msg = 0, dump_submsg = 0;

int decode = 0;		/* decode grib file flag */

char *item_deliminator = ":";
char *nl = "\n\t";

/* current grib location */
long int pos, len;
int submsg, msg_no;

/* lat lon array .. used by Latlon.c to store location information */
float *lat = NULL, *lon = NULL;
int latlon = 0;
int new_GDS = 0;
int old_GDS_size = 0;
int GDS_max_size = 100;
unsigned char *old_gds;
int nx, ny, npnts, res, scan;

/* may disappear in future versions  grib2 decoder variables */
gribfield *grib_data;

/*
 * wgrib2
 *
 * simple wgrib for GRIB2 files
 *
 */


int main(int argc, char **argv) {
    FILE *in;
    unsigned char *msg, *sec[9];
    long int last_pos;

    int file_arg, i, j, num_submsgs;
    int n_arg;
    unsigned int k, ndata;
    float *data, missing_c_val_1, missing_c_val_2;
    g2int *bitmap, has_bitmap;
    g2float *g2_data;

    struct ARGLIST arglist[N_ARGLIST];
    int narglist = 0;
    char *new_argv[N_ARGLIST];
    void *local[N_ARGLIST];
    int has_inv_option, last_submsg;
    int err;

    /* no arguments .. help screen */
    if (argc == 1) {
	f_help(-1,NULL,NULL,0,inv_out,local,"most");
	printf("%s\n", inv_out);
	exit(8);
    }

    /* copy argv */

    for (i = 0; i < argc; i++) {
	new_argv[i] = argv[i];
    }

    /* scan for "inv" and input file */
    has_inv_option = 0;
    file_arg = 0;
    for (i = 1; i < argc; i++) {
	if (new_argv[i][0] != '-') {
	    /* must be filename */
            file_arg = i;
            continue;
        }
	/* must be an option */
	for (j = 0; j < nfunctions; j++) {
	    if (strcmp(&(new_argv[i][1]),functions[j].name) == 0) {
	        if (functions[j].type == inv) has_inv_option = 1;
		i += functions[j].nargs;
                break;
            }
        }
    }

    /* if no inv option, use default inventory .. put in beginning */
    if (has_inv_option == 0) {
        for (i = argc-1; i > 0; i--) {
	    new_argv[i+1] = new_argv[i];
	}
        new_argv[1] = "-s";
	argc += 1;
    } 


    /* parse parameters */
    file_arg = 0;
    for (i = 1; i < argc; i++) {

	if (new_argv[i][0] != '-') {
	    /* must be filename */
	    if (file_arg == 0) {
		file_arg = i;
		continue;
	    } else {
		fatal_error("too many grib files .. 2nd=%s", new_argv[i]);
	    }
	}

	/* must be an option */

	for (j = 0; j < nfunctions; j++) {
	    if (strcmp(&(new_argv[i][1]),functions[j].name) == 0) {
#ifdef DEBUG
		fprintf(stderr,"match .. -%s %d args\n",  functions[j].name, functions[j].nargs);
#endif
                /* add to function argument list */
		arglist[narglist].fn = j;
		arglist[narglist].i_argc = i+1;
		i += functions[j].nargs;
		if (i >= argc) fatal_error("missing arguments option=%s",functions[j].name);
		narglist++;
		break;
	    }
	}

	if (j == nfunctions) {
	    fatal_error("unknown option %s", new_argv[i]);
	}
    }

    /* initialize options mode = -1 */

#ifdef DEBUG
    fprintf(stderr,"init options narglist %d\n",narglist);
#endif

    for (j = 0; j < narglist; j++) {
	inv_out[0] = 0;
	n_arg = functions[arglist[j].fn].nargs;
        err = 0;
        if (n_arg == 0) err = functions[arglist[j].fn].fn(-1,NULL,NULL,0, inv_out,local+j);
	else if (n_arg == 1) err = functions[arglist[j].fn].fn(-1,NULL,NULL,0, inv_out,local+j,
		new_argv[arglist[j].i_argc]);
	else if (n_arg == 2) err = functions[arglist[j].fn].fn(-1,NULL,NULL,0, inv_out,local+j,
		new_argv[arglist[j].i_argc],new_argv[arglist[j].i_argc+1]);
	else if (n_arg == 3) err = functions[arglist[j].fn].fn(-1,NULL,NULL,0, inv_out,local+j,
		new_argv[arglist[j].i_argc],new_argv[arglist[j].i_argc+1],new_argv[arglist[j].i_argc+2]);
	else if (n_arg == 4) err = functions[arglist[j].fn].fn(-1,NULL,NULL,0, inv_out,local+j,
		new_argv[arglist[j].i_argc],new_argv[arglist[j].i_argc+1],
		new_argv[arglist[j].i_argc+2], new_argv[arglist[j].i_argc+3]);


        if(inv_out[0] != 0)  printf("%s", inv_out);
        if (err) exit(8);
    }

    if (file_arg == 0 && argc > 1) fatal_error("no input file", "");
    if (latlon == 1 && output_order != wesn) 
           fatal_error("latitude-longitude information is only available with -order we:sn","");


    /* open input file */
    if ((in = fopen(new_argv[file_arg],"rb")) == NULL) {
        fatal_error("could not open file: %s", new_argv[file_arg]);
    }

    ndata = 0;
    data = NULL;
    msg_no = 1;
    len = pos = 0;
    submsg = 0;
    msg = NULL;

    if ((old_gds = (unsigned char *) malloc(GDS_max_size * sizeof(char)) ) == NULL) {
	fatal_error("memory allocation problem old_gds in wgrib2.main","");
    }
    
    last_pos = -1;
    last_submsg = -1;

    /* if dump mode .. position io stream */

    if (input == dump_mode) {
        while (msg_no < dump_msg) {
            if ((msg = rd_grib2_msg(in, &pos, &len,&num_submsgs)) == NULL) {
                fatal_error_i("record %d not found", dump_msg);
            }
            last_pos = pos;
            pos += len;
            msg_no++;
        }
#ifdef DEBUG
        printf("dump mode msg=%d\n", msg_no);
#endif
    }

    /* 
     * submsg = 0 .. beginning of unread record
     * submsg = i .. start at ith submsg
     * num_submsgs = number of submessages in grib message
     */


    /* inventory loop */ 

    for (;last_message == 0;) {

        /* need position and submessage number of message */
        if (input == inv_mode || input == dump_mode) {
            if (input == inv_mode) {
                if (rd_inventory(&msg_no,&submsg, &pos)) break;
            }
            else if (input == dump_mode) {
                if (dump_msg == -1) break;
                submsg = dump_submsg;
                dump_msg = -1;
	    }

            if (pos != last_pos) {
	        if ((msg = rd_grib2_msg(in, &pos, &len, &num_submsgs)) == NULL) {
                    fatal_error_i("grib message #%d not found", msg_no);
                    break;
                }
                last_pos = pos;
            }

            if (pos == last_pos && submsg == last_submsg + 1) {
                /* read previous submessage */
		if (parse_next_msg(sec) != 0) {
                    fprintf(stderr,"\n*** grib message #%d.%d not found ***\n\n", msg_no, submsg);
                    break;
		}
            }
            else {
                /* need to get desired submessage into sec */
		if (parse_1st_msg(sec) != 0) {
                    fprintf(stderr,"\n*** grib message #%d.1 not found ***\n\n", msg_no);
                    break;
		}
                for (i = 2; i <= submsg; i++) {
		    if (parse_next_msg(sec) != 0) {
                        fprintf(stderr,"\n*** grib message #%d.%d not found ***\n\n", msg_no, i);
                        break;
                    }
		}
	    }
            last_submsg = submsg;
	}
        else if (input == all_mode) {
	    if (submsg == 0) {
	        if ((msg = rd_grib2_msg(in, &pos, &len, &num_submsgs)) == NULL) break;
                submsg = 1;
	    }
	    else if (submsg > num_submsgs) {
		pos += len;
                msg_no++;
	        if ((msg = rd_grib2_msg(in, &pos, &len, &num_submsgs)) == NULL) break;
                submsg = 1;
	    }
            if (submsg == 1) {
		if (parse_1st_msg(sec) != 0) {
		    fprintf(stderr,"illegal format: parsing 1st submessage\n");
		}
            }
            else {
		if (parse_next_msg(sec) != 0) {
                    fprintf(stderr,"illegal format: parsing submessages\n");
                }
	    }
	}
        if (only_submsg > 0 && only_submsg != submsg) {
	    submsg++;
	    continue;
	}

        /* only option */
#ifdef USE_REGEX
        if (match) {
	   inv_out[0] = 0;
	   if (num_submsgs > 1) {
	       sprintf(inv_out,"%d.%d:", msg_no, submsg);
	   }
           else {
	       sprintf(inv_out,"%d:", msg_no);
	   }

           f_s(0, sec, NULL, 0, inv_out+strlen(inv_out), NULL);
           if (is_match(inv_out) != 0) {
              submsg++;
              continue;
           }
        }
#endif

        /* see if new GDS */

	if ((i = GB2_Sec3_size(sec)) != old_GDS_size) {
	    new_GDS = 1;
	}
	else {
	    new_GDS = 0;
	    for (j = 0; j < i; j++) {
		if (old_gds[j] != sec[3][j]) new_GDS = 1;
	    }
	}
	if (new_GDS) {
	    if (i > GDS_max_size) {
		free(old_gds);
		GDS_max_size = i;
    		if ((old_gds = (unsigned char *) malloc(GDS_max_size) ) == NULL) {
			fatal_error("memory allocation problem old_gds in wgrib2.main","");
		}
	    }
	    for (j = 0; j < i; j++) {
		old_gds[j] = sec[3][j];
            }
	    old_GDS_size = i;
	    /* update grid information */
            get_nxny(sec, &nx, &ny, &npnts, &res, &scan);	 /* get nx, ny, and scan mode of grid */
            if (latlon) get_latlon(sec);			 /* get lat lon of grid points */
	}

	if (decode) {
//	    j = code_table_5_0(sec);
//		printf("codetable 5.0 is %d\n",j);

	    /* unpack grib field */
	    err = g2_getfld(msg,submsg,1,1,&grib_data);

            if (err != 0) {
                fprintf(stderr,"Fatal decode err=%d msg %d.%d\n",err, msg_no, submsg);
                exit(8);
            }

	    /* allocate data */
	    if (GB2_Sec3_npts(sec) != ndata) {
		if (ndata) free(data);
	        ndata = GB2_Sec3_npts(sec);
		if (ndata) {
		    data = (float *) malloc(ndata * sizeof(float));
		    if (data == NULL) fatal_error_i("memory allocation failed nbytes=%d",(int)ndata*sizeof(float));
		}
		else { data = NULL; }
	    }

	    has_bitmap = grib_data->ibmap;
	    g2_data = &(grib_data->fld[0]);
	    if (has_bitmap == 0 || has_bitmap == 254) {
		bitmap = grib_data->bmap;
		for (k = 0; k < ndata; k++) {
		     data[k] = (bitmap[k] == 0) ? UNDEFINED : g2_data[k];
		}
	    }
	    else {
		for (k = 0; k < ndata; k++) {
		    data[k] = *g2_data++;
		}
	    }

	    /* complex packing uses special values for undefined */
	    i = code_table_5_5(sec);
	    if (i == 1) {
		missing_c_val_1 = ieee2flt(sec[5]+23);
		for (k = 0; k < ndata; k++) {
		    if (data[k] == missing_c_val_1) data[k] = UNDEFINED;
		}
	    }
  	    else if (i == 2) {
		missing_c_val_1 = ieee2flt(sec[5]+23);
		missing_c_val_2 = ieee2flt(sec[5]+27);
		for (k = 0; k < ndata; k++) {
		    if (data[k] == missing_c_val_1) data[k] = UNDEFINED;
		    if (data[k] == missing_c_val_2) data[k] = UNDEFINED;
		}
	    }

	    /* convert to standard output order we:sn */
	    if (output_order == wesn) to_we_sn_scan(data);
	    else if (output_order == wens) to_we_ns_scan(data);
	}
        else {
	    if (ndata) free(data);
            ndata = 0;
            data = NULL;
        }

	if (num_submsgs > 1) {
	    printf("%d.%d%s%ld", msg_no, submsg, ":", pos);
	}
        else {
	    printf("%d%s%ld", msg_no, ":", pos);
	}

	for (j = 0; j < narglist; j++) {
            if (functions[arglist[j].fn].type == inv) printf("%s", item_deliminator);
            if (functions[arglist[j].fn].type != setup) {
		inv_out[0] = 0;
	        n_arg = functions[arglist[j].fn].nargs;
		if (n_arg == 0) 
                    functions[arglist[j].fn].fn(mode, sec, data, ndata, inv_out, local+j);
		else if (n_arg == 1)
		    functions[arglist[j].fn].fn(mode, sec, data, ndata, inv_out, local+j,
			 new_argv[arglist[j].i_argc]);
		else if (n_arg == 2)
		    functions[arglist[j].fn].fn(mode, sec, data, ndata, inv_out, local+j,
			new_argv[arglist[j].i_argc], new_argv[arglist[j].i_argc+1]);
		else if (n_arg == 3)
		    functions[arglist[j].fn].fn(mode, sec, data, ndata, inv_out, local+j,
			new_argv[arglist[j].i_argc], new_argv[arglist[j].i_argc+1],
			new_argv[arglist[j].i_argc+2]);
		else if (n_arg == 4)
		    functions[arglist[j].fn].fn(mode, sec, data, ndata, inv_out, local+j,
			new_argv[arglist[j].i_argc], new_argv[arglist[j].i_argc+1],
			new_argv[arglist[j].i_argc+2], new_argv[arglist[j].i_argc+3]);

                if (functions[arglist[j].fn].type == inv) printf("%s", inv_out);
           }
	}

	submsg++;

	if (decode) g2_free(grib_data);

	printf("\n");
	if (flush_mode) fflush(stdout);
	if (dump_msg > 0) break;
    }

    /* finalize all functions, call with mode = -2 */

    for (j = 0; j < narglist; j++) {
        if (functions[arglist[j].fn].type != setup) {
	    n_arg = functions[arglist[j].fn].nargs;
	    if (n_arg == 0) 
                functions[arglist[j].fn].fn(-2, NULL, NULL, 0, inv_out, local+j);
	    else if (n_arg == 1)
		functions[arglist[j].fn].fn(-2, NULL, NULL, 0, inv_out, local+j,
			new_argv[arglist[j].i_argc]);
	    else if (n_arg == 2)
		functions[arglist[j].fn].fn(-2, NULL, NULL, 0, inv_out, local+j,
			new_argv[arglist[j].i_argc], new_argv[arglist[j].i_argc+1]);
	    else if (n_arg == 3)
		functions[arglist[j].fn].fn(-2, NULL, NULL, 0, inv_out, local+j,
			new_argv[arglist[j].i_argc], new_argv[arglist[j].i_argc+1],
			new_argv[arglist[j].i_argc+2]);
	    else if (n_arg == 4)
		functions[arglist[j].fn].fn(-2, NULL, NULL, 0, inv_out, local+j,
			new_argv[arglist[j].i_argc], new_argv[arglist[j].i_argc+1],
			new_argv[arglist[j].i_argc+2], new_argv[arglist[j].i_argc+3]);
        }
    }
    exit(0);
}

/*
 * reads inventory and pulls out the address of the record
 * and submessage number
 */

int rd_inventory(int *rec_num, int *submsg, long int *pos) {

    long int tmp;
    int c, i;

    /* format: [digits].[submessage]:[pos]: or digits]:[pos]: */

    i = 0;

    c=getchar();
    while (c == ' ') c = getchar();
    if (c == EOF) return 1;
    if (!isdigit(c)) {
	fatal_error_i("bad inventory on line: %d",*rec_num);
    }

    /* record number */
    while (isdigit(c) ) {
	i = 10*i + c - '0';
        c=getchar();
    }

    *rec_num = i;

    if (c == '.') {
        i = 0;
	while (isdigit(c = getchar()) ) {
	    i = 10*i + c - '0';
	}
	*submsg = i;
    }
    else {
        *submsg = 1;
    }

    if (c != ':') fatal_error_i("bad inventory on line: %d",*rec_num);

    tmp = 0;
    while (isdigit(c = getchar()) ) {
        tmp = 10*tmp + c - '0';
    }

    /* if (c != ':') fatal_error_i("bad inventory on line: %d",*rec_num); */

    /* read the rest of the line */
    while (c != EOF && c != '\n') {
	c = getchar();
    }
    *pos = tmp;
    return 0;
}

void set_mode(int new_mode) {
	mode = new_mode;
}
