#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "grib2.h"
#include "grb2.h"
#include "wgrib2.h"
#include "fnlist.h"

/*
 * Grib_out
 *
 * routines to encode data into grib2
 *
 * 12/2007 Public Domain by Wesley Ebisuzaki
 *
 */

unsigned char *mk_bms(float *data, int *ndata);
int mk_sec5_7(float *data, int n, 
	unsigned char **sec5, unsigned char **sec7);
int grib_out(unsigned char *sec0, unsigned char *sec1, unsigned char *sec2,
     unsigned char *sec3, unsigned char *sec4, float *data, int ndata, FILE *out);


extern int decode, nx, ny, scan;
extern int flush_mode;
extern enum output_order_type output_order;

/*
 * HEADER:100:grib_out:output:1:writes decoded/modified data in grib-2 format to file X
 */

int f_grib_out(ARG1) {

    int i, j, k;
    float *data_tmp;

    if (mode == -1) {
        decode = 1;
        if ((*local = (void *) fopen(arg1, "wb")) == NULL) 
		fatal_error("Could not open %s", arg1);
    }
    else if (mode == -2) {
	if (fclose((FILE *)*local) != 0) fatal_error("Failed in closing file: %s", arg1);
    }
    else if (mode >= 0) {
	if ((data_tmp = (float *) malloc(ndata * sizeof(float))) == NULL)
	 	fatal_error("memory allocation - data_tmp","");
	// copy data to data_tmp in original scan order!
	if (output_order == raw) {
	    for (i = 0; i < ndata; i++) data_tmp[i] = data[i];
	}
	else if (output_order == wesn) {
	    for (j = 0; j < ny; j++) {
	        for (i = 0; i < nx; i++) {
		    k = ij2p(i, j, scan);
		    data_tmp[k] = data[i+j*nx];
	        }
	    }
	}
	else fatal_error("grib_out does not work on this ordering","");

        grib_out(sec[0], sec[1], sec[2], sec[3], sec[4], data_tmp, ndata, (FILE *)*local);
        if (flush_mode) fflush((FILE *) *local);
	free(data_tmp);
    }
    return 0;
}

/*
 * write a grib-2 file 
 *
 * sec0..sec4 predefined sections 0 to 4
 * data[] = values to encode into grib
 * ndata = size of data
 * out = output file
 *
 * currently this program writes in simple packing 12 bits per data point
 */

int grib_out(unsigned char *sec0, unsigned char *sec1, unsigned char *sec2, 
     unsigned char *sec3, unsigned char *sec4, float *data, int ndata, FILE *out) {

    long int size;
    int n_defined;
    unsigned char s[8];
    unsigned char *sec5, *sec6, *sec7;


    /* make a new data section (5) */
    n_defined = ndata;

    sec6 = mk_bms(data, &n_defined);			// make bitmap section
    mk_sec5_7(data, n_defined, &sec5, &sec7);		// make sec 5 and 7

    size = (long int) GB2_Sec0_size + GB2_Sec8_size +
         (sec1 ? uint4(sec1) : 0) +
         (sec2 ? uint4(sec2) : 0) +
         (sec3 ? uint4(sec3) : 0) +
         (sec4 ? uint4(sec4) : 0) +
         (sec5 ? uint4(sec5) : 0) +
         (sec6 ? uint4(sec6) : 0) +
         (sec7 ? uint4(sec7) : 0);

    /* section 0 */
    fwrite((void *) sec0, sizeof(char), 8, out);
    uint8_char(size, s);
    fwrite((void *) s, sizeof(char), 8, out);

    if (sec1) fwrite((void *)sec1, sizeof(char), uint4(sec1), out);
    if (sec2) fwrite((void *)sec2, sizeof(char), uint4(sec2), out);
    if (sec3) fwrite((void *)sec3, sizeof(char), uint4(sec3), out);
    if (sec4) fwrite((void *)sec4, sizeof(char), uint4(sec4), out);
    if (sec5) fwrite((void *)sec5, sizeof(char), uint4(sec5), out);
    if (sec6) fwrite((void *)sec6, sizeof(char), uint4(sec6), out);
    if (sec7) fwrite((void *)sec7, sizeof(char), uint4(sec7), out);
    
    fwrite((void *) "7777", sizeof(char), 4, out);

    free(sec5);
    free(sec6);
    free(sec7);
    return 0;
}

/*
 * based on mk_BMS (from gribw)
 *
 * public domain 12/2007 Wesley Ebisuzaki
 */

unsigned char *mk_bms(float *data, int *ndata) {

    int nn, bms_size, c, start, imask, i;
    unsigned char *bms, *cbits;

    nn = *ndata;
    bms_size = 6 + (nn+7) / 8;
    bms = (unsigned char *) malloc(bms_size);

    uint_char(bms_size, bms);		// length of section 6
    bms[4] = 6;				// section 6
    bms[5] = 0;				// has bitmap

    cbits = bms + 6;
    c = start = 0;
    imask = 128;
    for (i = 0; i < nn; i++) {
	if (data[i] < UNDEFINED_LOW || data[i] > UNDEFINED_HIGH) {
	    c += imask;
	    data[start++] = data[i];
	}
	if ((imask >>= 1) == 0) {
	    *cbits++ = c;
	    c = 0;
	    imask = 128;
	}
    }
    if (imask != 128) *cbits = c;
    *ndata = start;
    if (nn == start) {
	free (bms);
	bms = (unsigned char *) malloc(6);
	uint_char(6, bms);		// length of section 6
	bms[4] = 6;			// section 5
	bms[5] = 255;			// no bitmap
    }
    return bms;
}


static int max_bits = 12;

/*
 * HEADER:100:grib_bits:misc:1:X set number of binary bits for grib_out packing
 */

int f_grib_bits(ARG1) {
    int i;
    char *s;

    s = arg1;
    i = 0;
    while (*s >= '0' && *s <= '9') {
	i = i*10 + *s++ - '0';
    }
    if (i != 0) max_bits = i;
    return 0;
}

/*
 * make sec 5 and 7 using simple packing
 */


int mk_sec5_7(float *data, int n, unsigned char **sec5, unsigned char **sec7) {

    float min_val, max_val;
    int i, nbits, bin_scale;
    double frange, scale, fmin;
    unsigned long int sec5_size, sec7_size;
    unsigned char *p;

    for (max_val = min_val = data[0], i = 1; i < n; i++) {
	if (min_val > data[i]) min_val = data[i];
        if (max_val < data[i]) max_val = data[i];
    }
    /* ecmwf style */
    fmin = min_val;
    frange = max_val - fmin;
    if (frange != 0.0) {
	frexp(frange, &i);
	bin_scale = i - max_bits;
	nbits = max_bits;
	scale = ldexp(1.0, -bin_scale);
	frange = floor((max_val-fmin)*scale + 0.5);
	frexp(frange, &i);
	if (i != nbits) bin_scale++;
    }
    else {
	bin_scale = nbits = 0;
    }
    if (bin_scale) {
        scale = ldexp(1.0, -bin_scale);
	for (i = 0; i < n; i++) {
	    data[i] = (data[i] - fmin)*scale + 0.5;
	}
    }
    else {
	for (i = 0; i < n; i++) {
	    data[i] = data[i] - fmin + 0.5;
	}
    }

    sec5_size = 21;
    sec7_size = 5 + (nbits * (n / 8)) + (nbits * (n % 8) + 7) / 8;

    // section 7
    *sec7 = p = (unsigned char *) malloc(sec7_size);
    uint_char(sec7_size, p);
    p[4] = 7;
    flist2bitstream(data,p + 5,n,nbits);

    // section 5
    *sec5 = p = (unsigned char *) malloc(sec5_size);
    uint_char(sec5_size, p);		// length of section 5
    p[4] = 5;				// section 5
    uint_char(n, p+5);			// number of defined points
    uint2_char(0,p+9);			// template 5.0
    flt2ieee(fmin,p+11);		// ieee reference value
    int2_char(bin_scale,p+15);
    int2_char(0,p+17);
    p[19] = nbits;
    p[20] = 0;				// template 5.1 - set to floating

    return 0;
}



/*
 * grib: convert linear list of ints to a bitstream
 *
 * in public domain : 2007 Wesley Ebisuzaki
 */


static unsigned int mask[] = {0,1,3,7,15,31,63,127,255};

void flist2bitstream(float *list, unsigned char *bitstream, int ndata, int nbits) 
{

    int cbits, jbits;
    unsigned int j, c;

    if (nbits == 0) {
	return;
    }
    if (nbits < 0) {
	fprintf(stderr,"nbits < 0!  nbits = %d\n", nbits);
	exit(0);
    }

    cbits = 8;
    c = 0;
    while (ndata-- > 0) {
	/* note float -> unsigned int .. truncate */
        j = (unsigned int) *list++;
	jbits = nbits;
	while (cbits <= jbits) {
	    if (cbits == 8) {
	        jbits -= 8;
	        *bitstream++ = (j >> jbits) & 255;
	    }
	    else {
	        jbits -= cbits;
	        *bitstream++ = (c << cbits) + ((j >> jbits) & mask[cbits]);
		cbits = 8;
	        c = 0;
	    }
	}
	/* now jbits < cbits */
	if (jbits) {
	    c = (c << jbits) + (j & mask[jbits]);
	    cbits -= jbits;
	}
    }
    if (cbits != 8) *bitstream++ = c << cbits;
}

