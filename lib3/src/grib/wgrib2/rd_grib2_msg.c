#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include <float.h>

#include "grb2.h"
#include "wgrib2.h"

/*
 * rd_grib2_msg.c *                              Wesley Ebisuzaki
 *
 * unsigned char rd_grib2_msg(FILE *input, long int pos, *int len)
 *
 * to read grib2 message
 *
 *    msg = rd_grib_msg(input, position, &len)
 *
 *    input is the file
 *    position is the byte location (should be set to zero for first record 
 *    msg = location of message
 *
 *    to get next record: *  position = position + len;
 *
 * rd_grib_msg allocates its own buffer which can be seen by the
 * parsing routines
 *
 * 1/2007 cleanup M. Schwarb
 */

#define BUFF_ALLOC0	(1024*64)
#define MSEEK 2048

static unsigned char *buffer = NULL, *Msg = NULL, *Sec[9], *Sec6_bitmap;
static int buffer_size = 0;

unsigned char *rd_grib2_msg(FILE *input, long int *pos, long int *len, int *num_submsgs){

    long int len_grib, position, i;
    unsigned long int tmp;
    long int n_bytes;
    unsigned char *p;

    /* setup grib buffer */
    if (buffer == NULL) {
        if ((buffer = (unsigned char *) malloc(BUFF_ALLOC0)) == NULL) {
	    fatal_error("not enough memory: rd_grib2_msg","");
	}
        buffer_size = BUFF_ALLOC0;
    }

    /* find Msg and length of message */

    position = *pos;
    Msg = seek_grib2(input, &position, &len_grib, buffer, MSEEK, &n_bytes);

    if (Msg == NULL) {
        *len = 0;
	return NULL;
    }

    /* read all whole grib record .. to save I/O time, add to end of buffer */

    if (len_grib + Msg - buffer > buffer_size) {
	tmp = Msg - buffer;
        buffer_size = len_grib + Msg - buffer + 5000;
        buffer = (unsigned char *) realloc((void *) buffer, buffer_size);
        if (buffer == NULL) fatal_error("rd_grib2_msg: ran out of memory","");
	Msg = buffer + tmp;
    }
    if (fseek(input, *pos+n_bytes, SEEK_SET) == -1) 
           fatal_error("rd_grib2_msg seek, outside the file, bad grib file","");

    i = len_grib + Msg - buffer - n_bytes; 	/* no. of bytes need to read */
    if (i > 0 && fread(buffer+n_bytes, sizeof (unsigned char), i, input) != i) 
        fatal_error("rd_grib2_msg, read outside of file, bad grib file","");

    Sec[8] = Msg + len_grib - 4;
    if (Sec[8][0] != '7' || Sec[8][1] != '7' || Sec[8][2] != '7' || Sec[8][3] != '7') {
        fatal_error("rd_grib2_msg, missing end section ('7777')","");
    }
    Sec[0] = Msg;

    /* scan message for number of submessages and perhaps for errors */
    p = Msg +  GB2_Sec0_size;
    i = 0;
    while (p < Sec[8]) {
	if (p[4] == 7) i++;
	p += uint4(p);
    }
    if (p != Sec[8]) {
	fprintf(stderr, "rd_grib2_msg: illegal format, end section expected\n");
	exit(8);
    }
    *num_submsgs = i;

/*    *len = len_grib + (position-*pos); */

    *pos = position;
    *len = len_grib;
    return Msg;
}

/*
 * with grib 1, a message = 1 field
 * with grib 2, a message can have more than one field
 *
 * this routine parses a grib2 message that has already been read into buffer
 *
 * parse_1st_msg .. returns 1st message starting at Msg
 */ 

int parse_1st_msg(unsigned char **sec) {

	unsigned char *p;
	int i;
 
	if (Msg == NULL) fatal_error("parse_1st_msg .. Msg == NULL","");

	Sec[0] = Msg;
	Sec[1] = Sec[2] = Sec[3] = Sec[4] = Sec[5] = Sec[6] = Sec[7] = 
	Sec6_bitmap = NULL;

	p = Msg + 16;

	while (Sec[8] - p > 0) {
	    if (p[4] > 8) fatal_error_i("parse_1st_msg illegal section %d", (int) p[4]);
	    Sec[p[4]] = p;

	    /* bitmap section */
	    if (p[4] == 6) {
		if (p[5] == 0) {
		    Sec6_bitmap = p;
		}
		else if (p[5] == 254) {
	            fatal_error("parse_1st_msg: illegal grib msg, missing bitmap","");
		}
	    }

	    /* last section */
	    if (p[4] == 7) {
		for (i = 0; i < 9; i++) {
		    sec[i] = Sec[i];
		}
		return 0;
	    }
	    p += uint4(p);
	}
	fatal_error("parse_1st_msg illegally format grib","");
	return 1;
}

int parse_next_msg(unsigned char **sec) {

	unsigned char *p;
	int i;
 
	p = sec[7];
	if (p[4] != 7) {
            fatal_error("parse_next_msg: parsing error","");
	}
	p += uint4(p);

	while (p < Sec[8]) {
	    Sec[p[4]] = p;
	    if (p[4] == 6 && p[5] == 254) {
		Sec[p[4]] = Sec6_bitmap;
		if (Sec6_bitmap == NULL) fatal_error("parse_next_msg: illegal grib msg, missing bitmap","");
	    }
	    if (p[4] == 7) {
		for (i = 0; i < 9; i++) {
		    sec[i] = Sec[i];
		}
		return 0;
	    }
	    p += uint4(p);
	}
	return 1;
}
