#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "grb2.h"
#include "wgrib2.h"

/*
 * various conversion routines
 *
 * 2006: Public Domain Wesley Ebisuzaki
 * 1/2007: uint8 fix Wesley Ebisuzaki
 */

/* routines to return various sized integers from GRIB file */

unsigned int uint2(unsigned char *p) {
	return (p[0] << 8) + p[1];
}

unsigned int uint4(unsigned char *p) {
    return ((p[0] << 24) + (p[1] << 16) + (p[2] << 8) + p[3]);
}

int uint4_missing(unsigned char *p) {
    int t;

    t = p[0] & 128 ? -1 : 0;
    t = t << 8 | p[0];
    t = t << 8 | p[1];
    t = t << 8 | p[2];
    t = t << 8 | p[3];

    return t;
}

unsigned long int uint8(unsigned char *p) {

#if (ULONG_MAX == 4294967295UL) 
	if (p[0] || p[1] || p[2] || p[3]) {
		fprintf(stderr,"unsigned value (8 byte integer) too large for machine\n");
		fprintf(stderr,"fatal error .. run on 64-bit machine\n");
		exit(8);
	}
	return  ((unsigned long int)p[4] << 24) + ((unsigned long int)p[5] << 16) + 
                ((unsigned long int)p[6] << 8) + (unsigned long int)p[7];
#else
	return  ((unsigned long int)p[0] << 56) + ((unsigned long int)p[1] << 48) + 
                ((unsigned long int)p[2] << 40) + ((unsigned long int)p[3] << 32) + 
                ((unsigned long int)p[4] << 24) + ((unsigned long int)p[5] << 16) +
		((unsigned long int)p[6] << 8) + (unsigned long int)p[7];
#endif
}

int int2(unsigned char *p) {
	int i;
	if (p[0] & 0x80) {
		i = -(((p[0] & 0x7f) << 8) + p[1]);
	}
	else {
		i = (p[0] << 8) + p[1];
	}
	return i;
}
int int4(unsigned char *p) {
	int i;
	if (p[0] & 0x80) {
		i = -(((p[0] & 0x7f) << 24) + (p[1] << 16) + (p[2] << 8) + p[3]);
	}
	else {
		i = (p[0] << 24) + (p[1] << 16) + (p[2] << 8) + p[3];
	}
	return i;
}

float scaled2flt(int scale_factor, int scale_value) {
   if (scale_factor == 0) return (float) scale_value;
   /* old return scale_value * Int_Power(0.1, scale_factor); */ 
   if (scale_factor < 0) return scale_value * Int_Power(10.0, -scale_factor);
   return scale_value / Int_Power(10.0, scale_factor);
}

void uint8_char(unsigned long int i, unsigned char *p) {
    int j;
    for (j = 0; j < 8; j++) {
	p[7-j] = i & 255;
        i = i >> 8;
    }
}

void uint_char(unsigned int i, unsigned char *p) {
    p[0] = (i >> 24) & 255;
    p[1] = (i >> 16) & 255;
    p[2] = (i >>  8) & 255;
    p[3] = (i      ) & 255;
}

void int_char(int i, unsigned char *p) {
    int sign = 0;
    if (i < 0) {
	sign = 128;
	i = -i;
    }
    p[0] = ((i >> 24) & 127) | sign;
    p[1] = (i >> 16) & 255;
    p[2] = (i >>  8) & 255;
    p[3] = (i      ) & 255;
    return;
}

void uint2_char(unsigned int i, unsigned char *p) {
    p[0] = (i >>  8) & 255;
    p[1] = (i      ) & 255;
    return;
}

void int2_char(int i, unsigned char *p) {
    int sign = 0;
    if (i < 0) {
	sign = 128;
	i = -i;
    }
    p[0] = ((i >> 8) & 127) | sign;
    p[1] = i & 255;
    return;
}
