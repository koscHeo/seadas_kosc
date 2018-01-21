#include <stdio.h>
#include "wgrib2.h"
#include "grb2.h"



/*
 * writes the grib message sections into a grib message
 *
 * public domain: wesley ebisuzaki 1/2007
 */

int wrt_grib_secs(unsigned char *sec0, unsigned char *sec1, unsigned char *sec2, unsigned char *sec3,
	unsigned char *sec4, unsigned char *sec5, unsigned char *sec6, unsigned char *sec7, FILE *out) {

    unsigned long int size;
    unsigned char s[8];
    unsigned int size1, size2, size3, size4, size5, size6, size7;

    if (sec0 == NULL) fatal_error("wrt_grib_sec: sec0 is undefined","");
    size = (unsigned long int) GB2_Sec0_size;

    if (sec1 == NULL) fatal_error("wrt_grib_sec: sec1 is undefined","");
    size += (size1 = uint4(sec1));

    size2 = 0;
    if (sec2) size += (size2 = uint4(sec2));

    if (sec3 == NULL) fatal_error("wrt_grib_sec: sec3/gds is undefined","");
    size += (size3 = uint4(sec3));

    if (sec4 == NULL) fatal_error("wrt_grib_sec: sec4/pds is undefined","");
    size += (size4 = uint4(sec4));

    if (sec5 == NULL) fatal_error("wrt_grib_sec: sec5/drs is undefined","");
    size += (size5 = uint4(sec5));

    size6 = 0;
    if (sec6) size += (size6 = uint4(sec6));

    if (sec7 == NULL) fatal_error("wrt_grib_sec: sec7/ds is undefined","");
    size += (size7 = uint4(sec7));

    size += GB2_Sec8_size;

    /* section 0 */
    fwrite((void *) sec0, sizeof(char), 8, out);
    uint8_char(size, s);
    fwrite((void *) s, sizeof(char), 8, out);

    fwrite((void *) sec1, sizeof(char), size1, out);
    if (sec2) fwrite((void *) sec2, sizeof(char), size2, out);
    fwrite((void *) sec3, sizeof(char), size3, out);
    fwrite((void *) sec4, sizeof(char), size4, out);
    fwrite((void *) sec5, sizeof(char), size5, out);
    if (sec6) fwrite((void *) sec6, sizeof(char), size6, out);
    fwrite((void *) sec7, sizeof(char), size7, out);
    fwrite((void *) "7777", sizeof(char), 4, out);
    return ferror(out);
}

