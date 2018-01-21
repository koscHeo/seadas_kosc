#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "grib2.h"
#include "grb2.h"
#include "wgrib2.h"
#include "fnlist.h"


extern int decode, flush_mode;
extern int need_output_file;
extern int file_append;
extern enum output_order_type output_order;

extern float *lat, *lon;
extern int latlon;

extern int nx, ny, npnts, res, scan, new_GDS;

extern int header;
extern char *text_format;
extern int text_column;

/*
 * HEADER:100:lola:inv:4:lon-lat grid values X=lon0:nlon:dlon Y=lat0:nlat:dlat Z=file A=[bin|text|spread]
 */

int f_lola(ARG4) {

    int nx, ny, nxny, i, j, k;
    float latitude, longitude;
    double x0,dx, y0,dy;

    struct local_struct {
        int nlat, nlon;
        double lat0, lon0, dlat, dlon;
        FILE *out;
	/* int iptr[nxny]; */
	int *iptr;
    };
    struct local_struct *save;
    char open_mode[4];

    /* initialization phase */

    if (mode == -1) {
        decode = latlon = 1;
        if (sscanf(arg1,"%lf:%d:%lf", &x0, &nx, &dx) != 3) {
            fatal_error("lola parsing longitudes lon0:nx:dlon  %s", arg1);
        }
        if (sscanf(arg2,"%lf:%d:%lf", &y0,&ny,&dy) != 3) {
            fatal_error("lola parsing latitudes lat0:nx:dlat  %s", arg2);
        }

	if (strcmp(arg4,"spread") != 0 && strcmp(arg4,"text") != 0 && strcmp(arg4,"grib") != 0
		&& strcmp(arg4,"bin") != 0) fatal_error("lola bad write mode %s", arg4);

	strcpy(open_mode, file_append ? "a" : "w");
	if (strcmp(arg4,"bin") == 0 || strcmp(arg4,"grib") == 0) strcat(open_mode,"b");

        nxny = nx*ny;
        save = (struct local_struct *)malloc( sizeof(struct local_struct));
        if (save == NULL) fatal_error("lola memory allocation ","");
        *local = save;
	save->nlon = nx;
	save->lon0 = x0;
	save->dlon = dx;

	save->nlat = ny;
	save->lat0 = y0;
	save->dlat = dy;
	save->iptr = (int *) malloc(nx*ny * sizeof(int));

	save->out = fopen(arg3,open_mode);
	if (save->out == NULL) fatal_error("lola could not open file %s", arg3);
	return 0;
    }

    /* cleanup phase */

    if (mode == -2) return 0;

    /* processing phase */

    save = *local;
    nx = save->nlon;
    ny = save->nlat;
    nxny = nx*ny;

    if (new_GDS) {

        if (output_order != wesn) fatal_error("lola only works in we:sn order","");
        if (lat == NULL || lon == NULL || data == NULL) fatal_error("lola: no val","");

	/* find the nearest points for the grid */
	closest_init(sec);
        k = 0;
        for (j = 0; j < ny; j++) {
	    latitude = save->lat0 + j*save->dlat;
            for (i = 0; i < nx; i++) {
               longitude = save->lon0 + i*save->dlon;
               save->iptr[k++] = closest(sec, latitude, longitude);
            }
        }
    }

    if (strcmp(arg4,"spread") == 0) {
        k = 0;
	fprintf(save->out, "longitude, latitude, value,\n");
        for (j = 0; j < ny; j++) {
	    latitude = save->lat0 + j*save->dlat;
            for (i = 0; i < nx; i++) {
                longitude = save->lon0 + i*save->dlon;
	        fprintf(save->out, "%.2f, %.2f, %g,\n",longitude,latitude,data[save->iptr[k++]]);
	    }
        }
    }

    else if (strcmp(arg4,"bin") == 0) {
	i = nxny * sizeof(float);
        if (header) fwrite((void *) &i, sizeof(int), 1, save->out);
	for (i = 0; i < nxny; i++) {
	    fwrite(data + save->iptr[i], sizeof(float), 1, save->out);
	}
        if (header) fwrite((void *) &i, sizeof(int), 1, save->out);
    }
    else if (strcmp(arg4,"text") == 0) {
	/* text output */
        if (header == 1) {
            fprintf(save->out ,"%d %d\n", nx, ny);
        }
	for (i = 0; i < nxny; i++) {
	    fprintf(save->out, text_format, data[save->iptr[i]]);
	    fprintf(save->out, ((i+1) % text_column) ? " " : "\n");
	}
    }
    else if (strcmp(arg4,"grib") == 0) {
	printf("grib writer");

	for (i = 0; i < nxny; i++) {
	    fprintf(save->out, "%g\n", data[save->iptr[i]]);
	}
    }

    if (flush_mode) fflush(save->out);
    return 0;
}

/*
 * create the grib2 gds for a lat-lon grid
 */

unsigned char *gds_lola(int nx, float x0, float dx, int ny, float y0, float dy) {
	static unsigned char gds[100];

    double r;

    uint_char(72, gds);		/* 1-4 length of section */

    gds[4] = 3;			/* number of section */
    gds[5] = 0;

    uint_char(nx*ny, gds+6);	/* number of data points */
	
    gds[10] = 0;		/* number of octtets for optional list of number */
    gds[11] = 0;		/* code 3.11 */

    gds[12] = gds[13] = 0;	/* template 3.0 */

    /* template 3.0 values */

    gds[14] = 0;		/* shape of earth */
    gds[15] = 0;		/* scale factor of radius of spherical earth */
    uint_char(0, gds+16);	/* scaled value of radius of spherical earth */
    gds[20] = 0;		/* scale factor of major radious of oblate earth */
    uint_char(0, gds+21);	/* scaled value of major axis */
    gds[25] = 0;		/* scale factor of minor axis */
    uint_char(0, gds+26);	/* scaled value of minor axis */

    uint_char(nx, gds+30);
    uint_char(ny, gds+34);
    uint_char(0, gds+38);
    uint_char(0, gds+42);

    uint_char((int) (y0*1000000.0), gds+46);	/* lat1 */
    uint_char((int) (x0*1000000.0), gds+50);	/* lon1 */

    gds[54] = 0;				/* flag table 3.3 */

    r = y0 + (ny-1)*dy;
    uint_char((int) (r*1000000.0), gds+55);	/* lat1 */
    r = x0 + (nx-1)*dx;
    if (r >= 360.0) r -= 360.0;
    uint_char((int) (r*1000000.0), gds+59);	/* lat1 */

    int_char((int) (dx*1000000.0), gds+63);	/* Di i direction increment */
    int_char((int) (dy*1000000.0), gds+67);	/* Dj i direction increment */
    gds[71] = 0x40;				/* scanning mode */

	return gds;
}	
