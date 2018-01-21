#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "wgrib2.h"
#include "fnlist.h"

extern int nx, ny, npnts, scan;
extern int *variable_dim, *raw_variable_dim;
extern enum output_order_type output_order;

/*
 * undo the scan mode madness
 *
 * i, j is grid coordinate (0,0) is south west corner
 * nx, ny are grid dimensions
 * scan_mode is grib2 scan mode parameter
 *
 * returns integer 0 .. nx*ny-1 which is the location of the i, jth data
 *         -1 for error
 */

int ij2p(int i, int j, int scan_mode) {

    if (i < 0 || j < 0) return -1;

    /* regular grid */
    if (nx > 0 || ny > 0) {
        if (i > nx || j > ny) return -1;

        j = (scan_mode & 64) ? j : ny-1 - j;
        i = ((scan_mode & 16) && (j % 2 == 1)) ?  nx - 1 - i : i;
        i = (scan_mode & 128) ? nx-1 - i : i;

        return (scan_mode & 32) ?  j + i*ny : i + nx*j;
    }
    /* thinned longitudes */
    if (nx == -1 && ny > 0) {
	return -1;
    }
   return -1;
}

/*
 * to_we_sn_scan
 *    this routine converts scanning order to standard we:sn
 *    default for binary and text output
 */

int to_we_sn_scan(float *data) {

    float *data2;
    int ix, iy, dx;
    float *p0, *p1, *p2;
    int lscan;

    if (nx == 1 || ny == 1) return 0;
    lscan = scan >> 4;

    if (lscan == 4) return 0; 			/* already we:sn order */

    if ((data2 = (float *) malloc(npnts * sizeof(float))) == NULL) 
	fatal_error("allocation of memory error","");

    if (lscan == 0 && nx > 0 && ny > 0) {	/* regular grid: convert from we:ns to we:sn */
	p0 = data;
	p1 = data2 + npnts;
        for (iy = 0; iy < ny; iy++) {
	    p1 -= nx;
	    memcpy(p1, p0, nx * sizeof(float));
	    p0 += nx;
	}
        memcpy(data, data2, npnts * sizeof(float));
	free(data2);
	return 0;
    }

    if (lscan == 0 && nx == -1 && ny > 0) { /* quasi-regular grid: convert from we:ns to we:sn */
	p0 = data;
	p1 = data2 + npnts;
        for (iy = 0; iy < ny; iy++) {
	    dx = raw_variable_dim[iy];
	    p1 -= dx;
	    memcpy(p1,p0,dx * sizeof(float));
	    p0 += dx;
	}
        memcpy(data, data2, npnts * sizeof(float));
	free(data2);
	return 0;
    }

    if (nx == -1 || ny == -1) fatal_error("not handled by to_we_sn_scan","");
    p0 = data2;
    for (iy = 0; iy < ny; iy++) {
	p1 = data + ij2p(0,iy,scan);
	p2 = data + ij2p(1,iy,scan);
	dx = p2 - p1;
	for (ix = 0; ix < nx; ix++) {
	    *p0++ = *p1;
	    p1 += dx;
	}
    }
    memcpy(data, data2, npnts * sizeof(float));
    free(data2);
    return 0;
}

/*
 * to_we_ns_scan
 *    this routine converts scanning order to standard we:ns
 */

int to_we_ns_scan(float *data) {

    float *data2;
    int ix, iy, dx;
    float *p0, *p1, *p2;
    int lscan;

    if (nx == 1 || ny == 1) return 0;
    lscan = scan >> 4;

    if (lscan == 0) return 0;                   /* already we:ns order */

    if ((data2 = (float *) malloc(npnts * sizeof(float))) == NULL)
        fatal_error("allocation of memory error","");

    if (lscan == 0 && nx > 0 && ny > 0) {       /* regular grid: convert from we:sn to we:ns */
        p0 = data;
        p1 = data2 + npnts;
        for (iy = 0; iy < ny; iy++) {
            p1 -= nx;
            memcpy(p1, p0, nx * sizeof(float));
            p0 += nx;
        }
        memcpy(data, data2, npnts * sizeof(float));
        free(data2);
        return 0;
    }

    if (lscan == 0 && nx == -1 && ny > 0) { /* quasi-regular grid: convert from we:sn to we:ns */
        p0 = data;
        p1 = data2 + npnts;
        for (iy = 0; iy < ny; iy++) {
            dx = raw_variable_dim[iy];
            p1 -= dx;
            memcpy(p1,p0,dx * sizeof(float));
            p0 += dx;
        }
        memcpy(data, data2, npnts * sizeof(float));
        free(data2);
        return 0;
    }

    /* uncommon scan order .. use general routine */

    if (nx == -1 || ny == -1) fatal_error("not handled by to_we_ns_scan","");
    p0 = data2;
    for (iy = ny-1; iy >= 0; iy--) {
        p1 = data + ij2p(0,iy,scan);
        p2 = data + ij2p(1,iy,scan);
        dx = p2 - p1;
        for (ix = 0; ix < nx; ix++) {
            *p0++ = *p1;
            p1 += dx;
        }
    }
    memcpy(data, data2, npnts * sizeof(float));
    free(data2);
    return 0;
}

/*
 * HEADER:200:order:setup:1:decoded data in X (raw|we:sn|we:ns) order, we:sn is default
 */


int f_order(ARG1) {
    if (mode == -1) {
	if (strcmp(arg1,"raw") == 0) output_order = raw;
	else if (strcmp(arg1,"we:sn") == 0) output_order = wesn;
	else if (strcmp(arg1,"we:ns") == 0) output_order = wens;
	else {
	    fatal_error("order: arg=%s expecting raw|we:sn|we:ns", arg1);
	}
    }
    return 0;
}
