/*
 * aviris.h
 *
 *  Created on: May 18, 2015
 *      Author: rhealy
 */

#ifndef SRC_L2GEN_AVIRIS_H_
#define SRC_L2GEN_AVIRIS_H_
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sort_double.h>

typedef struct aviris_l1b_t {

    // info

    int     npixels;                      /**< number of pixels in AVIRIS              */
    int     nscans;                       /**< number of scans in AVIRIS              */
    int     nbands;                       /**< number of visible bands in AVIRIS      */


} aviris_l1b_t;

typedef struct aviris_struct {
    int32_t year,day,month, doy, msec;
    int32_t npix,nscan,wgs_nscan,wgs_npix;
    double **sena, **senz, **sola, **solz, **utc;
//    double *lat, *lon, *elev;
    double *elev,lat0,lon0;
    double *gain;
    double *wave,*fwhm;
    projPJ *pj_ortho,*pj_latlong;
    double easting,northing;
    double pixelSize;
    int    utmZone, numBands;
    int    interleave, eastbyscan;
    FILE   *av_fp;
    gsl_spline *spline;
    gsl_interp_accel *spl_acc;
} aviris_t;

void readWavInfo(FILE* fp, char* tag, char* val);
void aviris_proj4_convert(aviris_t * data, int numPoints, double *x, double *y);

#endif /* SRC_L2GEN_AVIRIS_H_ */
