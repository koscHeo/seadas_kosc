/*
 * prism.h
 *
 *  Created on: June 2015
 *      Author: rhealy
 */

#ifndef SRC_L2GEN_PRISM_H_
#define SRC_L2GEN_PRISM_H_
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sort_double.h>

typedef struct prism_l1b_t {

    // info

    int     npixels;                      /**< number of pixels in prism              */
    int     nscans;                       /**< number of scans in prism              */
    int     nbands;                       /**< number of visible bands in prism      */


} prism_l1b_t;

typedef struct prism_struct {
    int32_t year,day,month, doy, msec,eyear,edoy,emsec;
    int32_t npix,nscan,wgs_nscan,wgs_npix;
    double **sena, **senz, **sola, **solz, **utc;
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
    float alt;
} prism_t;

void readNextLine_av(FILE* fp, char* tag, char* val);
char* getinbasename(char *file);
void prism_proj4_convert(prism_t * data, int numPoints, double *x, double *y);
char* checkTagLine_av(char *line, char* tag);
double getValidAngle(double *ang, int32_t npix, int32_t skip);
int readBinScanLine_float(float *Lt, int32_t recnum, int32_t npix, double *gain, int nbands, int numBands, int interleave, int swap, FILE *ptr);
prism_t* createPrivateData_pr (int numBands, int32_t nscan, int32_t npix);
void freePrivateData_pr(prism_t* data);
void get_zenaz(float *pos, float lon, float lat, float *senz, float *sena);
extern void l_sun_(int32_t *, int32_t *, double *, float *, float *);

#endif /* SRC_L2GEN_PRISM_H_ */
