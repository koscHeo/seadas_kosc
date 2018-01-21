/*
 * olci.h
 *
 *  Created on: Jul 29, 2015
 *      Author: rhealy
 */

#ifndef SRC_L2GEN_OLCI_H_
#define SRC_L2GEN_OLCI_H_

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sort_double.h>

#define MAXOLCI_RADFILES 21

typedef struct olci_struct {
    int32_t year,day,month, doy, msec;
    int32_t npix,nscan;
    double **sena, **senz, **sola, **solz, **utc;
    double *lat, *lon, *elev;
    double lat0,lon0;
    double *gain;
    double *wave,*fwhm;
    double pixelSize;
    int    numBands, numRadFiles;
    char  *olci_radfiles[MAXOLCI_RADFILES+1];
    char  *geoCoordinatesFile,*tieGeoCoordinatesFile, *tieGeometriesFile,*instrumentFile, *time_coordinatesFile, *tieMeteoFile;
    gsl_spline *spline;
    gsl_interp_accel *spl_acc;
    char olci_varname[MAXOLCI_RADFILES][FILENAME_MAX];
} olci_t;


#endif /* SRC_L2GEN_OLCI_H_ */
