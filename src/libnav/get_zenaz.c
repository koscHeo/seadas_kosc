#include <libnav.h>

#include <proj_api.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>


/*
 * get_zenaz
 *
 *  Created on: Aug 23, 2013
 *      Author: dshea
 *
 *  Calculate the zenith and  azimuth for input lat/lon values
 *  given the satellite position.
 *
 *   pos(3)   Orbit Position Vector (km)
 *   lat      Pixel geodetic latitudes
 *   lon      Pixel geodetic longitudes
 *   zenith   Pixel sensor zenith angle
 *   azimuth  Pixel sensor azimuth angle
 *
 */
void get_zenaz(float *pos, float lon, float lat, float *senz, float *sena) {
    int i;
    double gv[3];
    float ea[3];
    float no[3];
    float up[3];
    double rh[3];
    double rl[3];

    double re = 6378.137;
    double f = 1 / 298.257;
    double omf2 = (1.0 - f) * (1.0 - f);
    double xlat = lat * DEG_TO_RAD;
    double xlon = lon * DEG_TO_RAD;

    // Compute the local vertical, East and North unit vectors
    up[0] = cos(xlat) * cos(xlon);
    up[1] = cos(xlat) * sin(xlon);
    up[2] = sin(xlat);
    double upxy = sqrt(up[0] * up[0] + up[1] * up[1]);
    ea[0] = -up[1] / upxy;
    ea[1] = up[0] / upxy;
    ea[2] = 0.0;
    crossp_(up, ea, no);

    // Compute geocentric position vector
    double xlatg = atan(tan(xlat) * omf2);
    gv[0] = cos(xlatg) * cos(xlon);
    gv[1] = cos(xlatg) * sin(xlon);
    gv[2] = sin(xlatg);
    double r = re * (1.0 - f) / sqrt(1.0 - (2.0 - f) * f * pow(cos(xlatg), 2));

    //Transform the pixel-to-spacecraft and Sun vectors into the local frame
    gsl_matrix* xlMatrix = gsl_matrix_alloc(3, 3);
    gsl_matrix* rhMatrix = gsl_matrix_alloc(3, 1);
    gsl_matrix* rlMatrix = gsl_matrix_alloc(3, 1);

    for (i = 0; i < 3; i++) {
        gsl_matrix_set(xlMatrix, 0, i, ea[i]);
        gsl_matrix_set(xlMatrix, 1, i, no[i]);
        gsl_matrix_set(xlMatrix, 2, i, up[i]);

        gsl_matrix_set(rhMatrix, i, 0, pos[i] - r * gv[i]);
        gsl_matrix_set(rlMatrix, i, 0, 0.0);
    }
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, xlMatrix, rhMatrix, 0.0,
            rlMatrix);

    for (i = 0; i < 3; i++) {
        rl[i] = gsl_matrix_get(rlMatrix, i, 0);
    }

    gsl_matrix_free(xlMatrix);
    gsl_matrix_free(rhMatrix);
    gsl_matrix_free(rlMatrix);

    // Compute the sensor zenith and azimuth
    *senz = RAD_TO_DEG * atan2(sqrt(rl[0] * rl[0] + rl[1] * rl[1]), rl[2]);

    // Check for zenith close to zero
    if (*senz > 0.05)
        *sena = RAD_TO_DEG * atan2(rl[0], rl[1]);
    else
        *sena = 0.0;

    if (*sena < 0.0)
        *sena += 360.0;

}
