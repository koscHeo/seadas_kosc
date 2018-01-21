#ifndef LIBNAV_H
#define LIBNAV_H

/*
 * Function prototypes for routines defined in src/libnav
 */

#ifdef  __cplusplus
extern "C" {
#endif

    int crossp_(float *v1, float *v2, float *v3);
    void sunangs_( int *year, int *day, float *gmt, float *lon, float *lat,
                   float *sunz, float *suna);

    void get_zenaz(float *pos, float lon, float lat, float *zenith, float *azimuth);
    
    void compute_alpha(float lon[], float lat[],
                       float senz[], float sena[],
                       double mnorm[3], int npix, float alpha[]);

#ifdef  __cplusplus
}
#endif


#endif /* LIBNAV_H */
