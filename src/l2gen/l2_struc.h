#ifndef _L2_STRUC_H
#define _L2_STRUC_H

#include "input_struc.h"
#include "target_struc.h"
#include "filehandle.h"
#include "hdf.h"

typedef struct l2_struct {
    int32_t   sensorID;
    int32_t   length;
    int32_t   npix;
    int32_t   nbands;
    int32_t   nbandsir;
    int32_t   nscans;
    int32_t  *bindx;
    int32_t   ndets;
    int32_t   iscan;
    int32_t   detnum;
    int32_t   mside;
    int32_t  *nobs;
    double fsol; // Earth-Sun distance correction factor
    float  tilt;
    char   *data;
    int32_t   *year;
    int32_t   *day;
    int32_t   *msec;
    float  *lon;
    float  *lat;
    float  *solz; //solar zenith angle
    float  *sola; //solar azimuth angle
    float  *senz; //sensor zenith angle
    float  *sena; //sensor azimuth angle
    float  *Lt; // Top of atmosphere radiance
    float  *Lt_unc;
    int32_t   *aermodmin;
    int32_t   *aermodmax;
    float  *aerratio2;
    int32_t   *aermodmin2;
    int32_t   *aermodmax2;
    float  *aerratio;
    float  *eps; // NIR aerosol reflectance ratio (single scattering)
    float  *taua; // aerosol optical thickness
    float  *TLg; // directly transmitted glint radiance
    float  *La; // aerosol radiance
    float  *Lw; // water-leaving radiance
    float  *nLw; // normalized water-leaving radiance
    float  *nLw_unc;
    float  *Rrs; //Remote sensing reflectance
    float  *Rrs_unc;
    float  *brdf; //bi-direction reflectance function
    float  *a; //absoprtion coefficient
    float  *bb; //backscattering coefficient
    float  *chl;
    int32_t   *num_iter;

    /* These are just pointers to data in L1 record */
    int32_t   *pixnum;
    unsigned char  *slot;
    float  *Ltir;
    float  *Bt;
    float  *rhof;
    float  *tLf; // white-cap radiance
    float  *Lr; // Rayliegh radiance
    float  *L_q; // Rayleigh polarization Q-component
    float  *L_u; // Rayleigh polarization U-component
    float  *polcor; // polarization correction
    float  *dpol; // degree of polarization
    float  *tg_sol; //gaseous transmittance in the solar path
    float  *tg_sen; //gaseous transmittance in the sensor path
    float  *t_sol; // transmittance in the solar path
    float  *t_sen; // transmittance in the sensor path
    float  *t_o2; // oxygen transmittance
    float  *t_h2o; // water vapor transmittance

    float  *ws; //wind speed
    float  *wd; //wind direction
    float  *mw; //meridonal wind vector
    float  *zw; //zonal wind vector
    float  *pr; //surface pressure
    float  *oz; //ozone
    float  *wv; //water vapor
    float  *rh; //relative humidity
    float  *no2_tropo; //tropospheric NO2 concentration
    float  *no2_strat; //stratospheric NO2 concentration
    float  *no2_frac; //
    float  *height;
    float *elev;
    float  alt; //altitude of sensor
    float  *glint_coef;
    float  *cloud_albedo;
    float  *aerindex;
    float  *sst;
    float  *sstref; // reference sea surface temperature
	short  *ssttype;
    float  *sssref;// reference sea surface salinity
    float  *sw_n; // seawater index of refraction
    float  *sw_a; // seawater absorption coefficient
    float  *sw_bb; //seawater backscattering coefficient
    float  *sw_a_avg; // band-averaged seawater absorption coefficient
    float  *sw_bb_avg; //band-averaged seawater backscattering coefficient
    float  *rhos; //surface reflectance (Rayleigh corrected)
    float  *rho_cirrus; //Cirrus cloud reflectance

    int32_t   *iwave;
    float  *fwave;
    float  *Fo; // Extraterrestrial irradiance
    float  *Fobar; // Mean Extraterrestrial irradiance
    float  *Fonom; // Nominal 11nm F0
    float  *Tau_r; //Rayleigh optical thickness
    float  *k_oz; // Ozone cross-sectional area
    float  *aw;
    float  *bbw;

    float  *delphi; //Relative azimuth
    float  *csolz; //cosine solar zenith
    float  *csenz; //cosine sensor zenith
    float  *alpha; //aerosol angstrom coefficient
    float  *scattang; //scattering angle

    char   *mask;
    int32_t   *flags;

    int32_t   *pixdet;
    float  *radcor;
    float **in_prods;

    instr  *input;
    tgstr  *tgrec;
    filehandle *fileInfo;
    
    float *Rrs_raman;

} l2str;

#endif



