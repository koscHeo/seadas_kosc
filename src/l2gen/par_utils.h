#ifndef _PARUTILS_H
#define _PARUTILS_H

#include <stdio.h>
#include <math.h>
#include "l2_struc.h"
#include "filehandle.h"
#include "l12_proto.h"

//typedef int logical;
// #define FALSE 0
// #define TRUE 1
#define NMODEL 12
#define NPHASE 75

float calc_par(
		l2str *l2rec,
		int ip,
		int nbands,
		float *Lt,
        float taua,
        float angstrom,
        float *wl,
        float *fo,
        float *ko3,
        float *taumolbar);

void GetAerPhase(l2str *l2rec, int ip, int32_t nbands, float angstrom, float *phasea,
		float *omegaa, float *modelAngstrom);

void read_aerosol_par(l2str *l2rec, int32_t nbands, float *tablewavelengths, float *tablephaseangles,
		float *tablealphas, float *tableomegas, float *tableaerphasefunc);

void* allocateMemoryPar(size_t numBytes, const char* name);

float EstimateDobson(int32_t year, int32_t month, int32_t day, float lat);

float EstimateWatVap(int32_t year, int32_t month, int32_t day, float lat);

float varsol(int32_t jday);
void triseset(int32_t *jday, float xlon, float xlat, float *trise, float *tset);
int Greg2Jul(int32_t year, int32_t month, int32_t day);

// return solar zenith - ignore solar azimuth angle
float get_solz(int jday, float time, float lon, float lat);

float interp_as_taulow(float csz, float tau);

float interp_as_tauhigh(float csz);

#endif
