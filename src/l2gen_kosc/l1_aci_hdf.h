#ifndef  _L1_ACI_HDF_H
#define  _L1_ACI_HDF_H

/* noaa avhrr satellites */                                                     
#define NO06  6    /* the 2nd character is a capital oh, the 3rd is a zero */   
#define NO07  7                                                                 
#define NO08  8                                                                 
#define NO09  9                                                                 
#define NO10 10                                                                 
#define NO11 11                                                                 
#define NO12 12                                                                 
#define NO14 14                                                                 
#define NO15 15                                                                 
#define NO16 16                                                                 
#define NO17 17                                                                 
#define NO18 18                                                                 
#define NO19 19                                                                 

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mfhdf.h"
#include "passthebuck.h"
#include "l1_struc.h"
#include "filehandle.h"

int closel1_aci_hdf (filehandle *l1file);
int openl1_aci_hdf  (filehandle *l1file);
int readl1_aci_hdf  (filehandle *l1file, int32 recnum, l1str *l1rec); 
const char* xsatid2name(int xsatid);
int satname2xsatid(const char* satname);

float etinvert_(int32 *ich, float *radi);
float etintegrate_(int32 *ich, float *etemp);
void etloadresp_(int32_t *lin, char *cal);
int32_t avconsh_(int32_t *lunin, int32_t *npix, int32_t *jday);
int32_t avlooph_(float *satz, float *solz, float *delphi, float *rayly,
        float *aersol, float *aglint);



#endif
