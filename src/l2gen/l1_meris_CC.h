#ifndef  _L1_MERIS_CC_H
#define  _L1_MERIS_CC_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <netcdf.h>
#include "mfhdf.h"
#include "passthebuck.h"
#include "l1_struc.h"
#include "filehandle.h"

int closel1_meris_CC(filehandle *l1file);
int openl1_meris_CC(filehandle *l1file);
int readl1_meris_CC(filehandle *l1file,int32 recnum,l1str *l1rec); 
int readl1_lonlat_meris_CC(filehandle *l1file,int32 recnum,l1str *l1rec); 

#endif
