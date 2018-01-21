#ifndef  _L1B_VIIRS_NC_H
#define  _L1B_VIIRS_NC_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <netcdf.h>

#include "l1_struc.h"
#include "filehandle.h"

int closel1b_viirs_nc();
int openl1b_viirs_nc  (filehandle *l1file);
int readl1b_viirs_nc  (filehandle *l1file, int32 recnum, l1str *l1rec); 
int readl1b_lonlat_viirs_nc  (filehandle *l1file, int32 recnum, l1str *l1rec); 

#endif
