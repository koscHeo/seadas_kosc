#ifndef  _L1_HDF_GENERIC_READ_H
#define  _L1_HDF_GENERIC_READ_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mfhdf.h"
#include "passthebuck.h"
#include "l1_struc.h"
#include "filehandle.h"

int closel1_hdf_g     (filehandle *l1file);
int openl1_read_hdf_g (filehandle *l1file);
int readl1_hdf_g      (filehandle *l1file, int32 recnum, l1str *l1rec); 

#endif
