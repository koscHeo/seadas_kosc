#ifndef  _L2_HDF_GENERIC_H
#define  _L2_HDF_GENERIC_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "l12_parms.h"
#include "l2_struc.h"
#include "hdf.h"
#include "mfhdf.h"
#include "passthebuck.h"
#include "filehandle.h"
#include "dfutils.h"
#include <timeutils.h>

int closel2_hdf     (filehandle *l2file);
int openl2_hdf      (filehandle *l2file);
int writel2_hdf     (filehandle *l2file, int32_t recnum, l2str *l2rec); 

#endif
