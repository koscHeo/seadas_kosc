#ifndef  _L1_MOS_HDF_H
#define  _L1_MOS_HDF_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "hdf.h"
#include "mfhdf.h"
#include "l1_struc.h"
#include "filehdr_struc.h"
#include "filehandle.h"
#include "l12_parms.h"
#include "l12_proto.h"

int openl1_read_mos_hdf(filehandle *l1file);
int readl1_mos_hdf(filehandle *l1file, int32_t recnum, l1str *l1rec);
int closel1_mos_hdf(filehandle *l1file);

#endif
