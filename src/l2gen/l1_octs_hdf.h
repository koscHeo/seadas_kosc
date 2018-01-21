#ifndef  _L1_OCTS_HDF_H
#define  _L1_OCTS_HDF_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "l1_struc.h"
#include "filehdr_struc.h"
#include "filehandle.h"
#include "l12_parms.h"
#include "l12_proto.h"
#include "mfhdf.h"


int openl1_read_octs_hdf(filehandle *l1file);
int readl1_octs_hdf(filehandle *l1file, int32_t recnum, l1str *l1rec);
int closel1_octs_hdf(filehandle *l1file);
int navigation(int32 fileID);
int CalcViewAngle(float32 lon1, float32 lat1, float32 pos[3], 
         float32 usun[]);
int LeapCheck(int yr);

#endif

