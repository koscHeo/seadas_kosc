#ifndef  _L1_PRISM_H
#define  _L1_PRISM_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "passthebuck.h"
#include "l1_struc.h"
#include "filehandle.h"

int closel1_prism(filehandle *l1file);
int openl1_prism(filehandle *l1file);
int readl1_prism(filehandle *l1file,int32 recnum,l1str *l1rec);

#endif
