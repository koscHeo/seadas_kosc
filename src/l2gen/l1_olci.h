/*
 * l1_olci.h
 */

#ifndef L1_OLCI_H_
#define L1_OLCI_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "passthebuck.h"
#include "l1_struc.h"
#include "filehandle.h"
#include "olci.h"

int closel1_olci (filehandle *file);
int openl1_olci  (filehandle *file);
int readl1_olci  (filehandle *file, int32_t recnum, l1str *l1rec);
olci_t* createPrivateData_olci (int numBands);

#endif /* L1_OLCI_H_ */

