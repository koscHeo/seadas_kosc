#ifndef L3MAPGEN_H /* avoid re-inclusion */
#define L3MAPGEN_H

#include <clo.h>
#include "version.h"

enum InterpType {
    Interp_Nearest,
    Interp_Bin,
    Interp_Linear,
    Interp_Area
};

extern clo_optionList_t* optionList;

int l3mapgen_init_options(clo_optionList_t* list, const char* softwareVersion);
int l3mapgen_read_options(clo_optionList_t* list, int argc, char* argv[]);

#endif /* L3MAPGEN_H */
