#ifndef _INPUT_STR_H
#define _INPUT_STR_H

#include <stdio.h>
#include "version.h"

#ifdef __cplusplus
extern "C" {
#endif

#include "netcdf.h"  // Needs to be first to define netcdf stuff JMG
#include "hdf.h"
#include <ctype.h>
#include "hdf5.h"
#include "clo.h"
#include "genutils.h"
#include <assert.h>


typedef struct input_struct {

    char ifile[FILENAME_MAX];
    char ofile[FILENAME_MAX];
    char oformat[20];
    char l3bprod[2048];
    float west;
    float east;
    float north;
    float south;
    float lat;
    float lon;
    int32_t bin_number;
    float radius;
    int verbose;

} instr;

void par_option_cb(clo_option_t *option);
int l3bindump_init_options(clo_optionList_t* list);
int l3bindump_read_options(clo_optionList_t* list, int argc, char* argv[]);

int l3bindump_load_input(clo_optionList_t *list, instr *input);
void l3bindump_input_init(instr *input);
int l3bindump_usage();


#ifdef __cplusplus
}
#endif

#endif
