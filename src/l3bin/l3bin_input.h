#ifndef _INPUT_STR_H
#define _INPUT_STR_H

#include <stdio.h>
#include <stdint.h>

typedef struct input_struct {

  char    infile [FILENAME_MAX];
  char    ofile  [FILENAME_MAX];
  char    pfile  [FILENAME_MAX];
  char    out_parm[1024];
  char    tflag;
  char    parms[16384];
  char    pversion[16];

  int32_t    syear;
  int32_t    sday;
  int32_t    eyear;
  int32_t    eday;

  int32_t    sorbit;
  int32_t    eorbit;

  int32_t    reduce_fac;

  int32_t    noext;

  char    merged[4096];

  float   loneast;
  float   lonwest;
  float   latnorth;
  float   latsouth;

  int32_t    verbose;
  int32_t    unit_wgt;
  int32_t    median;

  uint32_t deflate;
  char    oformat [20];

} instr;

//int l3bin_input(int argc, char **argv, instr *input); 

#endif



