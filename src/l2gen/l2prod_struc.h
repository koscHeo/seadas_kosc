#ifndef _L2PROD_STRUC_H
#define _L2PROD_STRUC_H

#include "mfhdf.h"

#define MAX_DIM           3  /* Max number of dimensions for an SDS  */
#define DIMNAMELEN       32  /* String length on dimension names     */
#define UNITLEN          32  /* String length on unit names          */
#define TITLELEN        256  /* String length for long names         */

#define PARAM_TYPE_NONE     0
#define PARAM_TYPE_VIS_WAVE 1
//#define PARAM_TYPE_IR_WAVE  2
#define PARAM_TYPE_ALL_WAVE 3
#define PARAM_TYPE_BAND     4
#define PARAM_TYPE_INT      5

typedef struct l2prod_index_struct {
    int     param_type;
    char    name_prefix[UNITLEN];
    char    name_suffix[UNITLEN];
    int     cat_ix;
    int     prod_ix;
    int32   datatype;
    float32 slope;
    float32 offset;
    double   min;
    double   max;
    int32   rank;
    int32   dim[MAX_DIM];
    char    dimname[MAX_DIM][DIMNAMELEN];
    char    title_format[TITLELEN];
    char    title[TITLELEN];
    char    units[UNITLEN];
    float   badData;
    char    product_id[UNITLEN];
    char    algorithm_id[UNITLEN];
    char    standard_name[TITLELEN];
} l2prodstr;

#endif
