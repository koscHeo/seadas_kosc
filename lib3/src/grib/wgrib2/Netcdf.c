/******************************************************************************************
 Copyright (C) 2007 Kristian Nilssen
 This file is part of wgrib2 and is distributed under terms of the GNU General Public License
 For details see, Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 Boston, MA  02110-1301  USA

 3/2007 1st release version, changes by Wesley Ebisuzaki

 Version 09.01.2008 (Sergey Varlamov):

  -with the help of Rich Signell added some support for export
   of Polar (grid_template=20) and Lambert (grid_template=30)
   projected data under the CF-1.0 convention.
   Mapping is done using "coordinates" attribute and directly
   defining in the netcdf file the lat-lon values for all points.
   Other possibility would be to specify in the netcdf file "grid_mapping"
   attributes (see CF convention) when mapping to lat-lon will be done
   by user processing software. It would extend range of eligible
   for conversion grid templates in future.
   CF convention NetCDF files are not directly compatible with GrADS v1.9b4 but
   potentially could be used with user-provided "GrADS data description file (catalog)".
  -changed behavior for -nc_pack option: with 'float' packing it does only the check
   of data to be in valid data range, with 'short' or 'byte' directives
   next the simple packing of data is going. Data outside of valid range
   are writted as _FillValue values (missing). Originally it was an error
   if non-missing data were found out of packing range.
  -updated help file.
  _FillValue=UNDEFINED for float non-packed data.

 Version 13.07.2007 (Sergey Varlamov):

  -help file was updated; removed dublicating instructions from this code

  -corrected warnings with some compilers (like ecc) when printing size_t X values
   usind the %d format: now apply %lu format and (unsigned long)X conversion
   in print operations.

  -added some attributes to the generated netcdf file that could help to make it more
   "self-describing" as it is the strong feature of netcdf. These include time attributes
   like reference_time/date if it is possible to fix, time_step etc.
   The global attributes like center and sub-center information etc could be added.
   Does somebody need it?

 Version 18.06.2007 (Sergey Varlamov):

  -in the user defined table of GRIB2 to NETCDF conversion parameters
   it was added support for next optional directives, mainly for
   the 'advanced' users, see help:
    $nlev 5               Equivalent of '-nc_nlev'; has precedence over the command line option
    $levs 1000 500.0 ...  Vertical level values exported to the netcdf file, in netcdf units
    $lev_type 100:pressure:pressure level:mb:0.01
                          User-defined type of vertical level for 4D data.
    $grads 1              Equivalent of  '-nc_grads' if 1 or '-no_nc_grads' if 0; take precedence

  -the total number of fields successfully added or updated in the netcdf file
   is counted. If the user defined g2nc_table is used for data export to netcdf
   and no field were added/updated - warning message is issued at the cleanup stage
   as it could be caused by user error in filling the table. Error example:
   specify $levs in [mb], but left scaling 1,0. As pressure levels in grib file
   are default in [Pa] - no one level will fit to user table in [mb] and data
   will be skipped from output to the netcdf. Hope that warning could help in such cases.

  -corrected some found bugs. These include memory deallocation for the
   such shared objects as user-defined table for GRIB2 to NETCDF convertion
   and some others that could impact updating of existing partially undefined
   netcdf file (undefined vertical levels). Modified some Netcdf.c
   program variable names to be more self-explanatory,
   like g2nc_smlt -> g2nc_4Dlt etc.

 Version 13.06.2007 - changes by Sergey Varlamov:

   -in the user defined table of GRIB2 to NETCDF conversion parameters added
    keyword 'ignore'. New records format could be:

    wgib2_name:wgrib2_level|*:nc_name|ignore[:ignore|no|float|short|byte[:min:max]]

    Min and max values are significant for the short or byte packing only.
    Both could be omitted; it means that automatic scaling for the short
    or byte packing will be estimated from the first entered wgib2_name data
    at the wgrib2_level or at first level in case of '*' as level value.

    If the keyword 'ignore' is found as a netcdf variable name or
    as a packing type value, the corresponding data are ingnored
    and do not written to the netcdf file.
    'Ignore' keyword is recommended if the data from the same
    grib2 file are exported in number of output files (netcdf or other)
    by the same wgrib2 process when the same decoded data could be passed
    to output in other file of any supported type.
    Doing export to single netcdf file using 'ignore' keyword is not recommended
    as corresponding data are first decoded and after that skipped from writing
    to the netcdf file.

   -compare existing variables definition in open netcdf file
    with new data attributes before updating or adding data.

 Version 11.05.2007 - Sergey Varlamov:

   -added posibility to export data to netcdf as 4D data. Now supported are
    4D data, 3D data, or mixed 3D and 4D data in one netcdf file.
    4D data are defined in {TIME,LEV,LAT,LON} space.
    Grib2 types of vertical levels eligible for export to the netcdf as 4D data
    are included in nc_4Dlt (this file) and now include types
    20,100,104,105,107, and 160 of GRIB2 code.
    First found eligible type is used, error is generetad if other
    eligible type data are met in the input stream.
    To activate 4D data output plese use '-nc_nlev' option followed by
    integer 'max_number_of_vertical_levels'.
    The 'max_number_of_vertical_levels' defines vertical dimension of
    4D data exported to the netcdf file and do not apply to the data
    at the 'non-eligible' vertical levels (like mean sea level etc)
    or at levels included in to the user-defined
    GRID2 to NETCDF conversion table explicitly (see description below).
    Such 'invariant' level data are treated as 3D data
    defined in {TIME,LAT,LON} space with possible level information
    included into the variable names as were coded by Wesley Ebisuzaki.
    If 'max_number_of_vertical_levels' is zero - all coming data are threated as 3D,
    with level information included in variable names.
    When existing netcdf file is updated in the '-append' mode (see description below)
    the value of 'max_number_of_vertical_levels' must not exceed initial value provided
    when file was first created (defined).
    First creating netcdf file, vertical level values are not fixed (are undefined)
    and these are defined one-by-one when data at new level are added
    to the netcdf file, up to the 'max_number_of_vertical_levels'.
    It is user responsibility to define monotonically changing sequence
    of vertical levels, as it is required by the COARDS convention.
    Error is generated if new vertical level is specified in non-regular
    order.
    Example1: sequence of data at: 1000 mb, 850 mb, 950 mb - generate error for 950 mb data
    Example2: sequence of data at: 1000 mb, 950 mb, 850 mb, 1000 mb, 700 mb...
    is appropriate as 1000 mb level was already defined and data would be placed
    correctly in to the netcdf file.

   -it is possible to pack data in netcdf using '-nc_pack' option or user defined table
    of GRIB2 to NETCDF conversion parameters (see below).
    With '-nc_pack' option folowed by X=min:max:byte|short] all NEW input variables
    would be packed in short or byte with corresponding offset and scale fitting
    given range in specified data type with possible loss of precision.
    packed=(unpacked-offset)/scale; default no packing    
    NEW means that if some variable was already defined in the netcdf file
    and now is appended to it (in -append mode, see below), the initially defined
    packing parameters are used.
    Both zero min and max values activate 'auto' packing when scale and offset are defined
    from first entered field. When packing, fitting is checked, error if data
    do not fit packing limits.
    Packing in 'byte' is not recognized by GrADS v1.9b4.

   -added possibility to provide user defined table of GRIB2 to NETCDF
    conversion parameters.  It would be read from the user file with name provided
    following the '-nc_table' option.

   -option '-nc_grads' is introduced. GrADS (version 1.9b4) do not support
    non-constant data time stepping and silently generate wrong time stamps for such
    netcdf or opendap data sets when creating nice graphics. With this option an error
    is raised if time step is changing when writing data to the netcdf file. It is possible
    for the sequence of forecasts as example. Packing to byte also is not supported
    by GrADS v1.9b4.

   -added support for the '-append' option; it works fine with netcdf files
    created by wgrib2 as it uses dimension, variable and attribute names
    from the netcdf file generated by this utility.
    The existing netcdf file could be expanded in time as well as new variables
    at invariant vertical levels or in the range of fixed
    vertical levels could be added. For the netcdf (as 'direct access' data set)
    '-append' option really is also 'overwrite if same' if new data come with the
    variable name, level and verftime being same as already recorded in the file.
    As for all other output types - do not forget to delete invalid netcdf file
    if it was left from some errorneous runs. Else you could get unpredictable
    result!

   -changed float "_fillValue" attribute name to data type "_FillValue";
    No problems with GrADS to recognize this missing value,
    although not all is clear for me as in COARDS:
    "In cases where the data variable is packed via the scale_value attribute
     this implies that the missing_value flag is likewise packed. The same holds
     for the _FillValue attribute."...

   -check error code on verftime as in other Wesley Ebisuzaki files
   -many other small changes...
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <ctype.h>

#include "grb2.h"
#include "wgrib2.h"
#include "fnlist.h"

#ifdef USE_NETCDF
#include <netcdf.h>

//#define DEBUG_NC

#ifndef _MAX_PATH
#define _MAX_PATH   256 /* max. length of full pathname */
#endif

#define G2NC_MAX_NLEV 1000
/*
   max number of entries in GRIB2 to NETCDF conversion table,
   security only to avoid infinite read-reallocate cycle
 */
#define G2NC_MAX_VARS 1000

/* Sergey Varlamov:
   list of standard level types that could be clearly written as 4D in netcdf.
   smlt - 'supported multy-level types' of vertical coordinate.
   sname is used as variable name so must be different from all
   other dimension and physical variables written into the netcdf file.
   Selected from Level.c
   Or use code_table_3_15? It includes more candidates that are 'reserved' in Levels.c:
    case 110: string="Geometric height [m]"; break;
    case 111: string="Eta coordinate  []"; break;
    case 112: string="Geopotential height [gpm]"; break;
 */

typedef struct {
    int type;
    char * sname;
    char * lname;
    char * units;
    float scale;      /* multile on it, example: 0.01 to change Pa to gPa or mb */
} g2nc_4Dlt;

#define G2NC_NUM_4DLT 6

g2nc_4Dlt nc_4Dlt[G2NC_NUM_4DLT] = {
  { 20,"klevel","K level","K",1.,},            //"%g K level",
  {100,"plevel","pressure level","mb",0.01,},    //"%g hPa",
  {104,"slevel","sigma level","level",1.,},    //"%g sigma level",
  {105,"hlevel","hybrid level","level",1.,},   //"%g hybrid level",
  {107,"kilevel","K isentropic level","K",1.,}, //"%g K isentropic level",
  {160,"depth","ocean depth","m",1.,},         //"%g m below sea level",
};

/* Sergey Varlamov:
   user defined table of GRIB2 to NETCDF conversion parameters
   It would be read from user file provided with '-nc_table' option
   Format:  wgib2_name:[wgrib2_level|*]:nc_name[:[short|byte]:min:max]
   wgib2_name, wgrib2_level are as returned by wgrib2 in inventory;
   '*' used as wgrib2_level will apply for all levels not given explicitly,
   packing information is optional but it overwrites common packing
   if given with the '-nc_pack' option. Absence of packing info means no packing.
   Both zero min and max values activate 'auto' packing when these values are defined
   from first entered field. When packing, fitting is checked, error if do not fit.
 */
typedef struct {
    char * wgrib2_name;
    char * wgrib2_level;
    char * nc_name;
    int    nc_pack;       /* One of 0, NC_BYTE, NC_SHORT or NC_FLOAT. In last case check the valid_range only */
    float  nc_offset;     /* Used to check valid_range of values and pack data to short or to byte, replace "bad" by _FillValue */
    float  nc_scale;      /* Used to check valid_range of values and pack data to short or to byte, replace "bad" by _FillValue */
    float  nc_valid_min;  /* Used to check valid_range of values for float data, replace "bad" by _FillValue */
    float  nc_valid_max;  /* Used to check valid_range of values for float data, replace "bad" by _FillValue */
    int    ignore;
} g2nc_conv;

typedef struct {     /* default */
  int         used;  /* 0, if was used and is linked to one or more 'local' descriptors, could re-initialized */
  int         grads; /* 0, does check netcdf structure for: 1) fixed time step; 2) no byte packing? */
  int         nlev;  /* -1, number of vertical levels if defined in the table */
  g2nc_4Dlt * lt;    /* NULL, 4D level type description if defined in the table */
  float     * lv;    /* NULL, 4D level values to use if defined in the table */
  int         nvc;   /* 0, counter of variable conversion entries 'vc' in this table */
  g2nc_conv * vc;    /* NULL, variables conversion parameters */
} g2nc_table;

g2nc_table * nc_table = NULL; /* table undefined */

static int nc_grads_compatible = 0;
/*
 * HEADER:100:nc_grads:setup:0:require netcdf file to be grads v1.9b4 compatible (fixed time step only)
 */

int f_nc_grads(ARG0) {
  if (mode == -1) nc_grads_compatible = 1;
  return 0;
}

/*
 * HEADER:100:no_nc_grads:setup:0:netcdf file may be not grads v1.9b4 compatible, variable time step
 */

int f_no_nc_grads(ARG0) {
  if (mode == -1) nc_grads_compatible = 0;
  return 0;
}

/*
 * HEADER:100:nc_pack:setup:1:pack/check limits of all NEW input variables, X=min:max[:byte|short|float]
 *
 * NEW means that if some variable was already defined in the netcdf file and now is appended to it (-append mode)
 * initially defined packing parameters are used. min and max are used to estimate appropriate offset and scaling,
 * if both are zero - automatic scaling goes.
 */

int nc_pack = 0;
float nc_pack_offset = 0.;
float nc_pack_scale = 1.;
float nc_valid_min = 0.;
float nc_valid_max = 0.;
//short sfill_value = pow(2,15)-1;      /*_FillValue for short packing, max*/
//signed char bfill_value = pow(2,7)-1; /*_FillValue for byte packing, max*/
short sfill_value = 32767;     /*_FillValue for short packing, max */
signed char bfill_value = 127; /*_FillValue for byte packing, max */
//float ffill_value = 1e20;      /*_FillValue for float values */
float ffill_value = UNDEFINED;      /*_FillValue for float values */

int f_nc_pack(ARG1) {
  char * pack_to=NULL;
  float range;
  int i;
  if (mode == -1) {

//fprintf(stderr,"nc_pack: get arg: %s\n", arg1);

    pack_to = (char*) malloc((strlen(arg1)+12)*sizeof(char));
    if (pack_to == NULL) fatal_error("nc_pack: error allocating tmp string","");
    i = sscanf(arg1,"%g:%g:%s", &nc_valid_min, &nc_valid_max, pack_to);
    if ( i < 2) {
      fatal_error("nc_pack: bad value, expect min:max[:byte|short|float(default)], found %s", arg1);
//    if (sscanf(arg1,"%g:%g:%s", &nc_pack_offset, &nc_pack_scale, pack_to) != 3) {
//      fatal_error("netcdf: bad nc_pack, expect offset:scale:[byte|short], found %s", arg1);
    }
    range = (nc_valid_max - nc_valid_min);    
    if ( i == 2 ) strcpy(pack_to,"float");

    if (strcmp(pack_to,"byte")==0 || strcmp(pack_to,"BYTE")==0 || strcmp(pack_to,"Byte")==0)
      nc_pack = NC_BYTE;
    else if (strcmp(pack_to,"short")==0 || strcmp(pack_to,"SHORT")==0 || strcmp(pack_to,"Short")==0)
      nc_pack = NC_SHORT;
    else if (strcmp(pack_to,"float")==0 || strcmp(pack_to,"FLOAT")==0 || strcmp(pack_to,"Float")==0)
      nc_pack = NC_FLOAT; /* will check valid_range */
    else
      fatal_error("nc_pack: bad value, expect min:max:byte|short|float, found %s", arg1);
    if ( nc_pack == NC_FLOAT && fabs(range) < 1e-20 ) fatal_error(
      "nc_pack: small valid_range specified, expected min:max[:float], min<max, found %s", arg1);

    /* Center near 0 as signed char or short are used. */
    if (nc_pack != NC_FLOAT) {
      nc_pack_offset = (float)(nc_valid_max+nc_valid_min)*0.5;
      if ( nc_pack == NC_BYTE )
        nc_pack_scale = (float) (range/(bfill_value - 2))*0.5;
      else if ( nc_pack == NC_SHORT )
        nc_pack_scale = (float) (range/(sfill_value - 2))*0.5;
    }
    free(pack_to);

//fprintf(stderr,"nc_pack: nc_pack=%d\n", nc_pack);

  }
  return 0;
}

/*
 * HEADER:100:no_nc_pack:setup:0:no packing in netcdf for NEW variables
 *
 */
int f_no_nc_pack(ARG0) {
  if (mode == -1) {
    nc_pack = 0;
    nc_pack_scale = 1.;
    nc_pack_offset = 0.;
    nc_valid_min = 0.;
    nc_valid_max = 0.;
  }
  return 0;
}

/*
 * HEADER:100:nc_nlev:setup:1:netcdf, X = max LEV dimension for {TIME,LEV,LAT,LON} data
 */
static int nc_nlev = 0;

int f_nc_nlev(ARG1) {
  if (mode == -1) {
    if (sscanf(arg1,"%d", &nc_nlev) != 1) {
      fatal_error("netcdf: bad nc_nlev %s", arg1);
    }
    if (nc_nlev > G2NC_MAX_NLEV) {
      fatal_error("nc_nlev: value exceeds G2NC_MAX_NLEV, %s", arg1);
    }
  }
  return 0;
}

/*
 * HEADER:100:nc_table:setup:1:X is conversion_to_netcdf_table file name
 */
int f_nc_table(ARG1) {
  char  input[_MAX_PATH], on[_MAX_PATH], lv[_MAX_PATH], nn[_MAX_PATH], pk[_MAX_PATH];
  char * prd;
  FILE * fl;
  int i, ir, ierr, itmp;
  g2nc_conv * old_nct;
  float min, max, range, ftmp;

  if (mode == -1) {
    /* avoid 'lost memory' */
    if ( nc_table )
      if ( !nc_table->used)
        fatal_error("nc_table: second table name entered when first was not used: %s", arg1);

    fl = fopen(arg1,"r");
    if(fl == NULL)
      fatal_error("nc_table: can not open file for reading: %s", arg1);

    /* Deconnect from allocated memory, it is mapped in some 'local' structure,
     * and allocate new instance
     */
    nc_table = (g2nc_table*)malloc(sizeof(g2nc_table));
    if ( nc_table == NULL )
      fatal_error("nc_table: error allocating new table: %s", arg1);

    nc_table->vc = NULL;
    nc_table->lt = NULL;
    nc_table->lv = NULL;
    nc_table->nvc = 0;
    nc_table->used = 0;
    nc_table->nlev = -1;
    nc_table->grads = -1;

    ierr = 0;

    while(!feof(fl))
    {/* read strings including \n symbol */
      prd = fgets(input, _MAX_PATH, fl);
      if(prd == NULL) continue;
      if(strlen(input) < 2)continue;
      /* format:

           wgib2_name:wgrib2_level|*:nc_name|ignore[:ignore|no|float|[short|byte:min:max]]

         or/and special instructions that overcome command-line options
         for the netcdf file they are applied:

           # 4D level type description: grib2_type:user_short_name:user:long_name:user_units:scale_to_user_units
           $lev_type 100:plevel:pressure level:mb:0.01   #JMA MSM model upper layer (p) data; scale: Pa->mb
           $nlev 3   #number of vertical levels for 4D variables
           $levs 1000 500 100  #vertical level values in user_units, at least nlev;
                 50 10

         Length of input strings must be less then _MAX_PATH or 256 symbols,
         input on multiple lines is supported.
       */
//
//    next do not work as expected, parse string 'manually', it's C...
//    ir = sscanf(input," %[^:]s:%[^:]s%[^:\n]s%[^:\n]s:%g:%g",
//                     on,    lv,    nn,      pk,    &min, &max);
      ir=0;
      min=max=0;
      prd = input;
      i = sscanf(prd," %[^:]s",on);
      if (i < 1) continue;
      if (on[0]=='#' || on[0]=='\n') continue;  /* comment or line of spaces only, pass */
      if (on[0]=='$'){
        /* specification of vertical level for 4D variables in grib2 to netcdf conversion */
        if (strncmp(on,"$nlev",5) == 0 || strncmp(nn,"$NLEV",5) == 0){
          itmp = -1;
          i = sscanf(prd," %*s %d",&itmp);
          if ( i < 1 || itmp < 0 || itmp > G2NC_MAX_NLEV) {
            fprintf(stderr,"nc_table: large value or error in $nlev definition string:\n%s\n", prd);
            ierr = -1;
            break;
          }
          if (nc_table->nlev >= 0) {
            fprintf(stderr,"nc_table: found dublicate $nlev definition, was %d, new:\n%s\n",
            nc_table->nlev,prd);
            ierr = -2;
            break;
          }
          nc_table->nlev=itmp;
        }
        else if (strncmp(on,"$grads",6) == 0 || strncmp(nn,"$GRADS",6) == 0){
          itmp = -1;
          i = sscanf(prd," %*s %d",&itmp);
          if ( i < 1 || itmp < 0 || itmp > 1) {
            fprintf(stderr,"nc_table: use 0 or 1 in $grads definition string:\n%s\n", prd);
            ierr = -1;
            break;
          }
          if (nc_table->grads >= 0) {
            fprintf(stderr,"nc_table: found dublicate $grads definition, was %d, new:\n%s\n",
            nc_table->grads,prd);
            ierr = -2;
            break;
          }
          nc_table->grads=itmp;
        }
        else if (strncmp(on,"$lev_type",9) == 0 || strncmp(nn,"$LEV_TYPE",9) == 0){
          /* parse input string */
          i = sscanf(prd," %*s %d:",&itmp); /* type */
          if ( i >= 1 ) {
            ir++;
            i = strcspn( prd, ":");
            if (i > 0 && i < strlen(prd) ) {
              prd += (i+1);
              i = sscanf(prd," %[^: ]s",nn); /* sname, no spaces */
              if (i >= 1 ) {
                ir++;
                i = strcspn( prd, ":");
                if (i > 0 && i < strlen(prd) ) {
                  prd += (i+1);
                  i = sscanf(prd," %[^:]s",lv); /* lname */
                  if (i >= 1 ) {
                    ir++;
                    i = strcspn( prd, ":");
                    if (i > 0 && i < strlen(prd) ) {
                      prd += (i+1);
                      i = sscanf(prd," %[^:#\n]s",pk); /* units could be last */
                      if (i >= 1 ) {
                        ir++;
                        i = strcspn( prd, ":");
                        if (i > 0 && i < strlen(prd) ) {
                          prd += (i+1);
                          i = sscanf(prd," %g",&ftmp); /* scale */
                          if (i >= 1 ) {
                            ir++;
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
          if ( ir < 4 ) {
            ierr = -3;
            fprintf(stderr,"nc_table: error %d in $lev_type definition string:\n%s\n", ir,on);
            break;
          }
          if ( nc_table->lt != NULL ) {
            fprintf(stderr,"nc_table: found dublicate $lev_type definition:\n%s\n",on);
            ierr = -4;
            break;
          }
          nc_table->lt = (g2nc_4Dlt*)malloc(sizeof(g2nc_4Dlt));
          if ( ir > 4 )nc_table->lt->scale = ftmp;
          else         nc_table->lt->scale = 1;
          nc_table->lt->type = itmp;
          nc_table->lt->sname = (char *)malloc((strlen(nn)+1)*sizeof(char));
          if (nc_table->lt->sname) strcpy(nc_table->lt->sname,nn);
          else ierr = 1;
          nc_table->lt->lname = (char *)malloc((strlen(lv)+1)*sizeof(char));
          if (nc_table->lt->lname) strcpy(nc_table->lt->lname,lv);
          else ierr = 2;
          nc_table->lt->units = (char *)malloc((strlen(pk)+1)*sizeof(char));
          if (nc_table->lt->units) strcpy(nc_table->lt->units,pk);
          else ierr = 3;
          if ( ierr ) break;
        }
        else if (strncmp(on,"$levs",5) == 0 || strncmp(nn,"$LEVS",5) == 0){
          if ( nc_table->lv ) {
            fprintf(stderr,"nc_table: found dublicate $levs definition:\n%s\n",on);
            ierr = -5;
            break;
          }
          if (nc_table->nlev > 0){
            nc_table->lv = (float *)malloc(nc_table->nlev*sizeof(float));
            if ( nc_table->lv == NULL ) {
              ierr = 4;
              break;
            }
          }
          else {
            fprintf(stderr,"nc_table: found $levs directive for undef $nlev:\n%s\n",prd);
            ierr = -6;
            break;
          }
          i = strcspn( prd, "$");
          prd += (i+4+1);

          for (i=0; i < nc_table->nlev; i++) {
            search_lev: ir = strcspn(prd,"+-.0123456789");
            itmp = strcspn(prd,"#"); /* inline comment */
            if ( ir < 0 || ir >= strlen(prd) || itmp < ir ) { /* read next line */
              prd = fgets(input, _MAX_PATH, fl);
              if (prd == NULL ) {
                fprintf(stderr,"nc_table: error reading multy line $levs \n");
                ierr = -7;
                break;
              }
              goto search_lev;
            }
            itmp = strcspn(prd,":"); /* second important fields separator before number */
            if (itmp <= ir ) {
              fprintf(stderr,"nc_table: error entering $levs, check that there are $nlev values defined\n");
              ierr = -8;
              break;
            }
            prd += ir;
            ir = sscanf(prd,"%g",&nc_table->lv[i]);
            if ( ir < 1 ) {
              fprintf(stderr,"nc_table: $levs formatted input error\n");
              ierr = -9;
              break;
            }
            ir = strcspn(prd," ,;:\n");
            if(ir <= strlen(prd) ) prd += ir;
            else {
              fprintf(stderr,"nc_table: $levs parsing error, do no found fields separator\n");
              ierr = -10;
              break;
            }
          }
          if ( ierr ) break;
        }
        else {
          fprintf(stderr,"nc_table: found unrecognized directive:\n%s\n", prd);
          ierr = -11;
          break;
        }
        continue;
      }
      /* grib2 to netcdf variable convertion rule string parsing */
      ir++;
      i = strcspn( prd, ":");
      if (i > 0 && i < strlen(prd) ) {
        prd += (i+1);
        i = sscanf(prd," %[^:]s",lv);
        if (i >= 1 ) {
          ir++;
          i = strcspn( prd, ":");
          if (i > 0 && i < strlen(prd) ) {
            prd += (i+1);
            i = sscanf(prd," %[^: #\n]s",nn); /* could be last, no spaces */
            if (i >= 1 ) {
              ir++;
              i = strcspn( prd, ":");
              if (i > 0 && i < strlen(prd) ) {
                prd += (i+1);
                i = sscanf(prd," %[^: #\n]s",pk);  /* could be last, no spaces */
                if (i >= 1 ) {
                  ir++;
                  i = strcspn( prd, ":");
                  if (i > 0 && i < strlen(prd) ) {
                    prd += (i+1);
//                    i = sscanf(prd," %g: %g",&offset, &scale);
                    i = sscanf(prd," %g: %g",&min, &max);
                    ir += i;
                  }
                }
              }
            }
          }
        }
      }
      if ( ir < 3 ) {
        /* not eligible string, issue warning and pass */
//        fprintf(stderr,"nc_table: badly formatted string, ignore: %s",input);
//        continue;
        fprintf(stderr,"nc_table: badly formatted string:\n %s",input);
        ierr = -21;
        break;
      }
      if ( nc_table->nvc > G2NC_MAX_VARS ) {
        fprintf(stderr,"nc_table: nvc exceed G2NC_MAX_VARS: %d %d\n",
        nc_table->nvc, G2NC_MAX_VARS);
        ierr = -22;
        break;
      }
      i = nc_table->nvc;
      nc_table->nvc++;
      old_nct = nc_table->vc;

      nc_table->vc = (g2nc_conv*) realloc((void*)old_nct, nc_table->nvc*sizeof(g2nc_conv));
      nc_table->vc[i].ignore = 0;
      nc_table->vc[i].nc_pack = 0;
      nc_table->vc[i].nc_offset = 0.;
      nc_table->vc[i].nc_scale = 1.;
      nc_table->vc[i].nc_valid_min = 0.;
      nc_table->vc[i].nc_valid_max = 0.;
      nc_table->vc[i].wgrib2_name = (char *)malloc((strlen(on)+1)*sizeof(char));
      if (nc_table->vc[i].wgrib2_name)
        strcpy(nc_table->vc[i].wgrib2_name,on);
      else {
        ierr = 21;
        break;
      }
      nc_table->vc[i].wgrib2_level = (char *)malloc((strlen(lv)+1)*sizeof(char));
      if (nc_table->vc[i].wgrib2_level)
        strcpy(nc_table->vc[i].wgrib2_level,lv);
      else {
        ierr = 22;
        break;
      }
      nc_table->vc[i].nc_name = (char *)malloc((strlen(nn)+1)*sizeof(char));
      if (nc_table->vc[i].nc_name){
        strcpy(nc_table->vc[i].nc_name,nn);
        if (strcmp(nn,"ignore")==0 || strcmp(nn,"IGNORE")==0 || strcmp(nn,"Ignore")==0)
          nc_table->vc[i].ignore = 1;
      }
      else {
        ierr = 23;
        break;
      }
      if (ir < 4)continue;

      range = (max - min);
      if (strcmp(pk,"byte")==0 || strcmp(pk,"BYTE")==0 || strcmp(pk,"Byte")==0) {
        nc_table->vc[i].nc_pack = NC_BYTE;
        nc_table->vc[i].nc_scale = (float) (range/(bfill_value - 2))*0.5;
      }
      else if (strcmp(pk,"short")==0 || strcmp(pk,"SHORT")==0 || strcmp(pk,"Short")==0) {
        nc_table->vc[i].nc_pack = NC_SHORT;
        nc_table->vc[i].nc_scale = (float) (range/(sfill_value - 2))*0.5;
      }
      else if (strcmp(pk,"float")==0 || strcmp(pk,"FLOAT")==0 || strcmp(pk,"Float")==0) {
        nc_table->vc[i].nc_pack = NC_FLOAT;
      }
      else if (strcmp(pk,"ignore")==0 || strcmp(pk,"IGNORE")==0 || strcmp(pk,"Ignore")==0) {
        nc_table->vc[i].ignore = 1;
        continue;
      }
      else if (strcmp(pk,"no") == 0 || strcmp(pk,"NO") == 0 ) {
        continue;
      }
      else {
        fprintf(stderr,"nc_table: unsupported packing type, ignored: %s in %s",pk,input);
        continue;  /* not eligible string, warning and pass */
      }
      if (nc_table->vc[i].nc_pack == NC_FLOAT)
      {
        if (fabs(range) < 1e-20) {
          fprintf(stderr,"nc_table: small valid_range specified, expected ...:float:{min}:{max}, min<max: %s",input);
          ierr = 23;
          break;
        }
      }         
      else  nc_table->vc[i].nc_offset = (float) (min+max)*0.5;
      nc_table->vc[i].nc_valid_min=min;
      nc_table->vc[i].nc_valid_max=max;
    }
    fclose(fl);
    if( ierr ) {
      if ( ierr > 0 ) fprintf(stderr,"nc_table: allocation error\n");
      fatal_error("nc_table: fatal error parsing file %s", arg1);
    }

#ifdef DEBUG_NC
fprintf(stderr, "nc_table: total found: %d entries\n",nc_table->nvc);
for (i=0; i < nc_table->nvc; i++) {
  fprintf(stderr, "%2d) %s:%s:%s:%d:%d:%f:%f:%f:%f\n",i,
    nc_table->vc[i].wgrib2_name,
    nc_table->vc[i].wgrib2_level,
    nc_table->vc[i].nc_name,
    nc_table->vc[i].ignore,
    nc_table->vc[i].nc_pack,
    nc_table->vc[i].nc_valid_min,
    nc_table->vc[i].nc_valid_max,
    nc_table->vc[i].nc_offset,
    nc_table->vc[i].nc_scale);
}
if (nc_table->nlev >= 0) {
  fprintf(stderr, "nc_table: $nlev=%d\n",nc_table->nlev);
  if (nc_table->lv)
    for (i=0; i < nc_table->nlev; i++)
      fprintf(stderr, "lev(%d)=%g\n",i,nc_table->lv[i]);
}
if( nc_table->lt ) {
  fprintf(stderr, "nc_table: found $lev_type directive:\n");
  fprintf(stderr, "$lev_type = %d:%s:%s:%s:%g\n",
    nc_table->lt->type,  nc_table->lt->sname,
    nc_table->lt->lname, nc_table->lt->units,
    nc_table->lt->scale);
}
#endif
  }
  return 0;
}
int free_nc_table( g2nc_table * nc_table ) {
  // cleanup
  int i;
  if ( nc_table == NULL ) return 0;
  if ( nc_table->nvc > 0 && nc_table->vc ){
    for (i=0; i < nc_table->nvc; i++) {
      free(nc_table->vc[i].wgrib2_name);
      free(nc_table->vc[i].wgrib2_level);
      free(nc_table->vc[i].nc_name);
    }
  }
  if (nc_table->vc) free(nc_table->vc);
  if (nc_table->lt) free(nc_table->lt);
  if (nc_table->lv) free(nc_table->lv);
  free(nc_table);
  nc_table = NULL;
#ifdef DEBUG_NC
fprintf(stderr,"nc_table: cleaned-up...\n");
#endif
  return 0;
}
/*
 * HEADER:100:no_nc_table:setup:0:disable previously defined conversion_to_netcdf_table
 */
int f_no_nc_table(ARG0) {
  if (mode == -1) {
    if (nc_table) {
      if (nc_table->used) {
       /* deconnect from allocated memory, it is mapped in some 'local' structure */
       nc_table = NULL;
      }
      else fatal_error("no_nc_table: disabling not used -nc_table, is it your intention?","");
    }
  }
  return 0;
}
//vsm: are not used:
//extern int header;
//extern char *filename;
//extern int npts;
extern int file_append;
extern int decode, latlon;
extern float *lat, *lon;
extern int nx, ny, npnts;
extern enum output_order_type output_order;
extern int mode, new_GDS;

int level2(int mode, int type1, float value1, int type2, float value2, int center, int subcenter,
           char *inv_out);

void delchars(char *s, int c);
void rep_chars(char *s, int old, int new);
void fix_units(char *s, int n);
int get_unixtime(int year, int month, int day, int hour, int minute, int second);
void get_nc_dims(int ncid, int * time_dim, int * time_var,
                 int * time_ind, int * verf_utime, int * time_step,
                 int * y_dim, int * x_dim,
                 int * nlev, int * lev_dim, int * lev_var, int * lev_type,
                 int * lev_ind, int * lev_step, int dim_latlon);
void create_nc_dims(int ncid, int * time_dim, int * time_var, int * time_ind,
                    int * y_dim, int * x_dim,
                    int * nlev, int * lev_dim, int * lev_var, int * lev_type,
                    int * lev_ind, int * lev_step, g2nc_table * nc_table,
                    int dim_latlon, unsigned char **sec);
int get_nc_time_ind(int ncid, int verf_utime, int time_var, int * time_ind,
                    int * last_verf_time, int * tm_step);
int update_nc_ref_time(int ncid, int verf_utime, unsigned char **sec,
                   int time_var, int tm_step);
int get_nc_lev_ind(int ncid, g2nc_4Dlt * lt_4D, int level_type1, float lev_val,
                   int nlev,  int lev_dim, int lev_var,
                   int * lev_ind, int * lev_step, int * lev_type,
                   g2nc_table * nc_table);
int get_nc_conv_table(const char * name, const char * level,
                      const g2nc_table * nc_table);
/*
 * HEADER:100:netcdf:output:1:write netcdf data to X
 *
 *  As it seems, now wgrib2 options are activated in order
 *  as they are found on the command line.
 *  In this module all options-dependant actions are done when mode > 0,
 *  it guaranty that all global variables were assigned
 *
 */
int f_netcdf(ARG1) {

  // define a struct to store local info across iterations...
  typedef struct {
    int initialized;
    int ncid;
    int time_dim;
    int time_var;
    int y_dim;
    int x_dim;
    int lev_dim;
    int lev_var;
    char* ncfile;      /* name of file */
    int time_ind;      /* max written to netcdf time step index, 0,1,...*/
    int verf_utime;    /* last written to file verftime value */
    int grads_compatible;
    int time_step;     /* used to check that step is const with '-nc_grads_compatible' option */
    int file_append;   /* is real 'appending' going or it is first call and file is created? */
    int lev_type;      /* only one type allowed to be 4D data */
    int lev_ind;       /* max written to netcdf z-dim index, 0,1,...*/
    int lev_step;      /* to check that levels are going monotonically, save sign only */
    int nc_nlev;       /* max value of z-dim as when file was created */
    int nc_pack;       /* using specified packing */
    float nc_pack_scale, nc_pack_offset; /* packing options */
    float nc_valid_min, nc_valid_max;    /* test and packing options, secondary */
    g2nc_table *nc_table; /* user-defined conversion table */
    int free_nc_table;       /* as same table could be shared by differnt instances of netcdf - does free memory here? */
    int nx;
    int ny;
    int nid;           /* total number of fields written or updated in the netcdf file */
    int dim_latlon;    /* lat and lon coordinates dimension, 1 for COARDS or 2 for CF-1.0 conversion */
  } local_struct;
  local_struct* save;

  g2nc_4Dlt * lt_4D = NULL;
  g2nc_conv * nc_vc = NULL; /* variable conversion parameters */
  int varid, i, ok;
  char varname_buf[_MAX_PATH];
  char level_buf[_MAX_PATH];
  int time_ind;      /* current time step index to write in netcdf */
  int level_ind=-1;   /* current z-dim index to write in netcdf */
  size_t start[4];   /* record start location in 4D = {TIME0, LEV0, LAT0, LON0}; */
  size_t count[4];   /* record size in 4D = {TIMES, LEVS, LATS, LONS}; */
  int dimids[4],test_dimids[4];   //lev_types[2];

  int year, month, day, hour, minute, second;
  char name[_MAX_PATH], desc[_MAX_PATH], unit[_MAX_PATH];
  int grid_template, verf_utime, tm_step, ndims, ind_nct;
  int level_type1 = -1, level_type2 = -1, center, subcenter;
  float level_val1=0, dz, dz1; //level_val2=0,lev_vals[2], test_val1, test_val2;
  float scale_factor, add_offset, valid_min, valid_max;
  double max=0, min=0, range=0;
  nc_type var_type;
  int var_pack;
  float* fdata;
  short* sdata;
  signed char* bdata;
  char *str;

/* vsm: are not used, but could be used later to add variables to the netcdf
        that will specify additional data attributes like ftime...

  char ftime_buf[_MAX_PATH];
  char vt_buf[11];
  char btime_buf[11];
  unsigned char *p;
*/

  if (mode == -1) {    //initialization
    save = (local_struct *) malloc(sizeof(local_struct));
    if (!save) fatal_error("netcdf: %s","error doing malloc of save");
    save->ncfile = (char*) malloc(strlen(arg1)+1);
    if (!(save->ncfile)) fatal_error("netcdf: %s","error doing malloc of ncfile");
    strcpy(save->ncfile, arg1);
    decode = latlon = 1;
    save->initialized = 0; /* not yet */
    save->nid = 0;
    save->time_step = 0;   /* undef */
    save->nc_pack = nc_pack;
    save->nc_pack_scale = nc_pack_scale;
    save->nc_pack_offset = nc_pack_offset;
    save->nc_valid_min = nc_valid_min;
    save->nc_valid_max = nc_valid_max;
    if ( nc_table ) {
      if ( nc_table->grads >= 0)  save->grads_compatible = nc_table->grads;
      else                        save->grads_compatible = nc_grads_compatible;
      if ( nc_table->nlev >= 0)   save->nc_nlev = nc_table->nlev;
      else                        save->nc_nlev = nc_nlev;
      if ( nc_table->used )       save->free_nc_table = 0;
      else                        save->free_nc_table = 1;

      nc_table->used = 1;
      save->nc_table = nc_table;

      if (save->grads_compatible && nc_table->nlev > 2 && nc_table->lv) {
        dz1 = nc_table->lv[1]-nc_table->lv[0];
        for(i = 2; i < nc_table->nlev; i++) {
          dz = nc_table->lv[i]-nc_table->lv[i-1];
          if ( dz*dz1 <= 0 ) {
            fatal_error(
            "netcdf: -nc_grads require monotonic order of $levs in -nc_table for %s",
            save->ncfile);
          }
        }
      }
    }
    else {
      save->grads_compatible = nc_grads_compatible;
      save->nc_nlev = nc_nlev;
      save->nc_table = NULL;
      save->free_nc_table = 0;
    }
    save->file_append = 0;
    if (file_append) { /* try to open existing file for modifications,
                         if open fails - create new file */
      ok = nc_open (save->ncfile, NC_WRITE|NC_SHARE, &(save->ncid));
      if (ok == NC_NOERR) save->file_append = 1;
    }
    if ( !save->file_append )
      netcdf_command( nc_create(save->ncfile, 0, &(save->ncid)) );

    *local = save;

#ifdef DEBUG_NC
fprintf(stderr,"netcdf: initialization, nc_file=%s\n",save->ncfile);
fprintf(stderr,"netcdf: nlev=%d, append=%d ncid=%d\n",
save->nc_nlev,save->file_append, save->ncid);
#endif
  }
  else if (mode == -2) {  // cleanup
    save = *local;
    if ( save->nc_table && !save->nid ) {
      fprintf(stderr,
      "netcdf: WARNING! No data satisfying your selection criterias for %s!\n",save->ncfile);
      fprintf(stderr,"Check grib2 file, grep selection criteria(s) and -nc_table file\n");
    }
    else {
#ifdef DEBUG_NC
fprintf(stderr,"netcdf: were added/updated %d fields to %s, clean-up...\n",save->nid,save->ncfile);
#endif
    }
    netcdf_command( nc_close(save->ncid));
    free(save->ncfile);
    if (save->free_nc_table)
      free_nc_table(save->nc_table);
    free(save);
  }
  else if (mode >= 0) {  // mode 0,1,2 for decoding a message.
    if (output_order != wesn && output_order != wens) fatal_error("netcdf: only works in we:sn order","");
    save = *local;
    if (new_GDS && save->initialized ) fatal_error("netcdf: only one grid type at a time","");
    // Specify dimension of lat-lon coordinates (1 - COARDS, 2 - CF-1.0)
    grid_template = code_table_3_1(sec);
    if (grid_template == 0 || grid_template == 10 || grid_template == 40)
      // COARDS convention only allows lat-lon, mercator and gaussian    
      save->dim_latlon = 1; //"COARDS", rectangular non-rotated grids
    else if (grid_template == 20 || grid_template == 30) {
      // CF convention allows much more...
      save->dim_latlon = 2; //"CF-1.0"; include here also rotated grids when becomes possible
    }
    else
      fatal_error_i("netcdf: doesn't support grid_template %d",grid_template);
      
    if (save->dim_latlon > 1 && save->grads_compatible)
      fatal_error_i("netcdf: non-COARDS netCDF is not GrADS v1.9b4 compatible,\n"
                    "        unsupported grid_template=%d,\n"
                    "        but possibly could be open by GrADS with valid ctl (data description) file",
      grid_template);
    
    // To avoid problems in future: used occasionally are ndata (main), npnts, nx and ny (extern).
    // For the eligible grid_templates ndata==npnts==nx*ny, does it needs to check it once more?

    // What variable, levels we got?
/*#vsm# Question:
f_lev for some templates return 1 and leaves level_buf undefined.
Does treat it as a special case? Or stop the netcdf?
*/
    f_lev(mode, sec, data, ndata, level_buf, local);
    if (strcmp(level_buf, "reserved")==0) return(0);

    getName(sec, mode, NULL, name, desc, unit);

#ifdef DEBUG_NC
fprintf(stderr,"netcdf: Start processing of %s at %s\n", name, level_buf);
#endif

    if ( save->nc_nlev ) ndims = 4;  /* suppose 4D shape */
    else                 ndims = 3;  /* 3D shape */

    /* First check the 'ignore' keyword in the user defined vc (var conversion directives).
     * Look for the "exact" matching of name and level
     * in the user defined conversion table, then for 'any' level (*) */
    ind_nct = get_nc_conv_table(name, level_buf, save->nc_table);
    if ( ind_nct >= 0 ) nc_vc = save->nc_table->vc + ind_nct;
    else                nc_vc = NULL;
    if ( nc_vc ) {
      if ( nc_vc->ignore ) return(0);
      /*
      Found exactly the given variable on givel level: special level,
      special name, remove from candidates for treating as 4D
      */
      ndims = 3;
    }
    else {  /* search this variable for 'any' level value '*' */
      i = get_nc_conv_table(name, "*", save->nc_table);
      if ( i >= 0 )
        if ( save->nc_table->vc[i].ignore ) return(0);
    }
#ifdef DEBUG_NC
fprintf(stderr,"netcdf: Variable do not rejected, ndims=%d(?)\n",ndims);
#endif
    center = GB2_Center(sec);
    subcenter = GB2_Subcenter(sec);

    if ( !save->initialized ) { /* Got first eligible data, query or define netcdf file */
      if (save->file_append) {
        /* locate dimensions and corresponding variables, initialize them;
         * check does existing netcdf file structure is consistent
         * with the new coming data field structure? */
        get_nc_dims(save->ncid,
          &(save->time_dim), &(save->time_var),
          &(save->time_ind), &(save->verf_utime),  &(save->time_step),
          &(save->y_dim), &(save->x_dim),
          &(save->nc_nlev), &(save->lev_dim), &(save->lev_var), &(save->lev_type),
          &(save->lev_ind), &(save->lev_step), save->dim_latlon);
      }
      else {
        // define NEW nc file, assume nx and ny of the 1st message hold for all messages...
        create_nc_dims(save->ncid,
          &(save->time_dim), &(save->time_var), &(save->time_ind),
          &(save->y_dim), &(save->x_dim),
          &(save->nc_nlev), &(save->lev_dim), &(save->lev_var),
          &(save->lev_type),&(save->lev_ind), &(save->lev_step),
          save->nc_table, save->dim_latlon, sec);
      }
      save->nx = nx;
      save->ny = ny;
      save->initialized = 1;
    }
#ifdef DEBUG_NC
fprintf(stderr,"netcdf: initialized, nlev=%d, append=%d, file=%s, ndims=%d(?)\n",
save->nc_nlev,save->file_append,save->ncfile,ndims);
#endif
    /* Do next check for each new field added to the netcdf file */
    if (nx != save->nx || ny != save->ny ){
      fprintf(stderr,"netcdf: error processing %s at %s\n", name, level_buf);
      fprintf(stderr,"netcdf: existing grid (nx*ny) %d x %d, new grid %d x %d\n",
       save->nx,save->ny,nx,ny);
      fatal_error("netcdf: grid size must be same for all variables","");
    }
    //netcdf is open or created and defined, go ahead...

    /* work around z-dim:
     * is this data eligible for treating as 4D?
     * which level it is etc. Used code from Level.c */
    if ( save->nc_table ) lt_4D = save->nc_table->lt;   /* still could be NULL */
    else lt_4D = NULL;

//    level_type1 = sec[4][22];
//    level_type2 = sec[4][28];
    level_type1 = code_table_4_5a(sec);
    level_type2 = code_table_4_5b(sec);
    level_val1 = scaled2flt(INT1(sec[4][23]), int4(sec[4] + 24));
//    level_val2 = scaled2flt(INT1(sec[4][29]), int4(sec[4] + 30));
#ifdef DEBUG_NC
fprintf(stderr,"netcdf: levt1=%d, levt2=%d level_val1=%g\n",level_type1,level_type2, level_val1);
#endif

    if ( ndims == 4 ) {
      /* search level_type1 index in nc_4Dlt; if find - treat data as 4D */
      if (level_type2 == 255) {//must be undefined, now only single-level data are supported as 4D
        if ( lt_4D ) {
          if (lt_4D->type == level_type1 ) level_val1 *= lt_4D->scale;
          else lt_4D = NULL;
        }
        else { /* look in global conversion table */
          for (i=0; i < G2NC_NUM_4DLT; i++) {
            if (nc_4Dlt[i].type == level_type1) {
              level_val1 *= nc_4Dlt[i].scale;
              lt_4D = (g2nc_4Dlt *) &(nc_4Dlt[i]);
              break;
            }
          }
        }
      }
    }
    if ( !lt_4D ) ndims = 3; /* process data as 3D */
    if ( ndims == 4 ) { /* user defined level type for 4D data*/
      level_ind = get_nc_lev_ind( save->ncid, lt_4D,
                  level_type1, level_val1,
                  save->nc_nlev, save->lev_dim, save->lev_var,
                  &(save->lev_ind), &(save->lev_step),
                  &(save->lev_type), save->nc_table );
#ifdef DEBUG_NC
fprintf(stderr,"netcdf: lev=%f ind=%d\n",level_val1, level_ind);
#endif
      if ( level_ind < 0 ) {
        if ( save->nc_table ) {
          if( save->nc_table->lv ) { /* ignore data not at user-defined levels */
#ifdef DEBUG_NC
fprintf(stderr,"netcdf: Skip this field as not at one of user-defined levels\n");
#endif
//          fatal_error("netcdf: 4D data are not on user-defined level","");
            return 0;
          }
        }
        fatal_error("netcdf: can not define or find 4D z-dim level","");
      }
      /* look for user-defined name for variables at 'any' level, only 4D data*/
      if ( !nc_vc ) {
        ind_nct = get_nc_conv_table(name,"*",save->nc_table);
        if ( ind_nct >= 0 ) nc_vc = save->nc_table->vc + ind_nct;
      }
    }

    if ( nc_vc ) {
      /* it is user responsibility to provide legal variable names */
      strcpy(varname_buf,nc_vc->nc_name);
    }
    else {
      /* code to make variable names legal */
      varname_buf[0] = 0;
      if (isdigit(name[0])) strcpy(varname_buf, "N");  // no names starting with digit
      strcat(varname_buf, name);

      if ( ndims == 3 ) {
        strcat(varname_buf, "_");
        strcat(varname_buf, level_buf);
      }
      delchars(varname_buf, ' ');
      rep_chars(varname_buf, '-','_');
      rep_chars(varname_buf, '.','_');
      rep_chars(varname_buf, '+','_');
    }

    if( verftime(sec, &year, &month, &day, &hour, &minute, &second) != 0)
      fatal_error("netcdf: can not get verftime","");

    verf_utime = get_unixtime(year, month, day, hour, minute, second);

    time_ind = get_nc_time_ind(save->ncid, verf_utime,
                    save->time_var, &(save->time_ind),
                    &(save->verf_utime), &tm_step);
    if(time_ind < 0) {
      fatal_error("netcdf: unresolved timing problem","");
    }
    if (save->grads_compatible) {
      if( tm_step ) {
        if( save->time_step ) {
          if (save->time_step != tm_step) {
            fatal_error("netcdf: variable time stepping is not GrADS v1.9b4 compatible","");
          }
        }
        else save->time_step = tm_step;
      }
    }
    i = update_nc_ref_time(save->ncid, verf_utime, sec,
                           save->time_var, tm_step);

//vsm: seems now next code is not used at all, even for debugging...
/*
    f_ftime(mode, sec, data, ndata, ftime_buf, local);
    sprintf(vt_buf,"%.4d%.2d%.2d%.2d", year,month,day,hour);
    p = sec[1];
    sprintf(btime_buf,"%.4d%.2d%.2d%.2d", (p[12]<<8)+p[13], p[14],p[15],p[16]);
*/
    dimids[0] = save->time_dim;
    if ( ndims == 4 ) {
      dimids[1] = save->lev_dim;
      dimids[2] = save->y_dim;
      dimids[3] = save->x_dim;
    }
    else {
      dimids[1] = save->y_dim;
      dimids[2] = save->x_dim;
    }
    ok = nc_inq_varid (save->ncid, varname_buf, &varid);
    var_pack = 0;
    scale_factor = 0.;
    if (ok != NC_NOERR) {
      varid = -1; /* new variable, not in file yet */
      if ( nc_vc ) {
        var_pack = nc_vc->nc_pack;
        scale_factor = nc_vc->nc_scale;
        add_offset = nc_vc->nc_offset;
        valid_min = nc_vc->nc_valid_min;
        valid_max = nc_vc->nc_valid_max;
      }
      else if (save->nc_pack) {
        var_pack = save->nc_pack;
        scale_factor = save->nc_pack_scale;
        add_offset = save->nc_pack_offset;
        valid_min = save->nc_valid_min;
        valid_max = save->nc_valid_max;
      }
    }
    else {  /* same variable found in open netcdf file */
      netcdf_command( nc_inq_vartype (save->ncid, varid, &var_type) );
      /* Check existing variable definition in netcdf with requested */
      if ( nc_vc ) var_pack = nc_vc->nc_pack;
      else if (nc_pack) var_pack = save->nc_pack;
//      if ( var_type == NC_FLOAT ) var_type = 0;  //Jan 2007, vsm
      if ( var_pack != 0 && var_type != var_pack ) {
        fatal_error("netcdf: variable %s exists with other packing size",varname_buf);
      }
      netcdf_command( nc_inq_varndims (save->ncid, varid, &i) );
      if ( i != ndims ) {
        fprintf(stderr,"netcdf: ndims=%d in file, new ndims=%d\n",i,ndims);
        fatal_error("netcdf: variable %s exists having different dimension",varname_buf);
      }
      test_dimids[0]=test_dimids[1]=test_dimids[2]=test_dimids[3]=-1;
      netcdf_command( nc_inq_vardimid (save->ncid, varid, test_dimids) );
      for ( i = 0; i < ndims; i++ ) {
        if (test_dimids[i] != dimids[i]) {
          fatal_error("netcdf: variable %s exists with different dimension shape",varname_buf);
        }
      }
      /* packing information is taken from the open netcdf file */
      if ( var_type == NC_BYTE || var_type == NC_SHORT ) {
        ok = nc_get_att_float(save->ncid, varid, "scale_factor", &scale_factor);
        if (ok == NC_NOERR) {
          ok = nc_get_att_float(save->ncid, varid, "add_offset", &add_offset);
          if (ok == NC_NOERR) var_pack = var_type;
        }
        if (var_pack == 0)
          fatal_error("netcdf: no packing info for SHORT or BYTE var=%s,\n"
                      "        invalid netCDF file",varname_buf);
      }
      else if (var_type == NC_FLOAT) {
        // next attributes COULD present 
        ok = nc_get_att_float(save->ncid, varid, "valid_min", &valid_min);
        if (ok == NC_NOERR) {
          ok = nc_get_att_float(save->ncid, varid, "valid_max", &valid_max);
          if (ok == NC_NOERR) var_pack = NC_FLOAT; // will check the valid_range
        }
      }
    }
    if (var_pack == NC_BYTE && save->grads_compatible)
      fatal_error("netcdf: byte packing is not GrADS v1.9b4 compatible, var=%s,\n"
                  "        but such netCDF file still could be open by GrADS\n"
                  "        with valid ctl (data description) file", varname_buf);

//    if ( var_pack ) {
      // Find max-min for checking of packing parameters
      // vsm: From Jan 2007 version values that do not fit to provided valid_range
      // are replaced by UNDEFINED
    if ( scale_factor == 0. && (var_pack == NC_BYTE || var_pack == NC_SHORT) ) {
      max = min = ok = 0;
      for (i = 0; i < ndata; i++) {
        if (!UNDEFINED_VAL(data[i])) {
          if (ok) {
            max = max > data[i] ? max : data[i];
            min = min < data[i] ? min : data[i];
          }
          else {
            ok = 1;
            max = min = data[i];
          }
        }
      }
    }
    if (varid < 0){ /* new variable, not in file yet */

#ifdef DEBUG_NC
fprintf(stderr,"netcdf: add new var=%s to output file %s\n",varname_buf,save->ncfile);
#endif
      netcdf_command( nc_redef(save->ncid) );

      /* try to apply provided packing to all new variables */
      if ( var_pack == NC_BYTE ) {
        netcdf_command( nc_def_var (save->ncid, varname_buf, NC_BYTE, ndims, dimids, &varid) );
        netcdf_command( nc_put_att_schar(save->ncid, varid, "_FillValue", NC_BYTE, 1, &bfill_value) );
      }
      else if ( var_pack == NC_SHORT ) {
        netcdf_command( nc_def_var (save->ncid, varname_buf, NC_SHORT, ndims, dimids, &varid) );
        netcdf_command( nc_put_att_short(save->ncid, varid, "_FillValue", NC_SHORT, 1, &sfill_value) );
      }
      else { /* default case of unpacked float variables */
        netcdf_command( nc_def_var (save->ncid, varname_buf, NC_FLOAT, ndims, dimids, &varid) );
        netcdf_command( nc_put_att_float(save->ncid, varid, "_FillValue", NC_FLOAT, 1, &ffill_value) );
      }

//vsm: which 'short name' is correct for the NetCDF COARDS compliant file? 'name' or 'varname_buf'?
//      netcdf_command( nc_put_att_text(save->ncid, varid, "short_name", strlen(name), name) );
      netcdf_command( nc_put_att_text(save->ncid, varid, "short_name", strlen(varname_buf), varname_buf) );
      netcdf_command( nc_put_att_text(save->ncid, varid, "long_name", strlen(desc), desc) );
/*
      lev_types[0]=level_type1;
      lev_types[1]=level_type2;
      netcdf_command( nc_put_att_int(save->ncid, varid, "grib2_lt", NC_INT, 2, lev_types) );
*/
      if ( ndims == 4 )
        netcdf_command( nc_put_att_text(save->ncid, varid, "level",
           strlen(lt_4D->lname), lt_4D->lname) );
      else
      {
        netcdf_command( nc_put_att_text(save->ncid, varid, "level", strlen(level_buf), level_buf) );
/*
As an option, could get extended description from f_lev...
        lev_vals[0]=level_val1;
        lev_vals[1]=level_val2;
        netcdf_command( nc_put_att_float(save->ncid, varid, "grib2_lv", NC_FLOAT, 2, lev_vals) );
        ok = level2(mode, level_type1, level_val1, level_type2, level_val2,
                center, subcenter, desc);
        if ( strcmp(desc, level_buf) )
          netcdf_command( nc_put_att_text(save->ncid, varid, "level_description",
                          strlen(desc), desc) );
*/
      }
      fix_units(unit,sizeof(unit));
      netcdf_command( nc_put_att_text(save->ncid, varid, "units", strlen(unit), unit) );
      
      if (save->dim_latlon > 1 ) {
        str="longitude latitude";
        netcdf_command( nc_put_att_text(save->ncid, varid, "coordinates", strlen(str), str) );
      }
      if ( var_pack == NC_BYTE || var_pack == NC_SHORT ) {
        if ( scale_factor == 0.) {
          /* Auto-packing, center near 0 as signed char or short are used.
             Try secure scaling for all next input fields: extend the range on +-10%... */
          range = (max-min)*1.2;
          add_offset = (float) (min+max)*0.5;
          if ( var_pack == NC_BYTE )
            scale_factor = (float) (range/(bfill_value - 2))*0.5;
          else if ( var_pack == NC_SHORT )
            scale_factor = (float) (range/(sfill_value - 2))*0.5;
#ifdef DEBUG_NC
fprintf(stderr,"netcdf: auto-packing, min=%lf, max=%lf, offset=%f, scale=%f, var=%s\n",
min, max, add_offset, scale_factor, varname_buf);
#endif
        }
        netcdf_command( nc_put_att_float(save->ncid, varid, "scale_factor", NC_FLOAT, 1, &scale_factor) );
        netcdf_command( nc_put_att_float(save->ncid, varid, "add_offset", NC_FLOAT, 1, &add_offset) );
      }
      else if ( var_pack == NC_FLOAT ) {
        netcdf_command( nc_put_att_float(save->ncid, varid, "valid_min", NC_FLOAT, 1, &valid_min) );
        netcdf_command( nc_put_att_float(save->ncid, varid, "valid_max", NC_FLOAT, 1, &valid_max) );
      }
      netcdf_command( nc_enddef(save->ncid) );
    }

    // specify where to write new data portion - record size and location in 4D or 3D
    start[0] = time_ind; count[0] = 1;
    if ( ndims == 4 ) {
      start[1] = level_ind; count[1] = 1;
      start[2] = 0;         count[2] = ny;
      start[3] = 0;         count[3] = nx;
    }
    else {
      start[1] = 0;  count[1] = ny;
      start[2] = 0;  count[2] = nx;
    }
    // pack data and write one time step
    if ( var_pack == NC_BYTE || var_pack == NC_SHORT ) { /* get packing attributes from file */
      netcdf_command( nc_get_att_float(save->ncid, varid, "scale_factor", &scale_factor) );
      netcdf_command( nc_get_att_float(save->ncid, varid, "add_offset", &add_offset) );
      if ( fabs(scale_factor) < 1e-20) fatal_error("netcdf: zero packing scale factor %s",varname_buf);
/* vsm: From version Jan 2007 replace data out of range by UNDEFINED
      // stop if data do not fit to with packing parameters:      
      test_val1 = (min - add_offset) / scale_factor;
      test_val2 = (max - add_offset) / scale_factor;
      if ( var_pack == NC_BYTE ) {
        if ( test_val1 <= -bfill_value || test_val2 >= bfill_value ) {
          fprintf(stderr,"values: %lf %lf %f %f %f %f %d %d\n",
          min,max,test_val1,test_val2,add_offset,scale_factor,(int)(-bfill_value),(int)bfill_value );
          fatal_error("netcdf: values do not fit byte packing: %s",varname_buf);
        }
*/
      if ( var_pack == NC_BYTE ) {
        valid_min = -(bfill_value-2)*scale_factor + add_offset;
        valid_max =  (bfill_value-2)*scale_factor + add_offset;
        bdata = (signed char*) malloc(ndata*sizeof(signed char));
        if (!bdata) fatal_error("netcdf: %s","error doing malloc of bdata");
        for (i = 0; i < ndata; i++) {
          bdata[i] = (UNDEFINED_VAL(data[i]) || data[i] < valid_min || data[i] > valid_max ) ?
            bfill_value : (data[i] - add_offset) / scale_factor;
        }
        netcdf_command( nc_put_vara_schar(save->ncid, varid, start, count, bdata) );
        free(bdata);
      }
      else if (var_pack == NC_SHORT ) {
/*      
        if ( test_val1 <= -sfill_value || test_val2 >= sfill_value ) {
          fprintf(stderr,"values: %lf %lf %f %f %f %f %d %d\n",
          min,max,test_val1,test_val2,add_offset,scale_factor,(int)(-sfill_value),(int)sfill_value );
          fatal_error("netcdf: values do not fit short packing: %s",varname_buf);
        }
*/        
        valid_min = -(sfill_value-2)*scale_factor + add_offset;
        valid_max =  (sfill_value-2)*scale_factor + add_offset;
        sdata = (short*) malloc(ndata*sizeof(short));
        if (!sdata) fatal_error("netcdf: %s","error doing malloc of sdata");
        for (i = 0; i < ndata; i++) {
          sdata[i] = (UNDEFINED_VAL(data[i]) || data[i] < valid_min || data[i] > valid_max ) ?
            sfill_value : (data[i] - add_offset) / scale_factor;
        }
        netcdf_command( nc_put_vara_short(save->ncid, varid, start, count, sdata) );
        free(sdata);
      }
    }
    else { /* not packed float values */
      if (var_pack == NC_FLOAT) { // check valid_range
        fdata = (float*) malloc(ndata*sizeof(float));
        if (!fdata) fatal_error("netcdf: %s","error doing malloc of fdata");
        for (i = 0; i < ndata; i++) {
          fdata[i] = (UNDEFINED_VAL(data[i]) || data[i] < valid_min || data[i] > valid_max ) ?
            ffill_value : data[i];
        }
        netcdf_command( nc_put_vara_float(save->ncid, varid, start, count, fdata) );
        free(fdata);
      }  
      else {
        netcdf_command( nc_put_vara_float(save->ncid, varid, start, count, data) );
      }
    }
#ifdef DEBUG_NC
fprintf(stderr,"netcdf: added var=%s to output file %s\n",varname_buf,save->ncfile);
fprintf(stderr,"netcdf: varid=%d, time_ind=%d, lev_ind=%d, pack=%d\n",
varid,time_ind,level_ind,var_pack);
#endif
    save->nid++;
  }
  return 0;
}

void netcdf_command(int status) {
  if (status != NC_NOERR) fatal_error("netcdf error %s", nc_strerror(status));
}

/*
* remove characters from string
*/

void delchars(char *s, int c) {
  char *tail;
  tail = s;
  while (*s) {
    if (*s != c) *tail++ = *s;
    s++;
  }
  *tail = 0;
}

/*
* change any "old" characters into "new"
*/

void rep_chars(char *s, int old, int new) {
  while (*s) {
    if (*s == old) *s = new;
    s++;
  }
}

/*
* makes units COARDS compliant: C -> Celsius g -> gram
*/

void fix_units(char *s, int n) {

  char tmp[n+2], *p;

  tmp[0] = '+';
  strncpy(tmp+1, s, n);
  p = tmp+1;
  *s = 0;

  while (*p) {
    if (toupper(*p) == 'C' && !isalpha(p[1]) && !isalpha(p[-1])) {
      *s++ = 'C';
      *s++ = 'e';
      *s++ = 'l';
      *s++ = 's';
      *s++ = 'i';
      *s++ = 'u';
      *s++ = 's';
      p++;
    }
    else if (toupper(*p) == 'G' && !isalpha(p[1]) && !isalpha(p[-1])) {
      *s++ = 'g';
      *s++ = 'r';
      *s++ = 'a';
      *s++ = 'm';
      p++;
    }
    else {
      *s++ = *p++;
    }
  }
  *s = 0;
}
/*
 * Find UTC "seconds since 1970-01-01 00:00:00.0 0:00"
 */
int get_unixtime(int year, int month, int day, int hour, int minute, int second){
  struct tm t, *gmt_tm;
  time_t local_t, gmt_t;

  t.tm_sec = second;
  t.tm_min = minute;
  t.tm_hour = hour;
  t.tm_mday = day;
  t.tm_mon = month - 1;
  t.tm_year = year - 1900;
  t.tm_isdst = 0;
  local_t = mktime(&t);
  gmt_tm = gmtime(&local_t);
  gmt_t = mktime(gmt_tm);
  return (local_t + (local_t-gmt_t));
}
/* Code from verftime, but define reference time. Does move it to VerfTime.c ? */

int reftime(unsigned char **sec, int *year, int *month, int *day, int *hour, int *minute, int *second) {

    unsigned char *p;

    p = sec[1];
    *year = (p[12] << 8) | p[13];
    *month = p[14];
    *day = p[15];
    *hour = p[16];
    *minute = p[17];
    *second = p[18];

    return 0;
}

void get_nc_dims( int ncid, int * time_dim, int * time_var,
                    int * time_ind, int * verf_utime, int * time_step,
                    int * y_dim, int * x_dim,
                    int * nlev, int * lev_dim, int * lev_var,
                    int * lev_type, int * lev_ind, int * lev_step, int dim_latlon ) {
 /*
  * Check does existing open netcdf file is consistent  with the new coming data?
  * If nlev <= 0 vertical dimension is not looked for and is left undefined.
  * Else its parameters are defined and nlev is changed to number of vertical
  * levels as it is defined when the netcdf file was first created.
  * Sergey Varlamov
  */
  int i, wtime;
  size_t nyy, nxx, ntt, nzz;
  int nco_err, ok;
  float * test_ll;
  float dd, fill_value;
  const float test_d=1e-6;
  char * lev_name;
  size_t start[2], count[2];
  int lat_var, lon_var;

#ifdef DEBUG_NC
fprintf(stderr,"netcdf: step in get_nc_dims\n");
#endif

  nco_err=0;
  nzz = 0;
  /* locate dimensions and corresponding variables with fixed names */
  netcdf_command( nc_inq_dimid (ncid, "time", time_dim) );
  netcdf_command( nc_inq_dimlen (ncid, *time_dim, &ntt) );
  netcdf_command( nc_inq_varid (ncid, "time", time_var) );
  *time_ind = -1; /* undefined */
  *time_step = 0; /* undefined */
  if (ntt > 0) { /* if file is not empty - save last time index, get time value and step */
    *time_ind = ntt-1;
    start[0] = *time_ind;
    netcdf_command( nc_get_var1_int(ncid, *time_var, start, verf_utime) );
    if (ntt > 1) {
      start[0] = ntt-2;
      netcdf_command( nc_get_var1_int(ncid, *time_var, start, &wtime) );
      *time_step = *verf_utime-wtime;
    }
#ifdef DEBUG_NC
fprintf(stderr,"netcdf:get_nc_dims: time=%d %d %lu %d %d\n",
*time_dim,*time_var,(unsigned long)ntt,*time_ind,*verf_utime);
#endif
  }
//  else nco_err+=1;

  if (dim_latlon == 1) i = nx > ny ? nx:ny;
  else i = nx*ny;
  test_ll = (float*) malloc(i*sizeof(float));
  if (!test_ll) fatal_error("netcdf:get_nc_dims: %s","error doing malloc of test_ll");

  if (dim_latlon == 1) {
    // expect lat and lon dimension variables
    netcdf_command( nc_inq_dimid (ncid, "latitude", y_dim) );
    netcdf_command( nc_inq_dimlen (ncid, *y_dim, &nyy) );
    netcdf_command( nc_inq_varid (ncid, "latitude", &lat_var) );
    if (nyy == ny) {
      start[0] = 0; count[0] = ny;
      netcdf_command( nc_get_vara_float(ncid, lat_var, start, count, test_ll) );
      /* use 1% of first step value as criteria */
      if (ny > 1) dd = 0.01*fabs(lat[nx]-lat[0]); /* 'nx' here is not error, wgrib2 specifics! */
      else        dd = test_d*fabs(lat[0]);      /* assumed tolerance */
      for (i=0; i<ny ; i++) {
        if(fabs(lat[i*nx]-test_ll[i]) > dd ){
          fprintf(stderr,"different grid (latitude) in existing netcdf!\n");
          nco_err+=4;
          break;
        }
      }
#ifdef DEBUG_NC
fprintf(stderr,"netcdf:get_nc_dims: latitude=%d %d %lu %f %f\n",
*y_dim,lat_var,(unsigned long)nyy,test_ll[0],test_ll[ny-1]);
#endif
    }
    else nco_err+=2;

    netcdf_command( nc_inq_dimid (ncid, "longitude", x_dim) );
    netcdf_command( nc_inq_dimlen (ncid, *x_dim, &nxx) );
    netcdf_command( nc_inq_varid (ncid, "longitude", &lon_var) );
    if (nxx == nx) {
      start[0] = 0; count[0] = nx;
      netcdf_command( nc_get_vara_float(ncid, lon_var, start, count, test_ll) );
      if (nx > 1) dd = 0.01*fabs(lon[1]-lon[0]);
      else        dd = test_d*fabs(lon[0]);
      for (i=0; i<nx ; i++) {
        if(fabs(lon[i]-test_ll[i]) > dd ){
          fprintf(stderr,"different grid (longitude) in existing netcdf!\n");
          nco_err+=16;
          break;
        }
      }
#ifdef DEBUG_NC
fprintf(stderr,"netcdf:get_nc_dims: longitude=%d %d %lu %f %f\n",
*x_dim,lon_var,(unsigned long)nxx,test_ll[0],test_ll[nx-1]);
#endif
    }
    else nco_err+=8;
  }
  else {
    // expect x and y dimension + latitude and longitude variables
    netcdf_command( nc_inq_dimid (ncid, "y", y_dim) );
    netcdf_command( nc_inq_dimlen (ncid, *y_dim, &nyy) );
    netcdf_command( nc_inq_varid (ncid, "latitude", &lat_var) );
    if (nyy == ny) {
      start[0] = 0; count[0] = ny;
      start[1] = 0; count[1] = nx;
      netcdf_command( nc_get_vara_float(ncid, lat_var, start, count, test_ll) );
      /* use 1% of first step value as criteria */
      if (ny > 1) dd = 0.01*fabs(lat[nx]-lat[0]); /* 'nx' here is not error, wgrib2 specifics! */
      else        dd = test_d*fabs(lat[0]);      /* assumed tolerance */
      for (i=0; i<nx*ny ; i++) {
        if(fabs(lat[i]-test_ll[i]) > dd ){
          fprintf(stderr,"different grid (latitude) in existing netcdf!\n");
          nco_err+=4;
          break;
        }
      }
#ifdef DEBUG_NC
fprintf(stderr,"netcdf:get_nc_dims: latitude=%d %d %lu %f %f\n",
*y_dim,lat_var,(unsigned long)nyy,test_ll[0],test_ll[ny-1]);
#endif
    }
    else nco_err+=2;
    netcdf_command( nc_inq_dimid (ncid, "x", x_dim) );
    netcdf_command( nc_inq_dimlen (ncid, *x_dim, &nxx) );
    netcdf_command( nc_inq_varid (ncid, "longitude", &lon_var) );
    if (nxx == nx) {
      start[0] = 0; count[0] = ny;
      start[1] = 0; count[1] = nx;
      netcdf_command( nc_get_vara_float(ncid, lon_var, start, count, test_ll) );
      if (nx > 1) dd = 0.01*fabs(lon[1]-lon[0]);
      else        dd = test_d*fabs(lon[0]);
      for (i=0; i<nx*ny ; i++) {
        if(fabs(lon[i]-test_ll[i]) > dd ){
          fprintf(stderr,"different grid (longitude) in existing netcdf!\n");
          nco_err+=16;
          break;
        }
      }
#ifdef DEBUG_NC
fprintf(stderr,"netcdf:get_nc_dims: longitude=%d %d %lu %f %f\n",
*x_dim,lon_var,(unsigned long)nxx,test_ll[0],test_ll[nx-1]);
#endif
    }
    else nco_err+=8;
  }
  free(test_ll);

  *lev_ind = -1; /* undefined */
  *lev_step = 0; /* undefined */
  if (*nlev > 0) {
    netcdf_command( nc_inq_ndims (ncid, &i) );
    if (i < 4) {
      fprintf(stderr,"requested 4D output to existing %dD netcdf file\n",i);
      nco_err+=32;
    }
    else { /* find which is z-dimension?*/
      *lev_dim=-1;
      for (i=0;i<4;i++){
        if (*time_dim == i)continue;
        if (*y_dim == i)continue;
        if (*x_dim == i)continue;
        *lev_dim = i;
        break;
      }
      if (*lev_dim >= 0) {
        netcdf_command( nc_inq_dimlen (ncid, *lev_dim, &nzz) );
        if (nzz >= *nlev) { /* get lev_name, lev_var and lev_type from file */
          lev_name = (char*) malloc(NC_MAX_NAME+1);
          if (!lev_name) fatal_error("netcdf:get_nc_dims: %s",
            "error doing malloc of lev_name");
          netcdf_command( nc_inq_dimname (ncid, *lev_dim, lev_name) );
          /* Look for the same-named variable.
             It could be still default undefined (lev_type==-1)
             -this seems to be normal  */
          netcdf_command( nc_inq_varid (ncid, lev_name, lev_var) );
          netcdf_command( nc_get_att_int(ncid, *lev_var, "lev_type", lev_type) );

          test_ll = (float*) malloc(nzz*sizeof(float));
          if (!test_ll) fatal_error("netcdf:get_nc_dims: %s",
            "error doing malloc of test_ll for z_lev");

          start[0] = 0; count[0] = nzz;
          netcdf_command( nc_get_vara_float(ncid, *lev_var, start, count, test_ll) );

          ok = nc_get_att_float(ncid, *lev_var, "_FillValue",  &fill_value);
          if (ok != NC_NOERR) fill_value = NC_FILL_FLOAT;

          dd = fabs(fill_value)*test_d; /* differnce at least in 6th valid digit */

          *lev_ind = nzz-1;  /* suppose all are defined */
          for (i = 0; i < nzz; i++){
            if ( fabs(test_ll[i]-fill_value) <= dd ) {
              *lev_ind = i-1;
              break;
            }
          }
          if ( *lev_ind > 0){
            if (test_ll[1] > test_ll[0] ) *lev_step = 1; /*level increases */
            else *lev_step = -1;
          }
          free(test_ll);

#ifdef DEBUG_NC
fprintf(stderr,"netcdf:get_nc_dims: z-dim=%d %d %lu %d %s\n",
*lev_dim,*lev_var,(unsigned long)nzz,*lev_type,lev_name);
#endif
          /* Increse nlev to max value defined for existing file */
          *nlev = nzz;
          free(lev_name);
        }
        else {
          fprintf(stderr,"requested nz=%d exceeds this in existing netcdf file: %lu\n",
          *nlev,(unsigned long)nzz);
          nco_err+=64;
        }
      }
      else {
        fprintf(stderr,"z-dim get ID error\n");
        nco_err+=128;
      }
    }
  }
  if (nco_err) {
    fprintf(stderr,"*** Error appending data to existing netcdf file ***\n");
    fprintf(stderr,"*** Existing vise input dimensions [t,z,y,x]: [%lu, %lu, %lu, %lu]:[inf, <=%d, %d, %d] ***\n",
    (unsigned long)ntt,(unsigned long)nzz,(unsigned long)nyy,(unsigned long)nxx, *nlev,ny,nx);
    fatal_error_i("netcdf:get_nc_dims: dimension error: %d",nco_err);
  }
}

/* additional definition do not included into the grb2.h */
// #define GDS_Lambert_Lad(gds) (int4(gds+47)*0.000001)

void create_nc_dims(int ncid, int * time_dim, int * time_var, int * time_ind,
                  int * y_dim, int * x_dim,
                  int * nlev, int * lev_dim, int * lev_var, int * lev_type,
                  int * lev_ind, int * lev_step,
                  g2nc_table * nc_table, int dim_latlon, unsigned char **sec){
 /*
  * Define dimensions and dimension variables for new created netcdf file.
  * If nlev <= 0 vertical dimension is not looked for and is left undefined.
  * Else its parameters are defined and nlev is changed to number of vertical
  * levels as it is defined when the netcdf file was first created.
  * Modif. by Sergey Varlamov: time varialbe is left undefined;
  *   added vertical dimension (with tmp name and undefined values).
  *   May 7, 2007
  */
  int i, ref_time[2];
  float * test_ll;
  size_t start[2], count[2];
  int dimids[2];
  char * str, * name, * lname, * units;
  g2nc_4Dlt * g2nc_lt;
  float * g2nc_lv, dx=1, dy=1;
  int lat_var, lon_var, y_var, x_var;
  int grid_template;
  unsigned char *gds;

#ifdef DEBUG_NC
fprintf(stderr,"netcdf: in create_nc_dims\n");
#endif
  if (nc_table) {
    g2nc_lt = nc_table->lt;
    g2nc_lv = nc_table->lv;
  }
  else {
    g2nc_lt = NULL;
    g2nc_lv = NULL;
  }

  // define NEW nc file, assume nx and ny of the 1st message hold for all messages...
  str="";
  if (dim_latlon == 1)
    str = "COARDS";
  else if (dim_latlon == 2)
    str = "CF-1.0";
  else
    fatal_error("netcdf:create_nc_dims: %s","unsupported lat-lon dimension");
   
  netcdf_command( nc_put_att_text(ncid, NC_GLOBAL, "Conventions", strlen(str), str) );
  str = "created by wgrib2";
  netcdf_command( nc_put_att_text(ncid, NC_GLOBAL, "History", strlen(str), str) );

  // create coordinate variables...
  if (dim_latlon == 1) {
    netcdf_command( nc_def_dim (ncid, "latitude", ny, y_dim) );
    netcdf_command( nc_def_dim (ncid, "longitude", nx, x_dim) );
    dimids[0] = *y_dim;
    netcdf_command( nc_def_var (ncid, "latitude", NC_FLOAT, 1, dimids, &lat_var) );
    dimids[0] = *x_dim;
    netcdf_command( nc_def_var (ncid, "longitude", NC_FLOAT, 1, dimids, &lon_var) );
#ifdef DEBUG_NC
fprintf(stderr,"netcdf:create_nc_dims: latitude=%d %d %d\n",*y_dim,lat_var,ny);
fprintf(stderr,"netcdf:create_nc_dims: longitude=%d %d %d\n",*x_dim,lon_var,nx);
#endif
  }
  else {
    // get some info for grid_mapping
    grid_template = code_table_3_1(sec);
    gds = sec[3];
    if (grid_template == 20) {
      dy = fabs(GDS_Polar_dy(gds));
      dx = fabs(GDS_Polar_dx(gds));
//      lov = GDS_Polar_lov(gds);
//      lad = GDS_Polar_lad(gds);
    }
    else if (grid_template == 30) {
      dy = fabs(GDS_Lambert_dy(gds));
      dx = fabs(GDS_Lambert_dx(gds));
//      lov = GDS_Lambert_Lov(gds);
//      lad = GDS_Lambert_Lad(gds); //etc.
    }
    else
      fatal_error_i("netcdf:create_nc_dims: unsupported grid_template %d",grid_template);
  
    netcdf_command( nc_def_dim (ncid, "y", ny, y_dim) );
    netcdf_command( nc_def_dim (ncid, "x", nx, x_dim) );
    dimids[0] = *y_dim;
    netcdf_command( nc_def_var (ncid, "y", NC_FLOAT, 1, dimids, &y_var) );
//    str = "y-axis index";
    str = "y coordinate of projection";
//    str = "projection_y_coordinate";    //it is standard_name,
                                        //but could require grid_mapping definition that is
                                        //projection-dependant; could do it later;
                                        //do not need lat-lon values in this case as these would be
                                        //computed in processing software (CF convention)  
    netcdf_command( nc_put_att_text(ncid, y_var, "long_name", strlen(str), str) );
    dimids[0] = *x_dim;
    netcdf_command( nc_def_var (ncid, "x", NC_FLOAT, 1, dimids, &x_var) );
//    str = "x-axis index";
    str = "x coordinate of projection";
//    str = "projection_x_coordinate";
    netcdf_command( nc_put_att_text(ncid, x_var, "long_name", strlen(str), str) );
    str = "m";
    netcdf_command( nc_put_att_text(ncid, y_var, "units", strlen(str), str) );
    netcdf_command( nc_put_att_text(ncid, x_var, "units", strlen(str), str) );
    netcdf_command( nc_put_att_float(ncid, x_var, "grid_spacing", NC_FLOAT, 1, &dx) );
    netcdf_command( nc_put_att_float(ncid, y_var, "grid_spacing", NC_FLOAT, 1, &dy) );
    dimids[0] = *y_dim;
    dimids[1] = *x_dim;
    netcdf_command( nc_def_var (ncid, "latitude", NC_FLOAT, 2, dimids, &lat_var) );
    netcdf_command( nc_def_var (ncid, "longitude", NC_FLOAT, 2, dimids, &lon_var) );
#ifdef DEBUG_NC
fprintf(stderr,"netcdf:create_nc_dims: y=%d %d %d\n",*y_dim,*y_var,ny);
fprintf(stderr,"netcdf:create_nc_dims: x=%d %d %d\n",*x_dim,*x_var,nx);
fprintf(stderr,"netcdf:create_nc_dims: latitude=%d %d\n",*y_dim * *x_dim,lat_var);
fprintf(stderr,"netcdf:create_nc_dims: longitude=%d %d\n",*y_dim * *x_dim,lon_var);
#endif
  }
  str = "degrees_north";
  netcdf_command( nc_put_att_text(ncid, lat_var, "units", strlen(str), str) );
  str = "degrees_east";
  netcdf_command( nc_put_att_text(ncid, lon_var, "units", strlen(str), str) );
  str = "latitude";
  netcdf_command( nc_put_att_text(ncid, lat_var, "long_name", strlen(str), str) );
  str = "longitude";
  netcdf_command( nc_put_att_text(ncid, lon_var, "long_name", strlen(str), str) );
  
  netcdf_command( nc_def_dim (ncid, "time", NC_UNLIMITED, time_dim) );
  dimids[0] = *time_dim;
  netcdf_command( nc_def_var (ncid, "time", NC_INT, 1, dimids, time_var) );
  str = "seconds since 1970-01-01 00:00:00.0 0:00";
  netcdf_command( nc_put_att_text(ncid, *time_var, "units", strlen(str), str) );
  str = "verification time generated by wgrib2 function verftime()";
  netcdf_command( nc_put_att_text(ncid, *time_var, "long_name", strlen(str), str) );

  /* add more descriptions that could help to recognize what are in the file */
  ref_time[0] = 0; /* value - undefined */
  ref_time[1] = 0; /* type - undefined */
  netcdf_command( nc_put_att_int(ncid, *time_var, "reference_time", NC_INT, 2, ref_time) );
  str = "undefined";
  netcdf_command( nc_put_att_text(ncid, *time_var, "reference_date", strlen(str), str) );
  str = "no";
  netcdf_command( nc_put_att_text(ncid, *time_var, "description", strlen(str), str) );
  /* write time step attribute as zero - good for first or single time step, then update */
  netcdf_command( nc_put_att_int(ncid, *time_var, "time_step", NC_INT, 1, ref_time) );
  /* mark time variale as undefined */
  *time_ind = -1;

#ifdef DEBUG_NC
fprintf(stderr,"netcdf:create_nc_dims: time=%d %d unlimited\n",*time_dim,*time_var);
#endif

  if (*nlev > 0) { /* use tmp name, rename when eligible data will arrive... */
    if (g2nc_lt) name = g2nc_lt->sname;
    else         name = "level";
    netcdf_command( nc_def_dim (ncid, name, *nlev, lev_dim) );
  }
  if (*nlev > 0) { /* use tmp name, rename when eligible data will arrive... */
    dimids[0] = *lev_dim;
    if (g2nc_lt) {
      name = g2nc_lt->sname;
      *lev_type = g2nc_lt->type;
      units = g2nc_lt->units;
      lname = g2nc_lt->lname;
    }
    else {
      name ="level";
      *lev_type=-1;
      units = "undef";
      lname = units;
    }
    netcdf_command( nc_def_var (ncid, name, NC_FLOAT, 1, dimids, lev_var) );
    netcdf_command( nc_put_att_int(ncid, *lev_var, "lev_type", NC_INT, 1, lev_type) );
    netcdf_command( nc_put_att_text(ncid, *lev_var, "units", strlen(units), units) );
    netcdf_command( nc_put_att_text(ncid, *lev_var, "long_name", strlen(lname), lname) );
    // next force the return of "_FillValue" for non-written z-dim request:
    netcdf_command( nc_put_att_float(ncid, *lev_var, "_FillValue", NC_FLOAT, 1, &ffill_value) );

#ifdef DEBUG_NC
fprintf(stderr,"netcdf:create_nc_dims: z-dim=%d %d %d %d\n",
*lev_dim,*lev_var,*nlev,*lev_type);
#endif
  }
  netcdf_command( nc_enddef(ncid) );

  // populate coordinate variables...
  i = nx > ny ? nx:ny;
  test_ll = (float*) malloc(i*sizeof(float));
  if (!test_ll) fatal_error("netcdf:create_nc_dims: %s",
    "error doing malloc of test_ll");
  if (dim_latlon == 1) {
    for (i=0; i<ny ; i++) test_ll[i] = lat[i*nx];
    start[0] = 0; count[0] = ny;
    netcdf_command( nc_put_vara_float(ncid, lat_var, start, count, test_ll) );

    for (i=0; i<nx; i++) test_ll[i] = lon[i];
    start[0] = 0; count[0] = nx;
    netcdf_command( nc_put_vara_float(ncid, lon_var, start, count, test_ll) );
  }
  else {
//    for (i=0; i<ny ; i++) test_ll[i] = i+1.;
    for (i=0; i<ny ; i++) test_ll[i] = dy*i;
    start[0] = 0; count[0] = ny;
    netcdf_command( nc_put_vara_float(ncid, y_var, start, count, test_ll) );

//    for (i=0; i<nx; i++) test_ll[i] = i+1.;
    for (i=0; i<nx; i++) test_ll[i] = dx*i;
    start[0] = 0; count[0] = nx;
    netcdf_command( nc_put_vara_float(ncid, x_var, start, count, test_ll) );
    
    start[0] = 0;  count[0] = ny;
    start[1] = 0;  count[1] = nx;
    netcdf_command( nc_put_vara_float(ncid, lat_var, start, count, lat) );
    netcdf_command( nc_put_vara_float(ncid, lon_var, start, count, lon) );
  }
  free(test_ll);
  *lev_ind = -1;
  *lev_step = 0;

  if (*nlev > 0 && g2nc_lv) {
    *lev_ind = *nlev-1;
    start[0] = 0; count[0] = *nlev;
    netcdf_command( nc_put_vara_float(ncid, *lev_var, start, count, g2nc_lv) );
    if ( *lev_ind > 0){
      if (g2nc_lv[1] > g2nc_lv[0] ) *lev_step = 1; /*level increases */
      else *lev_step = -1;
    }
  }
}

int update_nc_ref_time(int ncid, int verf_utime, unsigned char **sec,
                       int time_var, int tm_step)
{
 /*
  * Update time attributes if necessary and possible
  * Sergey Varlamov
  */
  int ok, ref_utime, update_rt, nc_ref_time[2], nc_tm_step, update_ts;
  int year, month, day, hour, minute, second;
  char * str, ref_date[_MAX_PATH];

#ifdef DEBUG_NC
fprintf(stderr,"netcdf: step in update_nc_ref_time\n");
#endif
  str = "";
  ok = nc_get_att_int(ncid, time_var, "reference_time", nc_ref_time);
  if (ok != NC_NOERR)
  { /* attribute do not exists, could be with update of file from older wgrib2 version */
    nc_ref_time[0] = -1;
    nc_ref_time[1] = -1;
  }
  ok = nc_get_att_int(ncid, time_var, "time_step", &nc_tm_step);
  if (ok != NC_NOERR)
  { /* attribute do not exists, could be with update of file from older wgrib2 version */
    nc_tm_step = -1;
  }
  if ( nc_tm_step < 0 && nc_ref_time[1] < 0) return 1;

  ok = reftime(sec, &year, &month, &day, &hour, &minute, &second);
  ref_utime = get_unixtime(year, month, day, hour, minute, second);
  update_rt = 0;
  update_ts = 0;
  /* nc_ref_time[0] is reference time value,
     nc_ref_time[1] is the reference time type:
      0 - undefined
      1 - analyses, all for the same reference date, could be succeded by forecasts
      2 - analyses, but for different reference date/time
      3 - forecasts from the same reference date/time
     For the type 0 or 2 nc_ref_time[0] keeps first field reference date/time
  */
  if (nc_ref_time[0] == 0 && nc_ref_time[1] == 0) /* attribute exists, but value is undefined */
  { /* define for the first time */
    if (ref_utime == verf_utime)
    {
      nc_ref_time[1] = 1;
      str = "reference date is fixed, analyses";
    }
    else
    {
      nc_ref_time[1] = 3;
      str = "reference date is fixed, forecast or accumulated";
    }
    nc_ref_time[0] = ref_utime;
    sprintf(ref_date,"%.4d.%.2d.%.2d %.2d:%.2d:%.2d UTC",
    year,month,day,hour,minute,second);
    update_rt = 2;
  }
  else if (nc_ref_time[1] > 0)
  { /* exists, value is defined - check it against new and,
       if it is not consistent - modify type */
    switch (nc_ref_time[1]){
    case 1: /* 1 - analyses, all for the same date, could be succeded by forecasts */
      if (nc_ref_time[0] != ref_utime && ref_utime == verf_utime)
      { /* change to analyses for different ref_time */
        nc_ref_time[1] = 2;
        str = "reference date is first field date, analyses";
        update_rt = 1;
      }
      else if (nc_ref_time[0] == ref_utime && ref_utime != verf_utime)
      { /* possibly it is first forecast after analyses */
        nc_ref_time[1] = 3;
        str = "reference date fixed, forecast or accumulated"; /*including analyses*/
        update_rt = 1;
      }
      break;
    case 2: /*2 - analyses, but for different reference date/time */
      if (ref_utime != verf_utime)
      { /* mix, change type to undefined */
        nc_ref_time[1] = 0;
        str = "reference date is first field reference date";
        update_rt = 1;
      }
      break;
    case 3: /* 3 - forecasts from the same reference date/time */
      if (nc_ref_time[0] != ref_utime)
      { /* mix, change type to undefined */
        nc_ref_time[1] = 0;
        str = "reference date is first field reference date";
        update_rt = 1;
      }
      break;
    }
  }
  if (nc_tm_step == 0 && tm_step > 0) /* attribute exists, but value is undefined */
  {
    nc_tm_step = tm_step;
    update_ts = 1;
  }
  else if (nc_tm_step > 0 && tm_step > 0 && nc_tm_step !=  tm_step)
  {
    nc_tm_step = -1; /* variable time step, mark it as undefined */
    update_ts = 1;
  }
  if (update_rt || update_ts)
  {
#ifdef DEBUG_NC
fprintf(stderr,"netcdf: update_nc_ref_time, update_rt=%d\n",update_rt);
#endif
    netcdf_command( nc_redef(ncid) );
    if (update_rt)
    {
      netcdf_command( nc_put_att_int(ncid, time_var, "reference_time", NC_INT, 2, nc_ref_time) );
      netcdf_command( nc_put_att_text(ncid, time_var, "description", strlen(str), str) );
      if (update_rt > 1)
        netcdf_command( nc_put_att_text(ncid, time_var, "reference_date", strlen(ref_date), ref_date) );
    }
    if (update_ts)
      netcdf_command( nc_put_att_int(ncid, time_var, "time_step", NC_INT, 1, &nc_tm_step) );
    netcdf_command( nc_enddef(ncid) );
  }
  return 0;
}

int get_nc_time_ind(int ncid, int verf_utime, int time_var,
                    int * time_ind, int * last_verf_time, int * tm_step){
 /*
  * Compare provided verf_utime with this in open netcdf file and return
  * corresponding index for time dimension. If new value exceeds existing -
  * extend UNLIMITED time dimension writing new entry point.
  * Negative value is for error.
  * Sergey Varlamov
  */
  int t_ind, i, * wtime;
  size_t index[1], start[1], count[1];

#ifdef DEBUG_NC
fprintf(stderr,"netcdf: step in get_nc_time_ind\n");
#endif
  *tm_step = 0;
  if (*time_ind < 0) { /*time dim is undefined*/
    *last_verf_time = verf_utime-1;
    *time_ind = -1;
  }
  if(verf_utime == *last_verf_time) t_ind = *time_ind;
  else if(verf_utime > *last_verf_time) {/* write new time entry point */
    (*time_ind)++;
    t_ind = *time_ind;
    if (t_ind > 0) *tm_step = verf_utime - *last_verf_time;
    index[0] = t_ind;
    *last_verf_time = verf_utime;
    netcdf_command( nc_put_var1_int(ncid, time_var, index, last_verf_time) );
  }
  else { /* verf_utime < *last_verf_time
    Even if grib records are in 'normal' order, accumulated fields like precipitation
    are referenced by start time of accumulation (at least JMA MSM model output)
    - write them if possible. Other choice is to stop here with error
    and ask user to process these variables  separately.
    In this case time_ind and last_verf_time are not modified.
    */
    t_ind = -1;
    if(*time_ind > 0){
      wtime = (int *) malloc((*time_ind + 1)*sizeof(int));
      if (!wtime) fatal_error("netcdf:get_nc_time_ind: %s",
        "error doing malloc of wtime");

      start[0] = 0; count[0] = *time_ind + 1;
      netcdf_command( nc_get_vara_int(ncid, time_var, start, count, wtime) );

      for (i = *time_ind; i >= 0; i--){
        if (wtime[i] == verf_utime){
          t_ind = i;
          break;
        }
      }
      free(wtime);
    }
  }
  return t_ind;
}

int check_nc_timing(int ncid, int time_dim, int time_var) {
 /*
  * Check time values in open netcdf file and detect as error
  * the case with variable time step as GrADS (version 1.9b4)
  * do not support variable data time stepping and silently generates
  * wrong time stamps for such netcdf or opendap files when creating nice graphics.
  * Sergey Varlamov
  */
  int time0,time1, dt0, dt1, i;
  size_t index[1], ntt;

  netcdf_command( nc_inq_dimlen (ncid, time_dim, &ntt) );
  if (ntt > 2) { /* save last time index and get time value */
    index[0] = 0;
    netcdf_command( nc_get_var1_int(ncid, time_var, index, &time0) );
    index[0] = 1;
    netcdf_command( nc_get_var1_int(ncid, time_var, index, &time1) );
    dt1 = time1 - time0;
    for (i = 2; i< ntt; i++) {
      index[0] = i;
      time0 = time1;
      dt0 = dt1;
      netcdf_command( nc_get_var1_int(ncid, time_var, index, &time1) );
      dt1 = time1 - time0;
      if (dt0 != dt1) return 1;
    }
  }
  return 0;
}

int get_nc_lev_ind( int ncid, g2nc_4Dlt * lt_4D,
                    int level_type1, float lev_val,
                    int nlev,  int lev_dim, int lev_var,
                    int * lev_ind, int * lev_step, int * lev_type,
                    g2nc_table * nc_table){
 /*
  * Return index of defined z-dim or -1
  * If still undefined - redefine z_lev axes parameters - name and attributes
  * Sergey Varlamov
  */
  size_t index[1], start[1], count[1];
  int level_ind, i;
  const float test_z=1e-3;  /* const for z-dim test*/
  float dd, * test_lev;
  float * g2nc_lv;

#ifdef DEBUG_NC
fprintf(stderr,"netcdf: step in get_nc_lev_ind\n");
#endif

  if ( !lt_4D ) return -1;
  if (nlev <= 0) return -1;
  if ( nc_table ) g2nc_lv =  nc_table->lv;
  else            g2nc_lv =  NULL;

  /* level value test criteria dd */
  if (fabs(lev_val) < test_z) dd = test_z;
  else dd = fabs(lev_val)*test_z;

  if ( g2nc_lv ) { /* user provided list of levels, check first there */
    level_ind = -1;
    for (i = 0; i < nlev; i++){
      if ( fabs(g2nc_lv[i]-lev_val) < dd ) {
        level_ind = i;
        break;
      }
    }
    if ( level_ind < 0 ) return -1;  /* do not found in the list */
  }
  if (*lev_type < 0) { /* was not fixed yet in netcdf */
    *lev_type = level_type1;
    netcdf_command( nc_redef(ncid) );
    netcdf_command( nc_rename_dim (ncid, lev_dim, lt_4D->sname ) );
    netcdf_command( nc_rename_var (ncid, lev_var, lt_4D->sname ) );
    netcdf_command( nc_put_att_int(ncid, lev_var, "lev_type",NC_INT, 1, lev_type) );
    netcdf_command( nc_put_att_text(ncid, lev_var, "units",
            strlen(lt_4D->units), lt_4D->units) );
    netcdf_command( nc_put_att_text(ncid, lev_var, "long_name",
            strlen(lt_4D->lname), lt_4D->lname) );
    netcdf_command( nc_enddef(ncid) );
#ifdef DEBUG_NC
fprintf(stderr,"netcdf:get_nc_lev_ind: first defined 4D, lev_type=%d lev=%f\n",
level_type1,lev_val);
#endif
  }
  /* look z-dim info in the netcdf */
  level_ind = -1;

  test_lev = (float*) malloc(nlev*sizeof(float));
  if (!test_lev) fatal_error("netcdf:get_nc_lev_ind: %s",
    "error doing malloc of test_lev");

  start[0] = 0; count[0] = nlev;
  netcdf_command( nc_get_vara_float(ncid, lev_var, start, count, test_lev) );

  for (i = 0; i < nlev; i++){
    /* undefined return test_lev = _FillValue */
    if ( fabs(test_lev[i]-lev_val) < dd ) {
      level_ind = i;
      break;
    }
  }
  if ( level_ind < 0 ) {
   /* do not found this level in netcdf;
    * check that levels change monotonically and add new level value
    */
    level_ind = *lev_ind;  /*index of previous defined value */
    if(level_ind >= nlev-1) {
      fatal_error_i(
      "netcdf:get_nc_lev_ind: more vertical levels found then were defined, max=%d",
      nlev);
    }
    *lev_ind += 1;
    if (level_ind >= 0) {
      if ( *lev_step ) {
        if ( (*lev_step)*(lev_val - test_lev[level_ind]) < 0. )
          fatal_error("netcdf:get_nc_lev_ind: %s",
          "non-monotonic levels order in grib2, use -nc_table $levs to fix");
      }
      else {
        if (lev_val > test_lev[level_ind] ) *lev_step = 1; /*level increases */
        else *lev_step = -1;
      }
    }
    level_ind = *lev_ind;   /*index of current defined value, return it */
    index[0] = level_ind;
    netcdf_command( nc_put_var1_float(ncid, lev_var, index, &lev_val) );
  }
  else if ( level_ind > 0 && *lev_step == 0) { /* save levels direction info*/
    if (test_lev[1] > test_lev[0] ) *lev_step = 1; /* increasing order */
    else *lev_step = -1;
  }
  free(test_lev);
  return level_ind;
}

int get_nc_conv_table(const char * name, const char * level,
                      const g2nc_table * nc_table) {
  int i;

  if ( nc_table )
    if ( nc_table->nvc )
      for (i=0; i < nc_table->nvc; i++)
        if ( strcmp(nc_table->vc[i].wgrib2_name, name) == 0 &&
             strcmp(nc_table->vc[i].wgrib2_level, level) == 0 ) return i;
  return -1;
}
#else

int f_netcdf(ARG1) {
  if (mode == -1) {fprintf(stderr,"netcdf package not installed\n"); return 1;}
  return 0;
}
int f_nc_nlev(ARG1) {
  if (mode == -1) {fprintf(stderr,"netcdf package not installed\n"); return 1;}
  return 0;
}
int f_nc_pack(ARG1) {
  if (mode == -1) {fprintf(stderr,"netcdf package not installed\n"); return 1;}
  return 0;
}
int f_nc_table(ARG1) {
  if (mode == -1) {fprintf(stderr,"netcdf package not installed\n"); return 1;}
  return 0;
}
int f_no_nc_table(ARG0) {
  if (mode == -1) {fprintf(stderr,"netcdf package not installed\n"); return 1;}
  return 0;
}
int f_no_nc_pack(ARG0) {
  if (mode == -1) {fprintf(stderr,"netcdf package not installed\n"); return 1;}
  return 0;
}
int f_nc_grads(ARG0) {
  if (mode == -1) {fprintf(stderr,"netcdf package not installed\n"); return 1;}
  return 0;
}
int f_no_nc_grads(ARG0) {
  if (mode == -1) {fprintf(stderr,"netcdf package not installed\n"); return 1;}
  return 0;
}



#endif


