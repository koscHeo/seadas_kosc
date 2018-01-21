/*
   Copyright (C) 2004-2007,2009,2010,2012 Remik Ziemlinski @ noaa gov

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2, or (at your option)
   any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; see the file COPYING.
   If not, write to the Free Software Foundation,
   59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.  
*/

#include "nccmp.h"
#include "ncinfo.h"
#include "strlist.h"
#include <stdint.h>
#include <float.h>

varstruct vars1[(int)NC_MAX_VARS], vars2[(int)NC_MAX_VARS];
dimstruct dims1[(int)NC_MAX_DIMS], dims2[(int)NC_MAX_DIMS];
size_t nrec1, nrec2;
int nvars1, nvars2, ndims1, ndims2, recid1, recid2;

#define NCFORMATSTR(f)                                          \
    (f == NC_FORMAT_CLASSIC ? "NC_FORMAT_CLASSIC" :               \
      (f == NC_FORMAT_64BIT ? "NC_FORMAT_64BIT" :                 \
      (f == NC_FORMAT_NETCDF4 ? "NC_FORMAT_NETCDF4" :             \
     "NC_FORMAT_NETCDF4_CLASSIC")))                                 \

/* *********************************************************** */
/* Returns formatted string of dimension indices.
  Intended to print out locations of value differences. 
  'out' should be pre-allocated large enough for output. */
void getidxstr(varstruct* var, size_t* start, int curidx, char* out)
{
  int i;
  char tmp[8];
  memset(out,'\0',32);
  for(i=0; i < var->ndims-1; ++i) {
    sprintf(tmp, "%d ", (int)start[i]);
    strcat(out, tmp);
  }
  sprintf(tmp, "%d", curidx);
  strcat(out, tmp);
}
/* *********************************************************** */
/* Same as above but with fortran style indices, which means they're
  1-based (i.e. first index is 1, not 0) and fast-varying dimension is
  printed first (reverse order compared with C-style indices) */
void getidxstr_fortran(varstruct* var, size_t* start, int curidx, char* out)
{
  int i;
  char tmp[8];
  memset(out,'\0',32);
  sprintf(tmp, "%d", curidx + 1);
  strcat(out, tmp);
  
  for(i = var->ndims-2; i >= 0; --i) {
    sprintf(tmp, " %d", (int)start[i]+1);
    strcat(out, tmp);
  }
}
/* *********************************************************** */
/* Creates a string representation of a value into pre-allocated char array. */
void doubleToString(double v, char* out, char* formatprec) { \
  sprintf(out, formatprec, v);
}
/* *********************************************************** */
void floatToString(float v, char* out, char* formatprec) { \
  sprintf(out, formatprec, (double)v);
}
/* *********************************************************** */
void intToString(int v, char* out, char* ignore) { \
  sprintf(out, "%d", v);
}
/* *********************************************************** */
void shortToString(short v, char* out, char* ignore) { \
  sprintf(out, "%d", (int)v);
}
/* *********************************************************** */
void byteToString(unsigned char v, char* out, char* ignore) { \
  sprintf(out, "%c", v);
}
/* *********************************************************** */
void textToString(char v, char* out, char* ignore) { \
  sprintf(out, "%c", v);
}
/* *********************************************************** */
/* Creates a hex representation of a value into pre-allocated char array. */
void doubleToHex(double v, char* out) { \
  unsigned char *p = (unsigned char*)&v; 
  int i; 
  char tmp[3];
  
  strcpy(out, "0x");
  
  for(i = 0; i < sizeof(double); ++i) {
    sprintf(tmp, "%02X", p[i]);
    strcat(out, tmp);
  }
}
/* *********************************************************** */
void floatToHex(float v, char* out) { \
  unsigned char *p = (unsigned char*)&v; 
  int i; 
  char tmp[3];
  
  strcpy(out, "0x");
  
  for(i = 0; i < sizeof(float); ++i) {
    sprintf(tmp, "%02X", p[i]);
    strcat(out, tmp);
  }
}
/* *********************************************************** */
void intToHex(int v, char* out) { \
  unsigned char *p = (unsigned char*)&v; 
  int i; 
  char tmp[3];
  
  strcpy(out, "0x");
  
  for(i = 0; i < sizeof(int); ++i) {
    sprintf(tmp, "%02X", p[i]);
    strcat(out, tmp);
  }
}
/* *********************************************************** */
void shortToHex(short v, char* out) { \
  unsigned char *p = (unsigned char*)&v; 
  int i; 
  char tmp[3];
  
  strcpy(out, "0x");
  
  for(i = 0; i < sizeof(short); ++i) {
    sprintf(tmp, "%02X", p[i]);
    strcat(out, tmp);
  }
}
/* *********************************************************** */
void byteToHex(int8_t v, char* out) { \
  unsigned char *p = (unsigned char*)&v; 
  int i; 
  char tmp[3];
  
  strcpy(out, "0x");
  
  for(i = 0; i < sizeof(int8_t); ++i) {
    sprintf(tmp, "%02X", p[i]);
    strcat(out, tmp);
  }
}
/* *********************************************************** */
void textToHex(char v, char* out) { \
  unsigned char *p = (unsigned char*)&v; 
  int i; 
  char tmp[3];
  
  strcpy(out, "0x");
  
  for(i = 0; i < sizeof(char); ++i) {
    sprintf(tmp, "%02X", p[i]);
    strcat(out, tmp);
  }
}
/* *********************************************************** */
int 
excludevars(int ncid1, int ncid2, char** finallist,
                int nfinal, char** excludelist, int nexclude)
{
  int nvars;
  char** vars = NULL;
  int status = EXIT_SUCCESS;
  
  vars = NULL;
  
  if (    (nexclude == 0)
    )        
          return EXIT_SUCCESS;

  /* printf("%d: creating temporary var list array.\n", __LINE__); */
  /* get simple difference */
  if(newstringlist(&vars, &nvars, NC_MAX_VARS))
    status = EXIT_FATAL;
  
  /* printf("%d: getting all variables names from both input files.\n", __LINE__); */
  if(allvarnames(vars, nvars, ncid1, ncid2))
    status = EXIT_FATAL; 
  
  /*printf("vars=");
   printstrlist(vars, nvars, stdout);
  */
    
  if(strlistsd(  vars, excludelist, finallist,
                nvars, nexclude,   nfinal))
    status = EXIT_FATAL; 
  
  /*printf("excludelist=");
   printstrlist(excludelist, nexclude, stdout);
  */
  
  /*printf("finallist=");
   printstrlist(finallist, nfinal, stdout);
  */
  
  freestringlist(&vars, nvars);         

  return status;                               
}
/* *********************************************************** */
void handle_error(int status) {
     if (status != NC_NOERR) {
        fprintf(stderr, "%s\n", nc_strerror(status));
        exit(-1);
        }
     }

/* *********************************************************** */
/*
  Mimics incrementing odometer. 
  Returns 0 if rolled-over. 
  
  @param odo: the counters that are updated.
  @param limits: the maximum values allowed per counter.
  @param first: index of first counter to update.
  @param last: index of last counter to update.
  */
int odometer(size_t* odo, size_t* limits, int first, int last)
{
  int i = last;
  while (i >= first) {
    odo[i] += 1;
    if (odo[i] > limits[i]) 
      odo[i] = 0;
    else
      break;
    
    --i;
  }

#ifdef __DEBUG__
  printf("DEBUG : %d : odo = ", __LINE__);
  for(i=first; i <= last; ++i) {
    printf("%d ", odo[i]);
  }
  printf("\n");
#endif

  /* Test for roll over. */
  for(i=first; i <= last; ++i) {
    if (odo[i] != 0)
      return 1;
  }
    
  /* If we get here then rolled over. */
  return 0;
}
/* *********************************************************** */
/* Pretty prints attribute values into a string. 
  Assumes 'str' is pre-allocated large enough to hold output string.
*/
void prettyprintatt(int ncid, char* varname, int varid, char* name, char* str) 
{
  int status, i;
  nc_type type;
  size_t len;
  char* pc;
  int8_t* puc;
  short* ps;
  int* pi;
  float* pf;
  double* pd;
  char tmpstr[32];
        
  status = nc_inq_att(ncid, varid, name, &type, &len);
  if (status != NC_NOERR) {
    if (varid == NC_GLOBAL)
      fprintf(stderr, "ERROR : QUERYING GLOBAL ATTRIBUTE \"%s\"\n", name);
    else
      fprintf(stderr, "ERROR : QUERYING ATTRIBUTE \"%s\" FOR VARIABLE \"%s\"\n", name, varname);
    return;          
  }
  
  str[0] = '\0';
  if (len < 1) {
    return;
  }
  
  switch(type) 
  {
  case NC_BYTE:
    puc = XMALLOC(int8_t,len);
    
    status = nc_get_att(ncid, varid, name, puc);
    if (status != NC_NOERR)
    {
      XFREE(puc); return;
    }
    
    for(i=0; i < (int)len; ++i) {
      byteToString(puc[i], str+2*i, NULL);
      str[2*i+1] = ',';
    }
      
    XFREE(puc); 
    str[2*(int)len-1] = '\0';
    break;
  
  case NC_CHAR:
    pc = XMALLOC(char,len);
    status = nc_get_att_text(ncid, varid, name, pc);
    if (status != NC_NOERR)
    {
      XFREE(pc); return;
    }

    for(i=0; i < (int)len; ++i)
      textToString(pc[i], str+i, NULL);

    XFREE(pc); 
    str[(int)len] = '\0';
    break;
    
  case NC_SHORT:
    ps = XMALLOC(short,len);
    status = nc_get_att_short(ncid, varid, name, ps);
    if (status != NC_NOERR)
    {
      XFREE(ps); return;
    }

    for(i=0; i < (int)len; ++i) {
      sprintf(tmpstr, "%d,", ps[i]);
      strcat(str, tmpstr);
    }
    
    XFREE(ps); 
    str[strlen(str)-1] = '\0'; // Remove last comma.
    break;        

  case NC_INT:
    pi = XMALLOC(int,len);
    status = nc_get_att_int(ncid, varid, name, pi);
    if (status != NC_NOERR)
    {
      XFREE(pi); return;
    }

    for(i=0; i < (int)len; ++i) {
      sprintf(tmpstr, "%d,", pi[i]);
      strcat(str, tmpstr);
    }
    
    XFREE(pi); 
    str[strlen(str)-1] = '\0'; // Remove last comma.
    break;        

  case NC_FLOAT:
    pf = XMALLOC(float,len);
    status = nc_get_att_float(ncid, varid, name, pf);
    if (status != NC_NOERR)
    {
      XFREE(pf); return;
    }

    for(i=0; i < (int)len; ++i) {
      sprintf(tmpstr, "%.9g,", pf[i]);
      strcat(str, tmpstr);
    }
    
    XFREE(pf); 
    str[strlen(str)-1] = '\0'; // Remove last comma.
    break;        

  case NC_DOUBLE:
    pd = XMALLOC(double,len);
    status = nc_get_att_double(ncid, varid, name, pd);
    if (status != NC_NOERR)
    {
      XFREE(pd); return;
    }

    for(i=0; i < (int)len; ++i) {
      sprintf(tmpstr, "%.17g,", pd[i]);
      strcat(str, tmpstr);
    }
    
    XFREE(pd); 
    str[strlen(str)-1] = '\0'; // Remove last comma.
    break;        
  }
}
/* *********************************************************** */
int cmpatt(int ncid1, int ncid2, int varid1, int varid2,
           char* name, char* varname, nccmpopts* opts)
{
  int ncstatus, status, attid1, attid2;
  nc_type type1, type2;
  size_t lenp1, lenp2;
  char typestr1[256];
  char typestr2[256];

  status = EXIT_SUCCESS;
  ncstatus = nc_inq_att(ncid1, varid1, name, &type1, &lenp1);
  if (ncstatus != NC_NOERR) {
    fprintf(stderr, "DIFFER : VARIABLE \"%s\" IS MISSING ATTRIBUTE WITH NAME \"%s\" IN FILE \"%s\"\n", varname, name, opts->file1);
    
    if (!opts->warn[NCCMP_W_ALL]) 
      status = EXIT_DIFFER;

    return status;          
  }
  
  ncstatus = nc_inq_att(ncid2, varid2, name, &type2, &lenp2);
  if (ncstatus != NC_NOERR) {
    fprintf(stderr, "DIFFER : VARIABLE \"%s\" IS MISSING ATTRIBUTE WITH NAME \"%s\" IN FILE \"%s\"\n", varname, name, opts->file2);
    
    if (!opts->warn[NCCMP_W_ALL]) 
      status = EXIT_DIFFER;

    return status;          
  }

  if (type1 != type2) {
    type2string(type1, typestr1);
    type2string(type2, typestr2);                                                
    fprintf(stderr, "DIFFER : TYPES : ATTRIBUTE : %s : VARIABLE : %s : %s <> %s\n", name, varname, typestr1, typestr2);
            
    if (!opts->warn[NCCMP_W_ALL]) 
      status = EXIT_DIFFER;

    return status;
  }

  if (lenp1 != lenp2) {
    prettyprintatt(ncid1, varname, varid1, name, typestr1);
    prettyprintatt(ncid2, varname, varid2, name, typestr2);
    
    fprintf(stderr, "DIFFER : LENGTHS : ATTRIBUTE : %s : VARIABLE : %s : %lu <> %lu : VALUES : ", name, varname, (unsigned long)lenp1, (unsigned long)lenp2);
    
    switch(type1) {
      case NC_CHAR:
        /* Quote strings. */
        fprintf(stderr, "\"%s\" : \"%s\"\n", typestr1, typestr2);
        if (strcmp(typestr1,typestr2) == 0) {
          /* Same text, but invisible trailing nulls because lengths differ. */
          if (opts->warn[NCCMP_W_EOS] || opts->warn[NCCMP_W_ALL]) {
            /* Pass */
          } else {
            status = EXIT_DIFFER;
            return status;
          }
        }
        break;
      default:
        /* Unquoted. */
        fprintf(stderr, "%s : %s\n", typestr1, typestr2);
        if (!opts->warn[NCCMP_W_ALL]) {
          status = EXIT_DIFFER;
          return status;
        }
        break;
    } 
  }

  if (cmpattval(ncid1, ncid2, varid1, varid2, name, lenp1, type1) != NC_NOERR) {
    prettyprintatt(ncid1, varname, varid1, name, typestr1);
    prettyprintatt(ncid2, varname, varid2, name, typestr2);
    fprintf(stderr, "DIFFER : VARIABLE : %s : ATTRIBUTE : %s : VALUES : ", varname, name);

    switch(type1) {
      case NC_CHAR:
        /* Quote strings. */
        fprintf(stderr, "\"%s\" <> \"%s\"\n", typestr1, typestr2);
        break;
      default:
        /* Unquoted. */
        fprintf(stderr, "%s <> %s\n", typestr1, typestr2);
        break;
    }
    
    if (!opts->warn[NCCMP_W_ALL]) {
      status = EXIT_DIFFER;
      return status;
    }
  }
  
  return EXIT_SUCCESS;
}


/* *********************************************************** */
/* Assumes that both attributes have same length.
 */
int cmpattval(int nc1, int nc2, int varid1, int varid2, char* name, int len, nc_type type)
{
  char* c1;
  int8_t* uc1;
  short* s1;
  int* i1;
  float* f1;
  double* d1;
  char* c2;
  int8_t* uc2;
  short* s2;
  int* i2;
  float* f2;
  double* d2;
  int status;
  int i;
 
  if (name == NULL) return NC_EINVAL;
 
  switch(type)
  {
  case NC_BYTE:
          uc1 = XMALLOC(int8_t,len);
          uc2 = XMALLOC(int8_t,len);
          status = nc_get_att(nc1, varid1, name, uc1);
          if (status != NC_NOERR)
          {
                  XFREE(uc1); XFREE(uc2); return NC_EINVAL;
          }
          status = nc_get_att(nc2, varid2, name, uc2);
          if (status != NC_NOERR)
          {
                  XFREE(uc1); XFREE(uc2); return NC_EINVAL;
          }
          for(i = 0; i < len; ++i)
                    {
                            if (uc1[i] != uc2[i])
                            {
                                    XFREE(uc1); XFREE(uc2); return EXIT_DIFFER;
                            }
                    }
         XFREE(uc1); 
         XFREE(uc2); 
  break;
  case NC_CHAR:
          c1 = XMALLOC(char,len);
          c2 = XMALLOC(char,len);                
          status = nc_get_att_text(nc1, varid1, name, c1);
          if (status != NC_NOERR)
          {
                  XFREE(c1); XFREE(c2); return NC_EINVAL;
          }
          status = nc_get_att_text(nc2, varid2, name, c2);
          if (status != NC_NOERR)
          {
                  XFREE(c1); XFREE(c2); return NC_EINVAL;
          }
          if (strncmp(c1,c2,len) != 0)
          {
                  XFREE(c1); XFREE(c2); return EXIT_DIFFER;
          }
         XFREE(c1); 
         XFREE(c2); 
  break;
  case NC_SHORT:
          s1 = XMALLOC(short,len);
          s2 = XMALLOC(short,len);                
          status = nc_get_att_short(nc1, varid1, name, s1);
          if (status != NC_NOERR)
          {
                  XFREE(s1); XFREE(s2); return NC_EINVAL;
          }
          status = nc_get_att_short(nc2, varid2, name, s2);
          if (status != NC_NOERR)
          {
                  XFREE(s1); XFREE(s2); return NC_EINVAL;
          }
          for(i = 0; i < len; ++i)
          {
                  if (s1[i] != s2[i])
                  {
                          XFREE(s1); XFREE(s2); return EXIT_DIFFER;
                  }
          }
         XFREE(s1); 
         XFREE(s2); 
  break;        
  case NC_INT:
          i1 = XMALLOC(int,len);
          i2 = XMALLOC(int,len);                
          status = nc_get_att_int(nc1, varid1, name, i1);
          if (status != NC_NOERR)
          {
                  XFREE(i1); XFREE(i2); return NC_EINVAL;
          }
          status = nc_get_att_int(nc2, varid2, name, i2);
          if (status != NC_NOERR)
          {
                  XFREE(i1); XFREE(i2); return NC_EINVAL;
          }
          for(i = 0; i < len; ++i)
          {
                  if (i1[i] != i2[i])
                  {
                          XFREE(i1); XFREE(i2); return EXIT_DIFFER;
                  }
          }
         XFREE(i1); 
         XFREE(i2); 
  break;        
  case NC_FLOAT:
          f1 = XMALLOC(float,len);
          f2 = XMALLOC(float,len);                
          status = nc_get_att_float(nc1, varid1, name, f1);
          if (status != NC_NOERR)
          {
                  XFREE(f1); XFREE(f2); return NC_EINVAL;
          }
          status = nc_get_att_float(nc2, varid2, name, f2);
          if (status != NC_NOERR)
          {
                  XFREE(f1); XFREE(f2); return NC_EINVAL;
          }
          for(i = 0; i < len; ++i)
          {
                  if (f1[i] != f2[i])
                  {
                          XFREE(f1); XFREE(f2); return EXIT_DIFFER;
                  }
          }
         XFREE(f1); 
         XFREE(f2); 
  break;        
  case NC_DOUBLE:
          d1 = XMALLOC(double,len);
          d2 = XMALLOC(double,len);                
          status = nc_get_att_double(nc1, varid1, name, d1);
          if (status != NC_NOERR)
          {
                  XFREE(d1); XFREE(d2); return NC_EINVAL;
          }
          status = nc_get_att_double(nc2, varid2, name, d2);
          if (status != NC_NOERR)
          {
                  XFREE(d1); XFREE(d2); return NC_EINVAL;
          }
          for(i = 0; i < len; ++i)
          {
                  if (d1[i] != d2[i])
                  {
                          XFREE(d1); XFREE(d2); return EXIT_DIFFER;
                  }
          }
         XFREE(d1); 
         XFREE(d2); 
  break;        
  }
  
  return EXIT_SUCCESS;
}
/* *********************************************************** */
void type2string(nc_type type, char* str)
{
        switch(type)
        {
        case NC_BYTE:
          strcpy(str, "BYTE");
          break;
        case NC_CHAR:
          strcpy(str, "CHAR");
          break;
        case NC_SHORT:      
          strcpy(str, "SHORT");
          break;        
        case NC_INT:  
          strcpy(str, "INT");  
          break;        
        case NC_FLOAT:     
          strcpy(str, "FLOAT");
          break;        
        case NC_DOUBLE:
          strcpy(str, "DOUBLE");
          break;
        default:
          strcpy(str, "");
          break;
        }                 
} 
/* *********************************************************** */
int openfiles(nccmpopts* opts, int *ncid1, int *ncid2)
{
  int status;

  status = nc_open(opts->file1, NC_NOWRITE, ncid1);
  handle_error(status);

  status = nc_open(opts->file2, NC_NOWRITE, ncid2);
  handle_error(status);

  return 0;
}
/* *********************************************************** */
/* Compares record names and lengths. */
int
nccmprecinfo(nccmpopts* opts, int ncid1, int ncid2)
{
  char name1[256], name2[256];
  int status;

  status = EXIT_SUCCESS;
  
  if (opts->verbose)
    printf("INFO: Comparing record information.\n");

  status = nc_inq_unlimdim(ncid1, &recid1);
  handle_error(status);

  if (recid1 != -1) {
    status = nc_inq_dim(ncid1, recid1, name1, &nrec1);
    handle_error(status);
  } else {
    nrec1 = 0;
  }

  status = nc_inq_unlimdim(ncid2, &recid2);
  handle_error(status);

  if (recid2 != -1) {
    status = nc_inq_dim(ncid2, recid2, name2, &nrec2);
    handle_error(status);
  } else {
    nrec2 = 0;
  }

  if (instringlist(opts->excludelist, name1, opts->nexclude) ||
      instringlist(opts->excludelist, name2, opts->nexclude) ||
      !instringlist(opts->variablelist, name1, opts->nvariable) ||
      !instringlist(opts->variablelist, name2, opts->nvariable))
    return EXIT_SUCCESS;

  if (strcmp(name1, name2)) {
    fprintf(stderr, "DIFFER : NAMES OF RECORDS : %s <> %s\n", name1, name2);
    if (!opts->warn[NCCMP_W_ALL]) 
        status = EXIT_DIFFER;

    if(!opts->force) return status;
  }

  if (nrec1 != nrec2) {
    fprintf(stderr, "DIFFER : LENGTHS OF RECORDS : %s (%d) <> %s (%d)\n", name1, (int)nrec1, name2, (int)nrec2);
    if (!opts->warn[NCCMP_W_ALL]) 
        status = EXIT_DIFFER;

    if(!opts->force) return status;
  }

  return status;
}

int getgroupinfo (int ncid, int numgroups,GROUP_NODE *groups) {
    int *gids;
    int res, i;

    gids = malloc(sizeof(int) * numgroups);
    res = nc_inq_grps(ncid, NULL, gids);

    for (i=0; i<numgroups; i++){
        groups[i].groupID = gids[i];
        nc_inq_grpname(gids[i],groups[i].groupName);
    }
    return numgroups;
}

/* *********************************************************** */
/* Get dim info for file. */
void getdiminfo(int ncid, dimstruct* dims, int* ndims)
{
  int status, i;

  status = nc_inq_ndims(ncid, ndims);
  handle_error(status);

  // Query all dimids, which may not be from 0..N-1 in a HDF file.
  int dimids[(int)NC_MAX_DIMS];
  int include_parents = 1;
  status = nc_inq_dimids(ncid, ndims, dimids, include_parents);
  handle_error(status);

  for(i=0; i < *ndims; ++i) {
    dims[i].dimid = dimids[i];
    status = nc_inq_dim(ncid, dimids[i], dims[i].name, &dims[i].len);
    handle_error(status);
  }
}
/* *********************************************************** */
/* Copy attribute type to same type as var, just in case of mismatch. */
void broadcast_missing(nc_type var_type, nc_type att_type, missing_struct *values) {
  #define BROADCAST_MISSING(T) { \
    switch(att_type) { \
      case NC_CHAR: values->T = values->c; break; \
      case NC_BYTE: values->T = values->b; break; \
      case NC_SHORT: values->T = values->s; break; \
      case NC_INT: values->T = values->i; break; \
      case NC_FLOAT: values->T = values->f; break; \
      case NC_DOUBLE: values->T = values->d; break; \
    } \
  }
  
  switch(var_type) {
    case NC_CHAR: BROADCAST_MISSING(c); break;
    case NC_BYTE: BROADCAST_MISSING(b); break;
    case NC_SHORT: BROADCAST_MISSING(s); break;
    case NC_INT: BROADCAST_MISSING(i); break;
    case NC_FLOAT: BROADCAST_MISSING(f); break;
    case NC_DOUBLE: BROADCAST_MISSING(d); break;
  }
}
/* *********************************************************** */
char get_missing(int ncid, varstruct * var, const char* attname) {
  nc_type att_type;
  int status;
  
  status = nc_inq_atttype(ncid, var->varid, attname, &att_type);
  if (status != NC_NOERR) return 0;
  
  var->missing.type = att_type;
  
  switch(att_type) {
    case NC_CHAR:
      status = nc_get_att_text(ncid, var->varid, attname, &var->missing.c);
      if (status != NC_NOERR) return 0;
      break;
    case NC_BYTE: 
      status = nc_get_att(ncid, var->varid, attname, &var->missing.b);
      if (status != NC_NOERR) return 0;
      break;
    case NC_SHORT:
      status = nc_get_att_short(ncid, var->varid, attname, &var->missing.s);
      if (status != NC_NOERR) return 0;
      break;
    case NC_INT:
      status = nc_get_att_int(ncid, var->varid, attname, &var->missing.i);
      if (status != NC_NOERR) return 0;
      break;
    case NC_FLOAT:
      status = nc_get_att_float(ncid, var->varid, attname, &var->missing.f);
      if (status != NC_NOERR) return 0;
      break;
    case NC_DOUBLE:
      status = nc_get_att_double(ncid, var->varid, attname, &var->missing.d);
      if (status != NC_NOERR) return 0;
      break;
    default: return 0;
  }
  
  var->hasmissing = 1;
  broadcast_missing(var->type, att_type, &var->missing);
  
  return 1;
}
/* *****************************s****************************** */
/* Read all the file's metadata for variables. */
void getvarinfo(int ncid, varstruct* vars, int* nvars, int verbose)
{
  int status, i, j, recid;
  char name[NC_MAX_NAME];
  nc_type att_type;
  
  status = nc_inq_nvars(ncid, nvars);
  handle_error(status);

  status = nc_inq_unlimdim(ncid, &recid);
  handle_error(status);

  for(i=0; i < *nvars; ++i) {
    vars[i].varid = i;

    status = nc_inq_var(ncid, i, vars[i].name, &vars[i].type,
                        &vars[i].ndims, vars[i].dimids, &vars[i].natts);
    handle_error(status);

    vars[i].len = 1;
    for(j=0; j < vars[i].ndims; ++j) {
      status = nc_inq_dimlen(ncid, vars[i].dimids[j], &vars[i].dimlens[j]);
      handle_error(status);
#ifdef __DEBUG__
      printf("DEBUG : %d : %s dimid %d, len %d\n", __LINE__, vars[i].name, j, vars[i].dimlens[j]);
#endif
      vars[i].len *= vars[i].dimlens[j];
    }

    vars[i].hasrec = (vars[i].dimids[0] == recid);
    
    /* Get missing_value or _FillValue. */
    if (get_missing(ncid, &vars[i], "missing_value") == 0)
      get_missing(ncid, &vars[i], "_FillValue");
    
    if (verbose) {
      if (vars[i].hasmissing) {
        printf("INFO: \"%s\" missing value: ", vars[i].name);
        switch(vars[i].missing.type) {
          case NC_BYTE: printf("%d (byte)\n", vars[i].missing.b); break;
          case NC_CHAR: printf("%d (char)\n", vars[i].missing.c); break;
          case NC_SHORT: printf("%d (short)\n", vars[i].missing.s); break;
          case NC_INT: printf("%d (int)\n", vars[i].missing.i); break;
          case NC_FLOAT: printf("%g (float)\n", vars[i].missing.f); break;
          case NC_DOUBLE: printf("%g (double)\n", vars[i].missing.d); break;
        }
      }
    }
  }
}
/* *********************************************************** */
/* Returns index to varstruct in list, otherwise -1. */
int isinvarstructlist(char* name, varstruct* vars, int nvars) 
{
  int i;

  for(i=0; i < nvars; ++i) {
    if (strcmp(name, vars[i].name) == 0) 
      return i;
  }

  return -1;
}
/* *********************************************************** */
/* Get vars to use to do cmp based on input and exclusion lists. */
int makecmpvarlist(nccmpopts* opts, int ncid1, int ncid2)
{
  int status, i;

  newstringlist(&opts->cmpvarlist, &opts->ncmpvarlist, (int)NC_MAX_VARS);
  if(opts->variable) {
    if (opts->verbose)
      printf("INFO: Using variables provided in list.\n");

    status = strlistu(opts->variablelist, opts->cmpvarlist, opts->cmpvarlist, 
                      opts->nvariable, opts->ncmpvarlist, opts->ncmpvarlist);
  } else if (opts->exclude) {
    if (opts->verbose)
      printf("INFO: Excluding variables in provided list.\n");

    status = excludevars(ncid1, ncid2, opts->cmpvarlist, opts->ncmpvarlist, opts->excludelist, opts->nexclude);
  } else {
    if (opts->verbose)
      printf("INFO: Using all variables.\n");

    status = allvarnames(opts->cmpvarlist, opts->ncmpvarlist, ncid1, ncid2);
  }

  opts->ncmpvarlist = getnumstrlist(opts->cmpvarlist, (int)NC_MAX_VARS);
  if (opts->verbose) {
    printf("INFO: Variables to compare (%d):\n", opts->ncmpvarlist);
    for(i=0; i < opts->ncmpvarlist-1; ++i)
      printf("%s, ", opts->cmpvarlist[i]);
    
    if (opts->ncmpvarlist)
      printf("%s\n", opts->cmpvarlist[i]);
  }
  
  return status;
}
/* *********************************************************** */
/* Gets list of all variable names in both input files. */
int 
allvarnames(char** list, int nvars, int ncid1, int ncid2)
{
  char** tmplist = NULL;
  int ntmp;

  newstringlist(&tmplist, &ntmp, NC_MAX_VARS);

  if( ncallvars(ncid1, tmplist, ntmp) )
      return 1;
  
  /* printf("%d: ncallvars returned.\n", __LINE__); */
  
  if( strlistu(tmplist, list, list, ntmp, nvars, nvars) )
      return 1;

  /* printf("%d: Variable names from file 1 collected.\n", __LINE__); */
  clearstringlist(tmplist, NC_MAX_VARS);
      
  if( ncallvars(ncid2, tmplist, ntmp) )
      return 1;

  if( strlistu(tmplist, list, list, ntmp, nvars, nvars) ) 
      return 1;

  /* printf("%d: Variable names from file 2 collected.\n", __LINE__); */
  freestringlist(&tmplist, ntmp);        

  return EXIT_SUCCESS;
}
/* *********************************************************** */

int 
nccmpformats(nccmpopts* opts, int ncid1, int ncid2)
{
  int status, fmt1, fmt2;

  status = nc_inq_format(ncid1, &fmt1);
  handle_error(status);

  status = nc_inq_format(ncid2, &fmt2);
  handle_error(status);
  
  if (fmt1 != fmt2) {
    fprintf(stderr, "DIFFER : FILE FORMATS : %s <> %s\n", NCFORMATSTR(fmt1), NCFORMATSTR(fmt2));
    
    if (!opts->warn[NCCMP_W_ALL] &&
        !opts->warn[NCCMP_W_FORMAT])
      return EXIT_DIFFER;
  }
  
  return EXIT_SUCCESS;
}
/* *********************************************************** */

int 
nccmpglobalatts(nccmpopts* opts, int ncid1, int ncid2)
{
  int ngatts1, ngatts2, nattsex1, nattsex2, i, status, status2, attid1, attid2, nprocessedatts;
  nc_type type1, type2;
  size_t len1, len2;
  char name1[NC_MAX_NAME], name2[NC_MAX_NAME];
  char** processedatts = NULL;
  char typestr1[1024], typestr2[1024];
  
  status = EXIT_SUCCESS;
  status2 = status;
  if (!opts->global)
    return status;

  if (opts->history == 0) 
    appendstringtolist(&opts->globalexclude, "history", &opts->nglobalexclude);
  
  /*  
  printf("globalexclude =");
  printstrlist(opts->globalexclude, opts->nglobalexclude, stdout);
  */
    
  /* Number of global atts to compare with exclusion taken into account. */
  nattsex1 = 0;
  nattsex2 = 0;
  
  status = nc_inq_natts(ncid1, &ngatts1);
  handle_error(status);
  
  status = nc_inq_natts(ncid2, &ngatts2);
  handle_error(status);
  
  for(i = 0; i < ngatts1; ++i) {
    attid1 = i;
    status = nc_inq_attname(ncid1, NC_GLOBAL, attid1, name1);
    handle_error(status);
    
    if (!instringlist(opts->globalexclude, name1, opts->nglobalexclude)) {
      ++nattsex1;
    }
  }
  
  for(i = 0; i < ngatts2; ++i) {
    attid2 = i;
    status = nc_inq_attname(ncid2, NC_GLOBAL, attid2, name2);
    handle_error(status);
    
    if (!instringlist(opts->globalexclude, name2, opts->nglobalexclude)) {
      ++nattsex2;
    }
  }
  
  if(nattsex1 != nattsex2) {
    fprintf(stderr, "DIFFER : NUMBER OF GLOBAL ATTRIBUTES : %d <> %d\n", nattsex1, nattsex2);
    
    if (!opts->warn[NCCMP_W_ALL]) 
      status2 = EXIT_DIFFER;
      
    if(!opts->force) return status2;
  }
    
  if(newstringlist(&processedatts, &i, NC_MAX_VARS)) {
    fprintf(stderr, "ERROR: Failed to allocated string list for comparing  global attributes.\n");
    return EXIT_FATAL;
  }
  
  for(i = 0; i < ngatts1; ++i) {
    attid1 = i;
    status = nc_inq_attname(ncid1, NC_GLOBAL, attid1, name1);
    handle_error(status);
    
    /* Log that this gatt was processed. */
    addstringtolist(processedatts, name1, NC_MAX_VARS);
    
    if (instringlist(opts->globalexclude, name1, opts->nglobalexclude)) 
      continue;
      
    status = nc_inq_att(ncid1, NC_GLOBAL, name1, &type1, &len1);
    if (status != NC_NOERR) {
      fprintf(stderr, "Query failed on global attribute in %s\n", opts->file1);
      if (!opts->warn[NCCMP_W_ALL]) 
        status2 = EXIT_DIFFER;
        
      if(opts->force) continue; else return status2;
    }
    
    status = nc_inq_att(ncid2, NC_GLOBAL, name1, &type2, &len2);
    if (status != NC_NOERR) {
      fprintf(stderr, "DIFFER : NAME OF GLOBAL ATTRIBUTE : %s : GLOBAL ATTRIBUTE DOESN'T EXIST IN \"%s\"\n", name1, opts->file2);
      if (!opts->warn[NCCMP_W_ALL]) 
        status2 = EXIT_DIFFER;
        
      if(opts->force) continue; else return status2;
    }
    
    if (type1 != type2)  {
      type2string(type1, typestr1);
      type2string(type2, typestr2);     
      fprintf(stderr, "DIFFER : GLOBAL ATTRIBUTE TYPES : %s : %s <> %s\n", name1, typestr1, typestr2);
      if (!opts->warn[NCCMP_W_ALL]) 
        status2 = EXIT_DIFFER;
        
      if(opts->force) continue; else return status2;
    }
    
    if (len1 != len2)   {
      prettyprintatt(ncid1, NULL, NC_GLOBAL, name1, typestr1);
      prettyprintatt(ncid2, NULL, NC_GLOBAL, name1, typestr2);
      fprintf(stderr, "DIFFER : LENGTHS OF GLOBAL ATTRIBUTE : %s : %lu <> %lu : VALUES : %s <> %s\n", name1, 
              (unsigned long)len1, (unsigned long)len2, typestr1, typestr2);
      if (!opts->warn[NCCMP_W_ALL]) 
        status2 = EXIT_DIFFER;
        
      if(opts->force) continue; else return status2;
    }
    
    if (cmpattval(ncid1,ncid2,NC_GLOBAL,NC_GLOBAL,name1,len1,type1) != NC_NOERR) {
      /* Pretty print values. */
      prettyprintatt(ncid1, NULL, NC_GLOBAL, name1, typestr1);
      prettyprintatt(ncid2, NULL, NC_GLOBAL, name1, typestr2);
      fprintf(stderr, "DIFFER : VALUES OF GLOBAL ATTRIBUTE : %s : %s <> %s\n", name1, typestr1, typestr2);
      if (!opts->warn[NCCMP_W_ALL]) 
        status2 = EXIT_DIFFER;
        
      if(opts->force) continue; else return status2;
    }
  }
  
  for(i = 0; i < ngatts2; ++i) {
    attid2 = i;
    status = nc_inq_attname(ncid2, NC_GLOBAL, attid2, name2);
    if (status != NC_NOERR) {
      fprintf(stderr, "Query failed for global attribute name in %s\n", opts->file2);
      if (!opts->warn[NCCMP_W_ALL]) 
        status2 = EXIT_DIFFER;
        
      if(opts->force) continue; else return status2;
    }
    
    /* Skip if already processed (or excluded). */
    if (instringlist(processedatts, name2, NC_MAX_VARS))
      continue;
    
    /* Log that this att was processed. */
    addstringtolist(processedatts, name2, NC_MAX_VARS);
    
    status = nc_inq_att(ncid2, NC_GLOBAL, name2, &type2, &len2);
    if (status != NC_NOERR) {
      fprintf(stderr, "Query failed on global attribute in %s\n", opts->file2);
      if (!opts->warn[NCCMP_W_ALL]) 
        status2 = EXIT_DIFFER;
        
      if(opts->force) continue; else return status2;
    }
    
    status = nc_inq_att(ncid1, NC_GLOBAL, name2, &type1, &len1);
    if (status != NC_NOERR) {
      fprintf(stderr, "DIFFER : NAME OF GLOBAL ATTRIBUTE : %s : GLOBAL ATTRIBUTE DOESN'T EXIST IN %s\n", name2, opts->file1);
      if (!opts->warn[NCCMP_W_ALL]) 
        status2 = EXIT_DIFFER;
        
      if(opts->force) continue; else return status2;
    }
    
    if (type1 != type2) {
      type2string(type1, typestr1);
      type2string(type2, typestr2);
      fprintf(stderr, "DIFFER : GLOBAL ATTRIBUTE TYPE : %s : %s <> %s\n", name1, typestr1, typestr2);
      if (!opts->warn[NCCMP_W_ALL]) 
        status2 = EXIT_DIFFER;
        
      if(opts->force) continue; else return status2;
    }
    
    if (len1 != len2)  {
      prettyprintatt(ncid1, NULL, NC_GLOBAL, name1, typestr1);
      prettyprintatt(ncid2, NULL, NC_GLOBAL, name1, typestr2);
      
      fprintf(stderr, "DIFFER : LENGTHS OF GLOBAL ATTRIBUTE : %s : %lu <> %lu : VALUES : ", name1, (unsigned long)len1, (unsigned long)len2);

      switch(type1) {
        case NC_CHAR:
          /* Quote strings. */
          fprintf(stderr, "\"%s\" : \"%s\"\n", typestr1, typestr2);
          if (strcmp(typestr1,typestr2) == 0) {
            /* Same text, but invisible trailing nulls because lengths differ. */
            if (opts->warn[NCCMP_W_EOS] || opts->warn[NCCMP_W_ALL]) {
              /* Pass */
            } else {
              status2 = EXIT_DIFFER;
              if(opts->force) continue; else return status2;
            }
          }
          break;
        default:
          /* No quotes. */
          fprintf(stderr, "%s : %s\n", typestr1, typestr2);
          if (!opts->warn[NCCMP_W_ALL]) {
            status2 = EXIT_DIFFER;
            if(opts->force) continue; else return status2;
          }
          break;
      } 
    }
    
    if (cmpattval(ncid1,ncid2,NC_GLOBAL,NC_GLOBAL,name1,len1,type1) != NC_NOERR) {
      /* Pretty print values. */
      prettyprintatt(ncid1, NULL, NC_GLOBAL, name1, typestr1);
      prettyprintatt(ncid2, NULL, NC_GLOBAL, name1, typestr2);
      fprintf(stderr, "DIFFER : VALUES OF GLOBAL ATTRIBUTE : %s : %s <> %s\n", name1, typestr1, typestr2);
      if (!opts->warn[NCCMP_W_ALL]) 
        status2 = EXIT_DIFFER;
        
      if(opts->force) continue; else return status2;
    }
  }
    
  /* Clear the list. */
  freestringlist(&processedatts, NC_MAX_VARS);
  processedatts = NULL;

  return status2;
}
/* *********************************************************** */
int 
nccmpmetadata(nccmpopts *opts, int ncid1, int ncid2)
{
  int i, j, j1, j2, status, ncstatus, dimid1, dimid2, tmp1, tmp2, attid1, attid2, natts1, natts2;
  size_t len1, len2;
  char name1[NC_MAX_NAME], name2[NC_MAX_NAME], recname1[NC_MAX_NAME], recname2[NC_MAX_NAME], typestr1[1024], typestr2[1024];
  char** processedatts = NULL;

  status = EXIT_SUCCESS;

  if (opts->verbose)
    printf("INFO: Comparing metadata.\n");

  if (opts->verbose)
    printf("INFO: Comparing number of dimensions.\n");

  if (ndims1 != ndims2) {
    fprintf(stderr, "DIFFER : NUMBER OF DIMENSIONS IN FILES : %d <> %d\n", ndims1, ndims2);
    if (!opts->warn[NCCMP_W_ALL]) 
      status = EXIT_DIFFER;

    if (!opts->force)
      return status;
  }

  if (opts->verbose)
    printf("INFO: Getting record dimension names, if they exist.\n");

  if (recid1 != -1) {
    ncstatus = nc_inq_dimname(ncid1, recid1, recname1);
    handle_error(ncstatus);
  } else
    strcpy(recname1, "");
  
  if (recid2 != -1) {
    ncstatus = nc_inq_dimname(ncid2, recid2, recname2);
    handle_error(ncstatus);
  } else
    strcpy(recname1, "");
  
  /* Dimensions */
  if (opts->verbose)
    printf("INFO: Comparing dimension lengths.\n");

  for(i = 0; i < ndims1; ++i) {
    dimid1 = dims1[i].dimid;
    ncstatus = nc_inq_dim(ncid1, dimid1, name1, &len1);
    if (ncstatus != NC_NOERR) {
      if (!opts->warn[NCCMP_W_ALL]) 
        status = EXIT_DIFFER;
        
      fprintf(stderr, "Failed to query dimension id %d in file %s.\n", dimid1, opts->file1);
      if(opts->force) continue; else return status;
    }

    ncstatus = nc_inq_dimid(ncid2, name1, &dimid2);
    if (ncstatus != NC_NOERR) {
      fprintf(stderr, "DIFFER : NAME : DIMENSION : %s : DIMENSION DOESN'T EXIST IN \"%s\"\n", name1, opts->file2);
      
      if (!opts->warn[NCCMP_W_ALL]) 
        status = EXIT_DIFFER;
      
      if(opts->force) continue; else return status;
    }

    ncstatus = nc_inq_dim(ncid2, dimid2, name2, &len2);
    if (ncstatus != NC_NOERR) {
      if (!opts->warn[NCCMP_W_ALL]) 
        status = EXIT_DIFFER;
        
      fprintf(stderr, "Failed to query dimension \"%s\" in file \"%s\".\n", name1, opts->file2);
      if(opts->force) continue; else return status;
    }

    if (len1 != len2) {
      fprintf(stderr, "DIFFER : LENGTHS : DIMENSION : %s : %lu <> %lu\n", name1, 
              (unsigned long)len1, (unsigned long)len2);
              
      if (!opts->warn[NCCMP_W_ALL]) 
        status = EXIT_DIFFER;
      
      if(opts->force) continue; else return status;
    }
  }

  /* Variables */
  if (opts->verbose)
    printf("INFO: Comparing variable datatypes and rank.\n");

  for(i=0; i < opts->ncmpvarlist; ++i) {
    j1 = findvar(opts->cmpvarlist[i], vars1);
    if (j1 == -1) {
      fprintf(stderr, "DIFFER : NAME : VARIABLE : %s : VARIABLE DOESN'T EXIST IN \"%s\"\n",   opts->cmpvarlist[i], opts->file1);
      
      if (!opts->warn[NCCMP_W_ALL]) 
        status = EXIT_DIFFER;                                         
        
      if(opts->force) continue; else goto recover;
    }

    j2 = findvar(opts->cmpvarlist[i], vars2);
    if (j2 == -1) {
      fprintf(stderr, "DIFFER : NAME : VARIABLE : %s : VARIABLE DOESN'T EXIST IN \"%s\"\n",   opts->cmpvarlist[i], opts->file2);
      
      if (!opts->warn[NCCMP_W_ALL]) 
        status = EXIT_DIFFER;                                         
      
      if(opts->force) continue; else goto recover;
    }

    if(vars1[j1].type != vars2[j2].type) {
      type2string(vars1[j1].type, typestr1);
      type2string(vars2[j2].type, typestr2);                        
      fprintf(stderr, "DIFFER : TYPES : VARIABLE : %s : %s <> %s\n", opts->cmpvarlist[i], typestr1, typestr2);
      
      if (!opts->warn[NCCMP_W_ALL]) 
        status = EXIT_DIFFER;             
                                  
      if(opts->force) continue; else goto recover;
    }

    if(vars1[j1].ndims != vars2[j2].ndims) {
      fprintf(stderr, "DIFFER : NUMBER : DIMENSIONS : VARIABLE : %s : %d <> %d\n", opts->cmpvarlist[i], ndims1, ndims2);
      
      if (!opts->warn[NCCMP_W_ALL]) 
        status = EXIT_DIFFER;                                         
        
      if(opts->force) continue; else goto recover;
    }
  }
  
  /*printf("DEBUG : %d : \n", __LINE__);*/
  if (opts->verbose)
    printf("INFO: Comparing variables' dimension names.\n");

  for(i=0; i < opts->ncmpvarlist; ++i) {
    j1 = findvar(opts->cmpvarlist[i], vars1);
    j2 = findvar(opts->cmpvarlist[i], vars2);
    
    if ((j1 == -1) || (j2 == -1))
      continue;
    
    /* dimensions */
    for(j = 0; j < vars1[j1].ndims; ++j) {
      dimid1 = vars1[j1].dimids[j];
      
      if (j < vars2[j2].ndims)
        dimid2 = vars2[j2].dimids[j];
      else
        break;

      /*printf("DEBUG : %d : %s,  %s, %s\n", __LINE__, opts->cmpvarlist[i], dims1[dimid1].name, dims2[dimid2].name);*/
      
      if (strcmp(dims1[dimid1].name, dims2[dimid2].name) != 0) {
        fprintf(stderr, "DIFFER : DIMENSION NAMES FOR VARIABLE %s : %s <> %s\n", opts->cmpvarlist[i], dims1[dimid1].name, dims2[dimid2].name);
        
        if (!opts->warn[NCCMP_W_ALL]) 
          status = EXIT_DIFFER;
        
        if(opts->force) continue; else goto recover;        
      }

      tmp1 = strcmp(dims1[dimid1].name, recname1);
      tmp2 = strcmp(dims2[dimid2].name, recname2);

      if ( (tmp1 == 0) && (tmp2 != 0) ) {
        fprintf(stderr, "DIFFER : VARIABLE : %s : DIMENSION %s IS RECORD IN FILE \"%s\" BUT NOT IN \"%s\"\n", vars1[j1].name, dims1[dimid1].name, opts->file1, opts->file2);
        
        if (!opts->warn[NCCMP_W_ALL]) 
          status = EXIT_DIFFER;                                         
          
        if(opts->force) continue; else goto recover;                     
      } else if ( (tmp1 != 0) && (tmp2 == 0) ) {
        fprintf(stderr, "DIFFER : VARIABLE : %s : DIMENSION %s IS RECORD IN FILE \"%s\" BUT NOT IN \"%s\"\n", vars1[j1].name, dims2[dimid2].name, opts->file2, opts->file1);
        
        if (!opts->warn[NCCMP_W_ALL]) 
          status = EXIT_DIFFER;                                         
          
        if(opts->force) continue; else goto recover;                     
      }
    }
  }
  
  if (opts->verbose)
    printf("INFO: Comparing variables' attributes.\n");

  /*printf("DEBUG : %d : \n", __LINE__);*/
  /* Pass in 'i' as dummy. */
  if (newstringlist(&processedatts, &i, NC_MAX_VARS)) {
    fprintf(stderr, "ERROR: Failed to allocate string list for comparing attributes.\n");
    return EXIT_FATAL;
  }
  /*printf("DEBUG : %d : \n", __LINE__); fflush(stdout);*/
  
  for(i=0; i < opts->ncmpvarlist; ++i) {
    /*printf("DEBUG : %d : i = %d\n", __LINE__, i); fflush(stdout);*/
    j1 = findvar(opts->cmpvarlist[i], vars1);
    j2 = findvar(opts->cmpvarlist[i], vars2);
    
    if ((j1 == -1) || (j2 == -1))
      continue;

    natts1 = natts2 = 0;
    clearstringlist(processedatts, NC_MAX_VARS);
    
    /*printf("DEBUG : %d : var=%s\n", __LINE__, opts->cmpvarlist[i]); fflush(stdout);*/
    
    /* Attributes */
    for(attid1=0; attid1 < vars1[j1].natts; ++attid1) {
      ncstatus = nc_inq_attname(ncid1, vars1[j1].varid, attid1, name1);
      if (ncstatus != NC_NOERR) {
        if (!opts->warn[NCCMP_W_ALL]) 
          status = EXIT_DIFFER;
          
        if(opts->force) continue; else goto recover;
      }
  
      if (instringlist(opts->excludeattlist, name1, opts->nexcludeatt) || instringlist(processedatts, name1, NC_MAX_VARS))
        continue;

      /* Log that this att was processed. */
      addstringtolist(processedatts, name1, NC_MAX_VARS);
      ++natts1;
      /*printf("natts1 %s, %d\n", name1, natts1);*/
      ncstatus = cmpatt(ncid1, ncid2, vars1[j1].varid, vars2[j2].varid, name1, vars1[j1].name, opts);
      if (ncstatus == EXIT_DIFFER) {
        status = ncstatus;
        if(opts->force) continue; else goto recover;
      }
    }
    
    /*printf("DEBUG : %d : \n", __LINE__);*/
    for(attid2=0; attid2 < vars2[j2].natts; ++attid2) {
      ncstatus = nc_inq_attname(ncid2, vars2[j2].varid, attid2, name2);
      if (ncstatus != NC_NOERR) {
        fprintf(stderr, "Failed to query variable %s attribute in file \"%s\"\n", vars2[j2].name, opts->file2);
        
        if (!opts->warn[NCCMP_W_ALL]) 
          status = EXIT_DIFFER;
          
        if(opts->force) continue; else goto recover;
      }

      if (instringlist(opts->excludeattlist, name2, opts->nexcludeatt))
        continue;
      
      /* Count non-excluded attribute. */
      ++natts2;
      /*printf("natts2 %s, %d\n", name2, natts2);*/
      if (instringlist(processedatts, name2, NC_MAX_VARS))
        continue;

      /* Log that this att was processed. */
      addstringtolist(processedatts, name2, NC_MAX_VARS);
      
      /* Do comparison. */
      ncstatus = cmpatt(ncid1, ncid2, vars1[j1].varid, vars2[j2].varid, name2, vars2[j2].name, opts);
      if (ncstatus == EXIT_DIFFER) {
        status = ncstatus;
        if(opts->force) continue; else goto recover;
      }
    }
    
    if(natts1 != natts2) {
      fprintf(stderr, "DIFFER : NUMBER OF ATTRIBUTES : VARIABLE : %s : %d <> %d\n", opts->cmpvarlist[i], natts1, natts2);
      if (!opts->warn[NCCMP_W_ALL]) 
        status = EXIT_DIFFER;                                         
        
      if(opts->force) continue; else goto recover;
    }
  }
  
 recover:
  freestringlist(&processedatts, NC_MAX_VARS);
  processedatts = NULL;
  
  return status;  
}
/* *********************************************************** */
/* Returns index into varstruct array if found using name otherwise -1. */
int findvar(char * name, varstruct *vars)
{
  int i;
  for(i=0; i < NC_MAX_VARS; ++i) {
    if (strcmp(name, vars[i].name) == 0)
      return i;
  }

  return -1;
}
/* *********************************************************** */
/* Element compare and return first offset that differs. 
   Arrays should have sentinals at end that are guaranteed to differ.
*/
off_t cmp_text(char* in1, char* in2) {
  char const *p1, *p2;
  for(p1 = in1, p2 = in2; *p1 == *p2; ++p1, ++p2)
    continue;
  return p1 - in1;
}
off_t cmp_text_missing(char* in1, char* in2, char m1, char m2) {
  char const *p1, *p2;
  for(p1 = in1, p2 = in2; (*p1 == *p2) || ((*p1 == m1) && (*p2 == m2)); ++p1, ++p2)
    continue;
  return p1 - in1;
}
off_t cmp_byte(int8_t* in1, int8_t* in2) {
  int8_t const *p1, *p2;
  for(p1 = in1, p2 = in2; *p1 == *p2; ++p1, ++p2)
    continue;
  return p1 - in1;
}
off_t cmp_byte_missing(int8_t* in1, int8_t* in2, int8_t m1, int8_t m2) {
  int8_t const *p1, *p2;
  for(p1 = in1, p2 = in2; (*p1 == *p2) || ((*p1 == m1) && (*p2 == m2)); ++p1, ++p2)
    continue;
  return p1 - in1;
}
off_t cmp_short(short* in1, short* in2) {
  short const *p1, *p2;
  for(p1 = in1, p2 = in2; *p1 == *p2; ++p1, ++p2)
    continue;
  return p1 - in1;
}
off_t cmp_short_missing(short* in1, short* in2, short m1, short m2) {
  short const *p1, *p2;
  for(p1 = in1, p2 = in2; (*p1 == *p2) || ((*p1 == m1) && (*p2 == m2)); ++p1, ++p2)
    continue;
  return p1 - in1;
}
off_t cmp_int(int* in1, int* in2) {
  int const *p1, *p2;
  for(p1 = in1, p2 = in2; *p1 == *p2; ++p1, ++p2)
    continue;
  return p1 - in1;
}
off_t cmp_int_missing(int* in1, int* in2, int m1, int m2) {
  int const *p1, *p2;
  for(p1 = in1, p2 = in2; (*p1 == *p2) || ((*p1 == m1) && (*p2 == m2)); ++p1, ++p2)
    continue;
  return p1 - in1;
}
off_t cmp_float(float* in1, float* in2) {
  float const *p1, *p2;
  for(p1 = in1, p2 = in2; fabs(*p1 - *p2)<  FLT_EPSILON; ++p1, ++p2)
    continue;
  return p1 - in1;
}
off_t cmp_float_missing(float* in1, float* in2, float m1, float m2) {
  float const *p1, *p2;
  for(p1 = in1, p2 = in2; (fabs(*p1 - *p2)<  FLT_EPSILON) || ((fabs(*p1 - m1)<  FLT_EPSILON) && (fabs(*p2 - m2)<  FLT_EPSILON)); ++p1, ++p2)
    continue;
  return p1 - in1;
}
off_t cmp_double(double* in1, double* in2) {
  double const *p1, *p2;
  for(p1 = in1, p2 = in2; fabs(*p1 - *p2)<  DBL_EPSILON; ++p1, ++p2)
    continue;
  return p1 - in1;
}
off_t cmp_double_missing(double* in1, double* in2, double m1, double m2) {
  double const *p1, *p2;
  for(p1 = in1, p2 = in2; (fabs(*p1 - *p2)<  DBL_EPSILON) || ((fabs(*p1 - m1)<  DBL_EPSILON) && (fabs(*p2 - m2)<  DBL_EPSILON)); ++p1, ++p2)
    continue;
  return p1 - in1;
}

/* *********************************************************** */
/* Do the comparision of variables.
   Record index (rec) is optional; -1 if not applicable.
   Returns comparison result success or failure.
*/
int cmpvar(char* name, int rec, nccmpopts* opts, int ncid1, int ncid2)
{
  int8_t *b1,*b2;
  char *c1,*c2;
  short *s1,*s2;
  int *i1,*i2;
  float *f1,*f2;
  double *d1,*d2;
  varstruct *v1,*v2;
  int idx1, idx2, status, i;
  size_t start[NC_MAX_DIMS], count[NC_MAX_DIMS];
  char tmpstr1[256], tmpstr2[256], idxstr[256];
  off_t nitems, diff, cmplen;
  size_t odomax[NC_MAX_DIMS];
  int diffstatus = EXIT_SUCCESS;
  char *message = "DIFFER : VARIABLE : %s : POSITION : %s : VALUES : %s <> %s\n";
  char value1str[32], value2str[32];
  double doublevalue1, doublevalue2;
  int intvalue1, intvalue2;
  short shortvalue1, shortvalue2;
  float floatvalue1, floatvalue2;
  char textvalue1, textvalue2;
  int8_t bytevalue1, bytevalue2;
  char printHex = strchr(opts->precision, 'X') || strchr(opts->precision, 'x');
  char do_missing = 0;
  
  if (opts->verbose) {
    if (rec != -1)
      printf("INFO: Comparing data for variable \"%s\" at record %d.\n", name, (int)rec);
    else
      printf("INFO: Comparing non-record data for variable \"%s\".\n", name);
  }
  
  idx1 = findvar(name, vars1);
  idx2 = findvar(name, vars2);

  v1 = &vars1[idx1];
  v2 = &vars2[idx2];

  /* Setup missing mode. */
  if (opts->missing && v1->hasmissing && v2->hasmissing) {
    do_missing = 1;
  }
  
  /*printf("DEBUG : %s len : %d <> %d\n", name, v1->len, v2->len); */
  
  if (v1->len != v2->len) {
    fprintf(stderr, "DIFFER : SIZE OF VARIABLE \"%s\" : %d <> %d\n", name, (int)v1->len, (int)v2->len);
    
    if (!opts->warn[NCCMP_W_ALL]) \
      return EXIT_DIFFER;
    else
      return EXIT_SUCCESS;
  }
  
  for(i=0; i < v1->ndims; ++i) {
    start[i] = 0;
    odomax[i] = v1->dimlens[i] - 1;
  }

#ifdef __DEBUG__
  printf("DEBUG : %d : odomax = ", __LINE__);
  for(i=0; i < v1->ndims; ++i) {
    printf("%d ", odomax[i]);
  }
  printf("\n");
#endif

  /* If has record dim. */
  if (v1->hasrec && (rec >= 0)) 
    start[0] = rec;
  
  /* Read in slab for last dimension at-a-time only. 
    We'll then loop over all the other outer dimensions. */
  for(i=0; i < v1->ndims-1; ++i) {
    count[i] = 1;
  }

    /* We'll always read in entire last dimension
       except if only dimension is record. */
  if ((v1->ndims == 1)  && (v1->hasrec)) {
        nitems = 1;
    } else
        nitems = v1->dimlens[v1->ndims-1];
  
    count[v1->ndims-1] = nitems;
    
  /*printf("DEBUG : %d : nitems = %d\n", __LINE__, nitems);\*/

#define CMP_VAR(TYPE, NCFUNTYPE, ncid1, ncid2, P1, P2, M) {\
    P1 = (TYPE*)malloc(sizeof(TYPE) * (nitems + 1));      \
    P2 = (TYPE*)malloc(sizeof(TYPE) * (nitems + 1));      \
     \
    do { \
      /* printf("start = %d %d %d %d, count = %d %d %d %d\n", (int)start[0], (int)start[1], (int)start[2], (int)start[3], (int)count[0], (int)count[1], (int)count[2], (int)count[3]); */ \
      status = nc_get_vara(ncid1, v1->varid, start, count, P1); \
      handle_error(status); \
      status = nc_get_vara(ncid2, v2->varid, start, count, P2); \
      handle_error(status); \
      /* Sentinels. */ \
      P1[nitems] = 0; \
      P2[nitems] = 1; \
      \
      cmplen = nitems; \
            /* for(i=0; i<nitems; ++i) { \
                printf("nitems = %d, rec = %d, P1[%d] = %g, P2[%d] = %g\n", nitems, rec, i, P1[i], i, P2[i]); \
            } */ \
      if (do_missing) \
        diff = cmp_##NCFUNTYPE##_missing(P1, P2, v1->missing.M, v2->missing.M); \
      else \
        diff = cmp_##NCFUNTYPE(P1, P2); \
      \
      while (diff < cmplen) { \
        if (!opts->warn[NCCMP_W_ALL]) \
          diffstatus = EXIT_DIFFER; \
        else \
          diffstatus = EXIT_SUCCESS; \
        \
        if (opts->fortran) \
          getidxstr_fortran(v1, start, nitems-cmplen+diff, idxstr); \
        else \
          getidxstr(v1, start, nitems-cmplen+diff, idxstr); \
        \
        NCFUNTYPE##value1 = P1[nitems-cmplen+diff]; \
        NCFUNTYPE##value2 = P2[nitems-cmplen+diff]; \
        if (printHex) {\
          NCFUNTYPE##ToHex(NCFUNTYPE##value1, value1str); \
          NCFUNTYPE##ToHex(NCFUNTYPE##value2, value2str); \
        } else { \
          NCFUNTYPE##ToString(NCFUNTYPE##value1, value1str, opts->precision); \
          NCFUNTYPE##ToString(NCFUNTYPE##value2, value2str, opts->precision); \
        } \
        fprintf(stderr, message, v1->name, idxstr, value1str, value2str); \
        /* printf("%d nitems=%d, cmplen=%d, diff=%d\n", __LINE__, nitems, cmplen, diff); */\
        if (opts->force) { \
          cmplen = cmplen - (diff+1); \
          diff = cmp_##NCFUNTYPE(P1+diff+1, P2+diff+1); \
        } else  \
          goto break_##NCFUNTYPE; \
      } \
      /* Increment all but first (if record) and last dimensions. */\
    } while(odometer(start, odomax, (rec >= 0), v1->ndims-2)); \
    break_##NCFUNTYPE:\
    free(P1); \
    free(P2); \
    }

  switch(v1->type) {
  case NC_BYTE: 
    CMP_VAR(int8_t,
                byte, ncid1, ncid2,
                b1, b2, b);
    break;
    
  case NC_CHAR:
    CMP_VAR( char, 
                text, ncid1, ncid2,
                c1, c2, c);
    break;
    
  case NC_SHORT:
    CMP_VAR(short, short, ncid1, ncid2,s1, s2, s)
    break;
    
  case NC_INT:
    CMP_VAR(int, int, ncid1, ncid2,i1, i2, i)
    break;
  
  case NC_FLOAT:
    CMP_VAR(float, float, ncid1, ncid2,f1, f2, f);
    break;
  
  case NC_DOUBLE:
    CMP_VAR(double, double,ncid1, ncid2, d1, d2, d);
    break;
  }

  return diffstatus;
}
/* *********************************************************** */
/* Do the comparision of variables.
   Record index (rec) is optional; -1 if not applicable.
   Returns comparison result success or failure.
*/
int cmpvartol(char* name, int rec, nccmpopts* opts, int ncid1, int ncid2)
{
  int8_t *b1,*b2;
  char *c1,*c2;
  short *s1,*s2;
  int *i1,*i2;
  float *f1,*f2;
  double *d1,*d2;
  varstruct *v1,*v2;
  int idx1, idx2, status, i;
  size_t start[NC_MAX_DIMS], count[NC_MAX_DIMS];
  char tmpstr1[256], tmpstr2[256], idxstr[256];
  off_t nitems, diff, cmplen;
  size_t odomax[NC_MAX_DIMS];
  int diffstatus = EXIT_SUCCESS;
  char *message = "DIFFER : VARIABLE : %s : POSITION : %s : VALUES : %s <> %s : PERCENT : %g\n";
  char value1str[32], value2str[32];
  double absdelta;
  double doublevalue1, doublevalue2;
  int intvalue1, intvalue2;
  short shortvalue1, shortvalue2;
  float floatvalue1, floatvalue2;
  char textvalue1, textvalue2;
  int8_t bytevalue1, bytevalue2;
  char printHex = strchr(opts->precision, 'X') || strchr(opts->precision, 'x');
  char do_missing = 0;
  
  if (opts->verbose) {
    if (rec != -1)
      printf("INFO: Comparing data for variable \"%s\" at record %d.\n", name, (int)rec);
    else
      printf("INFO: Comparing non-record data for variable \"%s\".\n", name);
  }

  idx1 = findvar(name, vars1);
  idx2 = findvar(name, vars2);

  if (idx1 < 0) {
    if (! opts->metadata) /* This gets reported in cmpmeta. */
      fprintf(stderr, "DIFFER : Failed to find variable \"%s\" in file \"%s\".\n", name, opts->file1);
    
    if (!opts->warn[NCCMP_W_ALL])
      return EXIT_DIFFER;
    else
      return EXIT_SUCCESS;
  }

  if (idx2 < 0) {
    if (! opts->metadata) /* This gets reported in cmpmeta. */
      fprintf(stderr, "DIFFER : Failed to find variable \"%s\" in file \"%s\".\n", name, opts->file2);
  
    if (!opts->warn[NCCMP_W_ALL])
      return EXIT_DIFFER;
    else
      return EXIT_SUCCESS;
  }

  v1 = &vars1[idx1];
  v2 = &vars2[idx2];
    
  /* Setup missing mode. */
  if (opts->missing && v1->hasmissing && v2->hasmissing) {
    do_missing = 1;
  }

  for(i=0; i < v1->ndims; ++i) {
    start[i] = 0;
    odomax[i] = v1->dimlens[i] - 1;
  }

  /* If has record dim. */
  if (v1->hasrec && (rec >= 0)) 
    start[0] = rec;
  
  /* Read in slab for last dimension at-a-time only. 
    We'll then loop over all the other outer dimensions. */
  for(i=0; i < v1->ndims-1; ++i) {
    count[i] = 1;
  }
  
    /* We'll always read in entire last dimension
       except if only dimension is record. */
  if ((v1->ndims == 1)  && (v1->hasrec)) {
        nitems = 1;
    } else
        nitems = v1->dimlens[v1->ndims-1];
  
    count[v1->ndims-1] = nitems;

    /* todo: make cmpvar and cmpvartol same function to re-use code immediately above; just use conditional to choose CMP_VAR or CMP_VARTOL macro. */
#define CMP_VARTOL(TYPE, NCFUNTYPE, ncid1, ncid2, P1, P2, M) {       \
    P1 = (TYPE*)malloc(sizeof(TYPE) * nitems);      \
    P2 = (TYPE*)malloc(sizeof(TYPE) * nitems);      \
     \
    do { \
            /* printf("start = %d %d %d %d, count = %d %d %d %d\n", (int)start[0], (int)start[1], (int)start[2], (int)start[3], (int)count[0], (int)count[1], (int)count[2], (int)count[3]); */ \
      status = nc_get_vara(ncid1, v1->varid, start, count, P1); \
      handle_error(status); \
      status = nc_get_vara(ncid2, v2->varid, start, count, P2); \
      handle_error(status); \
      \
            /* for(i=0; i<nitems; ++i) { \
                printf("nitems = %d, rec = %d, P1[%d] = %g, P2[%d] = %g\n", nitems, rec, i, P1[i], i, P2[i]); \
            } */ \
      for(i=0; i < nitems; ++i) \
      {                       \
        if (do_missing) { \
          if ((v1->missing.M == P1[i]) && (v2->missing.M == P2[i])) continue; \
        } \
        \
        absdelta = fabs((double)(P1[i]-P2[i]));     \
        if (opts->abstolerance ? (absdelta > opts->tolerance) : (double)absdelta*100./(fabs((double)P1[i]) > fabs((double)P2[i]) ? fabs((double)P1[i]) : fabs((double)P2[i])) > opts->tolerance)  \
        {                           \
          if (!opts->warn[NCCMP_W_ALL]) \
            diffstatus = EXIT_DIFFER; \
          else \
            diffstatus = EXIT_SUCCESS; \
            \
          if ((v1->ndims == 1)  && (v1->hasrec)) {\
            if (opts->fortran) \
              getidxstr_fortran(v1, start, rec, idxstr); \
            else \
              getidxstr(v1, start, rec, idxstr); \
          } else {\
            if (opts->fortran) \
              getidxstr_fortran(v1, start, i, idxstr); \
            else \
              getidxstr(v1, start, i, idxstr); \
          }\
          NCFUNTYPE##value1 = P1[i]; \
          NCFUNTYPE##value2 = P2[i]; \
          if (printHex) {\
            NCFUNTYPE##ToHex(NCFUNTYPE##value1, value1str); \
            NCFUNTYPE##ToHex(NCFUNTYPE##value2, value2str); \
          } else {\
            NCFUNTYPE##ToString(NCFUNTYPE##value1, value1str, opts->precision); \
            NCFUNTYPE##ToString(NCFUNTYPE##value2, value2str, opts->precision); \
          } \
          fprintf(stderr, message, v1->name, idxstr, value1str, value2str, (double)absdelta*100./(fabs((double)P1[i]) > fabs((double)P2[i]) ? fabs((double)P1[i]) : fabs((double)P2[i])) ); \
          if (! opts->force)                      \
          {                                       \
            goto break_##NCFUNTYPE;                   \
          }                                       \
        }                                         \
      }                                           \
    } while(odometer(start, odomax, (rec >= 0), v1->ndims-2)); \
    break_##NCFUNTYPE:\
    free(P1); \
    free(P2); \
}

  switch(v1->type) {
  case NC_BYTE: 
    CMP_VARTOL(int8_t,
                byte, ncid1, ncid2,
                b1, b2, b);
    break;
    
  case NC_CHAR:
    CMP_VARTOL( char,
                text, ncid1, ncid2,
                c1, c2, c);
    break;
    
  case NC_SHORT:
    CMP_VARTOL(short, short, ncid1, ncid2,s1, s2, s)
    break;
    
  case NC_INT:
    CMP_VARTOL(int, int, ncid1, ncid2,i1, i2, i)
    break;
  
  case NC_FLOAT:
    CMP_VARTOL(float, float,ncid1, ncid2, f1, f2, f);
    break;
  
  case NC_DOUBLE:
    CMP_VARTOL(double, double,ncid1, ncid2, d1, d2, d);
    break;
  }

  return diffstatus;
}
/* *********************************************************** */
int 
nccmpdatarecvartol(int ncid1, int ncid2, char* varname, nccmpopts* opts,
                size_t recstart, size_t recend)
{
    int cmpstatus;
    int status = EXIT_SUCCESS;
  
    for(; recstart <= recend; ++recstart)
    {
      cmpstatus = cmpvartol(varname, recstart, opts,ncid1, ncid2) || status;
      if (cmpstatus != EXIT_SUCCESS) {
        status = cmpstatus;
        if (opts->force)
          continue;
        else
          break;
      }
    }

    return status;
}
/* *********************************************************** */
int 
nccmpdatarecvar(int ncid1, int ncid2, char* varname, nccmpopts* opts,
                size_t recstart, size_t recend)
{
    int cmpstatus;
    int status = EXIT_SUCCESS;
  
    for(; recstart <= recend; ++recstart)
    {
      cmpstatus = cmpvar(varname, recstart, opts, ncid1, ncid2) || status;
      if (cmpstatus != EXIT_SUCCESS) {
        status = cmpstatus;
        if (opts->force)
          continue;
        else
          break;
      }
    }

    return status;
}
/* *********************************************************** */
int nccmpdatatol(int ncid1, int ncid2, nccmpopts* opts)
{
  int i, idx1, idx2;
  int status, cmpstatus, nprocessed;
  char** processed = NULL;

  status = EXIT_SUCCESS;

  if (opts->verbose)
    printf("INFO: Comparing data with tolerance.\n");
  
  if(newstringlist(&processed, &nprocessed, NC_MAX_VARS)) {
    fprintf(stderr, "ERROR: Failed to allocated string list for comparing data.\n");
    return EXIT_FATAL;
  }
  
  for(i=0; i < opts->ncmpvarlist; ++i) {
    if (instringlist(processed, opts->cmpvarlist[i], nprocessed))
      /* Skip varnames already processed. */
      continue;
    
    if (opts->verbose)
      printf("INFO: Comparing data for variable \"%s\".\n", opts->cmpvarlist[i]);

    addstringtolist(processed, opts->cmpvarlist[i], nprocessed);

    /* Has rec? */
    idx1 = findvar(opts->cmpvarlist[i], vars1);
    if (vars1[idx1].hasrec) { 

      /* Compare only if # recs are equal and not zero. */ 
      if ((nrec1 == nrec2) && (nrec1+nrec2)) {
        cmpstatus = nccmpdatarecvartol(ncid1, ncid2, opts->cmpvarlist[i], opts, 0, nrec1-1);
        if (cmpstatus) {
          status = cmpstatus;
          if (opts->force)
            continue;
          else
            break;
        }
      }
    } else {
      cmpstatus = cmpvartol(opts->cmpvarlist[i], -1, opts,ncid1, ncid2);
      if (cmpstatus) {
        status = cmpstatus;
        if (opts->force)
          continue;
        else
          break;
      }
    }
  }
  
  if (opts->verbose)
    printf("INFO: Finished comparing data.\n");

  return status;
}
/* *********************************************************** */
int nccmpdata(nccmpopts* opts, int ncid1, int ncid2)
{
  int i, idx1, idx2;
  int status, cmpstatus, nprocessed;
  char** processed = NULL;
    char str1[32], str2[32];
    
  status = EXIT_SUCCESS;

  if (opts->verbose)
    printf("INFO: Comparing data.\n");

  if (opts->tolerance != 0)
    return nccmpdatatol(ncid1, ncid2, opts);
  
  if(newstringlist(&processed, &nprocessed, NC_MAX_VARS)) {
    fprintf(stderr, "ERROR: Failed to allocated string list for comparing data.\n");
    return EXIT_FATAL;
  }
  
  for(i=0; i < opts->ncmpvarlist; ++i) {
    if (instringlist(processed, opts->cmpvarlist[i], nprocessed))
      /* Skip varnames already processed. */
      continue;
    
    if (opts->verbose)
      printf("INFO: Comparing data for variable \"%s\".\n", opts->cmpvarlist[i]);

    addstringtolist(processed, opts->cmpvarlist[i], nprocessed);
        
    idx1 = findvar(opts->cmpvarlist[i], vars1); 
    idx2 = findvar(opts->cmpvarlist[i], vars2); 

      if (idx1 < 0) {
        if (! opts->metadata) /* This gets reported in cmpmeta. */
          fprintf(stderr, "DIFFER : Failed to find variable \"%s\" in file \"%s\".\n", opts->cmpvarlist[i], opts->file1);
        
        if (!opts->warn[NCCMP_W_ALL])
          status = EXIT_DIFFER;
      if (opts->force)
          continue;
            else
          break;
      }

      if (idx2 < 0) {
        if (! opts->metadata) /* This gets reported in cmpmeta. */
          fprintf(stderr, "DIFFER : Failed to find variable \"%s\" in file \"%s\".\n", opts->cmpvarlist[i], opts->file2);

        if (!opts->warn[NCCMP_W_ALL])
          status = EXIT_DIFFER;
      if (opts->force)
          continue;
            else
          break;
      }
                
        if (vars1[idx1].len != vars2[idx2].len) {
        fprintf(stderr, "DIFFER : SIZE OF VARIABLE \"%s\" : %d <> %d\n", opts->cmpvarlist[i], (int)vars1[idx1].len, (int)vars2[idx2].len);

        if (!opts->warn[NCCMP_W_ALL])
          status = EXIT_DIFFER;
      if (opts->force)
          continue;
            else
          break;
      }
        
        if (vars1[idx1].type != vars2[idx2].type) {
            type2string(vars1[idx1].type, str1);
            type2string(vars2[idx2].type, str2);
        fprintf(stderr, "DIFFER : TYPE OF VARIABLE \"%s\" : %s <> %s\n", opts->cmpvarlist[i], str1, str2);
    
      if (!opts->warn[NCCMP_W_ALL])
        status = EXIT_DIFFER;

      if (opts->force)
          continue;
            else
          break;
      }
        
        /* Has rec? */
    if (vars1[idx1].hasrec) { 

      /* Compare only if # recs are equal and not zero. */ 
      if ((nrec1 == nrec2) && (nrec1+nrec2)) {
        /* TODO: Check if ignorem missing. */
        cmpstatus = nccmpdatarecvar(ncid1,ncid2, opts->cmpvarlist[i], opts, 0, nrec1-1);
        if (cmpstatus)
        {
          status = cmpstatus;
          if (opts->force)
            continue;
          else
            break;
        }
      }
    } else {
      /* TODO: Check if ignorem missing. */
      cmpstatus = cmpvar(opts->cmpvarlist[i], -1, opts, ncid1, ncid2);
      if (cmpstatus)
      {
        status = cmpstatus;
        if (opts->force)
          continue;
        else
          break;
      }
    }
  }
  
  if (opts->verbose)
    printf("INFO: Finished comparing data.\n");

  return status;
}
/* *********************************************************** */
int 
nccmp(nccmpopts* opts)
{
  int ncid1, ncid2;
  int numgrps1, numgrps2;
  int *gids1, *gids2;
  int i, status;
  
  if (opts->verbose)
    printf("INFO: Opening input files.\n");

  status = openfiles(opts,&ncid1, &ncid2);
  if (status)
    return status;

  if (opts->verbose)
    printf("INFO: Creating variable comparison list.\n");
  
  GROUP_NODE *groups1 = NULL;
  GROUP_NODE *groups2 = NULL;

  if (opts->verbose)
    printf("INFO: Comparing file formats.\n");

  status += nccmpformats(opts, ncid1, ncid2);
  if (status && !opts->force)
    return status;

  if (opts->verbose)
    printf("INFO: Comparing global attributes.\n");
  status += nccmpglobalatts(opts, ncid1, ncid2);
    
  if (status && !opts->force)
    return status;
  
  if (opts->verbose)
    printf("INFO: Collecting dimension information for first file.\n");

  getdiminfo(ncid1, dims1, &ndims1);

  if (opts->verbose)
    printf("INFO: Collecting dimension information for second file.\n");
  
  getdiminfo(ncid2, dims2, &ndims2);

  nc_inq_grps(ncid1, &numgrps1, NULL);
  if (numgrps1 > 0){
      groups1 = malloc(sizeof(GROUP_NODE)*numgrps1);
      getgroupinfo(ncid1, numgrps1, groups1);
  }
  nc_inq_grps(ncid2, &numgrps2, NULL);
    if (numgrps2 > 0){
        groups2 = malloc(sizeof(GROUP_NODE)*numgrps2);
        getgroupinfo(ncid2, numgrps2, groups2);
    }

    if (numgrps1 == numgrps2 && numgrps1 > 0)  {

        for (i = 0; i < numgrps1; i++){
            status += makecmpvarlist(opts,groups1[i].groupID, groups2[i].groupID);
            if (status && !opts->force)
                return status;

            if (opts->verbose)
                printf("INFO: Comparing group attributes [%s]:\n",groups1[i].groupName);
            status += nccmpglobalatts(opts, groups1[i].groupID, groups2[i].groupID);
            if (status && !opts->force)
                return status;

            if (opts->verbose)
                printf("INFO: Collecting variable information for first file.\n");

            getvarinfo(groups1[i].groupID, vars1 , &nvars1, opts->verbose);

            if (opts->verbose)
                printf("INFO: Collecting variable information for second file.\n");

            getvarinfo(groups2[i].groupID, vars2, &nvars2, opts->verbose);

            status += nccmprecinfo(opts,groups1[i].groupID,groups2[i].groupID);
            if (status && !opts->force)
                return status;

            if (opts->metadata) {
                status += nccmpmetadata(opts,groups1[i].groupID,groups2[i].groupID);
            }

            if (status && !opts->force)
                return status;

            if (opts->data) {
                status += nccmpdata(opts,groups1[i].groupID,groups2[i].groupID );
            }
            if (status && !opts->force)
                return status;

            clearstringlist(opts->cmpvarlist, opts->ncmpvarlist);
        }
    }
    // Check root group
    status += makecmpvarlist(opts,ncid1, ncid2);
    if (status && !opts->force)
        return status;

    if (opts->verbose)
        printf("INFO: Collecting variable information for first file.\n");

    getvarinfo(ncid1, vars1 , &nvars1, opts->verbose);

    if (opts->verbose)
        printf("INFO: Collecting variable information for second file.\n");

    getvarinfo(ncid2, vars2, &nvars2, opts->verbose);

    status += nccmprecinfo(opts,ncid1, ncid2);
    if (status && !opts->force)
        return status;

    if (opts->metadata) {
        status += nccmpmetadata(opts,ncid1, ncid2);
    }

    if (status && !opts->force)
        return status;

    if (opts->data) {
        status += nccmpdata(opts,ncid1, ncid2 );
    }
    if (status && !opts->force)
        return status;

  if (opts->verbose)
    printf("INFO: Comparisons complete. Freeing memory.\n");
  
  return status;
}

/* *********************************************************** */

int
main(int argc, char** argv)
{
  int status;
  nccmpopts opts;

  status = EXIT_SUCCESS;

  initnccmpopts(&opts);
  
  /* parse command-line args. & options */
  status = getnccmpopts(argc, argv, &opts);
  
  if (status != EXIT_SUCCESS) 
    goto end;
  
  if (opts.verbose)
    printf("INFO: Command-line options parsed.\n");

  status = nccmp(&opts);
  
end:
  if (opts.report_identical) {
    if ( !(opts.help || opts.version) && 
          (status == EXIT_SUCCESS) ) {
      printf("Files \"%s\" and \"%s\" are identical.\n",
        opts.file1, opts.file2);
    }
  }
  
  freenccmpopts(&opts);
  
  exit(status);
} 
