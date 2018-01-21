/*----------------------------------------------------------------------------|
|                                                                             |
|        Copyright (C) 1999  Emergent IT Inc.  and Raytheon Systems Company   |
|                                                                             |
|Permission to use, modify, and distribute this software and its documentation|
|for any purpose without fee is hereby granted, provided that the above       |
|copyright notice appear in all copies and that both that copyright notice and|
|this permission notice appear in supporting documentation.                   |
|                                                                             |
|-----------------------------------------------------------------------------|
|                                                                             |
| This is the main header file to be distributed with the HDF-EOS library.    |
|                                                                             |
| Last date updated: June 5, 2001                                             |
|                    Aug 23, 2001  A.M. Added thread-safe related blocks.     |
|                    May 29, 2002  S.Z  Added ZA interface                    |
|                    August, 2003  S.Z  Added szip compression methods.       |
|                    April,  2004  S.Z  Added a data type flag HE5T_CHARSTRING|
|                    August, 2004  S.Z  Added field number type in            |
|                                       HE5_CmpDTSinfo                        |
-----------------------------------------------------------------------------*/


#ifndef HE5_HDFEOSDEF_H_
#define HE5_HDFEOSDEF_H_

#include <hdf5.h>

#ifdef H5_USE_16_API
#include <H5DSpublic.h>
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef __cplusplus
#include  <cfortHdf.h>
#endif

#ifdef   _HDFEOS5_THREADSAFE
#include <pthread.h>
#endif


#define HE5_HDFEOSVERSION "1.14.01"

#define H5_USE_16_API 1

#ifndef FALSE
#   define FALSE 0
#endif
#ifndef TRUE
#   define TRUE (!FALSE)
#endif

#ifndef SUCCEED
#   define SUCCEED           0
#endif

#ifndef FAIL
#   define FAIL            (-1)
#endif

#ifndef HDstrcmp
#   define  HDstrcmp(X,Y)		strcmp(X,Y)
#endif

#ifndef  MAX
#    define  MAX(X,Y)           ((X)>(Y)?(X):(Y))
#endif

#ifndef  MIN
#    define  MIN(X,Y)           ((X)<(Y)?(X):(Y))
#endif

#define  CHECKPOINTER(p)    {                              \
     status = HE5_EHchkptr((void *)p,#p);                  \
     if (status == FAIL) goto COMPLETION;                  \
}

#define  CHECKNAME(p)    {                                 \
     status = HE5_EHchkname((char *)p,#p);                 \
     if (status == FAIL) goto COMPLETION;                  \
}

#ifndef HDFI_H
typedef unsigned uintn;
#endif

#ifdef WIN32
#define LONGLONG __int64
#else
#define LONGLONG long long
#endif

/*
 ----------------------------------------------
 |          HDF-EOS Defined Sizes             |
 ----------------------------------------------
 */

#define HE5_DTSETRANKMAX         8
#define HE5_FLDNUMBERMAX       500
#define HE5_OBJNAMELENMAX      256
#define HE5_BLKSIZE         640000
#define HE5_CHUNKSIZE         1000

#define HE5_HDFE_TYPESTRSIZE    80
#define HE5_HDFE_DIMBUFSIZE    256
#define HE5_HDFE_NAMBUFSIZE    256
#define HE5_HDFE_ERRBUFSIZE    256
#define HE5_HDFE_UTLBUFSIZE   1024

/*
 ----------------------------------------------
 |     HDF-EOS Global Data Structures         |
 ----------------------------------------------
 */

typedef struct
{
  hid_t  ID;     /* Field-associated dataset ID    */
  char   *name;  /*      HDF-EOS field name        */
}HE5_DTSinfo;    /* Field-associated dataset Info  */

/* 
 ----------------------------------------------
 |     Compound Dataset Information Data      |
 |               Structure                    |
 ----------------------------------------------
 */

typedef struct
{
  int            nfields;                  /* Number of data fields    */
  int            rank[HE5_FLDNUMBERMAX];   /* Fields rank array        */
  int            array[HE5_FLDNUMBERMAX];  /* Flag if field is an array*/
  char           *fieldname[HE5_FLDNUMBERMAX];/* Array of field names     */
                                              /* Array of dimension sizes */ 
  size_t         dims[HE5_FLDNUMBERMAX][HE5_DTSETRANKMAX];
  size_t         datasize;                 /* Size of data (bytes)     */
  size_t         offset[HE5_FLDNUMBERMAX]; /* Array of field offsets   */
  hid_t          dtype[HE5_FLDNUMBERMAX];  /* Array of field type IDs  */
  hid_t          numtype[HE5_FLDNUMBERMAX];/* Array of field number type IDs */
  H5T_class_t    dclass[HE5_FLDNUMBERMAX]; /* Array of field class IDs */
}HE5_CmpDTSinfo;

typedef  struct 
{
  long        count;     /* Object counter */
  long        strsize;   /* Operator data  */
  char        *name;     /* Object name    */ 
}HE5_OBJINFO;

/* File access flags for TOOLKIT */
/* ----------------------------- */
#ifndef HDF4_ACC_RDONLY
#   define  HDF4_ACC_RDONLY 1
#endif

#ifndef HDF5_ACC_RDONLY
#   define  HDF5_ACC_RDONLY 11
#endif

#ifndef HDF4_ACC_RDWR
#   define  HDF4_ACC_RDWR 3
#endif

#ifndef HDF5_ACC_RDWR
#   define  HDF5_ACC_RDWR 13
#endif

#ifndef HDF4_ACC_CREATE
#   define  HDF4_ACC_CREATE 4
#endif

#ifndef HDF5_ACC_CREATE
#   define  HDF5_ACC_CREATE 14
#endif

/* Data type flags for FORTRAN wrappers */
/* ------------------------------------ */
#ifndef    HE5F_ACC_RDWR
#define    HE5F_ACC_RDWR           100 
#endif
#ifndef    HE5F_ACC_RDONLY
#define    HE5F_ACC_RDONLY         101 
#endif
#ifndef    HE5F_ACC_TRUNC
#define    HE5F_ACC_TRUNC          102 
#endif

#define    HE5T_NATIVE_INT           0 
#define    HE5T_NATIVE_UINT          1 
#define    HE5T_NATIVE_SHORT         2 
#define    HE5T_NATIVE_USHORT        3 
#define    HE5T_NATIVE_SCHAR         4 
#define    HE5T_NATIVE_UCHAR         5 
#define    HE5T_NATIVE_LONG          6 
#define    HE5T_NATIVE_ULONG         7 
#define    HE5T_NATIVE_LLONG         8 
#define    HE5T_NATIVE_ULLONG        9 
#define    HE5T_NATIVE_FLOAT        10 
#define    HE5T_NATIVE_REAL         10
#define    HE5T_NATIVE_DOUBLE       11 
#define    HE5T_NATIVE_LDOUBLE      12 
#define    HE5T_NATIVE_INT8         13 
#define    HE5T_NATIVE_UINT8        14 
#define    HE5T_NATIVE_INT16        15 
#define    HE5T_NATIVE_UINT16       16 
#define    HE5T_NATIVE_INT32        17 
#define    HE5T_NATIVE_UINT32       18 
#define    HE5T_NATIVE_INT64        19 
#define    HE5T_NATIVE_UINT64       20 
#define    HE5T_NATIVE_B8           21 
#define    HE5T_NATIVE_B16          22 
#define    HE5T_NATIVE_B32          23 
#define    HE5T_NATIVE_B64          24 
#define    HE5T_NATIVE_HSIZE        25 
#define    HE5T_NATIVE_HERR         26 
#define    HE5T_NATIVE_HBOOL        27 
#define    HE5T_STD_I8BE            28 
#define    HE5T_STD_I8LE            29 
#define    HE5T_STD_I16BE           30 
#define    HE5T_STD_I16LE           31 
#define    HE5T_STD_I32BE           32 
#define    HE5T_STD_I32LE           33 
#define    HE5T_STD_I64BE           34 
#define    HE5T_STD_I64LE           35 
#define    HE5T_STD_U8BE            36 
#define    HE5T_STD_U8LE            37 
#define    HE5T_STD_U16BE           38 
#define    HE5T_STD_U16LE           39 
#define    HE5T_STD_U32BE           40 
#define    HE5T_STD_U32LE           41 
#define    HE5T_STD_U64BE           42 
#define    HE5T_STD_U64LE           43 
#define    HE5T_STD_B8BE            44 
#define    HE5T_STD_B8LE            45 
#define    HE5T_STD_B16BE           46 
#define    HE5T_STD_B16LE           47 
#define    HE5T_STD_B32BE           48 
#define    HE5T_STD_B32LE           49 
#define    HE5T_STD_B64BE           50 
#define    HE5T_STD_B64LE           51 
#define    HE5T_IEEE_F32BE          52                 
#define    HE5T_IEEE_F32LE          53                 
#define    HE5T_IEEE_F64BE          54                 
#define    HE5T_IEEE_F64LE          55                 
#define    HE5T_NATIVE_CHAR         56
#define    HE5T_CHARSTRING          57

#define    HE5S_UNLIMITED_F         -1 
#define    HE5S_UNLIMITED_F_64      -1 


/* Field Merge Flags   */
/* ------------------- */
#define HE5_HDFE_NOMERGE   0
#define HE5_HDFE_AUTOMERGE 1


/* XXentries Codes     */
/* ------------------- */
#define HE5_HDFE_NENTDIM   0
#define HE5_HDFE_NENTMAP   1
#define HE5_HDFE_NENTIMAP  2
#define HE5_HDFE_NENTGFLD  3
#define HE5_HDFE_NENTDFLD  4

/* Angle Conversion Codes */
/* ---------------------- */
#define HE5_HDFE_RAD_DEG      0
#define HE5_HDFE_DEG_RAD      1
#define HE5_HDFE_DMS_DEG      2
#define HE5_HDFE_DEG_DMS      3
#define HE5_HDFE_RAD_DMS      4
#define HE5_HDFE_DMS_RAD      5


/* Swath Subset Modes       */
/* ------------------------ */
#define HE5_HDFE_MIDPOINT       0
#define HE5_HDFE_ENDPOINT       1
#define HE5_HDFE_ANYPOINT       2
#define HE5_HDFE_INTERNAL       0
#define HE5_HDFE_EXTERNAL       1
#define HE5_HDFE_NOPREVSUB     -1


/* Grid Origin Codes       */
/* ----------------------- */
#define HE5_HDFE_GD_UL         0
#define HE5_HDFE_GD_UR         1
#define HE5_HDFE_GD_LL         2
#define HE5_HDFE_GD_LR         3


/* Pixel Registration Codes */
/* ------------------------ */
#define HE5_HDFE_CENTER         0
#define HE5_HDFE_CORNER         1


/* GCTP Projection Codes */
/* --------------------- */
#define HE5_GCTP_GEO         0
#define HE5_GCTP_UTM         1
#define HE5_GCTP_SPCS        2
#define HE5_GCTP_ALBERS      3
#define HE5_GCTP_LAMCC       4
#define HE5_GCTP_MERCAT      5
#define HE5_GCTP_PS          6
#define HE5_GCTP_POLYC       7
#define HE5_GCTP_EQUIDC      8
#define HE5_GCTP_TM          9
#define HE5_GCTP_STEREO     10
#define HE5_GCTP_LAMAZ      11
#define HE5_GCTP_AZMEQD     12
#define HE5_GCTP_GNOMON     13
#define HE5_GCTP_ORTHO      14
#define HE5_GCTP_GVNSP      15
#define HE5_GCTP_SNSOID     16
#define HE5_GCTP_EQRECT     17
#define HE5_GCTP_MILLER     18
#define HE5_GCTP_VGRINT     19
#define HE5_GCTP_HOM        20
#define HE5_GCTP_ROBIN      21
#define HE5_GCTP_SOM        22
#define HE5_GCTP_ALASKA     23
#define HE5_GCTP_GOOD       24
#define HE5_GCTP_MOLL       25
#define HE5_GCTP_IMOLL      26
#define HE5_GCTP_HAMMER     27
#define HE5_GCTP_WAGIV      28
#define HE5_GCTP_WAGVII     29
#define HE5_GCTP_OBLEQA     30
#define HE5_GCTP_CEA        97
#define HE5_GCTP_BCEA       98
#define HE5_GCTP_ISINUS     99


/* Tiling Codes            */
/* ----------------------- */
#define HE5_HDFE_NOTILE        0
#define HE5_HDFE_TILE          1


/* Compression Methods       */
/* ------------------------- */
#define HE5_HDFE_COMP_NONE                0
#define HE5_HDFE_COMP_RLE                 1
#define HE5_HDFE_COMP_NBIT                2
#define HE5_HDFE_COMP_SKPHUFF             3
#define HE5_HDFE_COMP_DEFLATE             4
#define HE5_HDFE_COMP_SZIP_CHIP           5
#define HE5_HDFE_COMP_SZIP_K13            6
#define HE5_HDFE_COMP_SZIP_EC             7
#define HE5_HDFE_COMP_SZIP_NN             8
#define HE5_HDFE_COMP_SZIP_K13orEC        9
#define HE5_HDFE_COMP_SZIP_K13orNN       10
#define HE5_HDFE_COMP_SHUF_DEFLATE       11 
#define HE5_HDFE_COMP_SHUF_SZIP_CHIP     12
#define HE5_HDFE_COMP_SHUF_SZIP_K13      13
#define HE5_HDFE_COMP_SHUF_SZIP_EC       14
#define HE5_HDFE_COMP_SHUF_SZIP_NN       15
#define HE5_HDFE_COMP_SHUF_SZIP_K13orEC  16
#define HE5_HDFE_COMP_SHUF_SZIP_K13orNN  17


/* HDFEOS Group Codes              */
/* ------------------------------- */
#define HE5_HDFE_GEOGROUP           0
#define HE5_HDFE_DATAGROUP          1
#define HE5_HDFE_ATTRGROUP          2
#define HE5_HDFE_GRPATTRGROUP       3
#define HE5_HDFE_LOCATTRGROUP       4
#define HE5_HDFE_PROFGROUP          5
#define HE5_HDFE_PROFGRPATTRGROUP   6
#define HE5_HDFE_GEOGRPATTRGROUP    7   


/*
 ----------------------------------------------
 |    For  HDF-EOS5 Thread Safe Library       |
 ----------------------------------------------
 */

#ifdef   _HDFEOS5_THREADSAFE

typedef struct
{
  pthread_mutex_t  Lock;
  pthread_t        *MasterThread;
  pthread_cond_t   CondVar;
  unsigned int     LockCount;
} HE5_HDFE_MutexStruct;


HE5_HDFE_MutexStruct   GlobalMutex;


/* Macro for first thread initialization */
/* ------------------------------------- */
#define HE5_FIRST_THREAD_INIT  {                                             \
   status = pthread_once(&HE5_HDFE_TS_FirstInit, HE5_TSinitfirst);           \
}

/* Macro for locking the mutex */
/* --------------------------- */
#define HE5_LOCK        {                                                    \
   status = HE5_TSmutexlock(&GlobalMutex);                                   \
     if (status != SUCCEED) goto COMPLETION;                                 \
}

/* Macro for unlocking the mutex */
/* ----------------------------- */
#define HE5_UNLOCK      {                                                    \
   status = HE5_TSmutexunlock(&GlobalMutex);                                 \
     if (status != SUCCEED) goto COMPLETION;                                 \
}

#else  /* -D_HDFEOS5_THREADSAFE */

typedef struct
{
  int dummyVar;
} HE5_HDFE_MutexStruct;

/* disable any first thread init mechanism */
/* --------------------------------------- */

#define HE5_FIRST_THREAD_INIT

/* disable locks */
/* ------------- */
#define HE5_LOCK
#define HE5_UNLOCK

#endif  /* -D_HDFEOS5_THREADSAFE */


#ifdef __cplusplus
extern "C" {
#endif

/* Thread-safe function prototypes */
/* ------------------------------- */
void    
HE5_TSinitfirst( void);
herr_t   
HE5_TSmutexlock(HE5_HDFE_MutexStruct *);
herr_t   
HE5_TSmutexunlock(HE5_HDFE_MutexStruct *);

/* GCTP function prototypes */
/* ------------------------ */
#include  "HE5_GctpFunc.h"

/* HDF5 Error handling function (from "H5private.h") */
/* ------------------------------------------------- */
herr_t   H5E_push(H5E_major_t maj_num, H5E_minor_t min_num, const char *func_name, const char *file_name, unsigned line, const char *desc);

/*
 **********************************************************************
 *         E R R O R   H A N D L I N G    I N T E R F A C E           *
 **********************************************************************
 */


/* File access/info routines */

hid_t    HE5_EHopen(const char *filename, unsigned flags,  hid_t access_id);
herr_t   HE5_EHclose(hid_t fid);
herr_t   HE5_EHgetversion(hid_t fid, char *version);
herr_t   HE5_EHchkfid(hid_t fid, const char *name, hid_t *HDFfid, hid_t *grpID, uintn *access);
herr_t   HE5_EHidinfo(hid_t fid, hid_t *HDFfid, hid_t *gid);
long     HE5_EHattrcat(hid_t fid,   char *grpname, char *objectlist, long *strbufsize);
long     HE5_EHinquire(const char *filename, char *grpname, char *objectlist, long *strbufsize);
hsize_t HE5_EHset_error_on(int flag, int err_level);

/* PROTOTYPES: functions to turn HDFEOS5 error printing off */
herr_t HE5_EHprint(char *errbuf, const char *file, unsigned line);


/* Metadata access/info routines */

char     *HE5_EHmetagroup(hid_t fid , char *structname, char *structcode, char *groupname, char *metaptrs[]);
herr_t   HE5_EHgetmetavalue(char *metaptrs[], char *parameter, char *retstr);
herr_t   HE5_EHinsertmeta(hid_t fid, const char *structname, char *structcode, long metacode, char *metastr, hsize_t metadata[]) ;
herr_t   HE5_EHmetalist(char *instring, char *outstring);
herr_t   HE5_EHupdatemeta(hid_t fid, const char *structname, char *structcode, long metacode, char *metastr, hsize_t metadata[]);

/* Object info routines */

long     HE5_EHcntOBJECT(char *metabuf[]);
long     HE5_EHcntGROUP(char *metabuf[]);
herr_t   HE5_EHattr( hid_t attgrpID, const char *attrname, hid_t ntype, hsize_t count[], char *wrcode, void *datbuf);
herr_t   HE5_EHattrinfo( hid_t attgrpID, const char *attrname, hid_t ntype[], hsize_t *count);
long     HE5_EHdatasetcat(hid_t fid,  char *grpname, char *objectlist, long  *strbufsize);
int      HE5_EHobj_info(hid_t loc_id, const char *name, void *opdata);

/* Utility routines */

long     HE5_EHparsestr(const char *instring, char delim, char *pntr[], size_t len[]);
long     HE5_EHstrwithin(char *target, char *search, char delim);
herr_t   HE5_EHloadliststr(char *ptr[], long nentries, char *liststr, char delim);
double   HE5_EHconvAng(double inAngle, int code);
herr_t   HE5_EHrevflds(char *dimlist, char *revdimlist);
herr_t   HE5_EHbisect(double(*func) (double[]), double funcParms[], long nParms, double limLft, double limRgt, double convCrit, double *root);
hid_t    HE5_EHdtype2mtype(hid_t dtypeID);
hid_t    HE5_EHdtype2numtype(hid_t dtype);
hid_t    HE5_EHconvdatatype(int fortdatatype);
herr_t   HE5_EHwriteglbattr(hid_t fid, const char *attrname, hid_t ntype, hsize_t  count[], void *datbuf);
herr_t   HE5_EHreadglbattr(hid_t fid, const char *attrname, void * datbuf);
herr_t   HE5_EHglbattrinfo(hid_t fid, const char *attrname, hid_t *ntype, hsize_t *count);
long     HE5_EHinqglbattrs(hid_t fid, char *attrnames, long *strbufsize);
herr_t   HE5_EHinqglbdatatype(hid_t fid, const char *attrname, hid_t *dtype, H5T_class_t *classid, H5T_order_t *order, size_t *size);

/* Data type conversion wrappers */    
                      
long                 HE5_EHhid2long(hid_t invalue);
long                 HE5_EHint2long(int invalue);
int                  HE5_EHhid2int(hid_t invalue);
hid_t                HE5_EHint2hid(int invalue);
hid_t                HE5_EHlong2hid(long invalue);
int                  HE5_EHlong2int(long invalue);

hsize_t              HE5_EHhid2hsize(hid_t invalue);
hsize_t              HE5_EHint2hsize(int invalue);
hsize_t              HE5_EHlong2hsize(long invalue);
hid_t                HE5_EHhsize2hid(hsize_t invalue);
long                 HE5_EHhsize2long(hsize_t invalue);
int                  HE5_EHhsize2int(hsize_t invalue);
hssize_t             HE5_EHhsize2hssize(hsize_t invalue);

hssize_t             HE5_EHhid2hssize(hid_t invalue);
hssize_t             HE5_EHint2hssize(int invalue);
hssize_t             HE5_EHlong2hssize(long invalue);
hid_t                HE5_EHhssize2hid(hssize_t invalue);
long                 HE5_EHhssize2long(hssize_t invalue);
int                  HE5_EHhssize2int(hssize_t invalue);
hsize_t              HE5_EHhssize2hsize(hssize_t invalue);

unsigned LONGLONG    HE5_EHint2ullong(int invalue);
long                 HE5_EHullong2long(unsigned LONGLONG invalue);

herr_t               HE5_EHchkptr(void *p, char *name);
herr_t               HE5_EHchkname(char *p, char *name);


/*
 **********************************************************************
 *                    S W A T H    I N T E R F A C E                  *
 **********************************************************************
 */
  


/* File/Swath access routines */

hid_t   HE5_SWopen(const char *filename, uintn flags);
hid_t   HE5_SWcreate(hid_t fid, const char *swathname);
hid_t   HE5_SWattach(hid_t fid, const char *swathname);
herr_t  HE5_SWdetach(hid_t swathID);
herr_t  HE5_SWclose(hid_t fid);


/* Definition routines */

herr_t  HE5_SWdefdim(hid_t swathID,  char *dimname, hsize_t dim);
herr_t  HE5_SWdefdimmap(hid_t swathID, char *geodim, char *datadim, hsize_t offset, hsize_t increment);
herr_t  HE5_SWdefidxmap(hid_t swathID, char *geodim, char *datadim, long index[]);
herr_t  HE5_SWdefgeofield(hid_t swathID, const char *fieldname, char *dimlist, char *maxdimlist, hid_t ntype , int merge);
herr_t  HE5_SWdefdatafield(hid_t swathID, const char *fieldname, char *dimlist, char *maxdimlist, hid_t ntype, int merge);
herr_t  HE5_SWdefchunk(hid_t swathID, int ndims, const hsize_t *dim);
herr_t  HE5_SWdefcomp(hid_t swathID, int compcode, int *compparm);
herr_t  HE5_SWdefcomchunk(hid_t swathID, int compcode, int *compparm, int ndims, const hsize_t *dim);
herr_t  HE5_SWsetfillvalue(hid_t swathID, char *fieldname, hid_t ntype, void *fillval);
herr_t  HE5_SWsetalias(hid_t swathID, char *fieldname, const char *aliaslist);
herr_t  HE5_SWdropalias(hid_t swathID, int fldgroup, const char *aliasname);
herr_t  HE5_SWfldrename(hid_t swathID, char *oldfieldname, const char *newfieldname);
herr_t  HE5_SWsetdimscale(hid_t swathID, char *fieldname, char *dimname,const hsize_t dimsize, hid_t numbertype, void * data);

/* I/O routines */

herr_t  HE5_SWwritedatameta(hid_t swathID, const char *fieldname, char *dimlist, hid_t mvalue);
herr_t  HE5_SWwriteattr(hid_t swathID, const char *attrname, hid_t ntype, hsize_t count[], void *datbuf);
herr_t  HE5_SWwritegrpattr(hid_t swathID, const char *attrname, hid_t ntype, hsize_t count[], void *datbuf);
herr_t  HE5_SWwritegeogrpattr(hid_t swathID, const char *attrname, hid_t ntype, hsize_t count[], void *datbuf);
herr_t  HE5_SWwritelocattr(hid_t swathID, const char *fieldname, const char *attrname, hid_t ntype, hsize_t count[], void *datbuf);
herr_t  HE5_SWreadattr(hid_t swathID, const char *attrname, void *datbuf);
herr_t  HE5_SWreadgrpattr(hid_t swathID, const char *attrname, void *datbuf);
herr_t  HE5_SWreadgeogrpattr(hid_t swathID, const char *attrname, void *datbuf);
herr_t  HE5_SWreadlocattr(hid_t swathID, const char *fieldname, const char *attrname, void *datbuf);
herr_t  HE5_SWwritefield(hid_t swathID, char *fieldname, const hssize_t start[], const hsize_t stride[], const hsize_t edge[],  void *data);
herr_t  HE5_SWreadfield(hid_t swathID, char *fieldname, const hssize_t start[], const hsize_t stride[], const hsize_t edge[],  void *data);
herr_t  HE5_SWwritegeometa(hid_t swathID, const char *fieldname, char *dimlist, hid_t mvalue);
herr_t  HE5_SWwritedscaleattr(hid_t swathID, const char *fieldname, const char *attrname, hid_t numtype, hsize_t  count[],  void *datbuf);

/* Inquiry routines */

herr_t  HE5_SWchunkinfo(hid_t swathID, char *fieldname, int *ndims, hsize_t dims[]);
hsize_t HE5_SWdiminfo(hid_t swathID, char *dimname);
herr_t  HE5_SWmapinfo(hid_t swathID, char *geodim, char *datadim, long *offset, long *increment);
hsize_t HE5_SWidxmapinfo(hid_t swathID, char *geodim, char *datadim, long index[]);
int     HE5_SWfldsrch(hid_t swathID, char *fieldname, hid_t *fieldID, int *rank,  hsize_t dims[], hid_t *typeID);
herr_t  HE5_SWfieldinfo(hid_t swathID, char *fieldname, int *rank, hsize_t dims[], hid_t ntype[], char *dimlist, char *maxdimlist);
herr_t  HE5_SWcompinfo(hid_t swathID, char *fieldname, int *compcode, int compparm[]);
herr_t  HE5_SWattrinfo(hid_t swathID, const char *attrname, hid_t *ntype, hsize_t *count);
herr_t  HE5_SWgrpattrinfo(hid_t swathID, const char *attrname, hid_t *ntype, hsize_t *count);
herr_t  HE5_SWgeogrpattrinfo(hid_t swathID, const char *attrname, hid_t *ntype, hsize_t *count);
herr_t  HE5_SWlocattrinfo(hid_t swathID, const char *fieldname, const char *attrname, hid_t *ntype, hsize_t *count);
herr_t  HE5_SWinqdatatype(hid_t swathID, const char *fieldname, const char *attrname, int group, hid_t *dtype, H5T_class_t *classid, H5T_order_t *order, size_t *size);
long    HE5_SWinqdims(hid_t swathID, char *dimnames, hsize_t dims[]);
long    HE5_SWinqmaps(hid_t swathID, char *dimmaps, long offset[], long increment[]);
long    HE5_SWinqidxmaps(hid_t swathID, char *idxmaps, hsize_t idxsizes[]);
long    HE5_SWinqgeofields(hid_t swathID, char *fieldlist, int rank[], hid_t ntype[]);
long    HE5_SWinqdatafields(hid_t swathID, char *fieldlist, int rank[], hid_t ntype[]);
long    HE5_SWinqattrs(hid_t swathID, char *attrnames, long *strbufsize);
long    HE5_SWinqgrpattrs(hid_t swathID, char *attrnames, long *strbufsize);
long    HE5_SWinqgeogrpattrs(hid_t swathID, char *attrnames, long *strbufsize);
long    HE5_SWinqlocattrs(hid_t swathID, const char *fieldname, char *attrnames, long *strbufsize);
long    HE5_SWnentries(hid_t swathID, int entrycode, long *strbufsize);
long    HE5_SWinqswath(const char *filename, char *swathlist, long *strbufsize);
herr_t  HE5_SWregioninfo(hid_t swathID, hid_t regionID, char *fieldname, hid_t *ntype, int *rank, hsize_t dims[], size_t *size);
herr_t  HE5_SWperiodinfo(hid_t swathID, hid_t periodID, char *fieldname, hid_t *ntype, int *rank, hsize_t dims[], size_t *size);
herr_t  HE5_SWgeomapinfo(hid_t swathID, char *geodim);
herr_t  HE5_SWgetfillvalue(hid_t swathID, char *fieldname, void *fillval);
herr_t  HE5_SWaliasinfo(hid_t swathID, int fldgroup, const char *aliasname, int *length, char *buffer);
long    HE5_SWinqdfldalias(hid_t swathID, char *fldalias, long *strbufsize);
long    HE5_SWinqgfldalias(hid_t swathID, char *fldalias, long *strbufsize);
long    HE5_SWgetaliaslist(hid_t swathID, int fldgroup, char *aliaslist, long *strbufsize);
long    HE5_SWgetdimscale(hid_t swathID, char *fieldname, char *dimname, hsize_t *dimsize, hid_t *numbertype, void * databuff);
herr_t  HE5_SWdscaleattrinfo(hid_t swathID, const char *fieldname, const char *attrname, hid_t *ntype, hsize_t *count);
herr_t  HE5_SWreaddscaleattr(hid_t swathID, const char *fieldname, const char *attrname, void *datbuf);
long    HE5_SWinqdscaleattrs(hid_t swathID, const char *fieldname, char *attrnames, long *strbufsize);

/* Subsetting/Retrieving routines */

hid_t   HE5_SWdefboxregion(hid_t swathID, double cornerlon[], double cornerlat[], int mode);
hid_t   HE5_SWdefvrtregion(hid_t swathID, hid_t regionID, char *vertObj, double range[]);
hid_t   HE5_SWregionindex(hid_t swathID, double cornerlon[], double cornerlat[], int mode, char *geodim, hsize_t idxrange[]);
hid_t   HE5_SWdupregion(hid_t oldregionID);
hid_t   HE5_SWdeftimeperiod(hid_t swathID, double starttime, double stoptime, int mode);
herr_t  HE5_SWextractregion(hid_t swathID, hid_t regionID, char *fieldname, int externalflag, void *buffer);
herr_t  HE5_SWextractperiod(hid_t swathID, hid_t periodID, char *fieldname, int externalflag, void *buffer);
long    HE5_SWupdateidxmap(hid_t swathID, hid_t regionID, long indexin[], long indexout[], long indicies[]);
herr_t  HE5_SWupdatescene(hid_t swathID, hid_t regionID);
herr_t  HE5_SWindexinfo(hid_t regionID, char *object, int *rank, char *dimlist, hsize_t *indices[HE5_DTSETRANKMAX]);

/*
 ********************************
 *     PROFILE INTERFACE        *
 ********************************
 */ 

herr_t  HE5_PRdefine(hid_t swathID, const char *profilename, char *dimlist, char *maxdimlist, hid_t datatype_id);
herr_t  HE5_PRwrite(hid_t swathID, const char *profilename, const hssize_t start[], const hsize_t stride[], const hsize_t edge[], size_t size, void *buffer);
herr_t  HE5_PRread(hid_t swathID, const char *profilename, const hssize_t start[], const hsize_t stride[], const hsize_t edge[], void *buffer);
herr_t  HE5_PRreclaimspace(hid_t swathID, const char *profilename, void *buffer);
long    HE5_PRinquire(hid_t swathID, char *profnames, int *rank, H5T_class_t *classID);
herr_t  HE5_PRinfo(hid_t swathID, const char *profname, int *rank, hsize_t dims[], hsize_t maxdims[], hid_t *ntype, char *dimlist, char *maxdimlist);
herr_t  HE5_PRwritegrpattr(hid_t swathID, const char *attrname, hid_t ntype, hsize_t count[], void *datbuf);
herr_t  HE5_PRreadgrpattr(hid_t swathID, const char *attrname, void *datbuf);
herr_t  HE5_PRgrpattrinfo(hid_t swathID, const char *attrname, hid_t *ntype, hsize_t *count);
long    HE5_PRinqgrpattrs(hid_t swathID, char *attrnames, long *strbufsize);


/*
 *******************************
 *   EXTERNAL DATA FILES       *
 *******************************
 */

herr_t  HE5_SWsetextdata(hid_t swathID, const char *filelist, off_t offset[], hsize_t size[]);
int     HE5_SWgetextdata(hid_t swathID, char *fieldname, size_t namelength, char *filelist, off_t offset[], hsize_t size[]);


/*
 *******************************
 *  MOUNTING EXTERNAL FILES    *
 *******************************
 */
  
hid_t   HE5_SWmountexternal(hid_t swathID, int fldgroup, const char *extfilename);
herr_t  HE5_SWunmount(hid_t swathID, int fldgroup, hid_t fileID);
herr_t  HE5_SWreadexternal(hid_t swathID, int fldgroup, const char *fieldname, void *buffer);


/*
 **********************************************************************
 *                      G R I D    I N T E R F A C E                  *
 **********************************************************************
 */


/* File/Grid access routines */

hid_t    HE5_GDopen(const char *filename, uintn flags);
hid_t    HE5_GDcreate(hid_t fid, const char *gridname, long xdimsize, long ydimsize, double upleftpt[], double lowrightpt[]);
hid_t    HE5_GDattach(hid_t fid, const char *gridname);
herr_t   HE5_GDdetach(hid_t gridID);
herr_t   HE5_GDclose(hid_t fid);


/* Definition routines */

herr_t   HE5_GDdefdim(hid_t gridID,  char *dimname, hsize_t dim);
herr_t   HE5_GDdefproj(hid_t gridID, int projcode, int zonecode, int spherecode, double projparm[]);
herr_t   HE5_GDdefcomp(hid_t gridID, int compcode, int compparm[]);
herr_t   HE5_GDdeftile(hid_t gridID, int tilecode, int tilerank, const hsize_t *tiledims);
herr_t   HE5_GDdefcomtile(hid_t gridID, int compcode, int compparm[], int tilerank, const hsize_t *tiledims);
herr_t   HE5_GDdeforigin(hid_t gridID, int origincode);
herr_t   HE5_GDdefpixreg(hid_t gridID, int pixregcode);
herr_t   HE5_GDdeffield(hid_t gridID, const char *fieldname, char *dimlist, char *maxdimlist, hid_t ntype, int merge);
herr_t   HE5_GDsetfillvalue(hid_t gridID, const char *fieldname, hid_t ntype, void *fillval);
herr_t   HE5_GDsetalias(hid_t gridID, char *fieldname, const char *aliaslist);
herr_t   HE5_GDdropalias(hid_t gridID, int fldgroup, const char *aliasname);
herr_t   HE5_GDsetdimscale(hid_t gridID, char *fieldname, char *dimname,const hsize_t dimsize, hid_t numbertype, void * data);


/* I/O routines */

herr_t   HE5_GDwritefieldmeta(hid_t gridID, const char *fieldname, char *dimlist, hid_t ntype);
herr_t   HE5_GDwritefield(hid_t gridID, const char *fieldname, const hssize_t start[], const hsize_t stride[], const hsize_t edge[], void *data);
herr_t   HE5_GDreadfield(hid_t gridID, const char *fieldname, const hssize_t start[], const hsize_t stride[], const hsize_t edge[], void * buffer);
herr_t   HE5_GDwriteattr(hid_t gridID, const char *attrname, hid_t ntype, hsize_t  count[], void *datbuf);
herr_t   HE5_GDwritegrpattr(hid_t gridID, const char *attrname, hid_t ntype, hsize_t count[], void *datbuf);
herr_t   HE5_GDwritelocattr(hid_t gridID, const char *fieldname, const char *attrname, hid_t ntype, hsize_t  count[], void *datbuf);
herr_t   HE5_GDreadattr(hid_t gridID, const char *attrname, void *datbuf);
herr_t   HE5_GDreadgrpattr(hid_t gridID, const char *attrname, void *datbuf);
herr_t   HE5_GDreadlocattr(hid_t gridID, const char *fieldname, const char *attrname, void *datbuf);
herr_t   HE5_GDblkSOMoffset(hid_t gridID, long offset[], hsize_t count[], char *code);


/* Inquiry routines */

long     HE5_GDinqgrid(const char *filename, char *gridlist, long *strbufsize);
hsize_t  HE5_GDdiminfo(hid_t gridID, char *dimname);
herr_t   HE5_GDgridinfo(hid_t gridID, long *xdimsize, long *ydimsize, double upleftpt[], double lowrightpt[]);
herr_t   HE5_GDprojinfo(hid_t gridID, int *projcode, int *zonecode, int *spherecode, double projparm[]);
herr_t   HE5_GDorigininfo(hid_t gridID, int *origincode);
herr_t   HE5_GDpixreginfo(hid_t gridID, int *pixregcode);
herr_t   HE5_GDcompinfo(hid_t gridID, const char *fieldname, int *compcode, int compparm[]);
herr_t   HE5_GDfieldinfo(hid_t gridID, const char *fieldname, int *rank, hsize_t dims[], hid_t ntype[], char *dimlist, char *maxdimlist);
herr_t   HE5_GDregioninfo(hid_t gridID, hid_t regionID, const char *fieldname, hid_t *ntype, int *rank, hsize_t dims[], long *size, double upleftpt[], double lowrightpt[]);
long     HE5_GDnentries(hid_t gridID, int entrycode, long *strbufsize);
int      HE5_GDinqdims(hid_t gridID, char *dimnames, hsize_t  dims[]);
herr_t   HE5_GDattrinfo(hid_t gridID, const char *attrname, hid_t *ntype, hsize_t *count);
herr_t   HE5_GDgrpattrinfo(hid_t gridID, const char *attrname, hid_t *ntype, hsize_t *count);
herr_t   HE5_GDlocattrinfo(hid_t gridID, const char *fieldname, const char *attrname, hid_t *ntype, hsize_t *count);
long     HE5_GDinqattrs(hid_t gridID, char *attrnames, long *strbufsize);
long     HE5_GDinqgrpattrs(hid_t gridID, char *attrnames, long *strbufsize);
long     HE5_GDinqlocattrs(hid_t gridID, const char *fieldname, char *attrnames, long *strbufsize);
int      HE5_GDinqfields(hid_t gridID, char *fieldlist, int rank[], hid_t ntype[]);
herr_t   HE5_GDinqdatatype(hid_t gridID, const char *fieldname, const char *attrname, int fieldgroup, hid_t *dtype, H5T_class_t *classid, H5T_order_t *order, size_t *size);
herr_t   HE5_GDgetfillvalue(hid_t gridID, const char *fieldname, void *fillval);
herr_t   HE5_GDtileinfo(hid_t gridID, char *fieldname, int *tilecode, int *tilerank, hsize_t tiledims[]);
herr_t   HE5_GDaliasinfo(hid_t gridID, int fldgroup, const char *aliasname, int *length, char *buffer);
long     HE5_GDinqfldalias(hid_t gridID, char *fldalias, long *strbufsize);
long     HE5_GDgetaliaslist(hid_t gridID, int fldgroup, char *aliaslist, long *strbufsize);
long     HE5_GDgetdimscale(hid_t gridID, char *fieldname, char *dimname, hsize_t *dimsize, hid_t *numbertype, void * databuff);
herr_t   HE5_GDwritedscaleattr(hid_t gridID, const char *fieldname, const char *attrname, hid_t numtype, hsize_t  count[],  void *datbuf);
herr_t   HE5_GDdscaleattrinfo(hid_t gridID, const char *fieldname, const char *attrname, hid_t *ntype, hsize_t *count);
herr_t   HE5_GDreaddscaleattr(hid_t gridID, const char *fieldname, const char *attrname, void *datbuf);
long     HE5_GDinqdscaleattrs(hid_t gridID, const char *fieldname, char *attrnames, long *strbufsize);

/* Subsetting/Retrieving routines */

hid_t    HE5_GDdefboxregion(hid_t gridID, double cornerlon[], double cornerlat[]);
hid_t    HE5_GDdefvrtregion(hid_t gridID, hid_t regionID, char *vertObj, double range[]);
herr_t   HE5_GDdeftimeperiod(hid_t gridID, hid_t periodID, double starttime, double stoptime);
herr_t   HE5_GDextractregion(hid_t gridID, hid_t regionID, const char *fieldname, void *buffer);
hid_t    HE5_GDdupregion(hid_t oldregionID);
herr_t   HE5_GDgetpixels(hid_t gridID, long nLonLat, double lonVal[], double latVal[], long pixRow[], long pixCol[]);
long     HE5_GDgetpixvalues(hid_t gridID, long nPixels, long pixRow[], long pixCol[], const char *fieldname, void * buffer);
long     HE5_GDinterpolate(hid_t gridID, long nValues, double lonVal[], double latVal[], const char *fieldname, double interpVal[]);


/* Utility routine */

herr_t   HE5_GDij2ll(int, int, double[], int, long, long, double[], double[], long, long[], long[], double[], double[], int, int);
herr_t   HE5_GDll2ij(int, int, double[], int, long, long, double[], double[], long, double[], double[], long[], long[], double[], double[]);
herr_t   HE5_GDrs2ll(int projcode, double projparm[], long xdimsize, long ydimsize, double upleft[], double lowright[], int npnts, double r[], double s[], double longitude[], double latitude[], int pixcen, int pixcnr);


/*
 *******************************
 *   EXTERNAL DATA FILES       *
 *******************************
 */

herr_t  HE5_GDsetextdata(hid_t gridID, const char *filelist, off_t offset[], hsize_t size[]);
int     HE5_GDgetextdata(hid_t gridID, char *fieldname, size_t namelength, char *filelist, off_t offset[], hsize_t size[]);



/*
 **********************************************************************
 *                    P O I N T    I N T E R F A C E                  *
 **********************************************************************
 */


/* File/Point access routine */

hid_t    HE5_PTopen(const char *filename, uintn flags);
hid_t    HE5_PTcreate(hid_t fid, const char *pointname);
hid_t    HE5_PTattach(hid_t fid, const char *pointname);
herr_t   HE5_PTdetach(hid_t pointID);
herr_t   HE5_PTclose(hid_t fid);


/* Definition routines */

herr_t   HE5_PTdeflevel(hid_t pointID, const char *levelname, HE5_CmpDTSinfo *levelinfo);
herr_t   HE5_PTdeflinkage(hid_t pointID, char *parent, char *child, char *linkfield);


/* I/O routines */

herr_t   HE5_PTwritelevel(hid_t pointID, int level, hsize_t count[], size_t *size, void *data);
herr_t   HE5_PTupdatelevel(hid_t pointID, int level, char *fieldlist, hsize_t nrec, hssize_t recs[], void *data);
herr_t   HE5_PTreadlevel(hid_t pointID, int level, HE5_CmpDTSinfo *inStruct, size_t *size, void *datbuf);
herr_t   HE5_PTwriteattr(hid_t pointID, const char *attrname, hid_t ntype, hsize_t count[], void * datbuf);
herr_t   HE5_PTwritegrpattr(hid_t pointID, const char *attrname, hid_t ntype, hsize_t count[], void * datbuf);
herr_t   HE5_PTwritelocattr(hid_t pointID, const char *levelname, const char *attrname, hid_t ntype, hsize_t count[], void * datbuf);
herr_t   HE5_PTreadattr(hid_t pointID, const char *attrname, void * datbuf);
herr_t   HE5_PTreadgrpattr(hid_t pointID, const char *attrname, void * datbuf);
herr_t   HE5_PTreadlocattr(hid_t pointID, const char *levelname, const char *attrname, void *datbuf);

/* Inquiry routines */

hsize_t  HE5_PTnrecs(hid_t pointID, int level);
int      HE5_PTnlevels(hid_t pointID);
int      HE5_PTnfields(hid_t pointID, int level, char *fieldlist, long *strbufsize);
int      HE5_PTlevelindx(hid_t pointID, const char *levelname);
herr_t   HE5_PTgetlevelname(hid_t pointID, int level, char *levelname, long *strbufsize);
herr_t   HE5_PTbcklinkinfo(hid_t pointID, int level, char *linkfield);
herr_t   HE5_PTfwdlinkinfo(hid_t pointID, int level, char *linkfield);
herr_t   HE5_PTlevelinfo(hid_t pointID, int level, HE5_CmpDTSinfo *info);
herr_t   HE5_PTinqdatatype(hid_t pointID, const char *levelname, const char *attrname, int fieldgroup, hid_t *dtype, H5T_class_t *classid, H5T_order_t *order, size_t *size);
int      HE5_PTinqpoint(const char *filename, char *pointlist, long *strbufsize);
herr_t   HE5_PTgetrecnums(hid_t pointID, int inlevel, int outlevel, hsize_t inNrec, hssize_t inRecs[], hsize_t * outNrec, hssize_t outRecs[]);
herr_t   HE5_PTattrinfo(hid_t pointID, const char *attrname, hid_t *ntype, hsize_t *count);
herr_t   HE5_PTgrpattrinfo(hid_t pointID, const char *attrname, hid_t *ntype, hsize_t *count);
herr_t   HE5_PTlocattrinfo(hid_t pointID, const char *levelname, const char *attrname, hid_t *ntype, hsize_t *count);
long     HE5_PTinqattrs(hid_t pointID, char *attrnames, long *strbufsize);
long     HE5_PTinqgrpattrs(hid_t pointID, char *attrnames, long *strbufsize);
long     HE5_PTinqlocattrs(hid_t pointID, const char *levelname, char *attrnames, long *strbufsize);


/*
 **********************************************************************
 *                    Z A    I N T E R F A C E                        *
 **********************************************************************
 */
 
 
/* File/ZA access routines */
 
hid_t   HE5_ZAopen(const char *filename, uintn flags);
hid_t   HE5_ZAcreate(hid_t fid, const char *zaname);
hid_t   HE5_ZAattach(hid_t fid, const char *zaname);
herr_t  HE5_ZAdetach(hid_t zaID);
herr_t  HE5_ZAclose(hid_t fid);
 
 
/* Definition routines */
 
herr_t  HE5_ZAdefdim(hid_t zaID,  char *dimname, hsize_t dim);
herr_t  HE5_ZAdefine(hid_t zaID, const char *za_name, char *dimlist, char *maxdimlist, hid_t dtype);
herr_t  HE5_ZAdefchunk(hid_t zaID, int ndims, const hsize_t *dim);
herr_t  HE5_ZAdefcomp(hid_t zaID, int compcode, int *compparm);
herr_t  HE5_ZAdefcomchunk(hid_t zaID, int compcode, int *compparm, int ndims, const hsize_t *dim);
herr_t  HE5_ZAsetfillvalue(hid_t zaID, char *fieldname, hid_t ntype, void *fillval);
herr_t  HE5_ZAsetalias(hid_t zaID, char *fieldname, const char *aliaslist);
herr_t  HE5_ZAdropalias(hid_t zaID, int fldgroup, const char *aliasname);
herr_t  HE5_ZAfldrename(hid_t zaID, char *oldfieldname, const char *newfieldname);
herr_t  HE5_ZAsetdimscale(hid_t zaID, char *fieldname, char *dimname, const hsize_t dimsize, hid_t numbertype, void * data);
long    HE5_ZAgetdimscale(hid_t zaID, char *fieldname, char *dimname, hsize_t *dimsize, hid_t *numbertype, void * databuff);

/* I/O routines */
 
herr_t  HE5_ZAwritedatameta(hid_t zaID, const char *fieldname, char *dimlist, hid_t mvalue);
herr_t  HE5_ZAwriteattr(hid_t zaID, const char *attrname, hid_t ntype, hsize_t count[], void *datbuf);
herr_t  HE5_ZAwritegrpattr(hid_t zaID, const char *attrname, hid_t ntype, hsize_t count[], void *datbuf);
herr_t  HE5_ZAwritelocattr(hid_t zaID, const char *fieldname, const char *attrname, hid_t ntype, hsize_t count[], void *datbuf);
herr_t  HE5_ZAreadattr(hid_t zaID, const char *attrname, void *datbuf);
herr_t  HE5_ZAreadgrpattr(hid_t zaID, const char *attrname, void *datbuf);
herr_t  HE5_ZAreadlocattr(hid_t zaID, const char *fieldname, const char *attrname, void *datbuf);
herr_t  HE5_ZAwrite(hid_t zaID, char *za_name, const hssize_t start[], const hsize_t stride[], const hsize_t count[],  void *datbuf);
herr_t  HE5_ZAread(hid_t zaID, char *za_name, const hssize_t start[], const hsize_t stride[], const hsize_t count[],  void *datbuf);
herr_t HE5_ZAreaddscaleattr(hid_t zaID, const char *fieldname, const char *attrname, void *datbuf);
herr_t HE5_ZAwritedscaleattr(hid_t zaID, const char *fieldname, const char *attrname, hid_t numtype, hsize_t  count[],  void *datbuf); 
 
/* Inquiry routines */
 
hsize_t HE5_ZAdiminfo(hid_t zaID, char *dimname);
int     HE5_ZAfldsrch(hid_t zaID, char *fieldname, hid_t *fieldID, int *rank,  hsize_t dims[], hid_t *typeID);
herr_t  HE5_ZAinfo(hid_t zaID, char *za_name, int *rank, hsize_t dims[], hid_t dtype[], char *dimlist, char *maxdimlist);
herr_t  HE5_ZAcompinfo(hid_t zaID, char *fieldname, int *compcode, int compparm[]);
herr_t  HE5_ZAattrinfo(hid_t zaID, const char *attrname, hid_t *ntype, hsize_t *count);
herr_t  HE5_ZAgrpattrinfo(hid_t zaID, const char *attrname, hid_t *ntype, hsize_t *count);
herr_t  HE5_ZAlocattrinfo(hid_t zaID, const char *fieldname, const char *attrname, hid_t *ntype, hsize_t *count);
herr_t  HE5_ZAinqdatatype(hid_t zaID, const char *fieldname, const char *attrname, int group, hid_t *dtype, H5T_class_t *classid, H5T_order_t *order, size_t *size);
long    HE5_ZAinqdims(hid_t zaID, char *dimnames, hsize_t dims[]);
long    HE5_ZAinquire(hid_t zaID, char *za_name_list, int rank[], hid_t dtype[]);
long    HE5_ZAinqattrs(hid_t zaID, char *attrnames, long *strbufsize);
long    HE5_ZAinqgrpattrs(hid_t zaID, char *attrnames, long *strbufsize);
long    HE5_ZAinqlocattrs(hid_t zaID, const char *fieldname, char *attrnames, long *strbufsize);
long    HE5_ZAnentries(hid_t zaID, int entrycode, long *strbufsize);
long    HE5_ZAinqza(const char *filename, char *zalist, long *strbufsize);
herr_t  HE5_ZAgetfillvalue(hid_t zaID, char *fieldname, void *fillval);
herr_t  HE5_ZAaliasinfo(hid_t zaID, int fldgroup, const char *aliasname, int *length, char *buffer);
long    HE5_ZAinqfldalias(hid_t zaID, char *fldalias, long *strbufsize);
herr_t  HE5_ZAchunkinfo(hid_t zaID, char *fieldname, int *ndims, hsize_t dims[]);
long    HE5_ZAgetaliaslist(hid_t zaID, int fldgroup, char *aliaslist, long *strbufsize);
herr_t HE5_ZAdscaleattrinfo(hid_t zaID, const char *fieldname, const char *attrname, hid_t *ntype, hsize_t *count);
long HE5_ZAinqdscaleattrs(hid_t zaID, const char *fieldname, char *attrnames, long *strbufsize);

/*
 *******************************
 *   EXTERNAL DATA FILES       *
 *******************************
 */
 
herr_t  HE5_ZAsetextdata(hid_t zaID, const char *filelist, off_t offset[], hsize_t size[]);
int     HE5_ZAgetextdata(hid_t zaID, char *fieldname, size_t namelength, char *filelist, off_t offset[], hsize_t size[]);
 
 
/*
 *******************************
 *  MOUNTING EXTERNAL FILES    *
 *******************************
 */
 
hid_t   HE5_ZAmountexternal(hid_t zaID, int fldgroup, const char *extfilename);
herr_t  HE5_ZAunmount(hid_t zaID, int fldgroup, hid_t fileID);
herr_t  HE5_ZAreadexternal(hid_t zaID, int fldgroup, const char *fieldname, void *buffer);



#ifdef __cplusplus
}
#endif


#endif  /* #ifndef HE5_HDFEOSDEF_H_ */












