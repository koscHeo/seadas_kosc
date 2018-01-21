/***************************************************************************
 *		Mapic.h M-API utilities compile header file		*
 *									*/
/**************************************************************************
 * !C-INC
 * 
 * Purpose:	Utilities header file used by all subroutines contained in the
 *		MODIS Applications Programming Interface (API).
 *
 * !Description: The Header file mapic.h is part of a larger software system 
 *               called the MODIS Applications Programming Interface (API)
 *               Utility, abbreviated M-API.  The M-API Utility consists of
 *               subroutines which allow MODIS Science Team-supplied software
 *               to read in Level 1B radiance bands and write out output
 *               products and metadata to HDF files.  The functionality of the
 *               M-API is defined in the MODIS API User's Guide, Version 2.3.
 *
 * !Input parameters:  None
 *
 * !Output parameters: None
 *
 * !Revision history:
 *  $Log: mapic.h,v $
 *  Revision 6.1  2010/07/13 21:44:30  kuyper
 *  Corrected to reflect the fact that parse_string does modify the string it is
 *    passed.
 *
 *  Revision 5.1  2005/04/11 15:16:13  vlin
 *  pointer arguments to a M-API function that points to a
 *  read-only array of data are declared 'const'.
 *  vlin@saicmodis.com
 *
 *      Revision 01.00 1995/05/05
 *      Joan Baden
 *      Original development.
 *
 * !Team-unique header:
 *
 *      This software is developed by the MODIS Science Data Support
 *      Team for the National Aeronautics and Space Administration,
 *      Goddard Space Flight Center, under contract NAS5-32373.
 *
 *      HDF portions developed at the National Center for Supercomputing
 *      Applications at the University of Illinois at Urbana-Champaign.
 *
 * !END
 *********************************************************************/

#ifndef MAPIC_H
#define MAPIC_H

#ifdef  DEC_ALPHA
   #include <rpc/types.h>
#endif

#define PGS_MET_COMPILE
#define NEW_VERSION     /* if you use HDF 4.0r1p1 or older HDF
                           versions comment out this macro */
#ifdef NEW_VERSION
#define         VARIABLE        _HDF_VARIABLE
#define         DIMENSION       _HDF_DIMENSION
#define         CDF             _HDF_CDF
#endif

#include "hdf.h"
#include "mfhdf.h"
#undef ATTRIBUTE
/*#include "local_nc.h" */
/*#include "PGS_MET.h"
#undef VOIDP*/
#include "mapi.h"

typedef struct sdsinfo {
	int32	*dimsizes;
	int32	rank;
	int32	ntype;
}SDSINFO;

typedef struct vdinfo {
	char	access;		/* 'r' for read only, 'w' for write but the
				the vdata also readable */
	char	*read_fields;	/* the fields set for read */
	int	size;		/* size (nbytes) for the read fields */
}VDINFO;

#define HDc2fstr    MHDc2fstr

#define DATA_LABEL  "NOT EMPTY"

#define	DDS_PER_BLOCK	0

#define	NULLstr(c) (((c)==NULL) || (*(c) == '\0'))

#define MODFILEstruNexist(file)   (((file)==NULL) || (*(file)==NULL))

#define NULLMODFIL(file)     (((file) == NULL) || ((file)->access == 0))

#define MAX_REC    32767

#define RW_OK      1

#define WRITE_OK  -1
 
#define NOT_OK     0

#define APPEND_VDATA    -1

#define NO_OBJECT       -2

#define FDATATYPELENMAX   14
#define FI8		  "INTEGER*1"
#define FI16		  "INTEGER*2"
#define FI32		  "INTEGER*4"
#define FI64		  "INTEGER*8"
#define FR32		  "REAL*4"
#define FR64		  "REAL*8"
#define FUI8		  "UINTEGER*1"
#define FTXT		  "CHARACTER*(*)"
#define FUI16		  "UINTEGER*2"
#define FUI32		  "UINTEGER*4"
#define FUI64		  "UINTEGER*8"
#define P_SDID            0
#define P_HDFID           1
#define P_ACCESS          2
#define P_ADDR            3
#define MODFILLEN         5

/* Function name redefined for c routines called by FORTRAN */
#define nccrmtbl   FNAME(ccrmtbl) 
#define ncgmflds   FNAME(cgmflds)
#define ncgmtbl    FNAME(cgmtbl)
#define ncpmtbl    FNAME(cpmtbl) 
#define ncopmfil   FNAME(copmfil)
#define nccrmgrp   FNAME(ccrmgrp)
#define ncclmfil   FNAME(cclmfil)
#define ncadmgrp   FNAME(cadmgrp)
#define ncsrmgrp   FNAME(csrmgrp)
#define nccrmar    FNAME(ccrmar)
#define ncgmarin   FNAME(cgmarin)
#define ncgmecin   FNAME(cgmecin)
#define ncpmarin   FNAME(cpmarin)
#define ncpmfin    FNAME(cpmfin)
#define ncgmfin    FNAME(cgmfin)
#define nccpmfil   FNAME(ccpmfil)
#define ncpmdnam   FNAME(cpmdnam)
#define ncemobj	   FNAME(cemobj)
#define ncgmdnam   FNAME(cgmdnam)
#define ncgmardm   FNAME(cgmardm)
#define ncpmar	   FNAME(cpmar)
#define ncgmar	   FNAME(cgmar)
#define ncgmdmin   FNAME(cgmdmin)
#define ncpmdmin   FNAME(cpmdmin)
#define ncgmhoid   FNAME(cgmhoid)
#define nccmfh     FNAME(ccmfh)
#define ncrmfh     FNAME(crmfh)

/*#define SYSPRINT    1 */
   
#ifndef  SYSPRINT
#include <PGS_SMF.h>
#include "PGS_MAPI_39603.h"
#define  MAPIERR(b,fn)      PGS_SMF_SetDynamicMsg(MAPI_E_ERR,b,fn)
#define  MAPIWARN(b,fn)     PGS_SMF_SetDynamicMsg(MAPI_W_ARN,b,fn)
#else
#define PGS_SMF_MAX_MSGBUF_SIZE        481      /* max. 480 chars, 1 char for '\0'  */
#define  MAPIERR(b,fn)      fprintf(stderr,"%s",b)
#define  MAPIWARN(b,fn)     fprintf(stderr,"%s",b)
#endif


   int VFdatatypes(     	int32 vdata_id,
				int32 *stringlen,
	  			char *data_type);

   short int emptyVdata(        MODFILE *file,
                                int32 vdata_ref);

   short int Vdatattr_name(     char const *attribute,
                                int32 vdata_ref,
                                char *attribute_name);

   int32 datatype_to_DFNT(	char const *datatype);
	
   int DFNT_to_datatype(	int32   dfnt,
				char	*datatype);

   short int SDS_footprintOK_s(	int32	 sds_id,
				int32 start[ ],
				int32 dimsizes[ ]); 

   short int SDS_footprintOK(	int32 const sds_dimsizes[],
				int32 sds_rank,
				int32 const start[ ],
				int32 const dimsizes[ ]);

   int32	VSnametoref(	int32 	hdf_id,
				char 	*tablename);

   int		parse_string(	char		target_string[],
				char const	delimiter[],
				int		n_tokens,
				char		*tokens[]);

   short int    set_Vhasdata(   MODFILE *file, 
                                int32   vdata_ref);

   void ncgmardm (intf *modfil, _fcd arrnm, intf *lar, _fcd grpnm, intf
		  *lgr, _fcd dtype, intf *ldt, intf *rank, intf dims[], intf *ret); 

   void ncgmarin (intf modfil[MODFILLEN], _fcd arrnm, intf *lar, _fcd 
		  grpnm, intf *lgr, _fcd attr, intf *lat, _fcd dtype,
		  intf *ldt, intf *nms, VOIDP value, intf *ret);

   void nccrmar (intf *modfil, _fcd arrnm, intf *lar, _fcd grpnm, intf
		 *lgr, _fcd dtype, intf *ldt, intf *rank, intf *dims, intf *ret);

   void ncgmdmin(intf modfil[], _fcd arrnm, intf *lar, _fcd grpnm, intf
		 *lgr, intf *dim, _fcd attrnm, intf *lat, _fcd dtype, intf
		 *ldt, intf *nms, void *value, intf *ret);

   void ncpmdmin(intf modfil[], _fcd arrnm, intf *lar, _fcd grpnm, intf
		 *lgr, intf *dim, _fcd attrnm, intf *lat, _fcd dtype, intf
		 *ldt, intf *nms, void *value, intf *ret);

   void ncgmtbl (intf *mfile,
                 _fcd tbname, intf *ltb,
                 _fcd group,  intf *lgr,
                 _fcd field,  intf *lfi,
                 intf *start, intf *recno,
                 intf *bsize, VOIDP data,
                 intf *ret);

   void cpmarin(intf modfil[MODFILLEN], _fcd arrnm, intf *lar, _fcd
		grpnm, int *lgr, _fcd attr, intf *lat, _fcd dtype,
		intf *ldt, intf *nms, VOIDP value, intf *ret);

   void ncpmtbl (intf *mfile,
                 _fcd tbname, intf *ltb, 
		 _fcd group,  intf *lgr,
		 intf *start, intf *recno,
		 VOIDP data,  intf *ret);

   void nccrmtbl (intf *mfile,
		  _fcd tbname, intf *ltb,
		  _fcd class,  intf *lcl,
		  _fcd group,  intf *lgr,
		  _fcd field,  intf *lfi,
		  _fcd dtype,  intf *ldt,
		  intf *ret);

   void ncgmflds (intf *mfile,
		  _fcd tbname, intf *ltb,
		  _fcd group,  intf *lgr,
		  intf *strln, intf *recno,
		  intf *fldno, _fcd fldnm,
		  intf *lfl,   _fcd dtype,
		  intf *ldt,   _fcd clss,
		  intf *lcl,   intf *ret);

   void ncpmar (intf modfil[], _fcd arrnm, intf *lar, _fcd grpnm, intf *lgr,
		intf start[], intf dims[], void *data, intf *ret);

   void ncgmar (intf modfil[], _fcd arrnm, intf *lar, _fcd grpnm, intf *lgr,
		intf start[], intf dims[], void *data, intf *ret);

   int MTYPEc2f(char *c_data_type,
		char *f_data_type,
		int32 *f_length);

   int MTYPEf2c(char *f_data_type,
		char *c_data_type,
		int32 *c_length);

   void ncopmfil(_fcd fname, intf *lfn,
		_fcd access,intf *lac,
		intf modfil[MODFILLEN], intf *ret);

   void ncclmfil(intf modfil[MODFILLEN], intf *ret);

   void nccrmgrp(intf modfil[MODFILLEN], _fcd groupname,
		intf *lgr, _fcd classname, intf *lcl,
		intf *ret);

   void ncgmecin(intf modfil[MODFILLEN], _fcd pvlname, intf *lpv,
		 _fcd pname, intf *lpn, _fcd dtype, intf *ldt,
		 intf *nms, VOIDP pvalue, intf *ret);

   void nccpmfil(intf modfil[MODFILLEN], _fcd mdHandles, intf *lhandles,
		 _fcd HDFattrnms, intf *lattrs, intf *NumHandles, intf *ret);

   void cgmdnam(intf *modfil, _fcd arrnm, intf *lar, _fcd grpnm, intf
		*lgr, intf *dim, _fcd dname, intf *ldn, intf *ret);

   int32 MVSfind(int32 fid, char const *vsname);

   int addMODISgroup(MODFILE *file, char const *groupname, 
                     char const *classname, int32 tag, int32 ref);		   

   AGGREGATE MPVL2ODL(MODFILE *file, char const *PVLAttrName, 
                      char const *aggName); 

   void ncpmdnam(intf *modfil, _fcd arrnm, intf *lar, _fcd grpnm, intf
        *lgr, intf *dim, _fcd dname, intf *ldn, intf *ret);
 
   DATAID *getMODISarrayid(MODFILE *file, char const *arrayname, 
                           char const *groupname);

   DATAID *getMODIStableid(MODFILE *file, char const *tablename,
                           char const *groupname, char const *access);

   int addid(MODFILE *file, char const *name, char const *group, int32 id, 
             int32 type, char const *access);

   DATAID *searchid(MODFILE *file, char const *name, char const *group,
		    int32 type, char const *access);


   int32 searchMODISgroup(MODFILE *file, char const *groupname, 
                          char const *classname, char const *objectname, 
                          char const *objectclass, int32 objecttype);

   void ncemobj(intf *modfil, _fcd name, intf *lna, _fcd group, intf *lgr,
		intf *type, intf *ret);

   intn MHDc2fstr(char *str, intn len);

   void ncgmhoid(intf *modfil, _fcd name, intf *lna, _fcd group, intf *lgr,
		intf *type, _fcd access, intf *lac, intf *ret);

   void ncrmfh(intf  modfil[MODFILLEN], intf *ret);
   
   void nccmfh(intf *SWfid, intf modfil[MODFILLEN], intf *ret);
    
#endif
