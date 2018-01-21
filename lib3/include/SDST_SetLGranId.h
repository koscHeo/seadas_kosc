/************************************************************************
*			SDST_SetLGranId.h                               *
*                       utilities header file		            	*
*                                                                       *    	
*************************************************************************

***********************************************************************
* !C-INC
*
* !Purpose:	This utilities header file contains the prototype
*		and constants for the local granule ID routine
*		SDST_SetLocalGranId.
*
* !Description: The Header file SDST_SetLGranId.h is part of a larger 
*		software system called the SDST utility toolkit 
*		abbreviated SDST_TK.
*               SDST_TK consists of in-house developed utility routines
*               which the MODIS Science Team can link to MODIS production
*		software.
*
* !Input parameters:
*     none.
*
* !Output parameters:
*     none.
*
*
* !Revision history:
* $Log: SDST_SetLGranId.h,v $
* Revision 2.4  2001/11/15 16:11:56  pliu
* Changed VERSIONID_DEFLT definition to integer from string.
*
* Revision 2.3  2001/07/12 14:57:22  pliu
* Defined default values used by SDST_SetLocalGranId.c here.
*
 * Revision 2.2  1999/09/28  00:01:22  solanki
 * Cleaned up and working on all platforms.
 *
 * Revision 2.1  1999/05/11  17:38:31  jjb
 * Removed extraneous constants.
 * Changed return type of SDST_SetLocalGranID to PGSt_SMF_status.
 *
* Revision 1.1  1998/09/29  16:24:49  solanki
* Initial revision
*
*               Frederick J. Shaw 1996/04/02
*
* !Team-unique header:
*
* !References and Credits
*      This software is developed by the MODIS Science Data Support
*      Team for the National Aeronautics and Space Administration,
*      Goddard Space Flight Center, under contract NAS5-32373.
*
*      HDF portions developed at the National Center for Supercomputing
*      Applications at the University of Illinois at Urbana-Champaign.
*
* !Design Notes
*
* !END
********************************************************************
*/

#ifndef SDST_SetLGranId_
#define SDST_SetLGranId_

#include "PGS_MET.h"

/* Define the following default values that are used by 
   SDST_SetLocalGranId.c
*/
#define SAT_MODE 'A'
#define VERSIONID_DEFLT 0
#define TILE_DEFLT "000"
#define MAX_BUFF_SIZE 256
#define PROC_TIME_SIZE 80
#define STARTDATADAY "StartDataDay"
#define DD_DEFLT ""

/*Define constants for L1A, L1B, L2, L2G, L3, GEO for prod_type in an 
external file for including.  
*/
enum PROD_TYPE {L1A, L1B, L2, L2G, GEO, OCEANS, ATMOS, LAND};

PGSt_SMF_status SDST_SetLocalGranId(const enum PROD_TYPE prod_type, 
				    PGSt_MET_all_handles mdhandles); 
  
#endif
