
/***********************************************************************
* !C-INC
*
* !Purpose:	This utilities header file contains SDST Toolkit prototypes
*		and constants for MODIS science software.
*
* !Description: The Header file SDST_TK is part of a larger software
*               system called the SDST utility toolkit abbreviated SDST_TK.
*               SDST_TK consists of in-house developed utility routines
*               which the MODIS Science Team can link to.
*
* !Input parameters:
*     none.
*
* !Output parameters:
*     none.
*
* !Revision history:
* $Log: SDST_GetSetPSA.h,v $
* Revision 5.1  2009/06/25 16:19:13  kuyper
* Improved const-safety.
*
* Revision 1.1  1999/05/07 20:05:28  jayshree
* Initial revision
*
 * Revision 1.1  1999/04/19  18:34:56  jayshree
 * Initial revision
 *
*               Jayshree Murthy   04/19/199
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
****************************************************************************
*/

#ifndef SDST_GetSetPSA_
#define SDST_GetSetPSA_

#include "PGS_MET.h"
#include "PGS_MODIS_39604.h"


PGSt_SMF_status SDST_GetSetPSA ( PGSt_MET_handle mdHandle, 
                                 const char* PSANameStr, void* PSAValue);
#endif                        
