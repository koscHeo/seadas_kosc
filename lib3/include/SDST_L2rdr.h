/************************************************************************
*			SDST_L2rdr.h utilities header file		
*************************************************************************
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
* $Log: SDST_L2rdr.h,v $
* Revision 5.1  2009/06/25 16:40:23  kuyper
* Improved const safety.
*
* Revision 1.4  2001/11/02 20:26:06  pliu
* Added to include mapi.h
*
* Revision 1.3  1999/09/28 00:01:22  solanki
* Cleaned up and working on all platforms.
*
 * Revision 1.2  1999/04/29  16:34:45  solanki
 * Modified filename in prolog section.
 *
 * Revision 1.1  1999/02/05  18:58:30  solanki
 * Initial revision
 *
 * Revision 1.1  1999/02/05  18:58:30  solanki
 * Initial checkin
 *
*               Jayshree Murthy   10/30/1998
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

#ifndef MODISL2RDR_
#define MODISL2RDR_

#include "mapi.h"
int modisL2reader (
	MODFILE *,
	const int[],
	const char *,
	const char *,
	const int32[],
	const int32[],
	int32,
	double *
);
#endif

