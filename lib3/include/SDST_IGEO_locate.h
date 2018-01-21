#ifndef SDST_IGEO_LOCATE_H
#define SDST_IGEO_LOCATE_H
#include "mapi.h"

/*
!PROLOG*************************************************************************
!Description:
	C header file for SDST_IGEO_locate.h.
!Input Parameters:
	None

!Output Parameters:
	None

!Revision History:
$Log: SDST_IGEO_locate.h,v $
Revision 1.1  1999/09/29 15:24:43  kuyper
Initial revision

kuyper@ltpmail.gsfc.nasa.gov

Requirements:
	CCR 468

!Team-unique Header:
	This software is developed by the MODIS Science Data Support
	Team for the National Aeronautics and Space Administration,
	Goddard Space Flight Center, under contract NAS5-32373.

!END****************************************************************************
*/

PGSt_SMF_status SDST_IGEO_locate_line_frame
(
	MODFILE * const,
	int32 const,
	PGSt_double const[3],
	float32	* const,
	float32	* const
);

#endif
