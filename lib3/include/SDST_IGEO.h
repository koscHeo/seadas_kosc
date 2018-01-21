#ifndef SDST_IGEO_H
#define SDST_IGEO_H
/*
!PROLOG*************************************************************************
!Description:
	Header for using the C interface to the inverse geolocation routine.

!Input Parameters:	N/A

!Output Parameters:	N/A

!Revision History:
$Log: SDST_IGEO.h,v $
Revision 5.1  2009/06/26 17:32:22  kuyper
Improved const-safety.

James Kuyper	James.R.Kuyper@nasa.gov

Revision 1.2  2000/02/02 22:30:48  kuyper
Corrected sets to int16

 * Revision 1.1  1999/09/30  20:18:34  kuyper
 * Initial revision
 *
kuyper@ltpmail.gsfc.nasa.gov

Requirements:
	CCR 468

!Team-unique Header:
	This software is developed by the MODIS Science Data Support
	Team for the National Aeronautics and Space Administration,
	Goddard Space Flight Center, under contract NAS5-32373.

	HDF portions developed at the National Center for Supercomputing
	Applications at the University of Illinois at Urbana-Champaign.

!END****************************************************************************
*/

#include "mapi.h"
#define SDST_IGEO_MAX_MATCHES 3

int SDST_IGEO(
	MODFILE	*,
	int,
	const float32 [],
	const float32 [],
	float32	[][SDST_IGEO_MAX_MATCHES],
	float32	[][SDST_IGEO_MAX_MATCHES],
	int16 []
);

#endif

