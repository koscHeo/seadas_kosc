#ifndef SDST_IGEO_GET_GRAN_H
#define SDST_IGEO_GET_GRAN_H
/*
 *!C-INC************************************************************************
 *
 *!Description
 * !Revision History:
 *	$Log: SDST_IGEO_get_gran.h,v $
 *	Revision 1.2  2000/03/16 22:31:45  kuyper
 *	Added double-inclusion guards. Corrected #include.
 *
 * Revision 1.1  1999/12/16  15:31:11  lma
 * Initial revision
 *
 *       
 *!Team-unique Header:
 *	This software is developed by the MODIS Science Data Support
 *	Team for the National Aeronautics and Space Administration,
 *	Goddard Space Flight Center, under contract NAS5-32373.
 *
 *	HDF portions developed at the National Center for Supercomputing
 *	Applications at the University of Illinois at Urbana-Champaign.
 *!END**************************************************************************
 */

#include "SDST_IGEO_estimate_scan.h"


PGSt_SMF_status SDST_IGEO_get_graninfo
(
	MODFILE			*geo_file,
	graninfo_struct         *graninfo
);
#endif
