#ifndef SDST_IGEO_CALC_SCAN_H
#define SDST_IGEO_CALC_SCAN_H

/*
!C-INC****************************************************************************
!Description
	Information describing the plane containing a given frame.
	
!Revision History:
	$Log: SDST_IGEO_calc_scan.h,v $
	Revision 1.7  1999/12/16 15:25:36  lma
	updated prolog

 * Revision 1.6  1999/09/24  19:57:22  lma
 * added #include of SDST_IGEO_get_scan.h
 *
 * Revision 2.1  1999/09/21  19:00:29  lma
 * Remove typedef of scan_info_struct
 * Remove unneeded #includes of SDST_IGEO_get_scan.h, hdfi.h
 *
 * Revision 1.4  1999/09/03  15:01:48  jayshree
 * corrections made after code walkthrough
 * changed ifndef statement
 *
 * Revision 1.3  1999/08/26  14:40:38  jayshree
 * *** empty log message ***
 *
 * Revision 1.2  1999/08/16  15:45:02  jayshree
 * Initial Revision
 *
 * Revision 1.1  1999/07/21  18:35:10  jayshree
 * Initial revision
        
!Team-unique Header:	
 *	This software is developed by the MODIS Science Data Support
 *	Team for the National Aeronautics and Space Administration,
 *	Goddard Space Flight Center, under contract NAS5-32373.
 *
 *	HDF portions developed at the National Center for Supercomputing
 *	Applications at the University of Illinois at Urbana-Champaign.

!END****************************************************************************
*/



#include "PGS_TYPES.h"
#include "SDST_IGEO_get_scan.h"


typedef struct {
	int  iframe;
	int  first_det;
	int  last_det;
	PGSt_double normal[3];
} frame_info_struct;

PGSt_SMF_status SDST_IGEO_calc_scan_coords(
                scan_info_struct        *scan_info,
                frame_info_struct       frame_info[]
);
#endif
	
	
	
