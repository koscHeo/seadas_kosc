#ifndef SDST_IGEO_LINE_FRAME_H
#define SDST_IGEO_LINE_FRAME_H

/*
 *!C-INC**************************************************************************
 *
 *!Description
	Main header file for Inverse Geolocation modules 
	
 *!Revision History:
 * 	$Log: SDST_IGEO_line_frame.h,v $
 * 	Revision 1.4  1999/11/03 21:19:32  lma
 * 	Corrected #include files and Standardized the prolog.
 *
 * Revision 1.3  1999/09/29  15:37:35  kuyper
 * Corrected #include files.
 *
 * Revision 1.2  1999/09/03  15:03:12  jayshree
 * corrections made after code walkthrough
 * changed ifndef statement
 *
 * Revision 1.1  1999/08/26  14:41:03  jayshree
 * Initial revision
 *
 *
 *
 *!Team-unique Header:	
 *	This software is developed by the MODIS Science Data Support
 *	Team for the National Aeronautics and Space Administration,
 *	Goddard Space Flight Center, under contract NAS5-32373.
 *
 *	HDF portions developed at the National Center for Supercomputing
 *	Applications at the University of Illinois at Urbana-Champaign.
 *
 *!END****************************************************************************
*/

#include "SDST_IGEO_calc_scan.h"

int SDST_IGEO_estimate_line_frame (
			const scan_info_struct *,
			const PGSt_double ecr_position[],
			double *dline,
			float32 *dframe );


PGSt_SMF_status SDST_IGEO_calibrate_scan_par(
                        scan_info_struct *,
                        const frame_info_struct *);


#endif
	
	
	
