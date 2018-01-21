#ifndef SDST_IGEO_ESTIMATE_SCAN_H
#define SDST_IGEO_ESTIMATE_SCAN_H
/*
!C-INC**************************************************************************
!Description:
	Header file contains definition of graninfo_struct and prototypes
        of SDST_IGEO_calibrate_gran_par and SDST_IGEO_estimate_scan_number.


!Revision History:
$Log: SDST_IGEO_estimate_scan.h,v $
Revision 1.4  2000/02/28 17:04:02  lma
modified after walkthrough

 * Revision 1.3  1999/12/16  15:29:09  lma
 * updated prototype definition
 *
 * Revision 1.2  1999/12/07  20:39:06  lma
 * added definition of MAX_SCANS
 * changed mirror_axis to 3 dimension array
 * updated prolog
 *
 * Revision 1.1  1999/09/30  20:20:40  kuyper
 * Initial revision
 *
kuyper@ltpmail.gsfc.nasa.gov

Requirements:

!Team-unique Header:
	This software is developed by the MODIS Science Data Support
	Team for the National Aeronautics and Space Administration,
	Goddard Space Flight Center, under contract NAS5-32373.

	HDF portions developed at the National Center for Supercomputing
	Applications at the University of Illinois at Urbana-Champaign.

!END****************************************************************************
*/

#include "SDST_IGEO_get_scan.h"

#define MAX_SCANS 220

typedef struct {
	char	*filename;
	PGSt_double	T_ecr2gran[3][3];
	PGSt_double	mirror_axis[3];
	double		mirror_trans2, mirror_perp, scans_per_radian;
	int32	scans;
	int	first_scan;
	int8	quality[MAX_SCANS];
} graninfo_struct;


PGSt_SMF_status SDST_IGEO_calibrate_gran_par
(
	scan_info_struct * const [2],
	graninfo_struct * const
);


PGSt_SMF_status SDST_IGEO_estimate_scan_number
(
	PGSt_double const [3],
	MODFILE * const,
	graninfo_struct * const,
	double * const
);

#endif
