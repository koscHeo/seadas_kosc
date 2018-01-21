#ifndef SDST_IGEO_GET_SCAN_H
#define SDST_IGEO_GET_SCAN_H

/*
 *!C-INC************************************************************************
 *
 *!Description
 *	Header file contains definition of scan_info_struct and prototype
 *      of SDST_IGEO_get_scan.
 *
 *!Revision History:
 *	$Log: SDST_IGEO_get_scan.h,v $
 *	Revision 2.5  2000/04/17 20:17:36  kuyper
 *	Moved DEG2RAD in from source code.
 *
 * Revision 2.4  1999/12/07  20:42:05  lma
 * modified description
 *
 * Revision 2.3  1999/11/03  21:20:52  lma
 * modified EV_frames data type.
 *
 * Revision 2.2  1999/09/17  14:31:49  kuyper
 * Corrected scan_info_struct to include iscan.
 * Standardized the prolog.
 *
 * Revision 2.1  1999/09/16  22:58:12  kuyper
 * Corrected to typedef scan_info_struct.
 * Added prototype for SDST_IGEO_get_scan().
 *
 * Revision 1.3  1999/09/03  15:02:42  jayshree
 * corrections made after code walkthrough
 * changed ifndef statement
 *
 * Revision 1.2  1999/08/16  15:45:17  jayshree
 * Initial Revision
 *
 * Revision 1.1  1999/07/21  18:35:31  jayshree
 * Initial revision
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


#include "mapi.h"

#define MAX_FRAMES      1354
#define MAX_DETS         10
#define INVALID_SENSOR_RANGE 8
#define DEG2RAD		(3.14159265358979323846/180.0)

typedef struct
{
	PGSt_double	position[MAX_DETS][MAX_FRAMES][3];
	PGSt_double	T_ecr2scan[3][3];
	PGSt_double	sc_position[3];
	double		center_frame;
	double		center_line;
	double		frames_per_radian;
	double		lines_per_radian;
	int32		iscan;
	int16		EV_frames;
	uint8		gflags[MAX_DETS][MAX_FRAMES];
}	scan_info_struct;


PGSt_SMF_status SDST_IGEO_get_scan
(
	MODFILE			*geo_file,
	int32			iscan,
	scan_info_struct	**scan_info
);
#endif
