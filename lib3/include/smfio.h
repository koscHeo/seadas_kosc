/* file: smfio.h */

#ifndef smfio_h
#define smfio_h
/*
!C-INC****************************************************************************

!Description:

  Header file needed to use modsmf().

!Revision History:
 
$Log: smfio.h,v $
Revision 5.1  2009/06/25 21:38:49  kuyper
Removed (A LOT OF) obsolete stuff unrelated to the SDST Library.
Improved const safety.

James Kuyper	James.R.Kuyper@nasa.gov

Revision 1.3  1999/09/28 00:01:22  solanki
Cleaned up and working on all platforms.

 * Revision 1.2  1997/11/12  21:23:38  solanki
 * Modified for fortran wrapper to localgranule id.
 *
 * Revision 1.2.1.1  1997/07/18  22:40:58  kuyper
 * Merged in out-of-sequence changes.
 *
 * Revision 1.2  1997/07/18  21:58:00  kuyper
 * Baselined Version 1
 *
 *Parallel development:
 * Revision 1.13  1997/04/18  17:26:13  fhliang
 * added file name (L.1);
 * fixed prolog.
 *
 * Revision 1.12  1997/03/26  19:12:49  fhliang
 * Initial revision of SDST delivery of smfio.h.
 *
Revision 1.11  1996/12/06 19:07:23  fshaw
Initial revision

Revision 1.10  1996/08/01 21:09:25  fisher
updated to conform to modis standards

 * Revision 1.9  1996/06/15  06:41:29  fisher
 * ,
 *
 * Revision 1.8  1996/06/07  18:27:57  fisher
 * *** empty log message ***
 *
 * Revision 1.7  1996/05/13  19:44:53  fisher
 * updated for concat_mpfn()
 *
 * Revision 1.6  1996/04/08  18:05:13  fisher
 * prologue revision
 *
 * Revision 1.5  1996/04/08  17:55:52  fisher
 * Added prototype for rmpath
 *
 * Revision 1.4  1996/04/03  00:14:08  fisher
 * Added prototype for current_time_a
 *
 * Revision 1.3  1996/03/08  20:21:50  fisher
 * Added memory definition for current_time_a (TIMECODEASIZE)
 *
 * Revision 1.2  1995/07/14  17:57:12  fisher
 * Prologue update
 *
 * Revision 1.1  1995/07/12  18:58:26  fisher
 * Initial revision
 *
$Id: smfio.h,v 5.1 2009/06/25 21:38:49 kuyper Exp $

!Team-unique Header:
  This software has been created by the MODIS Science Data Support
  Team for the National Aeronautics and Space Administration,
  Goddard Space Flight Center, under contract NAS5-32373.

!References and Credits:

    Written by Paul S. Fisher
    Research and Data Systems Corporation
    SAIC/GSC MODIS Support Office
    7501 Forbes Blvd
    Seabrook MD 20706  

!Design Notes:

    This header file should be included in all relevant, integrated
    pieces of source code.

!END*****************************************************************************/


/* SDP Toolkit Status Message Facility header file */
#include "PGS_SMF.h"

/* Prototypes */
void modsmf(
	PGSt_SMF_code	mnemonicstring,
	const char	*infostring,
	const char	*functionstring
);

#define TIMECODEASIZE 28	/* Space needed to hold UTC-A time strings. */
#define TIMECODEBSIZE 27	/* Space needed to hold UTC-B time strings. */
#endif

