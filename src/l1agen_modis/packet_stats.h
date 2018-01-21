#ifndef PACKET_STATS
#define PACKET_STATS
/*!C-INC************************************************************************

!Description:
	Header file for packet filtering statistics.

!Input Parameters:
	N/A

!Output Parameters:
	N/A

Return Values:
	N/A

Externally Defined:
	stats		"print_stats.c"

Called by:
	N/A

Routines Called:
	N/A

!Revision History:
$Log: packet_stats.h,v $
Revision 6.1  2010/08/25 19:19:17  kuyper
Initial revision.

James Kuyper Jr.	James.R.Kuyper@NASA.gov

!Team-unique Header:
	This software is developed by the MODIS Science Data Support
	Team for the National Aeronautics and Space Administration,
	Goddard Space Flight Center, under contract NAS5-32373.

References and Credits

Design Notes
!END**************************************************************************
*/
#include "hdfi.h"
extern struct packet_stats
{
    uint32 packets;
    /* Single field validity checks */
    uint32 version;		/* Must be 0.	*/
    uint32 type;		/* Must be 0 (normal).	*/
    uint32 seq_flag;		/* Must be non-zero.	*/
    uint32 sec_hdr_flag;	/* Must be 1 (secondary header). */
    uint32 apid;		/* 64 or 127(?) (MODIS).	*/
    uint32 time_tag;
    uint32 quick_look;		/* Must be 0 (not selected)	*/
    uint32 pkt_type;		/* 0, 1, 2, 4 */

    /* Single packet consistency checks */
    uint32 seq_type;		/* == 3 only iff pkt_type == night	*/
    uint32 length_type;		/* short length iff pkt_type == night	*/
    uint32 type_flag;		/* type_flag non-zero only if pkt_type == Day */
    uint32 frame_count;		/* Consistent with type_flag, calib_type. */
    uint32 cksum;		/* Matches checksum of data field.	*/
} stats;

void print_stats(void);
#endif
