#include <stdio.h>
#include "PGS_SMF.h"
#include "packet_stats.h"
struct packet_stats stats;

void print_stats(void)
/*!C****************************************************************************

!Description:   
	Writes the packet filtering statistics to the Report Log.

!Input Parameters:
	N/A

!Output Parameters:
	N/A

Return Values:
	None

Externally Defined:
	N/A

Called by:
	process_a_granule()		"L1A_prototypes.h"

Routines Called:
	PGS_SMF_GenerateStatusReport()	"PGS_SMF.h"

!Revision History:
James Kuyper Jr.	James.R.Kuyper@NASA.gov

!Team-unique Header:
	This software is developed by the MODIS Science Data Support
	Team for the National Aeronautics and Space Administration,
	Goddard Space Flight Center, under contract NAS5-32373.

References and Credits

Design Notes
!END**************************************************************************
*/
{
    char buffer[1024];

    sprintf(buffer, "Packet filtering statistics.\n"
     "\tpackets:		%9u\n"	     "\tNon-zero version:	%9u\n"
     "\tTest:			%9u\n"	     "\tNon-zero sec_hdr_flag:	%9u\n"
     "\tBad APID:		%9u\n"	     "\tInvalid time_tag:	%9u\n"
     "\tQuick Look:		%9u\n"	     "\tBad Packet type:	%9u\n"
     "\tSeq_flag/type mismatch: %9u\n"	     "\tLength/type mismatch:	%9u\n"
     "\tinvalid frame_count:	%9u\n",
            stats.packets, stats.version, 
            stats.type, stats.sec_hdr_flag, 
            stats.apid, stats.time_tag, 
            stats.quick_look, stats.pkt_type, 
            stats.seq_type, stats.length_type, 
            stats.frame_count);

    PGS_SMF_GenerateStatusReport(buffer);
    /* The only way this function can fail is if it can't open the log file.
     * Therefore, there's no way to report the failure of this function.
     */
}

