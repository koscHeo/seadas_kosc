/*******************************************************************

   l1io.h

   purpose: include file for the use of the level 1 direct I/O routines

   Parameters: 
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      unsigned char *   l1a_path        I       path for level 1A file

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       15-Feb-1993     Original development

*******************************************************************/

/*
 *  Note that hdf.h is needed for this and navigation defs
 */
#include "hdf.h"
#include <mfhdf.h>
#include "nav_l1io.h"
/*
 *  the l1info_struct structure has file ids...
 */
typedef struct l1info_struct_d 
   {
   int32 fid;
   int32 sdfid; 
   }  l1info_struct;

/*
 *  prototypes
 */
 int32 l1io_open( char *, l1info_struct*, int32 *, int32 * );
 int32 l1io_read( l1info_struct, int, int16 *, 
                  navblockType * );
 void l1io_close( l1info_struct );
