/*
    This file contains all subroutines for the seawifs level 3 binning
*/


#include <sys/types.h>
#include <math.h>
#include <stdlib.h>

#include <seaproto.h>

#define		PI	3.141592653589793

/* 
   global variables only used inside this file, the call to 'bin_init'
   may pass these variables to the calling routine

	NUMROWS :  total number of rows for binning 
	*NUMBIN	:  no. of bin in each row 
	*BASEBIN:  1st bin of each row 
	TOTBINS:   total bin numbers 
	*LATBIN:   center latitude of each row 
*/

int32 *NUMBIN, *BASEBIN;
float *LATBIN;
int32 TOTBINS;
int32 NUMROWS;
float SEAM_LON;

/*
    bin_init: given the total row number, this subroutine returns 
              above variables.

*/

void bin_init(int32 nrow, int32 **nbin, int32 **bbin, float **lbin, int32 *tbin)
{
  int i;
  double radfac;

  if (nrow == -1)   /* free the memory if not the first call */
  {
    free(NUMBIN);
    free(BASEBIN);
    free(LATBIN);
    return;
  }

  SEAM_LON = -180.0;      /*  this value should be passed in  */
  NUMROWS = nrow;

  NUMBIN = (int32 *) calloc(NUMROWS, sizeof(int32));
  BASEBIN = (int32 *) calloc(NUMROWS, sizeof(int32));
  LATBIN = (float *) calloc(NUMROWS, sizeof(float));

  radfac = PI / 180.0;

  for (i=0; i<NUMROWS; i++)
  {
    *(LATBIN+i) = (i + 0.5) * (180.0 / NUMROWS) - 90.0;
    *(NUMBIN+i) = (int32) (cos(*(LATBIN+i)*radfac) * (2.0 * NUMROWS) + 0.5);
  }

  *BASEBIN = 1;

  for (i=1; i<NUMROWS; i++)
    *(BASEBIN+i) = *(BASEBIN+i-1) + *(NUMBIN+i-1);

  TOTBINS = *(BASEBIN+NUMROWS-1) + *(NUMBIN+NUMROWS-1) - 1;

  *nbin = NUMBIN;      /* pass pointer back to the calling routine */
  *bbin = BASEBIN;
  *lbin = LATBIN;
  *tbin = TOTBINS;
} 




/* 
   given a bin number, return the center lat/lon of that bin number 
   heuristic and binary search algorithm is used
*/
void bin2ll(int32 bin, float *lat, float *lon)
{
  int32   row, rlow, rhi, rmid;
  static int32 old_row=0;    /* 1-relative */

  if (old_row > 0 && *(BASEBIN+old_row-1) <= bin && *(BASEBIN+old_row) > bin)
  {
    row = old_row;
  }
  else
  {
    if (bin < 1)
      bin = 1;             /* south pole */
    if (bin > TOTBINS)
      bin = TOTBINS;       /* north pole */
 
    /* binary search for row in range [1..NUMROWS] */
    rlow = 1;            /* 1-relative */
    rhi = NUMROWS;       /* 1-relative */
    while (1)
    {
      rmid = (rlow + rhi - 1) / 2;     /* 0-relative */
      if (*(BASEBIN+rmid) > bin)
        rhi = rmid;      
      else
        rlow = rmid + 1;

      if (rlow == rhi)
      {
        row = rlow;
        break;
      }
    }
    old_row = row;
  }

  *lat = *(LATBIN+row-1);
  *lon = 360.0 * (bin - *(BASEBIN+row-1) + 0.5) / *(NUMBIN+row-1);
  *lon = *lon + SEAM_LON;  /* note, *lon returned here may be in 0 to 360 */
}



/* 
   given the lat/lon, return the bin number 
   lon has to be in the range of -180.0 to 180.0
*/
void ll2bin(float lat, float lon, int32 *bin)
{
  int32 row, col;       /* 0-relative */

  row = (int32) ((90.0 + lat) * (float) NUMROWS / 180.0);  
  col = (int32) ((float) (*(NUMBIN+row)) * (lon - SEAM_LON) / 360.0);
  *bin = *(BASEBIN+row) + col;
}




/* 
   given the lat/lon, return the row, column 
   lon has to be in the range of -180.0 to 180.0

   6/96. Cast double in order to get same bin #'s for SeaWiFS smap9
         and SeaDAS smigen.
*/
void ll2rc(float lat, float lon, int32 *row, int32 *col)
{
  *row = (int32) ((90.0 + (double) lat) * (double) NUMROWS / 180.0);  
  *col = (int32) ((double) (*(NUMBIN + (*row))) * ((double) lon - SEAM_LON) / 360.0);
  *row = *row + 1;
  *col = *col + 1; 
}




/*
   given row/column, return lat/lon
*/
void rc2ll(int32 row, int32 col, float *lat, float *lon)
{
  *lat = *(LATBIN+row-1);
  *lon = SEAM_LON + (360.0 * (col - 0.5) / *(NUMBIN+row-1));
}




/* 
   given a row/column number, return the bin number (1-relative) 
*/
void rc2bin(int32 row, int32 col, int32 *bin)
{
  *bin = *(BASEBIN+row-1) + col - 1;
}



/* 
   given a bin number, return the row and column (both are 1-relative) 
   heuristic and binary search algorithm is used
*/
void bin2rc(int32 bin, int32 *row, int32 *col)
{
  int32   rlow, rhi, rmid;
  static int32 old_row=0;    /* 1-relative */

  if (old_row > 0 && *(BASEBIN+old_row-1) <= bin && *(BASEBIN+old_row) > bin)
  {
    *row = old_row;
  }
  else
  {
    if (bin < 1)
      bin = 1;             /* south pole */
    if (bin > TOTBINS)
      bin = TOTBINS;       /* north pole */
 
    /* binary search for row in range [1..NUMROWS] */
    rlow = 1;            /* 1-relative */
    rhi = NUMROWS;       /* 1-relative */
    while (1)
    {
      rmid = (rlow + rhi - 1) / 2;     /* 0-relative */
      if (*(BASEBIN+rmid) > bin)
        rhi = rmid;      
      else
        rlow = rmid + 1;

      if (rlow == rhi)
      {
        *row = rlow;
        break;
      }
    }
    old_row = *row;
  }

  *col = bin - *(BASEBIN + (*row) - 1) + 1;
}



/* 
   given a bin number, return the center lat/lon of that bin number 
   no heuristic or binary search algorithm is used
   this routine is very slow due to the array reference (I think)
*/

/*
void old_bin2ll(bin, lat, lon)
int32 bin;
float *lat, *lon;
{
  int32 row;

  row = NUMROWS-1;

  while (bin < BASEBIN[row])
    row--;

  *lat = LATBIN[row];
  *lon = 360.0 * (bin - BASEBIN[row] + 0.5) / NUMBIN[row];
}

*/


/*  update version, much faster than using array reference */

void old_bin2ll(int32 bin, float *lat, float *lon)
{
  int32 row;
  int32 *tmpptr;

  row = NUMROWS - 1;
  tmpptr = BASEBIN + NUMROWS - 1;

  while (bin < *tmpptr--)
    row--;

  *lat = LATBIN[row];
  *lon = 360.0 * (bin - BASEBIN[row] + 0.5) / NUMBIN[row];
}

