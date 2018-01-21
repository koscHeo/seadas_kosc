/*

gcc -I../../inc/swfinc -g -I/home/joel/hdf/include -c -o MODll2snpx.o MODll2snpx.c

gcc MODll2snpx.o /home/joel/hdf/lib/libmfhdf.a /home/joel/hdf/lib/libdf.a /home/joel/hdf/lib/libz.a -ljpeg -lz -lm -o MODll2snpx

*/

/* ========================================================================
 * MODll2snpx - Determines bounding scans and pixels for given lon/lat box
 * 
 * Synopsis:
 *
 *   MSll2snpx geo-filename SW-lon SW-lat NE-lon NE-lat
 *
 * Description:
 * 
 * Modification history:
 *
 *     Programmer     Organization      Date       Description of change
 *   --------------   ------------    --------     ---------------------
 *   Joel M. Gales    Futuretech      22 Apr 2004  Original Development
 *
 * ======================================================================== */

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "mfhdf.h"

#define INT32   int32_t 
#define FLOAT32 float
#define BYTE    unsigned char
#define CMD_ARGS    "f:"   /* Valid commandline options */

#define RADEG 57.2958

void usage (char *prog, char *ver)
{
    int i;

    printf("\nThis is version %s of %s.\n",ver,prog);
    printf("Usage (box): %s geo_filename SWlon SWlat NElon NElat\n",prog);
    printf("Usage (pix): %s geo_filename lon lat\n",prog);
    /*
    printf("\nOptional arguments:\n");
    printf("\t[-f first-scan-number for search] (default=0)\n");
    */
    exit(-1);
}
  

/* -------------------------------------------------------------------- *
 *                              main                                    *
 * -------------------------------------------------------------------- */
int main(int argc, char *argv[])
{

    int32_t   iscan    = 0;          /* input scan number                  */
    int32_t   spix     = 0;          /* start pixel for subscene process   */
    int32_t   epix     = -1;         /* end pixel for subscene process     */
    int32_t   sscan    = 0;          /* start scan for subscene process    */
    int32_t   escan    = -1;         /* end scan for subscene process      */
    int32_t   ip;                    /* input pixel number                 */

    int32_t   npix     = 0;          /* Number of output pixels per scan   */
    int32_t   nscan    = 0;          /* Number of output scans             */

    int32_t   status; 

    int32_t    ipix;
    int32_t    jpix;
    int32_t    idet;
    float32 cornerlon[2];
    float32 cornerlat[2];
    float32 pixll[2];
    float32 lon;
    float32 lat;
    float64 uv[3];
    float64 uvpix[3];
    float64 maxcos = -1;
    float64 dot;
    float32 pixlon;
    float32 pixlat;

    int16   bndflag;
    int16   lonTest;
    int16   latTest;
    int16   first = 1;
    int16   found = 0;
    int16   pixsrch;

    int32   minScan = 2147483647;
    int32   maxScan = -1;
    int32   minPix  = 2147483647;
    int32   maxPix  = -1;
    int32   pixsn = -1;
    int32   pixpx = -1;

    int32 sd_id;
    int32 indx=-1;
    int32 sds_id_ll[2];
    int32 start[2], edges[2], dims[2];
    int32 rank, dtype, nattrs;

    float32 latlon[2][2000];

    char buffer[1000];

    /* Parameters for getopt() */
    extern int      opterr;
    extern int      optind;
    extern char    *optarg;
    int             c;  
    char            geofile[FILENAME_MAX] = "";


    /*									*/
    /* Process required command-line arguments                          */
    /*									*/
    switch (argc-optind+1) {
      case 4: 
        strcpy(geofile,argv[optind+0]);
	pixll[0]= atof(argv[optind+1]);
	pixll[1]= atof(argv[optind+2]);
	pixsrch = 1;

        uvpix[0] = cos(pixll[1]/RADEG) * cos(pixll[0]/RADEG);
        uvpix[1] = cos(pixll[1]/RADEG) * sin(pixll[0]/RADEG);
        uvpix[2] = sin(pixll[1]/RADEG);

        break;
      case 6: 
        strcpy(geofile,argv[optind+0]);
	cornerlon[0]= atof(argv[optind+1]);
	cornerlat[0]= atof(argv[optind+2]);
	cornerlon[1]= atof(argv[optind+3]);
	cornerlat[1]= atof(argv[optind+4]);

	bndflag = (cornerlon[0] < cornerlon[1]) ? 0 : 1;
	pixsrch = 0;
        break;
      default:
        usage(argv[0],"0.93");
        break;
    }

    sd_id = SDstart(geofile, DFACC_RDONLY);
    if (sd_id == -1) {
        printf("-E- %s: Error opening %s for reading.\n",argv[0],geofile);
        exit(1);
    }

    sds_id_ll[0] = SDselect(sd_id, SDnametoindex(sd_id, "Latitude"));
    if (sds_id_ll[0] == -1) {
        printf("-E- %s: Error opening Latitude field.\n",argv[0]);
        exit(1);
    }
    status = SDgetinfo(sds_id_ll[0], buffer, &rank, dims, &dtype, &nattrs);

    escan = dims[0];
    npix = dims[1];

    sds_id_ll[1] = SDselect(sd_id, SDnametoindex(sd_id, "Longitude"));
    if (sds_id_ll[1] == -1) {
        printf("-E- %s: Error opening Longitude field.\n",argv[0]);
        exit(1);
    }

    /* Generate corner lon/lat if single pixel search */
    /* ---------------------------------------------- */
    if (pixsrch) {

      cornerlon[0] = pixll[0] - 1;
      cornerlon[1] = pixll[0] + 1;

      cornerlon[0] = pixll[0] - 0.1;
      cornerlon[1] = pixll[0] + 0.1;

      if (cornerlon[0] < -180) cornerlon[0] = 360 - cornerlon[0];
      if (cornerlon[1] > +180) cornerlon[1] = cornerlon[0] - 360;


      cornerlat[0] = pixll[1] - 1;
      cornerlat[1] = pixll[1] + 1;

      cornerlat[0] = pixll[1] - 0.1;
      cornerlat[1] = pixll[1] + 0.1;

      if (cornerlat[0] < -90) cornerlat[0] = -89.99;
      if (cornerlat[1] > +90) cornerlat[0] = +89.99;
    }


    /*					 			*/
    /* Read file scan by scan                                   */
    /*								*/
    for (iscan=sscan; iscan<escan; iscan++) {

        if (iscan % 500 == 0) printf("#Reading scan: %d\n", iscan);

	start[0] = iscan;
	start[1] = 0;

	edges[0] = 1;
	edges[1] = dims[1];

	status = SDreaddata(sds_id_ll[0], start, NULL, edges, (VOIDP) latlon[0]);
	status = SDreaddata(sds_id_ll[1], start, NULL, edges, (VOIDP) latlon[1]);

        for (ipix = 0; ipix < npix; ipix++) {

	  lat = latlon[0][ipix];
	  lon = latlon[1][ipix];

	  if (lon > 180) lon = lon - 360;

	  latTest = (lat >= cornerlat[0] && lat <= cornerlat[1]);

	  if (bndflag == 1) {
	    lonTest = (lon >= cornerlon[1] && lon <= cornerlon[0]);
	    lonTest = 1 - lonTest;
	  }
	  else {
	    lonTest = (lon >= cornerlon[0] && lon <= cornerlon[1]);
	  }

	  if (lonTest + latTest == 2) {
	    if (ipix < minPix) minPix = ipix;
	    if (iscan < minScan) minScan = iscan;
	    if (ipix > maxPix) maxPix = ipix;
	    if (iscan > maxScan) maxScan = iscan;
	  }

	}

	if (pixsrch && minPix != -1 &&  maxPix != -1) {

	  for (ipix = minPix; ipix <= maxPix; ipix++) {
	    lat = latlon[0][ipix] / RADEG;
	    lon = latlon[1][ipix] / RADEG;
	    uv[0] = cos(lat) * cos(lon);
	    uv[1] = cos(lat) * sin(lon);
	    uv[2] = sin(lat);

	    dot = uv[0] * uvpix[0] + uv[1] * uvpix[1] + uv[2] * uvpix[2];
	    if (dot > maxcos) {
	      maxcos = dot;
	      pixsn = iscan;
	      pixpx = ipix;
	      pixlon = lon * RADEG;
	      pixlat = lat * RADEG;
	    }

	  }

	}
    }


    if (pixsrch) {
      if (pixsn != -1) {
	printf("#Lon=%f\n",pixlon);
	printf("#Lat=%f\n",pixlat);
	printf("sline=%d\n",pixsn+1);
	printf("spixl=%d\n",pixpx+1);
	printf("eline=%d\n",pixsn+1);
	printf("epixl=%d\n",pixpx+1);
      }
    }
    else {
      printf("# %d %d %d %d\n", minPix+1,maxPix+1,minScan+1,maxScan+1); 
      printf("sline=%d\n", minScan+1);
      printf("eline=%d\n", maxScan+1);
      printf("spixl=%d\n", minPix+1);
      printf("epixl=%d\n", maxPix+1);
    }

    SDendaccess(sds_id_ll[0]);
    SDendaccess(sds_id_ll[1]);
    SDend(sd_id);

    if (minScan == -1 || maxScan == -1 || minPix == -1 || maxPix == -1)
      exit(1);
    else
      exit(0);
}

