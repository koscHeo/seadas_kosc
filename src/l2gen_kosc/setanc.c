#include "l12_proto.h"

/* ============================================================================ */
/* no2conc() - retrieve no2 concentration from ancillary file                   */
/*                                                                              */
/* Written By: B. Franz, NASA OBPG, June 2006.                                  */
/*                                                                              */
/* ============================================================================ */

#define NXNO2 1440
#define NYNO2 720
#define NXANC 180
#define NYANC 90

void no2_frac(float lon, float lat, float *no2_frac_200)
  {
  static int firstCall = 1;
  static int nx = NXANC;
  static int ny = NYANC;
  static float dx = 360.0 / NXANC;
  static float dy = 180.0 / NYANC;
  static float map_frac[NYANC + 2][NXANC + 2];

  int i, j;
  float xx, yy;
  float t, u;

  float frac;

  *no2_frac_200 = 0.0;

  if (firstCall)
    {

    int32 sd_id;
    int32 sds_id;
    int32 status;
    int32 sds_index;
    int32 rank;
    int32 nt;
    int32 dims[H4_MAX_VAR_DIMS];
    int32 nattrs;
    int32 start[2];
    int32 edges[2];
    char name[H4_MAX_NC_NAME];
    char sdsname[H4_MAX_NC_NAME];
    float **map;
    char title[255];
    int16 month;
    char *no2_frac_fil, no2_frac_file[300];

    firstCall = 0;

    if ((no2_frac_fil = getenv("OCDATAROOT")) == NULL)
      {
      printf("-E- %s:  Error looking up environmental variable OCDATAROOT\n",
        __FILE__);
      exit(1);
      }
    strcpy(no2_frac_file, no2_frac_fil);
    strcat(no2_frac_file, "/common/trop_f_no2_200m.hdf");

    map = (float **) alloc2d_float(NXANC, NYANC);
    if (map == NULL)
      {
      printf("-E- %s:  Error allocating NO2 frac space for %s.\n", __FILE__,
        no2_frac_file);
      exit(1);
      }

    sd_id = SDstart(no2_frac_file, DFACC_RDONLY);
    if (sd_id == -1)
      {
      printf("-E- %s:  Error opening NO2 frac file %s.\n", __FILE__,
        no2_frac_file);
      exit(1);
      }

    printf("\nOpening NO2 frac file %s\n\n", no2_frac_file);

    strcpy(sdsname, "f_no2_200m");
    sds_index = SDnametoindex(sd_id, sdsname);
    if (sds_index == -1)
      {
      printf("-E- %s:  Error seeking %s SDS from %s.\n", __FILE__, sdsname,
        no2_frac_file);
      exit(1);
      }
    sds_id = SDselect(sd_id, sds_index);

    status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
    if (status != 0)
      {
      printf("-E- %s:  Error reading SDS info for %s from %s.\n", __FILE__,
        sdsname, no2_frac_file);
      exit(1);
      }
    if (dims[0] != ny || dims[1] != nx)
      {
      printf("-E- %s:  Dimension mis-match on %s array from %s.\n", __FILE__,
        sdsname, no2_frac_file);
      printf("  Expecting %d x %d\n", nx, ny);
      printf("  Reading   %d x %d\n", dims[1], dims[0]);
      exit(1);
      }

    start[0] = 0;
    start[1] = 0;
    edges[0] = ny;
    edges[1] = nx;

    status = SDreaddata(sds_id, start, NULL, edges, (VOIDP) & map[0][0]);
    if (status != 0)
      {
      printf("-E- %s:  Error reading SDS %s from %s.\n", __FILE__, sdsname,
        no2_frac_file);
      exit(1);
      }

    for (j = 0; j < ny; j++)
      {
      for (i = 0; i < nx; i++)
        {
        map_frac[j + 1][i + 1] = map[ny - j - 1][i];
        }
      }

    /* add boarders to simplify interpolation */

    for (j = 0; j < ny; j++)
      {
      map_frac[j + 1][0] = map_frac[j + 1][nx];
      map_frac[j + 1][nx + 1] = map_frac[j + 1][1];
      }
    for (i = 0; i < nx + 2; i++)
      {
      map_frac[0][i] = map_frac[1][i];
      map_frac[ny + 1][i] = map_frac[ny][i];
      }

    free2d_float(map);
    }

  /* interpolate to pixel location */

  i = MAX(MIN((int) ((lon + 180.0 + dx / 2) / dx), nx + 1), 0);
  j = MAX(MIN((int) ((lat + 90.0 + dy / 2) / dy), ny + 1), 0);

  xx = i * dx - 180.0 - dx / 2;
  yy = j * dy - 90.0 - dy / 2;

  t = (lon - xx) / dx;
  u = (lat - yy) / dy;

  frac = (1 - t) * (1 - u) * map_frac[j][i] + t * (1 - u) * map_frac[j][i + 1]
    + t * u * map_frac[j + 1][i + 1] + (1 - t) * u * map_frac[j + 1][i];

  /* return components of stratospheric and tropospheric no2  */

  *no2_frac_200 = MAX(frac, 0.0);

  return;
  }

void no2conc(char *no2file, float lon, float lat, int32_t doy, float *no2_tropo,
  float *no2_strat)
  {
  static int firstCall = 1;
  static int nx = NXNO2;
  static int ny = NYNO2;
  static float dx = 360.0 / NXNO2;
  static float dy = 180.0 / NYNO2;
  static float map_total[NYNO2 + 2][NXNO2 + 2];
  static float map_tropo[NYNO2 + 2][NXNO2 + 2];

  int i, j;
  float xx, yy;
  float t, u;

  float total;
  float tropo;

  *no2_tropo = 0.0;
  *no2_strat = 0.0;

  if (firstCall)
    {

    int32 sd_id;
    int32 sds_id;
    int32 status;
    int32 sds_index;
    int32 rank;
    int32 nt;
    int32 dims[H4_MAX_VAR_DIMS];
    int32 nattrs;
    int32 start[2];
    int32 edges[2];
    char name[H4_MAX_NC_NAME];
    char sdsname[H4_MAX_NC_NAME];
    float **map;
    char title[255];
    int16 month;
    char mstr[4] = "";

    firstCall = 0;

    map = (float **) alloc2d_float(NXNO2, NYNO2);
    if (map == NULL)
      {
      printf("-E- %s:  Error allocating space for %s.\n", __FILE__, no2file);
      exit(1);
      }

    sd_id = SDstart(no2file, DFACC_RDONLY);
    if (sd_id == -1)
      {
      printf("-E- %s:  Error opening NO2 file %s.\n", __FILE__, no2file);
      exit(1);
      }

    printf("\nOpening NO2 file %s\n\n", no2file);

    if (SDreadattr(sd_id, SDfindattr(sd_id, "Title"), (VOIDP) title) == 0)
      {
      if (strstr(title, "NO2 Climatology") != NULL)
        {
        month = (int) doy / 31;
        sprintf(mstr, "_%2.2i", month + 1);
        }
      }

    strcpy(sdsname, "tot_no2");
    strcat(sdsname, mstr);
    sds_index = SDnametoindex(sd_id, sdsname);
    if (sds_index == -1)
      {
      printf("-E- %s:  Error seeking %s SDS from %s.\n", __FILE__, sdsname,
        no2file);
      exit(1);
      }
    sds_id = SDselect(sd_id, sds_index);

    status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
    if (status != 0)
      {
      printf("-E- %s:  Error reading SDS info for %s from %s.\n", __FILE__,
        sdsname, no2file);
      exit(1);
      }
    if (dims[0] != ny || dims[1] != nx)
      {
      printf("-E- %s:  Dimension mis-match on %s array from %s.\n", __FILE__,
        sdsname, no2file);
      printf("  Expecting %d x %d\n", nx, ny);
      printf("  Reading   %d x %d\n", dims[1], dims[0]);
      exit(1);
      }

    start[0] = 0;
    start[1] = 0;
    edges[0] = ny;
    edges[1] = nx;

    status = SDreaddata(sds_id, start, NULL, edges, (VOIDP) & map[0][0]);
    if (status != 0)
      {
      printf("-E- %s:  Error reading SDS %s from %s.\n", __FILE__, sdsname,
        no2file);
      exit(1);
      }

    for (j = 0; j < ny; j++)
      {
      for (i = 0; i < nx; i++)
        {
        map_total[j + 1][i + 1] = map[ny - j - 1][i];
        }
      }

    strcpy(sdsname, "trop_no2");
    strcat(sdsname, mstr);
    sds_index = SDnametoindex(sd_id, sdsname);
    if (sds_index == -1)
      {
      printf("-E- %s:  Error seeking %s SDS from %s.\n", __FILE__, sdsname,
        no2file);
      exit(1);
      }
    sds_id = SDselect(sd_id, sds_index);

    status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
    if (status != 0)
      {
      printf("-E- %s:  Error reading SDS info for %s from %s.\n", __FILE__,
        sdsname, no2file);
      exit(1);
      }
    if (dims[0] != ny || dims[1] != nx)
      {
      printf("-E- %s:  Dimension mis-match on %s array from %s.\n", __FILE__,
        sdsname, no2file);
      printf("  Expecting %d x %d\n", nx, ny);
      printf("  Reading   %d x %d\n", dims[1], dims[0]);
      exit(1);
      }

    start[0] = 0;
    start[1] = 0;
    edges[0] = ny;
    edges[1] = nx;

    status = SDreaddata(sds_id, start, NULL, edges, (VOIDP) & map[0][0]);
    if (status != 0)
      {
      printf("-E- %s:  Error reading SDS %s from %s.\n", __FILE__, sdsname,
        no2file);
      exit(1);
      }

    for (j = 0; j < ny; j++)
      {
      for (i = 0; i < nx; i++)
        {
        map_tropo[j + 1][i + 1] = map[ny - j - 1][i];
        }
      }

    /* add boarders to simplify interpolation */

    for (j = 0; j < ny; j++)
      {
      map_total[j + 1][0] = map_total[j + 1][nx];
      map_total[j + 1][nx + 1] = map_total[j + 1][1];
      map_tropo[j + 1][0] = map_tropo[j + 1][nx];
      map_tropo[j + 1][nx + 1] = map_tropo[j + 1][1];
      }
    for (i = 0; i < nx + 2; i++)
      {
      map_total[0][i] = map_total[1][i];
      map_total[ny + 1][i] = map_total[ny][i];
      map_tropo[0][i] = map_tropo[1][i];
      map_tropo[ny + 1][i] = map_tropo[ny][i];
      }

    free2d_float(map);
    }

  /* interpolate to pixel location */

  i = MAX(MIN((int) ((lon + 180.0 + dx / 2) / dx), nx + 1), 0);
  j = MAX(MIN((int) ((lat + 90.0 + dy / 2) / dy), ny + 1), 0);

  xx = i * dx - 180.0 - dx / 2;
  yy = j * dy - 90.0 - dy / 2;

  t = (lon - xx) / dx;
  u = (lat - yy) / dy;

  total = (1 - t) * (1 - u) * map_total[j][i]
    + t * (1 - u) * map_total[j][i + 1] + t * u * map_total[j + 1][i + 1]
    + (1 - t) * u * map_total[j + 1][i];

  tropo = (1 - t) * (1 - u) * map_tropo[j][i]
    + t * (1 - u) * map_tropo[j][i + 1] + t * u * map_tropo[j + 1][i + 1]
    + (1 - t) * u * map_tropo[j + 1][i];

  /* return components of stratospheric and tropospheric no2  */

  *no2_strat = MAX(total - tropo, 0.0);
  *no2_tropo = MAX(tropo, 0.0);

  return;
  }

/* ============================================================================ */
/* ozone_climatology() - retrieve ozone concentration from daily climatology    */
/*                                                                              */
/* Written By: B. Franz, NASA OBPG, April 2009.                                 */
/* W. Robinson, SAIC, 14 Feb 2014  generalize code for varying grid size        */
/*                                                                              */
/* ============================================================================ */
float ozone_climatology(char *file, int day, float lon, float lat)
  {
  static int firstCall = 1;
  static int nx, ny, mapx;
  static float dx, dy, *map;

  int i, j;
  float xx, yy, *tmp;
  float t, u, val_11, val_12, val_21, val_22;
  float ozone;

  if (firstCall)
    {

    int32 sd_id;
    int32 sds_id;
    int32 status;
    int32 sds_index;
    int32 rank;
    int32 nt;
    int32 dims[H4_MAX_VAR_DIMS];
    int32 nattrs;
    int32 start[2];
    int32 edges[2];
    char name[H4_MAX_NC_NAME];
    char sdsname[H4_MAX_NC_NAME];
    char title[255];
    int16 month;
    /*
    char *path;
    char file[FILENAME_MAX] = "";
    */

    firstCall = 0;

    if (day < 1 || day > 366)
      {
      printf("-E- %s:  Bogus day number for ozone look-up: %d\n", __FILE__,
        day);
      exit(1);
      }

    /*
    if ((path = getenv("OCDATAROOT")) == NULL)
      {
      printf("-E- %s:  Error looking up environmental variable OCDATAROOT\n",
        __FILE__);
      exit(1);
      }
    strcpy(file, path);
    strcat(file, "/common/ozone_climatology.hdf");
    */

    sd_id = SDstart(file, DFACC_RDONLY);
    if (sd_id == -1)
      {
      printf("-E- %s:  Error openin file %s.\n", __FILE__, file);
      exit(1);
      }

    printf("\nOpening ozone file %s\n\n", file);

    sprintf(sdsname, "ozone_mean_%03d", day);
    sds_index = SDnametoindex(sd_id, sdsname);
    if (sds_index == -1)
      {
      printf("-E- %s:  Error seeking %s SDS from %s.\n", __FILE__, sdsname,
        file);
      exit(1);
      }
    sds_id = SDselect(sd_id, sds_index);

    status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
    if (status != 0)
      {
      printf("-E- %s:  Error reading SDS info for %s from %s.\n", __FILE__,
        sdsname, file);
      exit(1);
      }
    ny = dims[0];
    nx = dims[1];
    mapx = nx + 2;
    if( nx < 0 || ny < 0 ) 
      {
      printf("-E- %s:  grid has bad dimensions, sds %s from %s.\n", __FILE__,
        sdsname, file);
      printf("  Reading   %d x %d\n", nx, ny );
      exit(1);
      }

    dx = 360. / nx;
    dy = 180. / ny;
   /*
    *  allocate the storage for initial read and final interplation grid
    */
    if( ( ( tmp = (float *) malloc( nx * ny * sizeof(float) ) ) == NULL ) ||
        ( ( map = (float *) malloc( mapx * ( ny + 2 ) * sizeof(float) ) ) 
      == NULL ) )
      {
      printf(
        "-E- %s, %d: Unable to allocate space for climatology grid storage\n", 
        __FILE__, __LINE__ );
      exit(1);
      }

    start[0] = 0;
    start[1] = 0;
    edges[0] = ny;
    edges[1] = nx;

    status = SDreaddata(sds_id, start, NULL, edges, (VOIDP) tmp );
    if (status != 0)
      {
      printf("-E- %s:  Error reading SDS %s from %s.\n", __FILE__, sdsname,
        file);
      exit(1);
      }
   /*  upright the latitude during transfer ( ny - j - 1 part)  */
    for (j = 0; j < ny; j++)
      {
      for (i = 0; i < nx; i++)
        {
        *( map + i + 1 + mapx * ( j + 1 ) ) = *( tmp + i + nx * j );
        *( map + i + 1 + mapx * ( j + 1 ) ) = 
          *( tmp + i + nx * ( ny - j - 1 ) );
        }
      }

    /* add boarders to simplify interpolation - on sides (lon),
       duplicate the column from the other side, and at top, bottom (N, S),
       duplicate nearest neighbor  */

    for (j = 0; j < ny; j++)
      {
      *( map + mapx * ( j + 1 ) ) = *( map + nx + mapx * ( j + 1 ) );
      *( map + ( nx + 1 ) + mapx * ( j + 1 ) ) =
         *( map + 1 + mapx * ( j + 1 ) );
      }
    for (i = 0; i < mapx; i++)
      {
      *( map + i ) = *( map + mapx + i );
      *( map + i + mapx * ( ny + 1 ) ) = *( map + i + mapx * ny );
      }

    free(tmp);
    }

  /* interpolate to pixel location */

  i = MAX(MIN((int) ((lon + 180.0 + dx / 2) / dx), nx + 1), 0);
  j = MAX(MIN((int) ((lat + 90.0 + dy / 2) / dy), ny + 1), 0);

  xx = i * dx - 180.0 - dx / 2;
  yy = j * dy - 90.0 - dy / 2;

  t = (lon - xx) / dx;
  u = (lat - yy) / dy;

  val_11 = *( map + i + mapx * j );  /* [i,j] */
  val_12 = *( map + i + 1 + mapx * j ); /* [i+1,j] */
  val_21 = *( map + i + mapx * ( j + 1 ) );  /* [i,j+1] */
  val_22 = *( map + i + 1 + mapx * ( j + 1 ) );  /* [i+1,j+1] */

  if( ( val_11 < 0 ) || ( val_12 < 0 ) || ( val_21 < 0 ) || ( val_22 < 0 ) )
    {
    ozone = BAD_FLT;
    printf( "I %s, %d: Attempt to use missing climatology points\n",
      __FILE__, __LINE__ );
    }
  else
    ozone = (1 - t) * (1 - u) * val_11 + t * (1 - u) * val_12
      + t * u * val_22 + (1 - t) * u * val_21;
  return (ozone);
  }

/* -----------------------------------------------------------------------
 function setanc
 
 Returns windspeed, pressure, rel. humidity, water vapor, and ozone
 concentration for a specified time and place.  Uses climatology
 or realtime ancillary data files in SeaWiFS hdf format.  Ancillary
 file names are retrieved from the ancfiles.dat file.
 
 Returns 1 on error.
 
 Inputs:
 
 Outputs:
 zw(npix)	Zonal Wind		m/s
 mw(npix)	Meridional Wind		m/s
 ws(npix)	Wind Speed		m/s
 pr(npix)	Surface Pressure	millibars
 rh(npix)	Relative Humidity	%
 pw(npix)	Water Vapor  	        g/cm^2
 oz(npix)	Ozone			atm-cm
 no2(npix)	NO2			molecules cm^-2
 
 
 Written By: BA Franz, GSC, 6/97
 Conversion to C: G Fu, GSC, 3/99
 
 ----------------------------------------------------------------------- */

int setanc(l1str *l1rec, instr *input)
  {
  static float r2d = 57.29577951;
  static int firstCall = 1;
  static int32_t anc_id[] =
    {
    -1, -1
    };
  static short *ancqc = NULL;
  static char *no2file = NULL;
  static short npix;            // need this for l3gen  which changes the size of a line


  int32_t i, retval, status, iproc;
  float *lat = (float *) l1rec->lat;
  float *lon = (float *) l1rec->lon;
  short year = (short) *l1rec->year;
  short jday = (short) *l1rec->day;
  int32_t msec = (int32_t) * l1rec->msec;
  short parmID;

  status = 0;

  if (firstCall)
    {
    firstCall = 0;
    npix = (short) l1rec->npix;
    printf("\nOpening meteorological files.\n");
    printf("  met1   = %s\n", input->met1);
    printf("  met2   = %s\n", input->met2);
    printf("  met3   = %s\n", input->met3);
    printf("  ozone1 = %s\n", input->ozone1);
    printf("  ozone2 = %s\n", input->ozone2);
    printf("  ozone3 = %s\n", input->ozone3);
    printf("  no2    = %s\n", input->no2file);
    printf("\n");
    if ((ancqc = calloc(npix, sizeof(short))) == NULL)
      {
      fprintf(stderr, "-E- %s %d: Unable to allocate buffer space.\n", __FILE__,
        __LINE__);
      exit(1);
      }
    if (strcmp(input->no2file, "") != 0)
      {
      no2file = input->no2file;
      }
    /*
     *  do the setup and identification of alternate anc files
     *  (currently, only the ECMWF source)
     */
    if (anc_acq_init(input, anc_id) == -1)
      return 1;
    } else if(l1rec->npix > npix) {
        npix = (short) l1rec->npix;
        free(ancqc);
        if ((ancqc = calloc(npix, sizeof(short))) == NULL)
          {
          fprintf(stderr, "-E- %s %d: Unable to allocate buffer space.\n", __FILE__,
            __LINE__);
          exit(1);
          }
    }


  /*
   *  get standard met if anc_id = 1 or get met from ECMWF if anc_id = 0
   */
  if (anc_id[0] == 0)
    {
    if (anc_acq_lin(0, l1rec) == -1) /* arg 1 = 0 to access the met data */
      return 1;
    }
  else
    {
    /* relative humidity */
    /* ----------------- */

    parmID = 5;
    retval = get_ancillary(lat, lon, npix, year, jday, jday, msec, input->met1,
      input->met2, input->met3, input->anc_cor_file, parmID, l1rec->rh, 
      l1rec->ancqc);
    if (retval != 0)
      {
      fprintf(stderr,
        "-E- %s %d: Error loading relative humidity ancillary data. %s\n",
        __FILE__, __LINE__, input->anc_cor_file);
      status = 1;
      }

    /* wind speed */
    /* ---------- */

    parmID = 0;
    retval = get_ancillary(lat, lon, npix, year, jday, jday, msec, input->met1,
      input->met2, input->met3, input->anc_cor_file, parmID, l1rec->zw, ancqc);
    if (retval != 0)
      {
      fprintf(stderr,
        "-E- %s %d: Error loading Zonal wind speed ancillary data.\n", __FILE__,
        __LINE__);
      status = 1;
      }

    parmID = 1;
    retval = get_ancillary(lat, lon, npix, year, jday, jday, msec, input->met1,
      input->met2, input->met3, input->anc_cor_file, parmID, l1rec->mw, ancqc);
    if (retval != 0)
      {
      fprintf(stderr,
        "-E- %s %d: Error loading Meridional wind speed ancillary data.\n",
        __FILE__, __LINE__);
      status = 1;
      }

    for (i = 0; i < npix; i++)
      {
      l1rec->ws[i] = sqrt(
        l1rec->zw[i] * l1rec->zw[i] + l1rec->mw[i] * l1rec->mw[i]);
      l1rec->wd[i] = atan2f(-l1rec->zw[i], l1rec->mw[i]) * r2d;
      l1rec->ancqc[i] |= ancqc[i];
      }

    /* surface pressure */
    /* ---------------- */

    parmID = 2;
    retval = get_ancillary(lat, lon, npix, year, jday, jday, msec, input->met1,
      input->met2, input->met3, input->anc_cor_file, parmID, l1rec->pr, ancqc);
    if (retval != 0)
      {
      fprintf(stderr,
        "-E- %s %d: Error loading surface pressure ancillary data.\n", __FILE__,
        __LINE__);
      status = 1;
      }

    for (i = 0; i < npix; i++)
      {

      if (l1rec->pr[i] <= 0.0 || isnan(l1rec->pr[i]))
        l1rec->pr[i] = 1013.25;
      else if (l1rec->pr[i] < 900.0)
        l1rec->pr[i] = 900.0;
      else if (l1rec->pr[i] > 1100.0)
        l1rec->pr[i] = 1100.0;

      l1rec->ancqc[i] |= ancqc[i];

      /* if processing land, adjust pressure for terrain height */
      if (input->proc_land && l1rec->height[i] != 0.0)
        {
        l1rec->pr[i] *= exp(-l1rec->height[i] / 8434);
        }

      }

    /* precipitable water (water vapor) */
    /* -------------------------------- */

    parmID = 3;
    retval = get_ancillary(lat, lon, npix, year, jday, jday, msec, input->met1,
      input->met2, input->met3, input->anc_cor_file, parmID, l1rec->wv, ancqc);
    if (retval != 0)
      {
      fprintf(stderr,
        "-E- %s %d: Error loading precipitable water ancillary data.\n",
        __FILE__, __LINE__);
      status = 1;
      }

    /* convert from kg/m^2 to g/cm^2 */
    for (i = 0; i < npix; i++)
      {
      l1rec->wv[i] = l1rec->wv[i] / 10.0;
      l1rec->ancqc[i] |= ancqc[i];
      }
    } /* end of met portion ingest */

  /* ozone */
  /* ----- */
  if (anc_id[1] == 0)
    {
    if( anc_acq_lin( 1, l1rec ) == -1 )  /* arg 1 = 1 to acess ozone data */
      return 1;
      }
    else
      {
      if (strstr(input->ozone1,"ozone_climatology") != NULL)
        {
        for (i=0; i<npix; i++)
          {
	  l1rec->oz[i] = ozone_climatology(input->ozone1,*(l1rec->day),l1rec->lon[i],l1rec->lat[i]);
          if (l1rec->oz[i] < 0.0)
          ancqc[i] = 1;
          else
          ancqc[i] = 0;
          }
        }
      else
        {
        parmID = 4;
        retval = get_ancillary(lat, lon, npix, year, jday, jday, msec,
          input->ozone1, input->ozone2, input->ozone3, input->anc_cor_file,
          parmID, l1rec->oz, l1rec->ancqc);
        if (retval != 0)
          {
          fprintf(stderr,"-E- %s %d: Error loading Ozone ancillary data.\n",
            __FILE__,__LINE__);
          status = 1;
          }
        }

      /* convert from Dobson units to atm-cm */
      for (i=0; i<npix; i++)
        {
        l1rec->oz[i] = l1rec->oz[i] / 1000.0;
        l1rec->ancqc[i] |= ancqc[i];
        }

      }
   /*
    *  where l1rec->ancqc is set, turn on ATMFAIL in the flags
    */
    for (i = 0; i < npix; i++)
      if( l1rec->ancqc[i] != 0 )
        l1rec->flags[i] |= ATMWARN; 

    /* no2 and fraction */
    /* ---------------- */

  if ((input->gas_opt & NO2_BIT) != 0)
    for (i = 0; i < npix; i++)
      {
      no2conc(no2file, l1rec->lon[i], l1rec->lat[i], *(l1rec->day),
        &l1rec->no2_tropo[i], &l1rec->no2_strat[i]);
      l1rec->no2_tropo[i] *= 1e15;
      l1rec->no2_strat[i] *= 1e15;
      no2_frac(l1rec->lon[i], l1rec->lat[i], &l1rec->no2_frac[i]);
      }

  if (status != 0)
    {
    fprintf(stderr, "-E- %s %d: Error loading ancillary data.\n", __FILE__,
      __LINE__);
    return (1);
    }

  return (0);
  }

