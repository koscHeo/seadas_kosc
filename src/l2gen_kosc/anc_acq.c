/* ========================================================================= */
/* module anc_acq.c - functions to read alternate ancillary data             */
/*     Intended to initially do the MET and OZ with ECMWF data               */
/*                                                                           */
/* Written By: W. Robinson, SAIC, Aug, 2013                                  */
/*                                                                           */
/* ========================================================================= */
#include "l12_proto.h"
#include "met_cvt.h"
   /* status of the 3 ancillary files: MET1, 2, 3
      or OZONE1, 2, 3 */
#define ANC_STAT_1T 0  /* one time used */
#define ANC_STAT_2T_END 1  /* 2 different anc times in end of list 
                              = files[0] = [1] != [2] */
#define ANC_STAT_2T_START 2  /*2 different anc times in start of list
                              = files[0] ! [1] = [2] */
#define ANC_STAT_3T 3  /* all 3 files with different times */
#define ANC_STAT_CLIM 4  /* a climatology in  use */
#define NPRM 7
#define ANCBAD -999.
#define OZ_KG_M2_TO_DU   1. / 2.1415e-5
#define USE_PMSL 1  /* choice to use surface pressure (0)  or use the MSL 
                       (mean sea level) pressure (1) with appropriate 
                       adjustment for land height */
struct met_sto_str_d {
  float s_lon;  /* start longitude for a grid */
  float lon_step;  /* longitude incriment */
  int nlon;   /* convenience - # longitude points  */
  float e_lon;   /* end longitude for a grid */
  float s_lat;  /* start latitude for a grid */
  float lat_step;  /* latitude incriment */
  int nlat;  /* convenience - # latitude points  */
  float e_lat;   /* end latitude for a grid */
  double data_time[3];  /* time (in Julian days and fraction)
                           for the data 1, 2, 3 */
  int anc_f_stat;  /* status of the met data */
  float *data[3];  /*  storage for MET1, 2, 3 */
  };

typedef struct met_sto_str_d met_sto_str;

static met_sto_str met_sto[NPRM];
static int proc_land;

int anc_acq_init( instr *input, int32_t *anc_id )
/*******************************************************************

 anc_acq_init

 purpose: Identify the incoming MET files and if netcdf, set up the
   ancillary data so it can be accessed by anc_acq_line
   For now, only ECMWF netcdf files can be processed which contain
   only 1 time and (at least) the parameters listed in prm_nm

 Returns type: int - 0 - good -1 any trouble checking input anc files

 Parameters: (in calling order)
 Type              Name            I/O     Description
 ----              ----            ---     -----------
 instr *           input            I      program inputs
 int32_t *         anc_id           O      size 2 ID of anc data being handled: 
                                           for [met, ozone]
                                           0 - ECMWF files, 1 - NCEP/TOMS
                                           files, -1 - bad

 Modification history:
 Programmer        Date            Description of change
 ----------        ----            ---------------------
 W. Robinson       12-Aug-2013     Original development

 *******************************************************************/
  {
  char *files[3];  /* the 3 MET file names */
  /* ECMWF parameter names in the 2 arrays below */
  /*  name   descrip   ends up as value in met_sto
      sp     sfc pressure           pressure at surface (not currently used)
      tcwv   precip water           precip water
      msl    MSL pressure           seal lvl pressure
      u10    10 m u wind            uwind meridional N-S
      v10    10 m v wind            vwind zonal E-W
      t2m    2m temp           -    Relative humidity
      d2m    2m dewpoint temp  /  
      tco3   total Ozone            ozone
   */   
  char *prm_nm_met[] = { "sp", "tcwv", "msl", "u10", "v10", "t2m", "d2m" };
  char *prm_nm_oz[] = { "tco3" };
  int32_t n_prm_met = 7, n_prm_oz = 1, sto_ix;
 /*
  *  First, identify the type of the MET ancillary data
  */
  proc_land = input->proc_land;  /* carry to the anc_acq_lin routine */
  files[0] = input->met1;
  files[1] = input->met2;
  files[2] = input->met3;
  sto_ix = 0;  /* location of met in storage struct */

  anc_id[0] = anc_acq_ck( files, prm_nm_met, n_prm_met, sto_ix );
 /*
  *  Same for the ozone data
  */
  files[0] = input->ozone1;
  files[1] = input->ozone2;
  files[2] = input->ozone3;
  sto_ix = 6;  /* location of oz in storage struct */

  anc_id[1] = anc_acq_ck( files, prm_nm_oz, n_prm_oz, sto_ix );

  if( ( anc_id[0] == -1 ) || ( anc_id[1] == -1 ) )
    return -1;
  else
    return 0;
  }

int32_t anc_acq_ck( char **files, char **prm_nm, int n_prm, int32_t sto_ix )
/*******************************************************************

 anc_acq_ck

 purpose: Identify the incoming MET or OZONE file set and if netcdf, 
   set up the ancillary data so it can be accessed by anc_acq_line.
   Otherwise, do nothing and have getanc.c process the NCEP/TOMS data.
   For now, only ECMWF netcdf files can be processed which contain
   only 1 time and (at least) the parameters listed in prm_nm

 Returns type: int ancillary data identification: 0 - ECMWF data, 
  1 non-ECMWF data, -1 any trouble checking input anc files

 Parameters: (in calling order)
 Type              Name            I/O     Description
 ----              ----            ---     -----------
 char **           files            I      3 file name set, either MET 
                                           or OZONE source
 char **           prm_nm           I      array of parameter names to read 
                                           from the ECMWF
 int               n_prm            I      # parameters in prm_nm.  Note that
                                             the met parms (n_prm >1) have 7
                                             fields to read but the t(2 m) and
                                             td(2 m) combine to make RH so the 
                                             storage for met is 1 less 
                                             (nprm_sto)
 int32_t           sto_ix           I      position in storage array to place
                                           result
 int32_t           anc_typ          O      type of anc data 0 - ECMWF, 
                                           1 - non ECMWF, -1 - problem

 Modification history:
 Programmer        Date            Description of change
 ----------        ----            ---------------------
 W. Robinson       24-Sep-2013     Original development

 *******************************************************************/
  {
  int ids, iprm, npix, nlin, ilin, ipix;
  int t_days, t_hrs, lon_gt_180, ird;
  int dstpix, npix0, nlin0, status;
  int npix_ext, nlin_ext, ntim, f12_mtch, f23_mtch, anc_f_stat;
  int ncid, iprm_sto, nprm_sto;
  float s_lon, lon_step, e_lon, s_lat, lat_step, e_lat, time;
  float prm_bad[] = { ANCBAD, ANCBAD, ANCBAD, ANCBAD, ANCBAD, ANCBAD,
    ANCBAD, ANCBAD };
  float *base_data, *lat, *lon, *comp1, *comp2;
  double data_time;
  int64_t jd1900;

  nprm_sto = ( n_prm > 1 ) ? n_prm - 1 : n_prm;
 /*
  *  identify the type of files set from their existance and name matching
  */
  if( ( files[1] == NULL ) || ( files[2] == NULL ) ||
      ( files[1][0] == 0 ) || ( files[2][0] == 0 ) )
    anc_f_stat = ANC_STAT_CLIM;
  else
    {
    f12_mtch = 0;
    if( strcmp( files[0], files[1] ) == 0 ) f12_mtch = 1;
    f23_mtch = 0;
    if( strcmp( files[1], files[2] ) == 0 ) f23_mtch = 1;
  
    if( ( strcmp( files[0], files[2] ) == 0 ) && ( f12_mtch == 0 ) )
      {
      printf( "%s, %d E: ANC1 and 3 match while ANC2 different\n",
        __FILE__, __LINE__ );
      return -1;
      }
    if( ( f12_mtch == 1 ) && ( f23_mtch == 1 ) ) anc_f_stat = ANC_STAT_1T;
    else if( ( f12_mtch == 1 ) && ( f23_mtch == 0 ) )
      anc_f_stat = ANC_STAT_2T_END;
    else if( ( f12_mtch == 0 ) && ( f23_mtch == 1 ) )
      anc_f_stat = ANC_STAT_2T_START;
    else anc_f_stat = ANC_STAT_3T;
    }
 /*  END of file set identification  */

  if( anc_f_stat == ANC_STAT_CLIM )
    {
/*
    printf( "%s, %d I: Assuming standard (not ECMWF) climatology\n",
      __FILE__, __LINE__ );
*/
    return 1;
    }
  jd1900 = jd4713bc_get_jd( 1900, 1, 1 );
  for( ids = 0; ids < 3; ids++ )
    {
   /*
    *  ECMWF is only anc data in netcdf format for now, so if it opens
    *  it should be ECMWF format
    */
    if( ( ids > 0 ) && ( anc_f_stat== ANC_STAT_1T ) )
      break;
    if( ( ids == 1 ) && ( f12_mtch == 1 ) )
      continue;
    if( ( ids == 2 ) && ( f23_mtch == 1 ) )
      break;
  
    status = nc_open( files[ids], 0, &ncid );
    if( status != NC_NOERR )
      {
      if( ids == 0 ) 
        {
/*
        printf( 
        "%s, %d I: file: %s is not ECMWF, assuming standard ancillary file\n", 
          __FILE__, __LINE__, files[ids] );
*/
        return 1;
        }
      else
        {
        printf( "%s, %d: nc_open failed on file: %s\n", __FILE__, __LINE__,
          files[ids] );
        printf( "       Mismatch or bad ECMWF file\n" );
        return -1;
        }
      }
  
   /*
    *  get the basic information for MET1
    */
    if( ( ( npix0 = ncio_dim_siz( ncid, "longitude" ) ) == -1 ) ||
        ( ( nlin0 = ncio_dim_siz( ncid, "latitude" ) ) == -1 ) ||
        ( ( ntim = ncio_dim_siz( ncid, "time" ) ) == -1 ) )
      {
      printf( 
     "%s, %d: ncio_dim_siz error reading longitude, latitude or time datasets\n", 
        __FILE__, __LINE__ );
      return -1;
      }
    if( ids == 0 )
      {
      npix = npix0;
      nlin = nlin0;
      }
    else
      {
      if( ( npix != npix0 ) || ( nlin != nlin0 ) )
        {
        printf( "%s, %d: mismatch in size of MET array[%d]\n", __FILE__, 
          __LINE__, ids );
        return -1;
        }
      }
   /*
    * for now, if more than 1 time, we can't proces it
    */
    if( ntim > 1 )
      {
      printf( "%s, %d: Number of times > 1, can't deal with at this time\n",
        __FILE__, __LINE__ );
      return -1;
      }
    npix_ext = npix + 2;
    nlin_ext = nlin + 2;
   /*
    *  allocate storage in the structures for the data
    */
    for( iprm = 0; iprm < nprm_sto; iprm++ )
      {
      if( ( met_sto[ iprm + sto_ix ].data[ids] = ( float * ) 
        malloc( npix_ext * nlin_ext * sizeof(float) ) )  == NULL )
        {
        printf( "%s, %d: malloc failed for data[%d] in met_sto %d\n", __FILE__, 
          __LINE__, ids, iprm );
        return -1;
        }
      }
   /*
    *  for 1st dataset, make array to read data into initially and 
    *  arrays for lat, lon
    */
    if( ids == 0 )
      {
      if( ( ( base_data = ( float * ) 
            malloc( npix * nlin * sizeof(float) ) )  == NULL ) ||
          ( ( comp1 = ( float * )
            malloc( npix * nlin * sizeof(float) ) )  == NULL ) ||
          ( ( comp2= ( float * )
            malloc( npix * nlin * sizeof(float) ) )  == NULL ) )
        {
        printf( "%s, %d: malloc failed for base_data or comp arrays\n", 
          __FILE__, __LINE__ );
        return -1;
        }
    
      if( ( lat = ( float * ) malloc( nlin * sizeof(float) ) )  == NULL )
        {
        printf( "%s, %d: malloc failed for latitude\n", __FILE__,
          __LINE__ );
        return -1;
        }
      if( ( lon = ( float * ) malloc( npix * sizeof(float) ) )  == NULL )
        {
        printf( "%s, %d: malloc failed for longitude\n", __FILE__,
          __LINE__ );
        return -1;
        }
      if( ncio_grab_f_ds( ncid, "latitude", lat ) != 0 )
        {
        printf( "%s, %d: ncio_grab_f_ds failed on latitude\n", 
          __FILE__, __LINE__ );
        return -1;
        }
      if( ncio_grab_f_ds( ncid, "longitude", lon ) != 0 )
        {
        printf( "%s, %d: ncio_grab_f_ds failed on longitude\n",
          __FILE__, __LINE__ );
        return -1;
        }
     /*
      *  from the latitude, longitude arrays, determine the nav properties
      *  ECMWF longitudes go 0 -> 360 and we do -180 -> 180
      */
      s_lat = lat[0];
      e_lat = lat[ nlin - 1 ];
      lat_step = lat[1] - lat[0];
    
      lon_gt_180 = -1;
      for( ipix = 0; ipix < npix; ipix++ )
        {
        if( lon[ipix] > 180. )
          {
          s_lon = lon[ipix] - 360.;
          lon_gt_180 = ipix;
          e_lon = lon[ ipix - 1 ]; /* if lon_gt_180 = 0, need to upgrade this */
          break;
          }
        }
      lon_step = lon[1] - lon[0];
      }
   /*
    *  get the time from the dataset and convert to julian days
    */
    if( ncio_grab_f_ds( ncid, "time", &time ) != 0 )
      {
      printf( "%s, %d: error reading the time in\n", __FILE__, __LINE__ );
      return -1;
      }
  
    t_days = time / 24;
    t_hrs = (int) time % 24;
    data_time = (double) ( jd1900 + t_days ) + (double) t_hrs / 24.;
   /*
    *  for all params, read the data, put lon in range -180 - 180
    *  and add extra layer to make interpolation easier
    *  The RH is made from T and Td in else below
    */
    ird = 0;
    for( iprm = 0; iprm < nprm_sto; iprm++ )
      {
      iprm_sto = iprm + sto_ix;
      if( ird != 5 )
        {
        if( ncio_grab_stdsclf_ds( ncid, prm_nm[ ird ], 
            prm_bad[ird], base_data ) != 0 )
          {
          printf( "%s, %d: ncio_grab_stdsclf_ds failed on %s\n",
            __FILE__, __LINE__, prm_nm[ird] );
          return -1;
          }
        ird++;
        }
      else
        {
        if( ( ncio_grab_stdsclf_ds( ncid, prm_nm[ ird ], 
                 prm_bad[ ird + sto_ix ], comp1 ) != 0 ) ||
             ( ncio_grab_stdsclf_ds( ncid, prm_nm[ird + 1], 
                 prm_bad[ ird + sto_ix + 1 ], comp2 ) != 0 ) )
          {
          printf( "%s, %d: ncio_grab_stdsclf_ds failed on %s or %s\n",
            __FILE__, __LINE__, prm_nm[ird], prm_nm[ird + 1] );
          return -1;
          }
        ird = ird + 2;
     /*
      *  for the td and t at 2 m, we need to make a RH 
      */
        if( met_cvt_ttd_to_rh( npix * nlin, comp1, MET_UNITS__T_K, comp2,
            MET_UNITS__T_K, base_data ) != 0 )
          {
          printf( "met_cvt_ttd_to_rh had an error\n" );
          printf( "%s, %d: met_cvt_ttd_to_rh failure\n",
            __FILE__, __LINE__ );
          return -1;
          }
        }
     /*  rotate to -180 -> 180 */
      for( ilin = 0; ilin < nlin; ilin++ )
        {
        for( ipix = 0; ipix < npix; ipix++ )
          {
          dstpix = ipix - lon_gt_180;  /* put in with lon -180 -> 180 */
          if( dstpix < 0 ) dstpix += npix;
          *( met_sto[ iprm_sto ].data[ids] 
            + dstpix + 1 + ( ilin + 1 ) * npix_ext ) =
            *( base_data + ipix + npix * ilin );
          }
        }
     /*  now, the extra boarder: lat, then lon */
     /* for lat, repeat the nearest value */
      for( ipix = 0; ipix < npix; ipix++ )
        {
        *( met_sto[iprm_sto].data[ids] + ipix + 1 ) = 
          *( met_sto[iprm_sto].data[ids] + ipix + 1 + npix_ext );
        *( met_sto[iprm_sto].data[ids] + ipix + 1 + ( nlin + 1 ) * npix_ext ) =
          *( met_sto[iprm_sto].data[ids] + ipix + 1 + nlin * npix_ext );
        }
     /* for lon, use the opposite side value */
      for( ilin = 0; ilin < nlin_ext; ilin++ )
        {
        *( met_sto[iprm_sto].data[ids] + ilin * npix_ext ) = 
          *( met_sto[iprm_sto].data[ids] + npix + ilin * npix_ext );
        *( met_sto[iprm_sto].data[ids] + npix + 1 + ilin * npix_ext ) =
          *( met_sto[iprm_sto].data[ids] + 1 + ilin * npix_ext );
        }
     /*  put in the controls found above  */
      met_sto[iprm_sto].s_lon = s_lon;
      met_sto[iprm_sto].lon_step = lon_step;
      met_sto[iprm_sto].nlon = npix;
      met_sto[iprm_sto].e_lon = e_lon;
      met_sto[iprm_sto].s_lat = s_lat;
      met_sto[iprm_sto].lat_step = lat_step;
      met_sto[iprm_sto].nlat = nlin;
      met_sto[iprm_sto].e_lat = e_lat;
      met_sto[iprm_sto].data_time[ids] = data_time;
      met_sto[iprm_sto].anc_f_stat = anc_f_stat;
      }
   /*
    *  close the dataset
    */
    nc_close( ncid );
    }
  return 0;
  }

int anc_acq_lin( int32_t anc_class, l1str *l1rec )
/*******************************************************************

   anc_acq_lin
	
   purpose: get proper ancillary parameters for a particular line of 
     points

   Returns type: int - 0 if good, else -1

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int32_t           anc_class        I      anc data class to access the 
                                                correct stored grids: 0 for 
                                                MET grids and 1 for ozone grid
      l1str *           l1rec           I/O     structure with information
                                                for the line
   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       16 Aug 2013     original development

*******************************************************************/
  {
  double l_time;
  float data_val, uwnd, vwnd, data_val1, data_val2, dx, dy, lon_frac, lat_frac;
  float trg_lon, trg_lat, wt_t1, wt_t2, s_lon, s_lat;
  float *data;
  int iprm, xbox_st, ybox_st, nx, ny, t_interp, data1_ix, data2_ix, anc_f_stat;
  int npix, ipix, nprm, sto_st, sto_en;
 /*
  *  find places in the met_sto structure to look in
  */
  if( anc_class == 0 )
    {
    sto_st = 0;
    sto_en = 5;
    }
  else
    {
    sto_st = 6;
    sto_en = 6;
    }
 /*
  *  get the time of the current line
  */
  l_time = (double) jd4713bc_get_jd( *l1rec->year, 1, *l1rec->day );
  l_time += (double) *l1rec->msec / ( 86400.0 * 1000.0);

  npix = l1rec->npix;
 /*
  *  for this line, decide which of the 3 anc files will be needed based on 
  *  the line's time and the ancillary times
  */
  /*   ***** In this set up, all grids are the same, so only one
             determination is needed
   */
  anc_f_stat = met_sto[sto_st].anc_f_stat;
  if( anc_f_stat == ANC_STAT_1T )
    {
    /* use data[0] only */
    t_interp = 0;
    data1_ix = 0;
    /*  further along, when interpolating, use met_sto[0].data[data1_ix]
        to access the data */
    }
  else if( anc_f_stat == ANC_STAT_2T_START )
    {
   /* 2 different times in the 1, 2 positions */
    if( l_time < met_sto[sto_st].data_time[0] )
      {
      printf( "%s, %d: data time is before the ancillary data start time\n",
        __FILE__, __LINE__ );
      return -1;
      }
    else if( l_time > met_sto[sto_st].data_time[1] )
      {
     /* use MET2 only */
      t_interp = 0;
      data1_ix = 1;
      }
    else
      {
     /* in-between MET 1, 2 use data 0, 1 and time interpolate */
      t_interp = 1;
      data1_ix = 0;
      data2_ix = 1;
      wt_t1 = ( met_sto[sto_st].data_time[1] - l_time ) /
              ( met_sto[sto_st].data_time[1] - met_sto[sto_st].data_time[0] );
      wt_t2 = ( l_time - met_sto[sto_st].data_time[0] ) /
              ( met_sto[sto_st].data_time[1] - met_sto[sto_st].data_time[0] );
      }
    }
  else if( anc_f_stat == ANC_STAT_2T_END )
    {
    if( l_time < met_sto[sto_st].data_time[0] )
      {
     /* outside on the low end, use data[0] */
      t_interp = 0;
      data1_ix = 0;
      }
    else if( l_time > met_sto[sto_st].data_time[2] )
      {
     /* beyond the high end, Can't use end time alone */
      printf( "%s, %d: data time is after the ancillary data end time\n",
        __FILE__, __LINE__ );
      return -1;
      }
    else
      {
     /* between the MET 1 and 3 */
      t_interp = 1;
      data1_ix = 0;
      data2_ix = 2;
      wt_t1 = ( met_sto[sto_st].data_time[2] - l_time ) /
              ( met_sto[sto_st].data_time[2] - met_sto[sto_st].data_time[0] );
      wt_t2 = ( l_time - met_sto[sto_st].data_time[0] ) /
              ( met_sto[sto_st].data_time[2] - met_sto[sto_st].data_time[0] );
      }
    }
  else if( anc_f_stat == ANC_STAT_3T )
    {
    if( l_time < met_sto[sto_st].data_time[0] )
      {
      printf( "%s, %d: data time is before the ancillary data start time\n",
        __FILE__, __LINE__ );
      return -1;
      }
    else if( l_time > met_sto[sto_st].data_time[2] )
      {
      printf( "%s, %d: data time is after the ancillary data end time\n",
        __FILE__, __LINE__ );
      return -1;
      }
    else if( l_time < met_sto[sto_st].data_time[1] )
      {
     /* between data 0 and 1 */
      t_interp = 1;
      data1_ix = 0;
      data2_ix = 1;
      wt_t1 = ( met_sto[sto_st].data_time[1] - l_time ) /
              ( met_sto[sto_st].data_time[1] - met_sto[sto_st].data_time[0] );
      wt_t2 = ( l_time - met_sto[sto_st].data_time[0] ) /
              ( met_sto[sto_st].data_time[1] - met_sto[sto_st].data_time[0] );
      }
    else
      {
     /* what's left: between data 1 and 2 */
      t_interp = 1;
      data1_ix = 1;
      data2_ix = 2;
      wt_t1 = ( met_sto[sto_st].data_time[2] - l_time ) /
              ( met_sto[sto_st].data_time[2] - met_sto[sto_st].data_time[1] );
      wt_t2 = ( l_time - met_sto[sto_st].data_time[1] ) /
              ( met_sto[sto_st].data_time[2] - met_sto[sto_st].data_time[1] );
      }
    }
  else
    {
   /* this should not happen at this time - a status that is either 
      climatology or not defined */
    printf( "%s, %d: Undefined anc_f_stat - should not happen\n", 
      __FILE__, __LINE__ );
    return -1;
    }
 /*
  *  this found if time interpolation is needed and the and grids to use.
  *  next, for each pixel, find the bounding grid box and weights
  *  for bi-linear interpolation to be applied to each parameter
  *
  *  AGAIN note that since all parameters are from the same size grid,
  *  we can compute this information once.
  */
  dx = met_sto[sto_st].lon_step;
  dy = met_sto[sto_st].lat_step;
  s_lon = met_sto[sto_st].s_lon;
  s_lat = met_sto[sto_st].s_lat;
  nx = met_sto[sto_st].nlon;
  ny = met_sto[sto_st].nlat;
  for( ipix = 0; ipix < npix; ipix++ )
    {
    trg_lat = l1rec->lat[ipix];
    trg_lon = l1rec->lon[ipix];

/*
    xbox_st = 
      MAX( MIN( (INT)( ( trg_lon - s_lon + dx / 2. ) / dx ), nx + 1 ), 0 );
    ybox_st = 
      MAX( MIN( (INT)( ( trg_lat - s_lat + dy / 2. ) / dy ), ny + 1 ), 0 );
    x_dist = xbox_st * dx + s_lon - dx / 2;
    y_dist = ybox_st * dy + s_lat - dy / 2;

    I think below is correct for data at the grid points
*/
    xbox_st =
      MAX( MIN( (int)( ( trg_lon - s_lon + dx ) / dx ), nx + 1 ), 0 );
    ybox_st =
      MAX( MIN( (int)( ( trg_lat - s_lat + dy ) / dy ), ny + 1 ), 0 );

    lon_frac = ( trg_lon - s_lon ) / dx - (float) ( xbox_st - 1 );
    lat_frac = ( trg_lat - s_lat ) / dy - (float) ( ybox_st - 1 );

    for( iprm = sto_st; iprm < ( sto_en + 1 ); iprm++ )
      {
      data = met_sto[iprm ].data[data1_ix];
      data_val1 = bilin_interp( data, xbox_st, ( nx + 2 ), ybox_st, lon_frac, 
        lat_frac );

      if( t_interp == 1 )
        {
        data= met_sto[iprm].data[data2_ix];
        data_val2 = bilin_interp( data, xbox_st, ( nx + 2 ), ybox_st, lon_frac, 
          lat_frac );

       /*
        *  do time interpolation
        */
        if( data_val1 < ANCBAD + 1 )
          {
          if( data_val2 < ANCBAD + 1 )
            data_val = ANCBAD;
          else
            data_val = data_val2;
          }
        else
          data_val = wt_t1 * data_val1 + wt_t2 * data_val2;
        }
      else
        data_val = data_val1;
     /*
      *  place this interpolated value in proper l1rec slot
      */
      switch (iprm) 
        {
        case 0:  /*  sfc press */
         /*  Currently, no use for this, but it may be better than what
             is done with MSL pressure to take it to a height above sea 
             level.  USE_PMSL of 0 will use this.  Note that pressure on Mt
             Everest is nominally 337 mb, so enlarge range accordingly */
          if( USE_PMSL == 1 )
            break;
          else
            {
            data_val = ( data_val < ANCBAD + 1 ) ? ANCBAD : data_val / 100.;
            if( data_val < 0 ) l1rec->pr[ipix] = 1013.25;
            else if( data_val < 250. ) l1rec->pr[ipix] = 250.;
            else if( data_val > 1100. ) l1rec->pr[ipix] = 1100.;
            else l1rec->pr[ipix] = data_val;
            }
          break;
        case 1:  /*  precip water  */
         /* need to make from kg m^-2 into g cm^-2 */
          l1rec->wv[ipix] = ( data_val  < ANCBAD + 1 ) ? 0. : data_val / 10.;
          break;
        case 2:  /*  sea level pressure  */
         /* need to make from pascals to hectopascals (millibars) */
          if( USE_PMSL == 0 )
            break;
          else
            {
            data_val = ( data_val < ANCBAD + 1 ) ? ANCBAD : data_val / 100.;
            if( data_val < 0 ) l1rec->pr[ipix] = 1013.25;
            else if( data_val < 900. ) l1rec->pr[ipix] = 900.;
            else if( data_val > 1100. ) l1rec->pr[ipix] = 1100.;
            else l1rec->pr[ipix] = data_val;

            /* if processing land, adjust pressure for terrain height */
            if( proc_land && l1rec->height[ipix] != 0.0 )
              l1rec->pr[ipix] *= exp( -l1rec->height[ipix] / 8434 );
            }
          break;
        case 3:  /*  u wind, zonal W-E */
          uwnd = ( data_val < ANCBAD + 1 ) ? 0. : data_val;
          l1rec->zw[ipix] = uwnd;
          break;
        case 4:  /*  v wind, meridional S-N */
          vwnd = ( data_val < ANCBAD + 1 ) ? 0. : data_val;
          l1rec->mw[ipix] = vwnd;
  
          l1rec->ws[ipix] = sqrt( pow( uwnd, 2. ) + pow( vwnd, 2. ) );
          break;
        case 5:  /*  rel humidity % */
          l1rec->rh[ipix] = ( data_val < ANCBAD + 1 ) ? 0. : data_val;
          break;
        case 6:  /*  ozone  */
         /*  convert from kg m^-2 to DU (units of 10 um STP OZ thickness)
             to cm of thickness */
          l1rec->oz[ipix] = ( data_val < ANCBAD + 1 ) ? 
            0. : data_val * OZ_KG_M2_TO_DU / 1000.;
          break;
        }
      }
    }
  return 0;
  }

float bilin_interp( float *data, int xbox_st, int nx, int ybox_st, 
  float xfrac, float yfrac )
/*******************************************************************

   bilin_interp

   purpose: quick bi-linear interpolation.

   Returns type: float of interpolated result

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      float *           data             I      2-d grid of data
      int               xbox_st          I      x (longitude) index of 
                                                grid box to interpolate
      int               nx               I      # pixels in x
      int               ybox_st          I      x (latitude) index of
                                                grid box to interpolate
      float             xfrac            I      fractional grid box distance 
                                                from xbox_st to the point
      float             yfrac            I      fractional grid box distance 
                                                from ybox_st to the point

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       16 Aug 2013     original development

*******************************************************************/
  {
  float data_val;

  if( ( *( data + xbox_st + nx * ybox_st ) < ( ANCBAD + 1 ) ) ||
      ( *( data + xbox_st + nx * ( ybox_st + 1 ) ) < ( ANCBAD + 1 ) ) ||
      ( *( data + ( xbox_st + 1 ) + nx * ybox_st ) < ( ANCBAD + 1 ) ) ||
      ( *( data + ( xbox_st + 1 ) + nx * ( ybox_st + 1 ) ) < ( ANCBAD + 1 ) ) )
    data_val = ANCBAD;
  else
    data_val = 
      ( 1 - xfrac ) * ( 1 - yfrac ) *
          *( data + xbox_st + nx * ybox_st ) +
      ( 1 - xfrac ) * yfrac *
          *( data + xbox_st + nx * ( ybox_st + 1 ) ) +
      xfrac * ( 1 - yfrac ) *
          *( data + ( xbox_st + 1 ) + nx * ybox_st ) +
      xfrac * yfrac *
          *( data + ( xbox_st + 1 ) + nx * ( ybox_st + 1 ) );
  return data_val;
  }

int64_t jd4713bc_get_jd( int32_t year, int32_t month, int32_t day )
/*******************************************************************

   jd4713bc_get_jd

   purpose: get the julian day (from 4713 BC) for the year,
     month and day of month
     taken from the idl jd routine

   Returns type: int64_t - the julian date

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int32_t           year             I      standard 4 digit year
      int32_t           month            I      month of the year
      int32_t           day              I      day of month

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       5 Aug 2013     original development

*******************************************************************/
  {
  int64_t lyear, lmonth, lday, jday;

  lyear = (int64_t) year;
  lmonth = (int64_t) month;
  lday = (int64_t) day;
  jday = ( 367 * lyear  - 7 * ( lyear + ( lmonth + 9 ) / 12 ) / 4 
    + 275 * lmonth / 9  + lday + 1721014 );
  /*
  *  this additional step is only needed if you expect to work on dates 
  *  outside March 1, 1900 to February 28, 2100
  */
  jday = jday + 15 - 3 * ( ( lyear + ( lmonth - 9 ) / 7 ) / 100 + 1 ) / 4;
  return jday;
  }

int jd4713bc_get_date( int64_t jd, int32_t *year, int32_t *month, int32_t *day )
/*******************************************************************

   jd4713bc_get_date

   purpose: get the year, month, day from julian date
     (from 4713 BC)
     taken from the idl jddate routine

   Returns type: int - no set value now

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int64_t           jd               I      julian date
      int32_t *         year             O      standard 4 digit year
      int32_t *         month            O      month
      int32_t *         day              O      day of month

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       5 Aug 2013     original development

*******************************************************************/
  {
  int64_t v1, v2, v3, v4;
  v1 = jd + 68569;
  v2 = 4 * v1 / 146097;
  v1 = v1 - ( 146097 * v2 + 3 ) / 4;
  v3 = 4000 * ( v1 + 1 ) / 1461001;
  v1 = v1 - 1461 * v3 / 4 + 31;
  v4 = 80 * v1 / 2447;
  *day = v1 - 2447 * v4 / 80;
  v1 = v4 / 11;
  *month = v4 + 2 - 12 * v1;
  *year = 100 * ( v2 - 49 ) + v3 + v1;
  return 0;
  }
