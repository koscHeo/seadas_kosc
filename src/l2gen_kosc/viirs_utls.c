/*---------------------------------------------------------------
viirs_utls.c   - some needed viirs utilities

W. Robinson, SAIC, 11 Feb 2009
-----------------------------------------------------------------*/
#include <stdint.h>
#include <timeutils.h>
#include <genutils.h>


struct leap_str_def
  {
  int ncut;  /* number of cutoff times in the file */
  int16_t *year;  /* date of the leapsecond change */
  int16_t *mon;  /* runs 0 - 11 */
  int16_t *day;  /* runs 1 = 31 */
  double *leap;  /* seconds to correct TAI - UTC */
  };
typedef struct leap_str_def leap_str;

void viirs_u58_yds( int64_t u58, short *year, short *day, double *dsec )
/* ------------------------------------------------------------------------
   viirs_u58_yds

   purpose: get the year, day of year and seconds of day from the viirs
            time format

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      long long         u58              I      time in microsecs past 
                                                1 jan 1958
      short *           year             O      Year
      short *           day              O      day of year
      double *          dsec             O      seconds of day 

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC, 11 Feb 2009    Original development
      W. Robinson, SAIC, 17 Apr 2012    Added a leap second correction
                                        NOTE - the 34 is for leap sec 
                       correction and needs to be replaced with something 
                       to change as more leap seconds are added.

 ------------------------------------------------------------------------*/
#define SEC_PER_DAY 86400
#define DEL_DAY 4383  /* difference in days between 1958 and 1970  */
  {
  double sec, leapoff;
  double get_leap_u58( int64_t );

  leapoff = get_leap_u58( u58 ); /* get leap sec from table */
  sec = ( (double) u58 ) / 1.e6 - leapoff - DEL_DAY * SEC_PER_DAY;
  unix2yds( sec, year, day, dsec );
  return;
  }

double get_leap_u58( int64_t u58 )
/* ------------------------------------------------------------------------
   get_leap_u58

   purpose: get the proper leap second corection for an incoming IET time

   returns leap seconds for the time

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      int64_t           u58              I      time in microsecs past
                                                1 jan 1958

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC, 15 May 2012    Original development

 ------------------------------------------------------------------------*/
  {
  static int lasti = 0;
  static int rd_leap_file = 0;
  static double *leap_lst;
  static int64_t *u58_cutoff;
  static leap_str leap_info;
  double sec_d = 0, leap;
  int rd_leapsec( leap_str * );
  int icut;
  double sec70;
 /*
  *  read in the leapsec.dat and make u58 versions of information needed
  */
  if( rd_leap_file == 0 )
    {
   /*
    *  read general leapsec.dat file information
    */
    if( rd_leapsec( &leap_info ) != 0 )
      {
      printf( "%s,%d: Unable to read the leapsec.dat file\n", 
        __FILE__, __LINE__ );
      rd_leap_file = 2;
      }
    else
      {
      rd_leap_file = 1;

      if( ( leap_lst = 
        ( double *) malloc( leap_info.ncut * sizeof( double ) ) ) == NULL )
        {
        printf( "%s,%d: Unable to allocate leap second storage, leap_lst\n",
          __FILE__, __LINE__ );
        exit(1);
        }
      if( ( u58_cutoff = 
        ( int64_t *) malloc( leap_info.ncut * sizeof( int64_t ) ) ) == NULL )
        {
        printf( "%s,%d: Unable to allocate leap second storage, u58_cutoff\n",
          __FILE__, __LINE__ );
        exit(1);
        }
     /*
      *  convert the general information into u58 time units
      */
      for( icut = 0; icut < leap_info.ncut; icut++ )
        {
        sec70 = ymds2unix( leap_info.year[icut], leap_info.mon[icut]+1,
          leap_info.day[icut], sec_d );
        *( u58_cutoff + icut ) = ( ( sec70 + DEL_DAY * SEC_PER_DAY ) + 
          leap_info.leap[icut] ) * 1.e6;
        *( leap_lst + icut ) = leap_info.leap[icut];
        }
      free( leap_info.year );
      free( leap_info.mon );
      free( leap_info.day );
      free( leap_info.leap );
      }
    }
 /*
  *  find the leap second in the proper range
  *  first, if the new time is below the last gotten time, either
  *  return 0. leap seconds (below the range) or reset the search start to 
  *  the beginning of the cutoff times
  */
  if( rd_leap_file == 1 )
    {
    leap = -1.;
    if( u58 < u58_cutoff[lasti] )
      {
      if( lasti == 0 )
        leap = 0.;
      else
        lasti = 0;
      }
    if( leap == 0 ) return leap;
   /*
    *  look in the u58 time cutoff list for te proper bounding time and get 
    *  the leap sec for that time
    */
    for( icut = lasti; icut < leap_info.ncut; icut++ )
      {
      if( ( icut == leap_info.ncut - 1 ) || ( u58 < u58_cutoff[ icut + 1 ] ) )
        {
        lasti = icut;
        return leap_lst[lasti];
        }
      }
    }
 /*
  * an error happened if it came here.  A stop gap estimate could be done,
  * but just error out
  */
  printf( 
    "%s,%d: An error occured in getting the leap second\n", 
    __FILE__, __LINE__ );
 /*
  printf( "     Resorting to limited internal information\n" );
  leap = 34.;
  if( u58 >= 1719792035000000 ) leap = 35.;
  return leap;
  */
  exit(1);
  }

int rd_leapsec( leap_str *leap_info )
/* ------------------------------------------------------------------------
   rd_leapsec

   purpose: read contents of the leapsec.dat file and put in a structure

   returns 0 if OK

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      leap_str *        leap_info        O      structure of the leap second 
                                                information

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC, 16 Feb 2012    Original development

 ------------------------------------------------------------------------*/
  {
  char *filedir, filename[FILENAME_MAX], line[101], mon_str[5];
  char *mon_lst[] = { "JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC" };
  FILE *fp;
  int ncut, icut, i;
  int16_t year, mon, day;
  double leap;

  if( ( filedir = getenv( "OCVARROOT" ) ) == NULL )
    {
    printf( "-E- %s: OCVARROOT env variable undefined.\n", __FILE__ );
    exit(1);
    }
  strcpy( filename, filedir );
  strcat( filename, "/modis/leapsec.dat" );
  if(want_verbose) {
        printf( "%s: Reading leap second information from file %s\n", 
        __FILE__, filename );
  }
  if ( (fp = fopen( filename, "r" ) ) == NULL )
    {
    fprintf( stderr, 
      "-E- %s line %d: unable to open leapsec file: %s for reading\n",
      __FILE__, __LINE__, filename );
    exit(1);
    }
 /*
  *  first, count the lines of good data
  */
  ncut = 0;
  while( fgets(line, 100, fp ) )
    {
    if( sscanf( line, "%4hd %3s %2hd %*s %*s %*s %11lf",
    &year, mon_str, &day, &leap ) == 4 ) 
      ncut++;
    }
  rewind( fp );
 /*
  *  allocate storage for the data
  */
  if( ( leap_info->year = ( int16_t *) malloc( ncut * sizeof( int16_t ) ) )
    == NULL )
    {
    printf( "%s,%d: unable to allocate year storage\n", __FILE__, __LINE__ );
    exit(1);
    }
  if( ( leap_info->mon = ( int16_t *) malloc( ncut * sizeof( int16_t ) ) )
    == NULL )
    {
    printf( "%s,%d: unable to allocate mon storage\n", __FILE__, __LINE__ );
    exit(1);
    }
  if( ( leap_info->day = ( int16_t *) malloc( ncut * sizeof( int16_t ) ) )
    == NULL )
    {
    printf( "%s,%d: unable to allocate day storage\n", __FILE__, __LINE__ );
    exit(1);
    }
  if( ( leap_info->leap = ( double *) malloc( ncut * sizeof( double ) ) )
    == NULL )
    {
    printf( "%s,%d: unable to allocate leap  storage\n", __FILE__, __LINE__ );
    exit(1);
    }
 /*
  *  read in the info, convert string month and put in the structure
  */
  leap_info->ncut = ncut;
  icut = 0;
  while( fgets(line, 100, fp ) )
    {
    if( sscanf( line, "%4hd %3s %2hd %*s %*s %*s %11lf",
    &year, mon_str, &day, &leap ) == 4 )
      {
      mon = -1;
      for( i = 0; i < 12; i++ )
        {
        if( strcmp( mon_str, mon_lst[i] ) == 0 )
          {
          mon = i;
          leap_info->mon[icut] = mon;
          break;
          }
        }
      if( mon == -1 )
        {
        printf( "%s,%d: month string format in leapsec.dat file is incorrect\n",
          __FILE__, __LINE__ );
        exit(1);
        }
      leap_info->year[icut] = year;
      leap_info->mon[icut] = mon;
      leap_info->day[icut] = day;
      leap_info->leap[icut] = leap;
      icut++;
      }
    }
  fclose( fp );
  return 0;
  }
