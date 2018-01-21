#include "hdf5utils.h"

// IDL> print,num_leap_sec('2008-12-31T23:59:59.000')
#define NUM_LEAP_SECONDS 14 // # of leapsecs from 1980 thru 12/31/2008

// IDL>  print,long(str2gps('2009-01-01T00:00:00.000'))
#define GPS_JAN0109   914803215 // GPS  time for 01/01/2009

// IDL> print,long(str2gps('2012-07-01T00:00:00.000'))
#define GPS_JUL0112  1025136016 // GPS  time for 07/01/2012

// IDL> print,long(str2gps('2015-07-01T00:00:00.000'))
#define GPS_JUL0115  1119744017 // GPS  time for 07/01/2015

// date -u -d "Jan 01 00:00:00 UTC 2009" +%s
#define UNIX_JAN0109 1230768000 // Unix time for 01/01/2009

// date -u -d "Jul 01 00:00:00 UTC 2012" +%s
#define UNIX_JUL0112 1341100800 // Unix time for 07/01/2012

// date -u -d "Jul 01 00:00:00 UTC 2015" +%s
#define UNIX_JUL0115 1435708800 // Unix time for 07/01/2015


// gpstai   is seconds from 01/06/1980 (with no leap seconds adjustment)
// unixtime is seconds from 01/01/1970 (with    leap seconds adjustment)

double gpstai2utc2000( double gpstai) {
  if ( gpstai < GPS_JAN0109)
    return gpstai - DIFF2000_JAN0680 - NUM_LEAP_SECONDS;
  else if ( gpstai < GPS_JUL0112) 
    return gpstai - DIFF2000_JAN0680 - NUM_LEAP_SECONDS - 1;
  else if ( gpstai < GPS_JUL0115) 
    return gpstai - DIFF2000_JAN0680 - NUM_LEAP_SECONDS - 2;
  else
    return gpstai - DIFF2000_JAN0680 - NUM_LEAP_SECONDS - 3;
}


double gpstai2unix( double gpstai) {
  if ( gpstai < GPS_JAN0109)
    return gpstai + DIFFJAN0680_1970 - NUM_LEAP_SECONDS;
  else if ( gpstai < GPS_JUL0112)
    return gpstai + DIFFJAN0680_1970 - NUM_LEAP_SECONDS - 1;
  else if ( gpstai < GPS_JUL0115)
    return gpstai + DIFFJAN0680_1970 - NUM_LEAP_SECONDS - 2;
  else
    return gpstai + DIFFJAN0680_1970 - NUM_LEAP_SECONDS - 3;
}


double unix2gpstai( double unixtime) {
  if ( unixtime < UNIX_JAN0109)
    return unixtime - DIFFJAN0680_1970 + NUM_LEAP_SECONDS;
  else if ( unixtime < UNIX_JUL0112)
    return unixtime - DIFFJAN0680_1970 + NUM_LEAP_SECONDS + 1;
  else if ( unixtime < UNIX_JUL0115)
    return unixtime - DIFFJAN0680_1970 + NUM_LEAP_SECONDS + 2;
  else
    return unixtime - DIFFJAN0680_1970 + NUM_LEAP_SECONDS + 3;
}


// TAI time (sec since 01/06/1980)
double get_tai( int32_t year, int32_t doy, double millisec) {

  struct tm tm0;    
  double tai;

  int dmn[24] = {31,28,31,30,31,30,31,31,30,31,30,31,
		 31,29,31,30,31,30,31,31,30,31,30,31};

  tm0.tm_year = year - 1900;

  int isleap = 0;
  if ( ((year % 400) == 0) || (((year % 4) == 0) && ((year % 100) != 0)) )
    isleap = 1;

  // Determine mn dy from doy
  for (size_t i=0; i<12; i++) {
    if ( (doy - dmn[isleap*12+i]) <= 0) {
      tm0.tm_mon = i;
      tm0.tm_mday = doy;
      tm0.tm_hour = 0;
      tm0.tm_min = 0;
      tm0.tm_sec = 0;
      tm0.tm_isdst = -1;
      break;
    } else {
      doy -= dmn[isleap*12+i];
    }
  }

  tai = unix2gpstai( mktime( &tm0) + (millisec/1000.) + tm0.tm_gmtoff);
    
  return tai;
}

double get_tai( char *orbString) {

  // orbString = YYYYDOYHHMMSSSSS

  string ydhmsf_str = orbString;

  istringstream istr;
  int32_t year, doy, itemp;

  istr.clear(); istr.str( ydhmsf_str.substr(0,4)); istr >> year;
  istr.clear(); istr.str( ydhmsf_str.substr(4,3)); istr >> doy;
  istr.clear(); istr.str( ydhmsf_str.substr(7,2)); istr >> itemp;
  int32_t millisec = itemp * 3600000;
  istr.clear(); istr.str( ydhmsf_str.substr(9,2)); istr >> itemp;
  millisec += itemp * 60000;
  istr.clear(); istr.str( ydhmsf_str.substr(11,5)); istr >> itemp;
  millisec += itemp;

  struct tm tm0;    
  double tai;
    
  int dmn[24] = {31,28,31,30,31,30,31,31,30,31,30,31,
		 31,29,31,30,31,30,31,31,30,31,30,31};

  tm0.tm_year = year - 1900;

  int isleap = 0;
  if ( ((year % 400) == 0) || (((year % 4) == 0) && ((year % 100) != 0)) )
    isleap = 1;
  
  // Determine mn dy from doy
  for (size_t i=0; i<12; i++) {
    if ( (doy - dmn[isleap*12+i]) <= 0) {
      tm0.tm_mon = i;
      tm0.tm_mday = doy;
      tm0.tm_hour = 0;
      tm0.tm_min = 0;
      tm0.tm_sec = 0;
      tm0.tm_isdst = -1;
      break;
    } else {
      doy -= dmn[isleap*12+i];
    }
  }

  tai = unix2gpstai( mktime( &tm0) + (millisec/1000.) + tm0.tm_gmtoff);

  return tai;
}
