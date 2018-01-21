#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <timeutils.h>

int yds2tai( int16_t iyr, int16_t idy, double sec, double *tai93) {

  // Program to convert year, day, seconds to TAI93 (seconds since 1/1/1993)

  // Convert year and day to AJD
  int32_t jd = jday( iyr, 1, idy);

  if (getenv("OCVARROOT") == NULL) {
    printf("-E- OCVARROOT environment variable is not defined.\n");
    exit (-1);
  }

  // Compute uncorrected TAI93
  *tai93 = 86400*(jd - 2448989) + sec;

  // Get number of leapseconds 
  int nleap = leapseconds_since_1993( *tai93);

  // Compute corrected TAI93
  *tai93 = 86400*(jd - 2448989) + sec + nleap;

  return 0;
}
