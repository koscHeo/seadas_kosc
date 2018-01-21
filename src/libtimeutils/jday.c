#include <stdio.h>
#include <stdint.h>

int32_t jday( int16_t i, int16_t j, int16_t k) {

  //    This function converts a calendar date to the corresponding Julian
  //    day starting at noon on the calendar date.  The algorithm used is
  //    from Van Flandern and Pulkkinen, Ap. J. Supplement Series 41, 
  //    November 1979, p. 400.

  //    Written by Frederick S. Patt, GSC, November 4, 1992

  return 367*i - 7*(i+(j+9)/12)/4 + 275*j/9 + k + 1721014;
}
