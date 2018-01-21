#include <stdio.h>
#include <stdint.h>

int jdate( int32_t julian, int32_t *year, int32_t *doy) {

  int i;
  int32_t month, day;

  int32_t ja, jb, jc, jd, je, jalpha;

  jalpha = (int32_t)(((julian - 1867216L) - 0.25) / 36524.25);
  ja = julian + 1 + jalpha - (int32_t)(0.25 * jalpha);

  jb = ja + 1524;
  jc = (int32_t)(6680.0 + ((jb - 2439870L) - 122.1) / 365.25);
  jd = 365 * jc + (int32_t)(0.25 * jc);
  je = (int32_t)((jb - jd) / 30.6001);

  day = jb - jd - (int32_t)(30.6001 * je);
  month = je - 1;
  if (month > 12)
    month = month - 12;
  *year = jc - 4715;

  if ( month > 2) *year = *year - 1;
  if ( *year <= 0) *year = *year - 1;

  int dmn[24] = {31,28,31,30,31,30,31,31,30,31,30,31,
		 31,29,31,30,31,30,31,31,30,31,30,31};

  int isleap = 0;
  if ( ((*year % 400) == 0) || 
       (((*year % 4) == 0) && ((*year % 100) != 0)) )
    isleap = 1;

  *doy = day;
  for (i=1; i<month; i++) *doy = *doy + dmn[isleap*12 + (i-1)];

  return 0;
}
