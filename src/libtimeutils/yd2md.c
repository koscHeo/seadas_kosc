#include <timeutils.h>

/* -------------------------------------------------------------- */
/* z16_z20() - converts 16-char zulu time to 20-char zulu time    */
/* -------------------------------------------------------------- */
void z16_z20( char *z_16, char *z_20) {
  char tmp_str[2048];
  char buf1[32], buf2[32];
  int16_t year, month, day, dom, hr, mn, sc;

  strcpy( tmp_str, z_16);
  strncpy( buf1, tmp_str, 4);
  buf1[4] = 0;
  year = atol(buf1);
  strncpy( buf1, &tmp_str[4], 3);
  buf1[3] = 0;
  day = atol(buf1);
  yd2md(year, day, &month, &dom);
  
  strncpy( buf2, &tmp_str[7], 2);
  buf2[2] = 0;
  hr = atol(buf2);
  strncpy( buf2, &tmp_str[9], 2);
  buf2[2] = 0;
  mn = atol(buf2);
  strncpy( buf2, &tmp_str[11], 2);
  buf2[2] = 0;
  sc = atol(buf2);
  sprintf(z_20, "%d-%02d-%02dT%02d:%02d:%02dZ", year, month, dom, hr, mn, sc);
}

/* -------------------------------------------------------------- */
/* yd2md() - converts year and day to month and day-of-month      */
/* -------------------------------------------------------------- */
void yd2md( int16_t year, int16_t doy, int16_t *month, int16_t *dom)
{
    double sec = 0.0;
    int16_t   yr;

    unix2ymds(yds2unix(year,doy,sec),&yr,month,dom,&sec);

    return;
}

