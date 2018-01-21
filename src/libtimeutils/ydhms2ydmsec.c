#include <timeutils.h>

void ymdhms2ydmsec(int yy,int mm,int dd,int hh, int mn, int sc,
                   int32_t *year, int32_t *day, int32_t *msec)
{
    static int startOfMonth[2][12] = 
        {{0,31,59,90,120,151,181,212,243,273,304,334},
  	 {0,31,60,91,121,152,182,213,244,274,305,335}};
    int leap;

    if(yy%400 == 0 || (yy%100 != 0 && yy%4 == 0))
        leap = 1;
    else
        leap = 0;
    

    *year = yy;
    *day  = startOfMonth[leap][mm-1] + dd;
    *msec = 1000 * (sc + 60 * (mn + 60 * hh));
}

