#include <timeutils.h>

void date2ydmsec(char *date, int32_t *year, int32_t *day, int32_t *msec)
{
    int yy, mm, dd, hh, mn, sc, cc;

    sscanf (date, "%4d%2d%2d%2d%2d%2d%2d", 
            &yy, &mm, &dd, &hh, &mn, &sc, &cc);

    ymdhms2ydmsec(yy,mm,dd,hh,mn,sc,year,day,msec);
    *msec = *msec + (int32_t) cc*10;
}

void isodate2ydmsec(char *date, int32_t *year, int32_t *day, int32_t *msec)
{
    int yy, mm, dd, hh, mn, sc;
    sscanf (date, "%4d-%2d-%2dT%2d:%2d:%2dZ",
            &yy, &mm, &dd, &hh, &mn, &sc);

    ymdhms2ydmsec(yy,mm,dd,hh,mn,sc,year,day,msec);
}

