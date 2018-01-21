#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include "wgrib2.h"
#include "CodeTable4_4.h"

#define  FEB29   (31+29)
static int monthjday[12] = {
        0,31,59,90,120,151,181,212,243,273,304,334};

static int leap(int year) {
	if (year % 4 != 0) return 0;
	if (year % 100 != 0) return 1;
	return (year % 400 == 0);
}


/*
    add_time:  adds a positive value to a time code
    public domain 2006: wesley ebisuzaki
    1/2007 cleanup M. Schwarb, W. Ebisuzaki
 */

int add_time(int *year, int *month, int *day, int *hour, int *minute, int *second, unsigned int dtime, 
		int unit) {

    int y, m, d, h, mm, s, jday, i;

    y = *year;
    m = *month;
    d = *day;
    h = *hour;
    mm = *minute;
    s = *second;

    if (unit == YEAR) {
	*year = y + dtime;
	return 0;
    }
    if (unit == DECADE) {
	*year =  y + (10 * dtime);
	return 0;
    }
    if (unit == CENTURY) {
	*year = y + (100 * dtime);
	return 0;
    }
    if (unit == NORMAL) {
	*year = y + (30 * dtime);
	return 0;
    }
    if (unit == MONTH) {
        if (dtime < 0) {
           i = (-dtime) / 12 + 1;
           y -= i;
           dtime += (i * 12);
        }
	dtime += (m - 1);
	*year = y + (dtime / 12);
	*month = 1 + (dtime % 12);
	return 0;
    }

    if (unit == SECOND) {
	s += dtime;
	dtime = floor(s/60.0);
	*second = s - dtime*60;
	if (dtime == 0) return 0;
	unit = MINUTE;
    }
    if (unit == MINUTE) {
	mm += dtime;
	dtime = floor(mm/60.0);
	*minute = mm - dtime*60;
	if (dtime == 0) return 0;
	unit = HOUR;
    }
    if (unit == HOUR3) {
	dtime *= 3;
	unit = HOUR;
    }
    if (unit == HOUR6) {
	dtime *= 6;
	unit = HOUR;
    }
    if (unit == HOUR12) {
	dtime *= 12;
	unit = HOUR;
    }
    if (unit == HOUR) {
	h += dtime;
	dtime = floor(h/24.0);
	*hour = h - dtime*24;
	if (dtime == 0) return 0;
	unit = DAY;
    }

    /* this is the hard part */

    if (unit == DAY) {
	/* set m and day to Jan 0, and readjust dtime */
	jday = d + monthjday[m-1];
	if (leap(y) && m > 2) jday++;
        dtime += jday;

        while (dtime < 1) {
            y--;
	    dtime += 365 + leap(y);
        }

	/* one year chunks */
	while (dtime > 365 + leap(y)) {
	    dtime -= (365 + leap(y));
	    y++;
	}

	/* calculate the month and day */

	if (leap(y) && dtime == FEB29) {
	    m = 2;
	    d = 29;
	}
	else {
	    if (leap(y) && dtime > FEB29) dtime--;
	    for (i = 11; monthjday[i] >= dtime; --i);
	    m = i + 1;
	    d = dtime - monthjday[i];
	}
	*year = y;
	*month = m;
	*day = d;
	return 0;
   }
   fprintf(stderr,"add_time: undefined time unit %d\n", unit);
   return 1;
}

