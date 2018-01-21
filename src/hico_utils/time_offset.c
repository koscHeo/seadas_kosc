/*
*******************************************************************************
* time_offset.c
*
* Read start date/time, offset date/time and output resulting date/time
*
* time_offset start_year start_month start_day start_hour start_minute start_second
*     offset_year offset_month offset_day offset_hour offset_minute offset_second
*
* output:
*     result_year result_month result_day result_hour result_minute result_second
*
* WARNING: NOT SETUP FOR YEAR, MONTH AND DAY OFFSET VALUES AT THIS TIME!
*          Just short on time to finish this!
*
* Originally written by Karen Patterson
* November 2009
*******************************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void usage_call(void);

int main (int argc, char *argv[])
{
  float syear, smonth, sday, shour, sminute, ssecond; /* start time */
  float oyear, omonth, oday, ohour, ominute, osecond; /* offset time */
  int ryear, rmonth, rday, rhour, rminute, rsecond; /* result time */
  int jdate, jyear, jmonth, jday; /* temporary values for gregorian date calculations */
  int i,j,k,l,n;

  /* only run if there is the proper number of arguments */
  if (argc == 13)
    {
    syear = atof(argv[1]);
    smonth = atof(argv[2]);
    sday = atof(argv[3]);
    shour = atof(argv[4]);
    sminute = atof(argv[5]);
    ssecond = atof(argv[6]);
    oyear = atof(argv[7]);
    omonth = atof(argv[8]);
    oday = atof(argv[9]);
    ohour = atof(argv[10]);
    ominute = atof(argv[11]);
    osecond = atof(argv[12]);
    }
  else
    usage_call();

  /* first cut at result time */
  ryear = syear + oyear;
  rmonth = smonth + omonth;
  rday = sday + oday;
  rhour = shour + ohour;
  rminute = sminute + ominute;
  rsecond = ssecond + osecond;

  /* deal with overflows */
  /* time first - it's easy */
  while(rsecond >= 60)
    {
    rsecond = rsecond - 60;
    rminute = rminute + 1;
    }
  while(rsecond < 0)
    {
    rsecond = rsecond + 60;
    rminute = rminute - 1;
    }
  while(rminute >= 60)
    {
    rminute = rminute - 60;
    rhour = rhour + 1;
    }
  while(rminute < 0)
    {
    rminute = rminute + 60;
    rhour = rhour - 1;
    }
  while(rhour >= 24)
    {
    rhour = rhour - 24;
    rday = rday + 1;
    }
  while(rhour < 0)
    {
    rhour = rhour + 24;
    rday = rday - 1;
    }
  /* date is tricky */
  /* convert to Julian Day */
  jyear=ryear;
  jmonth=rmonth;
  jday=rday;
  jyear+=8000;
  if(jmonth<3)
    {
    jyear--;
    jmonth+=12;
    }
  jdate=(jyear*365) + (jyear/4) - (jyear/100) + (jyear/400) - 1200820 + (jmonth*153+3)/5 - 92 + jday - 1;

  /* convert back to Gregorian Date */
  /* from http://aa.usno.navy.mil/faq/docs/JD_Formula.php */
  l=jdate+68569;
  n=4*l/146097;
  l=l-(146097*n+3)/4;
  i=4000*(l+1)/1461001;
  l=l-1461*i/4+31;
  j=80*l/2447;
  k=l-2447*j/80;
  l=j/11;
  j=j+2-12*l;
  i=100*(n-49)+i+l;

  ryear=i;
  rmonth=j;
  rday=k;

  /* print the resulting date/time */
  fprintf(stdout, "%d %d %d %d %d %d\n",ryear,rmonth,rday,rhour,rminute,rsecond);
  return(0);
}

/* =============================================================
 *    usage_call()
 * ============================================================= */
void usage_call(void)
{
  printf("\n");
  printf("Usage: time_offset <year> <month> <day> <hour> <minute> <second>\n");
  printf("           <oyear> <omonth> <oday> <ohour> <ominute> <osecond>\n");
  printf("\n");
  printf("<year> = starting year\n");
  printf("<month> = starting month\n");
  printf("<day> = starting day of month\n");
  printf("<hour> = starting hour\n");
  printf("<minute> = starting minute\n");
  printf("<second> = starting second\n");
  printf("<oyear> = offset year\n");
  printf("<omonth> = offset month\n");
  printf("<oday> = offset day of month\n");
  printf("<ohour> = offset hour\n");
  printf("<ominute> = offset minute\n");
  printf("<osecond> = offset second\n");
  printf("\n ONLY USE <ohour> <ominute> and <osecond> in less than a day quantities\n");
  printf(" at this time!  Ran out of time to finish multi-day/month/year calculations!\n");
  exit (1);
}
