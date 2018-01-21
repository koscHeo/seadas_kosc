/*
 * dataday.c
 *
 *  Created on: Jan 22, 2015
 *      Author: rhealy
 *      Functions to calculate dataday based on equatorial crossing
 *      and swath vertices
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
//#include "cdata.h"
#include "dataday.h"
#define PI  3.1415926535897932384626433832795029L

void get_datadays(
        /*
         *:geospatial_lon_max = 152.5445f ;
          :geospatial_lon_min = 148.602f ;
         *
         */
time_t   starttime,  /* (in)  swath start time
                     (seconds since 1-Jan-1970 00:00:00 GMT) */
float    eqxhour,    /* (in)  sensor's nominal local equator crossing time
                       (expressed in hours) */
int      isnight,    /* (in)  =0 for daytime outlines, =1 for nighttime */
int      dateline,   /* (in)  indicates whether swath crosses dateline */
float    west,       /* (in)  westernmost longitude of swath */
float    east,       /* (in)  easternmost longitude of swath */
int32_t   *dataday0,  /* (out) dataday of swath (days since 1-Jan-1970) */
int32_t   *dataday1  /* (out) 2nd dataday for datline-spanning swaths */
){
  time_t    reftime;
  int       refday;
  float     refhour;


  reftime = (time_t)(starttime + (12 - (double)eqxhour)*3600);
  refday  = reftime/86400;
  refhour = (reftime%86400)/3600.0;

  if(dateline == DATELINE_NOT_CROSSED){
    *dataday1 = *dataday0 = refday;
  }
  else if(refhour < 6){
              *dataday0 = refday - 1;
              *dataday1 = refday;
   }
  else if(refhour > 18){
              *dataday0 = refday;
              *dataday1 = refday + 1;
  }
  else if(dateline == DATELINE_NORTH_POLE || dateline == DATELINE_SOUTH_POLE){
      *dataday1 = *dataday0 = refday;
 }
  else{
    float   wod;    /* number of degrees west of dateline */
    float   eod;    /* number of degrees east of dateline */

/*    if(isnight){
       The 0-degree meridian is the dateline.
      wod = fmod(360 - west, 360);
      eod = fmod(360 + east, 360);
    }
    else{*/
      /* The 180-degree meridian is the dateline. */
      wod = 180 - west;
      eod = east + 180;
   // }
    if(wod > eod){
      *dataday0 = refday - 1;
      *dataday1 = refday;
    }
    else{
      *dataday0 = refday;
      *dataday1 = refday + 1;
   }
  }
}
void get_coord_extrema(
int   isnight,    /* (in) =0 for daytime outlines, =1 for nighttime */
int   n,          /* (in) number of elements in outline arrays */
float *olat,      /* (in) outline latitudes  in clockwise order */
float *olon,      /* (in) outline longitudes in clockwise order */
float *north,     /* (out) northernmost latitude covered by swath */
float *south,     /* (out) southernmost latitude covered by swath */
float *west,      /* (out) westernmost longitude covered by swath */
float *east,      /* (out) easternmost longitude covered by swath */
int *dateline     /* (out) pole and 180-deg. meridian crossing info */
){

  float   maxlat, minlat, maxlon, minlon;
  int     p;
  int       w180, e180; /* westward and eastward dateline crossings */
  int       w000, e000; /* the same for the 0-degree dateline */
  int       first,s0done;
  int       irreg=0;    /* flag to signal irregular swath outline */

  /*
  Traverse the swath outline looking for 180-degree meridian crossings
  and keeping track of coordinate extrema.
  */
  minlat =  INVALID_COORD;
  minlon =  INVALID_COORD;
  maxlat = -INVALID_COORD;
  maxlon = -INVALID_COORD;
  w180 = e180 = w000 = e000 = 0;
  for(p=1; p<=n; p++){
    float deltalon;
    int   pmodn;

    pmodn = p % n;

    /*
    Trap cases where the outline just touches the 180-degree meridian.
    I don't want those recorded as a dateline crossing unless the outline
    continues on the other side of the meridian.
    */
    if     (olon[pmodn] == -180 && olon[p-1] > 0) olon[pmodn] =  180;
    else if(olon[pmodn] ==  180 && olon[p-1] < 0) olon[pmodn] = -180;

    deltalon = olon[pmodn] - olon[p-1] + 360*(e180 - w180);
    if(deltalon < -180){
      /* 180-degree meridian crossed heading eastward */
      e180++;
    }
    else if(deltalon > 180){
      /* 180-degree meridian crossed heading westward */
      w180++;
    }
    else if( olon[pmodn] >= 0 && olon[p-1] - 360*(e180 - w180) < 0 ){
      /* 0-degree meridian crossed heading eastward */
      e000++;
    }
    else if( olon[pmodn] < 0 && olon[p-1] - 360*(e180 - w180) >= 0 ){
      /* 0-degree meridian crossed heading westward */
      w000++;
    }
    if(abs(e180-w180) > 1){
     // fprintf(stderr,"-W- %s line %d: ",__FILE__,__LINE__);
      //fprintf(stderr,"180-degree meridian crossed twice\n ");
      irreg = 1;
    }
    if(abs(e000-w000) > 1){
     // fprintf(stderr,"-W- %s line %d: ",__FILE__,__LINE__);
     // fprintf(stderr,"0-degree meridian crossed twice \n");
      irreg = 1;
    }
    olon[pmodn] += 360*(e180 - w180);

    if(olon[pmodn] > maxlon) maxlon = olon[pmodn];
    if(olon[pmodn] < minlon) minlon = olon[pmodn];
    if(olat[pmodn] > maxlat) maxlat = olat[pmodn];
    if(olat[pmodn] < minlat) minlat = olat[pmodn];
  }

  /*
  Set the dateline argument. For daytime outlines, the dateline
  is the 180-degree meridian. For nighttime outlines, the dateline
  is the prime meridian.
  */
  if(isnight){
    if((w000 + e000) == 0){
      /* The swath did not cross the  0-degree meridian. */
      *dateline = DATELINE_NOT_CROSSED;
    }
    else if((w000 + e000) & 1){
      /*
      An odd number of crossings of the 0-degree meridian means that
      one of the poles is included in the swath.  I am assuming that a
      swath never includes both poles.
      */
      if(90 - maxlat < minlat + 90){
        /* It's the north pole. */
        *dateline = DATELINE_NORTH_POLE;
        maxlat = 90;
      }
      else{
        /* It's the south pole. */
        *dateline = DATELINE_SOUTH_POLE;
        minlat = -90;
      }
      minlon = -180;
      maxlon =  180;
    }
    else{
      /* The swath crosses the 0-degree meridian but doesn't include a pole. */
      *dateline = DATELINE_CROSSED;
    }
  }
  else{
    /* Set the dateline argument. */
    if((w180 + e180) == 0){
      /* The swath did not cross the 180-degree meridian. */
      *dateline = DATELINE_NOT_CROSSED;
    }
    else if((w180 + e180) & 1){
      /*
      An odd number of crossings of the 180-degree meridian means that
      one of the poles is included in the swath.  Since this code is
      for sun-synchronous daytime sensors, I am assuming that a swath
      never includes both poles.
      */
      if(90 - maxlat < minlat + 90){
        /* It's the north pole. */
        *dateline = DATELINE_NORTH_POLE;
        maxlat = 90;
      }
      else{
        /* It's the south pole. */
        *dateline = DATELINE_SOUTH_POLE;
        minlat = -90;
      }
      minlon = -180;
      maxlon =  180;
    }
    else{
      /* The swath crosses the 180-deg. meridian but doesn't include a pole. */
      *dateline = DATELINE_CROSSED;
    }
  }
  if(irreg) *dateline = DATELINE_IRREGULAR;

  /* Put the longitudes in the range [-180,180]. */
  for(p=0; p<n; p++){
    while(olon[p] < -180) olon[p] += 360;
    while(olon[p] >= 180) olon[p] -= 360;
  }
  while(minlon < -180) minlon += 360;
  while(maxlon >  180) maxlon -= 360;

  /* Copy the values to the calling function. */
  *north = maxlat;
  *south = minlat;
  *west  = minlon;
  *east  = maxlon;
}


int daynight_outlines(
int32_t  *year,
int32_t  *dayOfYear,
int32_t  *msecondOfDay,
int32_t   wid,        /* (in) width  of lat/lon arrays */
int32_t   hgt,        /* (in) height of lat/lon arrays */
float **lat,      /* (in) pixel latitudes  */
float **lon,      /* (in) pixel longitudes */
int32_t   n[2],       /* (out) number of vertices per outline */
float *olat[2],   /* (out) outline latitudes  (olat[2][n[i]]) */
float *olon[2],   /* (out) outline longitudes (olon[2][n[i]]) */
int8_t    **dorn     /* (out) 2D array, dorn[hgt][wid], will hold the    */
            /*       day (=0) or night (=1) state of each pixel */
            /*       in the input arrays.                       */
){

  struct coord  *list[2] = {NULL,NULL};
  int s,fgs,lgs,night_tally=0;
  int8_t isnight = -1;
  double secondOfDay;
  int idayOfYear, iyear;
  int status = 0;
  
  // search for first good scan
  for(s=0; s<hgt; s++){
    if ((lat[s][0] <= 90 && lat[s][0] >= -90) &&
        (lat[s][wid/2] <= 90 && lat[s][wid/2] >= -90) &&
        (lat[s][wid] <= 90 && lat[s][wid] >= -90)){
        break;
    }
  }

  if(s == hgt){
    fprintf(stderr,"-W- %s line %d: ",__FILE__,__LINE__);
    fprintf(stderr,"No valid latitudes found\n ");
    status = 1;
    return status;
  }
  else{
    fgs = s;
  }
  // search for last good scan
  for(s=hgt-1; s>fgs; s--){
    if ((lat[s][0] <= 90 && lat[s][0] >= -90) &&
        (lat[s][wid/2] <= 90 && lat[s][wid/2] >= -90) &&
        (lat[s][wid] <= 90 && lat[s][wid] >= -90)){
        break;
    }
  }
  
  lgs = s;
  if(lgs == fgs){
    fprintf(stderr,"-W- %s line %d: ",__FILE__,__LINE__);
    fprintf(stderr,"Insufficient valid latitudes found\n ");
    status = 1;
    return status;
  } 
  /*
  Save some time by first checking the four corners of
  the scene to see if it is day, night, or mixed.
  */
  for(s=fgs; s<lgs; s+=lgs-fgs-1){

    float   sv[3], sundist;
    int     p;

    /* Estimate the time of this scan line. */
   /* utime2yds(
      est_scan_time(stime,lat[fgs][wid/2],lon[fgs][wid/2],
                          lat[ s ][wid/2],lon[ s ][wid/2] ),
      &year,&dayOfYear,&secondOfDay
    );*/

    /* Get the solar position vector. */
    secondOfDay = msecondOfDay[s]/1000.;
    idayOfYear = dayOfYear[s];
    iyear = year[s];

    l_sun_(&iyear, &idayOfYear, &secondOfDay, sv, &sundist);


    for(p=0; p<wid; p+=wid-1){

      float pv[3];
      double    coslat;

      /* Compute the position vector of this pixel. */
      coslat = cos(lat[s][p] * PI/180);
      pv[0] = coslat * cos(lon[s][p] * PI/180);
      pv[1] = coslat * sin(lon[s][p] * PI/180);
      pv[2] =          sin(lat[s][p] * PI/180);

      if(pv[0]*sv[0] + pv[1]*sv[1] + pv[2]*sv[2] <= 0){
        /*
        If the dot product of the sun vector and the pixel
        vector is less than or equal to zero (i.e. the angle between
        the two vectors is greater than or equal to 90 degrees) then the
        pixel is on the night side of the planet.
        */
        night_tally++;
      }
      /*
      I don't expect any 1-pixel-wide scenes, but, just in case,
      the following line will prevent infinite loops.
      */
      if(wid < 2) break;
    }
    /*
    The following line prevents single-scan-line segments
    from causing an infinite loop.
    */
    if(lgs-fgs < 2) break;
  }

  if(night_tally == 0 || night_tally == 4){
    int     dn, i, j;

    dn = night_tally/4;     /* dn==0 : daytime    dn==1 : nighttime */

//    memset((*dorn)[0],dn,wid*(lgs-fgs));

    /* Compute the number of vertices in the scene outline. */
//    n[dn] = 2*(wid + hgt - fgs - 2);
    n[dn] = 2*(wid + lgs - fgs - 2);

    /* Allocate space for that many latitudes and longitudes. */
    MALLOC(olat[dn],float,n[dn]);
    MALLOC(olon[dn],float,n[dn]);
    float tlat, tlon;
    /* ... and fill them in. */
    for (j = 0, i = 0; i < wid - 1; i++) {
        tlat = lat[fgs][i];
        tlon = lon[fgs][i];
        if ((tlat >= -90 && tlat <= 90) &&
                (tlon >= -180 && tlon <= 180)) {
            olat[dn][j] = tlat;
            olon[dn][j] = tlon;
            j++;
        }

    }
    for (i = fgs; i < lgs - 1; i++) {
        tlat = lat[i][wid - 1];
        tlon = lon[i][wid - 1];
        if ((tlat >= -90 && tlat <= 90) &&
                (tlon >= -180 && tlon <= 180)) {
            olat[dn][j] = tlat;
            olon[dn][j] = tlon;
            j++;
        }

    }
    for (i = wid - 1; i > 0; i--) {
        tlat = lat[lgs - 1][i];
        tlon = lon[lgs - 1][i];
        if ((tlat >= -90 && tlat <= 90) &&
                (tlon >= -180 && tlon <= 180)) {
            olat[dn][j] = tlat;
            olon[dn][j] = tlon;
            j++;
        }

    }
    for (i = lgs - 1; i > fgs; i--) {
        tlat = lat[i][0];
        tlon = lon[i][0];
        if ((tlat >= -90 && tlat <= 90) &&
                (tlon >= -180 && tlon <= 180)) {
            olat[dn][j] = tlat;
            olon[dn][j] = tlon;
            j++;
        }
    }
    // set the real number of valid indices
    n[dn] = j;
    /* The other outline is empty in a non-mixed scene. */
    n[!dn] = 0;
    olat[!dn] = NULL;
    olon[!dn] = NULL;

  }
  else{
    /*
    This is a mixed scene, so I need to determine
    the day/night status of every pixel.
    */

    for(s=fgs; s<lgs; s++){

      float sv[3], sundist, prevlat, prevlon;
      int       p;

      /* Get the solar position vector. */
      secondOfDay = msecondOfDay[s]/1000.;
      idayOfYear = dayOfYear[s];
      iyear = year[s];
      l_sun_(&iyear, &idayOfYear, &secondOfDay, sv, &sundist);

      for(p=0; p<wid; p++){

        float       pv[3];
        double      coslat;
        int     daynightchange, i;

        if (lat[s][p] < -90 || lat[s][p] > 90 || lon[s][p] < -180 || lon[s][p] > 180 ) continue;
        /* Scan the pixels boustrophedonically. */
        if(s & 1){      /* If this an odd scan line ...  */
          i = wid - 1 - p;  /* then scan from right to left; */
        }
        else{           /* otherwise, ...           */
          i = p;        /* scan from left to right. */
        }

        /* Compute the position vector of this pixel. */
        coslat = cos(lat[s][i] * PI/180);
        pv[0] = coslat * cos(lon[s][i] * PI/180);
        pv[1] = coslat * sin(lon[s][i] * PI/180);
        pv[2] =          sin(lat[s][i] * PI/180);

        if(pv[0]*sv[0] + pv[1]*sv[1] + pv[2]*sv[2] > 0){
          /*
          If the dot product of the sun vector and the pixel
          vector is greater than zero (i.e. the angle between
          the two vectors is less than 90 degrees) then the
          pixel is on the day side of the planet.
          */
          daynightchange = isnight != 0;
          isnight = 0;
        }
        else{
          daynightchange = isnight != 1;
          isnight = 1;
        }

        /*
        If the pixel is on the edge of the scene or on the edge of a
        day/night transition, then add its coordinates to the appropriate
        end (head or tail) of the appropriate outline (day or night).
        */
        if(s==0 || i==wid-1){
          push(&list[isnight],tail,lat[s][i],lon[s][i]);
        }
        else if(i==0){
          push(&list[isnight],head,lat[s][i],lon[s][i]);
        }
        else if(s==lgs-1){
          if(s & 1){    /* moving towards the left */
            push(&list[isnight],tail,lat[s][i],lon[s][i]);
          }
          else{     /* moving towards the right */
            push(&list[isnight],head,lat[s][i],lon[s][i]);
          }
        }
        else if(daynightchange){
          if(s & 1){      /* moving towards the left */
            push(&list[ isnight],tail,lat[s][i],lon[s][i]);
            push(&list[!isnight],head,prevlat  ,prevlon  );
          }
          else{           /* moving towards the right */
            push(&list[ isnight],head,lat[s][i],lon[s][i]);
            push(&list[!isnight],tail,prevlat  ,prevlon  );
          }
        }
        dorn[s][i] = isnight;
        prevlat = lat[s][i];
        prevlon = lon[s][i];
      }
    }
    for(s=0;s<2;s++){
      int       count=0;
      struct coord  *p;
      float     *lats,*lons,*latp,*lonp,prevlat=-999,prevlon=-999;

      /* Count the number of vertices in this outline. */
      for(p=list[s]; p!=NULL; p=p->next) count++;

      if(count > 0){
        /* Allocate space for that many latitudes and longitudes. */
        MALLOC(lats,float,count);
        MALLOC(lons,float,count);

        /*
        Copy the latitudes and longitudes from the linked list to
        the arrays omitting any duplicate coordinates in the process.
        */
        latp = lats;
        lonp = lons;
        for(p=list[s]; p!=NULL; p=p->next){
          if(p->lat != prevlat && p->lon != prevlon){
            prevlat = *latp++ = p->lat;
            prevlon = *lonp++ = p->lon;
          }
          else{
            count--;
          }
        }

        if(count <= 0){
          fprintf(stderr,"-E- %s line %d: ",__FILE__,__LINE__);
          fprintf(stderr,"Program logic error. Likely some bad navigation got in - fix the code!.\n");
          status = 2;
        }

        /* Deallocate the memory in the linked list. */
        if(list[s] != NULL){
          struct coord  *pn;
          for(p=list[s]; p!=NULL; p=pn){
            pn = p->next;
            free(p);
            p = NULL;
          }
        }

        /* Return the output arrays to the calling program. */
        n[s] = count;
        olat[s] = lats;
        olon[s] = lons;
      }
      else{
        n[s] = 0;
        olat[s] = NULL;
        olon[s] = NULL;
      }
    }
  }
  return status;
}
void push(
struct coord    **vrtx,
enum hort   hort,
float     lat,
float     lon
){
  struct coord  *newvrtx;

  MALLOC(newvrtx,struct coord,1);

  newvrtx->lat = lat;
  newvrtx->lon = lon;

  if     (hort == head){
    newvrtx->next = *vrtx;
    *vrtx = newvrtx;
  }
  else if(hort == tail){
    struct coord    *p;

    newvrtx->next = NULL;

    if(*vrtx == NULL){
      *vrtx = newvrtx;
    }
    else{
      /*
      Go to the end of the list. Note that
      the final semicolon is intentional.
      */
      for(p=*vrtx; p->next!=NULL; p=p->next);

      p->next = newvrtx;
    }
  }
}

