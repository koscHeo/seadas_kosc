#include <stdlib.h>
#include <stdio.h>

#include "filehandle.h"
#include "l1_struc.h"
#include "l12_proto.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include <proj_api.h>

#include "goci.h"
#include "l1_goci.h"
#include <libnav.h>

static goci_l1b_t *goci_l1b = NULL;

static int year, month, day, hour, minute, second, doy, base_msec;
static uint32_t *buf;
static double *lat, *lon;

int
openl1_goci( filehandle *file )
/*
 *  openl1_goci  prepares for reading the GOCI L1B file
 */
{
    int   status;
    int   i;
    char tim_str[80];

    status = goci_l1b_open( file->name, &goci_l1b );
    if ( status ) {
        fprintf(stderr,"-E- %s line %d: unable to open %s\n",
            __FILE__,__LINE__,file->name);
        return 1;
    }
    if(want_verbose)
        printf("GOCI Level-1B %s\n", file->name );

    file->npix   = goci_l1b->npixels;
    file->nscan  = goci_l1b->nscans;
    file->nbands = goci_l1b->nbands;
    strcpy(file->spatialResolution,"500 m");

    // get time.  If slot navigation is available, get the start time 
    // for future modification based on the slot.  Else, use the center time

    if( goci_l1b->slot_nav_avail )
      strcpy( tim_str, "Scene Start time" );
    else
      strcpy( tim_str, "Scene center time" );

    goci_l1b_get_date( goci_l1b, tim_str, &year, &month, &day);
    goci_l1b_get_time( goci_l1b, tim_str, &hour, &minute, &second );

    ymdhms2ydmsec(year, month, day, hour, minute, second,
                  &year, &doy, &base_msec);

    if(want_verbose) {
        printf("GOCI %s: %4d-%02d-%02d %03d %02d:%02d:%02d %d\n",
              tim_str, year, month, day, doy, hour, minute, second, base_msec);
        printf("GOCI file has %d bands, %d samples, %d lines\n",
              file->nbands, file->npix, file->nscan );
    }

    buf = (uint32_t *) malloc( file->npix*sizeof(uint32_t) );
    if ( buf == NULL )
        goto memory_error;

    lat = (double *)   malloc ( file->npix*sizeof(double) );
    lon = (double *)   malloc ( file->npix*sizeof(double) );
    if ( lat == NULL || lon == NULL )
        goto memory_error;

    if(goci_proj4_open(goci_l1b)) {
        fprintf(stderr,"-E- %s line %d: unable to init proj4 geolocation for %s\n",
            __FILE__,__LINE__,file->name);
        goto memory_release;
    }
    if(want_verbose)
        printf("GOCI using internal navigation\n");

    return 0;

memory_error:
    fprintf(stderr, "%s: %d - memory error\n", __FILE__, __LINE__ );
memory_release:
    if ( buf ) free(buf);
    if ( lat ) free(lat);
    if ( lon ) free(lon);

    return 1;
}

int
readl1_goci( filehandle *file, int recnum, l1str *l1rec, int lonlat)
/*
 *  fill standard record with L1B line of data
 */
{
    int   status, min_msec = 86401 * 1000, max_msec = -1;

    herr_t    herr;

    double elev, solz, sola;
    float rel_sec;

    int    npix, nbands, ip, ib, ipb, msec;

    npix   = file->npix;
    nbands = file->nbands;
    msec = base_msec;

    //  set information about data

    l1rec->npix = file->npix;
    l1rec->sensorID = file->sensorID;

    *(l1rec->year) = year;
    *(l1rec->day)  = doy;
    *(l1rec->msec) = msec;

    //  get lat-lon
    for (ip=0; ip<npix; ip++) {
        lon[ip] = ip * goci_l1b->deltaX + goci_l1b->startX;
        lat[ip] = recnum * goci_l1b->deltaY + goci_l1b->startY;
    }
    goci_proj4_convert(goci_l1b, npix, lon, lat);

    for (ip=0; ip<npix; ip++) {

        l1rec->pixnum[ip] = ip;

        if ( isnan(lat[ip]) ) lat[ip] = -999.0;
        if ( isnan(lon[ip]) ) lon[ip] = -999.0;
        l1rec->lat[ip]    = lat[ip];
        l1rec->lon[ip]    = lon[ip];

        if (l1rec->lon[ip] < -181.0 || l1rec->lon[ip] > 181.0 ||
            l1rec->lat[ip] <  -91.0 || l1rec->lat[ip] >  91.0 )
          l1rec->navfail[ip] = 1;

        if(!lonlat) {
           //  get the msec based on goci nav availability
            if( goci_l1b->slot_nav_avail ) {
                unsigned char slot;
                slot = goci_slot_time( ip, recnum, goci_l1b, &rel_sec );
                l1rec->slot[ip] = slot + 1;
                msec = base_msec + rel_sec * 1000;
                if( msec > max_msec ) max_msec = msec;
                if( msec < min_msec ) min_msec = msec;
            } else
              msec = base_msec;
            float gmt = (float)msec / (1000.0 * 3600.0);
            sunangs_(&year, &doy, &gmt, l1rec->lon+ip, l1rec->lat+ip,
                    l1rec->solz+ip, l1rec->sola+ip);
            *(l1rec->sola+ip) = ( *(l1rec->sola+ip) > 180. ) ? 
              *(l1rec->sola+ip) - 360. : *(l1rec->sola+ip);

            get_zenaz(goci_l1b->sat_pos, l1rec->lon[ip], l1rec->lat[ip],
                    l1rec->senz+ip, l1rec->sena+ip);
            *(l1rec->sena+ip) = ( *(l1rec->sena+ip) > 180. ) ?
              *(l1rec->sena+ip) - 360. :  *(l1rec->sena+ip);
        }
    if( goci_l1b->slot_nav_avail )
      *( l1rec->msec ) = ( min_msec + max_msec ) / 2;
    }

    if(lonlat) {
        return(LIFE_IS_GOOD);
    }


    for ( ib = 0; ib < nbands; ib++ ) {

        status = goci_l1b_get_band( goci_l1b, ib, recnum, buf );
        if ( status ) {
            fprintf(stderr, "-E- %s line %d: Failed to read Lt, band %d, recnum %d\n",
                    __FILE__,__LINE__, ib, recnum );
            return 1;
        }

        for (ip=0; ip<npix; ip++) {
            ipb = ip*nbands+ib;
            if ( buf[ip] == 0 ) {
                l1rec->slot[ip]    = 0;
            }
            l1rec->Lt[ipb] = buf[ip]*1.0E-07;
        }
    }

    return 0;
}

int
closel1_goci( filehandle *file )
{
    int   status;

    if ( buf ) free(buf);
    if ( lat ) free(lat);
    if ( lon ) free(lon);
    if ( goci_l1b->slot_asg ) free( goci_l1b->slot_asg );
    if ( goci_l1b  ) free( goci_l1b );
    return 0;
}
