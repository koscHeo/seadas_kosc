/*
 * l1_aviris.c
 *
 *  Created on: May 18, 2015
 *      Author: Rick Healy SAIC
 *              NASA-GSFC OBPG
 */
#include <stdlib.h>
#include <stdio.h>

#include "filehandle.h"
#include "l1_struc.h"
#include "l12_proto.h"
#include <proj_api.h>
#include "jplaeriallib.h"
#include "aviris.h"

#define SKIP -9999

static const int maxBands = 224;
static double *lat, *lon;

aviris_t* createPrivateData_av (int numBands) {

    aviris_t* data = (aviris_t*)calloc(1, sizeof(aviris_t));
    if(data == NULL) {
        fprintf(stderr,"-E- %s line %d: unable to allocate private data for aviris\n",
                __FILE__,__LINE__);
        exit(1);
    }

    data->wave = (double *) malloc(numBands*sizeof(double) );
    data->fwhm = (double *) malloc(numBands*sizeof(double) );
    if(data->wave==NULL || data->fwhm==NULL) {
        fprintf(stderr,"-E- %s line %d: unable to allocate scale/offset data for aviris\n",
                __FILE__,__LINE__);
        exit(1);
    }


    return data;
}

void freePrivateData_av(aviris_t* data) {
    int k;
    free(data->wave);
    free(data->fwhm);
    free(data->gain);
    for(k=0;k<data->nscan;k++) {
        free(data->sena[k]);
        free(data->senz[k]);
        free(data->sola[k]);
        free(data->solz[k]);
        free(data->utc[k]);
    }
    free(data->sena);
    free(data->senz);
    free(data->sola);
    free(data->solz);
    free(data->utc);
    gsl_interp_accel_free(data->spl_acc);

}

int  openl1_aviris( filehandle *file ) {

    int16_t buffer[10];
    FILE *ptr;
    char tag[itemSize];
    char *val0;
    char val[itemSize];
    char *inbasename, *infile;
    int i, j, k, status,pos;
    double *indata, *elev;
    float  *indataf;
    float rotation;
    int numBands, num;
    char* result;
    char line[itemSize];
    int count;
    int year, month, day, hour, minute, second, doy;

    int isec = 0;
    double sec;

    char projStr[1024];

    aviris_t* data = file->private_data = createPrivateData_av(maxBands);

    inbasename = getinbasename_av(file->name);
    pos = strlen(inbasename);
    if (pos<=0) {
        fprintf(stderr,"-E- %s line %d: Not a avalid AVIRIS file %s\n",
                __FILE__,__LINE__,file->name);
        return 1;
    }

    sscanf(inbasename, "f%2d%2d%2d", &year, &month, &day);

    if (year >= 92) year = year + 1900; else year = year + 2000;

    sec = 0;
    hour = 0;
    minute = 0;
    isec = (int)sec;

    data->month = month;
    data->day   = day;

    ymdhms2ydmsec(year, month, day, hour, minute, isec,
                       &data->year, &data->doy, &data->msec);

    sec -= isec;
    data->msec += sec * 1000;

    printf("Date of AVIRIS flight: Y-%d M-%d D-%d\n",year,month,day);

    if ( (ptr = fopen(file->name,"r")) == NULL) {
        fprintf(stderr,"-E- %s line %d: unable to open %s\n",
                __FILE__,__LINE__,file->name);
        return 1;
    }


    int numLinesNeeded   = 1;
    int numSamplesNeeded = 1;
    int bandsNeeded      = 1;
    int waveLengthNeeded = 1;
    int fwhmNeeded       = 1;
    int utmZoneNeeded    = 1;
    int eastingNeeded    = 1;
    int northingNeeded   = 1;
    int pixelSizeNeeded  = 1;
    int interleaveNeeded = 1;
    int rotationAngNeeded= 1;

    // loop metadata
    while( numLinesNeeded ||
            numSamplesNeeded ||
            bandsNeeded      ||
            fwhmNeeded      ||
            pixelSizeNeeded ||
            utmZoneNeeded ||
            eastingNeeded ||
            northingNeeded ||
            waveLengthNeeded ||
            rotationAngNeeded ||
            interleaveNeeded) {

        result = fgets(line, itemSize, ptr);
         if(result == NULL) {
            fprintf(stderr,"-E- %s line %d: unable to read all of the required metadata from AVIRIS file\n",
                    __FILE__,__LINE__);
            exit(1);
        }
        trimBlanks(line);

        if((val0 = checkTagLine(line,"lines"))) {
            numLinesNeeded = 0;
            file->nscan = atoi(val0);
            data->nscan = file->nscan;
        }
        if((val0=checkTagLine(line,"samples"))) {
            numSamplesNeeded = 0;
            file->npix = atoi(val0);
            data->npix = file->npix;
       }
        if((val0=checkTagLine(line, "bands"))) {
            bandsNeeded = 0;
            numBands = atoi(val0);
            data->numBands = numBands;
            if (numBands > maxBands) {
                fprintf(stderr,"-E- %s line %d: number of bands (%d) from AVIRIS file > maxBands (%d)\n",
                        __FILE__,__LINE__,numBands,maxBands);
                exit(1);
            }
        }
        if((val0=checkTagLine(line,"interleave"))) {
            interleaveNeeded = 0;
            if (strstr(val0,"bip")) {
                data->interleave = BIP;
            } else if (strstr(val0,"bil")) {
                data->interleave = BIL;
            } else {
                fprintf(stderr,"Interleave = %s is not supported\n",val0);
                exit(1);
            }
       }

        if((val0=checkTagLine(line,"rotation angle"))) {
            rotationAngNeeded = 0;
            rotation = atof(val0);
            if (rotation > 45)
                data->eastbyscan = -1;
             else if (rotation < -45)
                 data->eastbyscan = 1;
             else
                 data->eastbyscan = 0;
        }

        if((val0=checkTagLine(line,"pixel size"))) {
            pixelSizeNeeded = 0;
            data->pixelSize = atof(val0);
        }
        if((val0=checkTagLine(line,"Northing"))) {
            northingNeeded = 0;
            data->northing = atof(val0);
        }
        if((val0=checkTagLine(line,"Easting"))) {
            eastingNeeded = 0;
            data->easting = atof(val0);
        }
        if((val0=checkTagLine(line,"UTM zone"))) {
            utmZoneNeeded = 0;
            data->utmZone = atoi(val0);
        }

        if((val0=checkTagLine(line,"wavelength"))) {
            waveLengthNeeded = 0;
            i = 0;
            readWavInfo_jpl(ptr, tag, val);
            while (i < maxBands && strcmp(tag,"}")) {
                data->wave[i] = atof(val);
                readWavInfo_jpl(ptr, tag, val);
                i++;
            }

            if (!strcmp(tag,"}") && i <= maxBands ) {
                data->wave[i] = atof(val);
                i++;
            } else { // if (i> maxBands) {

                fprintf(stderr,"-E- %s line %d: number of bands (%d) from AVIRIS file > maxBands (%d)\n",
                        __FILE__,__LINE__,file->nbands,maxBands);
                exit(1);
            }
        }

        numBands = i;

        if((val0=checkTagLine(line, "fwhm"))) {
            fwhmNeeded = 0;
            i = 0;
            readWavInfo_jpl(ptr, tag, val);
            while (i < maxBands && strcmp(tag,"}")) {
                data->fwhm[i] = atof(val);
                readWavInfo_jpl(ptr, tag, val);
                i++;
            }
            if (!strcmp(tag,"}") && i < maxBands ) {
                data->fwhm[i] = atof(val);
                i++;
            } else {
                fprintf(stderr,"-E- %s line %d: number of bands (%d) from AVIRIS file > maxBands (%d)\n",
                        __FILE__,__LINE__,file->nbands,maxBands);
                exit(1);
            }
        }
    }

    fclose(ptr);

    // Get info about the WGS84 data
    // Get the lat/lon/elev
    numLinesNeeded   = 1;
    numSamplesNeeded = 1;

    infile = malloc((pos+strlen("_obs.hdr"))*sizeof(char));
    strcpy(infile, inbasename);
    strcat(infile, "_obs.hdr");

    if ( (ptr = fopen(infile,"r")) == NULL) {
        fprintf(stderr,"-E- %s line %d: unable to open %s\n",
                __FILE__,__LINE__,infile);
        return 1;

    }
    // loop metadata
    while( numLinesNeeded ||
            numSamplesNeeded ) {

        readNextLine_jpl(ptr, tag, val);
        // skip blank lines
        if(tag[0] == 0)
            continue;

        // get date
        if(!strcmp(tag, "lines")) {
           numLinesNeeded = 0;
            data->wgs_nscan = atoi(val);
        }
        if(!strcmp(tag, "samples")) {
            numSamplesNeeded = 0;
            data->wgs_npix = atoi(val);
        }
    }

    fclose(ptr);


     //Get the lat/lon/elev

    infile = malloc((pos+strlen("_lonlat_eph"))*sizeof(char));
    strcpy(infile, inbasename);
    strcat(infile, "_lonlat_eph");

    printf("Reading lon/lat/elev information from file %s\n", infile);

    // allocate lat and lon storage
    num=6;
//    data->lat  = (double *) malloc(data->wgs_nscan*sizeof(double) );
//    data->lon  = (double *) malloc(data->wgs_nscan*sizeof(double) );
    elev   = (double *) malloc(data->wgs_nscan*sizeof(double) );
    indata = (double *) malloc(data->wgs_nscan*num*sizeof(double) );

//    if(data->lat==NULL || data->lon==NULL || data->elev==NULL || indata==NULL) {
        if(elev==NULL || indata==NULL) {
        fprintf(stderr,"-E- %s line %d: unable to allocate lat/lon data for aviris\n",
                __FILE__,__LINE__);
        exit(1);
    }

    if ( (ptr = fopen(infile,"rb")) == NULL) {
        fprintf(stderr,"-E- %s line %d: unable to open %s\n",
                __FILE__,__LINE__,infile);
        return 1;

    }
    status = fread(indata,sizeof(double),data->wgs_nscan*num,ptr);
    if (status != data->wgs_nscan*num) {
        printf("Wrong data read size: want %d got %d in file %s\n",data->wgs_nscan*num,status, infile);
        exit(1);
    }

    i = 0;

    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, data->wgs_nscan);

    double *dist = (double *) calloc( data->wgs_nscan, sizeof(double));
    double *xlon = (double *) calloc( data->wgs_nscan, sizeof(double));
    double *xlat = (double *) calloc( data->wgs_nscan, sizeof(double));

    data->spl_acc = gsl_interp_accel_alloc();
    data->lon0 = indata[0];
    data->lat0 = indata[1];
    while (i<data->wgs_nscan*num) {
        j=i/num;
        xlon[j]  = indata[i];
        xlat[j]  = indata[i+1];
        elev[j]  = indata[i+2];
        dist[j] = pow((xlon[j]-data->lon0),2.0) + pow((xlat[j]-data->lat0),2.0) ;
        //printf("RJH: elev[%d]=%lf dist=%lf lat=%f lon=%f lat0=%lf lon0=%lf\n",j,elev[j],dist[j],xlat[j],xlon[j], data->lat0, data->lon0);
        i+=num;
    }
   // Sort distances and corresponding elevation values
    gsl_sort2( dist, 1, elev, 1, data->wgs_nscan);

    // Initiate spline
    gsl_spline_init(spline, dist, elev, data->wgs_nscan);

    data->spline = spline;

    free(indata);
    fclose(ptr);


    // Get the sensor and solar data

    num = 10;
    indataf = (float *) malloc(file->npix*sizeof(float)*num);

    data->sena = (double **) malloc(file->nscan*sizeof(double*));
    data->senz = (double **) malloc(file->nscan*sizeof(double*));
    data->solz = (double **) malloc(file->nscan*sizeof(double*));
    data->sola = (double **) malloc(file->nscan*sizeof(double*));
    data->utc  = (double **) malloc(file->nscan*sizeof(double*));
    if(data->sena==NULL || data->senz==NULL || data->sola==NULL || data->solz==NULL || indataf == NULL) {
        fprintf(stderr,"-E- %s line %d: unable to allocate sensor and solar angle data for AVIRIS\n",
                __FILE__,__LINE__);
        exit(1);
    }

    for (i=0;i<file->nscan; i++) {
        data->sena[i] = (double *) malloc(file->npix*sizeof(double));
        data->senz[i] = (double *) malloc(file->npix*sizeof(double));
        data->sola[i] = (double *) malloc(file->npix*sizeof(double));
        data->solz[i] = (double *) malloc(file->npix*sizeof(double));
        data->utc[i]  = (double *) malloc(file->npix*sizeof(double));
        if(data->sena[i]==NULL || data->senz[i]==NULL || data->sola[i]==NULL || data->solz[i]==NULL) {
            fprintf(stderr,"-E- %s line %d: unable to allocate sensor and solar angle data for AVIRIS\n",
                    __FILE__,__LINE__);
            exit(1);
        }
    }

    //free(infile);
    infile = malloc((pos+strlen("_obs_ort"))*sizeof(char));
    strcpy(infile, inbasename);
    strcat(infile, "_obs_ort");

    printf("Reading sensor and solar angles information from file %s\n", infile);

    if ( (ptr = fopen(infile,"rb")) == NULL) {
        fprintf(stderr,"-E- %s line %d: unable to open %s\n",
                __FILE__,__LINE__,infile);
        return 1;

    }
    for (j=0;j<file->nscan;j++) {
        status = fread(indataf,sizeof(float),num*file->npix,ptr);
        if (status != num*file->npix) {
            fprintf(stderr,"-E- %s line %d: AVIRIS Wrong sensor and solar data read size: want %d got %d\n",
                __FILE__,__LINE__,file->npix*num,status);
            exit(1);
        }

        for(k=0;k<file->npix;k++) {
            data->sena[j][k] = indataf[1*file->npix+k];
            data->senz[j][k] = indataf[2*file->npix+k];
            data->sola[j][k] = indataf[3*file->npix+k];
            data->solz[j][k] = indataf[4*file->npix+k];
            data->utc [j][k] = indataf[9*file->npix+k];
            //printf("scan=%d utc[%d]=%f\n",j,k,data->utc[j][k]);
           //if (data->sena[j][k] > -999) printf("scan=%d sena[%d]=%lf senz[%d]=%lf sola[%d]=%lf solz[%d]=%lf solznight=%f\n",j,k,data->sena[j][k],k,data->senz[j][k],k,data->sola[j][k],k,data->solz[j][k],SOLZNIGHT);
        }
    }

    free(indataf);
    fclose(ptr);


    // Get the gain data
    free(infile);
    infile = malloc((pos+strlen(".gain"))*sizeof(char));
    strcpy(infile, inbasename);
    strcat(infile, ".gain");

    printf("Reading gain information from file %s\n", infile);

    if ( (ptr = fopen(infile,"r")) == NULL) {
        fprintf(stderr,"-E- %s line %d: unable to open %s\n",
                __FILE__,__LINE__,infile);
        return 1;

    }

    data->gain = (double *) malloc(numBands*sizeof(double));
    if(data->gain==NULL) {
        fprintf(stderr,"-E- %s line %d: unable to allocate gain data for AVIRIS\n",
                __FILE__,__LINE__);
        exit(1);
    }
    for (i=0;i<numBands && fscanf(ptr,"%lf %d",(data->gain+i),&k);i++) {
//        printf("gain: %lf %d\n",*(data->gain+i),k);
    }

    fclose(ptr);

    free(infile);
    infile = malloc((pos+strlen("_sc01_ort_img"))*sizeof(char));
    strcpy(infile, inbasename);
    strcat(infile, "_sc01_ort_img");
    printf("Opening AVIRIS image file %s\n", infile);

    if ( (data->av_fp = fopen(infile,"rb")) == NULL) {
        fprintf(stderr,"-E- %s line %d: unable to open %s\n",
                __FILE__,__LINE__,infile);
        return 1;

    }

    sprintf(projStr, "+proj=utm +zone=%d +ellps=WGS84 +datum=WGS84 +units=m +no_defs ",
            data->utmZone);
    //printf("RJH: projStr=%s\n",projStr);
    if (!(data->pj_ortho = pj_init_plus(projStr))) {
        printf("Error - AVIRIS first projection failed to init\n");
        exit(1);
    }
    if (!(data->pj_latlong = pj_latlong_from_proj(data->pj_ortho))) {
        fprintf(stderr,"-E- %s line %d: AVIRIS latlon projection failed to init\n",
            __FILE__,__LINE__);
        exit(1);
    }

    lat = (double *)   malloc ( file->npix*sizeof(double) );
    lon = (double *)   malloc ( file->npix*sizeof(double) );

    return(0);

}

int readl1_aviris( filehandle *file, int recnum, l1str *l1rec, int lonlat)
/*
 *  fill standard record with L1B line of data
 */
{
    int   status, min_msec = 86401 * 1000, max_msec = -1;

    static double last_good_hour=18;
    static int firstCall=1;
    double elev, solz, sola;
    float rel_sec;
    double sec;
    int year, month, day,hour, minute ;
    double dist;
    int isec;
    int    npix=file->npix, ip, ib, ipb, msec;
    static int swap;

    aviris_t* data = (aviris_t*)file->private_data;

    if (firstCall) {
        if (want_verbose)
            printf("file->nbands = %d, l1rec->nbands = %d\n",
                    (int) file->nbands, (int) l1rec->nbands);
        firstCall = 0;

        for (ip = 0; ip < npix; ip++) {
            l1rec->pixnum[ip] = ip;
            l1rec->flags[ip] = 0;
        }

        if ( endianess() == 1 )
            swap = 1;
        else
        	swap = 0;

    }

    //  set information about data
    l1rec->npix = file->npix;
    l1rec->sensorID = file->sensorID;

    int k=0;
    while ((hour = (int)(data->utc[recnum][k])) <0 && k<npix) k++;

    if (hour<0) hour = last_good_hour;

    last_good_hour = hour;

    minute = (data->utc[recnum][k] - hour)*60;
    sec = ((data->utc[recnum][k] - hour)*60 - minute)*60;
//    printf("Date:utc=%f sec=%f minute=%d (utc-hour)*60=%f ",data->utc[recnum][k],sec,minute,(data->utc[recnum][k] - hour)*60 );
    isec = (int)sec;


    ymdhms2ydmsec(data->year, data->month, data->day, hour, minute, isec,
                       &data->year, &data->doy, &data->msec);

    data->msec = sec * 1000;

    *(l1rec->year) = data->year;
    *(l1rec->day)  = data->doy;
    *(l1rec->msec) = data->msec;

//    printf("Date=%4d %d %d\n",*(l1rec->year),*(l1rec->day),*(l1rec->msec));

    npix   = file->npix;

    //  set information about data

    l1rec->npix = file->npix;
    l1rec->sensorID = file->sensorID;

    //  get lat-lon
    for (ip=0; ip<npix; ip++) {
        if (data->eastbyscan != 0) {
            lon[ip] = data->easting  + data->eastbyscan * recnum * data->pixelSize;     // rotated
            lat[ip] = data->northing + data->eastbyscan * ip     * data->pixelSize ;
        }else {
            lon[ip] = data->easting + ip * data->pixelSize;     // starts in upper left corner
            lat[ip] = data->northing - recnum*data->pixelSize ;
        }
    }

    aviris_proj4_convert(data, npix, lon, lat);

    for (ip=0; ip<npix; ip++) {

        l1rec->pixnum[ip] = ip;

        if ( isnan(lat[ip]) ) lat[ip] = -999.0;
        if ( isnan(lon[ip]) ) lon[ip] = -999.0;
        l1rec->lat[ip]    = lat[ip];
        l1rec->lon[ip]    = lon[ip];


        if (l1rec->lon[ip] < -181.0 || l1rec->lon[ip] > 181.0 ||
            l1rec->lat[ip] <  -91.0 || l1rec->lat[ip] >  91.0 )
          l1rec->navfail[ip] = 1;

        if (data->senz[recnum][ip] > SKIP)
            l1rec->senz[ip] = data->senz[recnum][ip];
        else
            l1rec->senz[ip] = getValidAngle(data->senz[recnum], npix, SKIP);

        if (data->sena[recnum][ip] > SKIP)
            l1rec->sena[ip] = data->sena[recnum][ip];
        else
            l1rec->sena[ip] = getValidAngle(data->sena[recnum], npix, SKIP);

        if (data->solz[recnum][ip] > SKIP)
            l1rec->solz[ip] = data->solz[recnum][ip];
        else
            l1rec->solz[ip] = getValidAngle(data->solz[recnum], npix, SKIP);

        if (data->sola[recnum][ip] > SKIP)
            l1rec->sola[ip] = data->sola[recnum][ip];
        else
            l1rec->sola[ip] = getValidAngle(data->sola[recnum], npix, SKIP);

    }

    // find interpolated elevation from wgs-84 lat/lon
    ip = npix/2;
        dist = (l1rec->lon[ip]-data->lon0)*(l1rec->lon[ip]-data->lon0) + (l1rec->lat[ip]-data->lat0)*(l1rec->lat[ip]-data->lat0);
        //printf("RJH: dist=%f lat=%f lon=%f lat0=%lf lon0=%lf\n",dist,l1rec->lat[ip],l1rec->lon[ip],data->lat0,data->lon0);
         l1rec->alt = (float)gsl_spline_eval( data->spline, dist, data->spl_acc);
//        printf("RJH: alt[%d]=%lf \n",ip,l1rec->alt);
//        ip=0;
//        printf("RJH:IPMIN: dist=%f lat=%f lon=%f lat0=%lf lon0=%lf\n",dist,l1rec->lat[ip],l1rec->lon[ip],data->lat0,data->lon0);
//        ip=npix-1;
//        printf("RJH:IPMAX: dist=%f lat=%f lon=%f lat0=%lf lon0=%lf\n",dist,l1rec->lat[ip],l1rec->lon[ip],data->lat0,data->lon0);

    readBinScanLine_int2(l1rec->Lt, recnum, l1rec->npix, data->gain, l1rec->nbands, data->numBands, data->interleave, swap, data->av_fp);

    return(LIFE_IS_GOOD);
}

void aviris_proj4_convert(aviris_t * data, int numPoints, double *x, double *y)
{
    int i;
    if(pj_transform(data->pj_ortho, data->pj_latlong, numPoints, 1, x, y, NULL)) {
        fprintf(stderr,"-E- %s line %d: AVIRIS proj4 transformation blew up\n",
            __FILE__,__LINE__);
        exit(1);
    }
    for(i=0; i<numPoints; i++) {
        x[i] *= RAD_TO_DEG;
        y[i] *= RAD_TO_DEG;
    }
}

int closel1_aviris(filehandle *file) {
    int ib;
    aviris_t* data = (aviris_t*) file->private_data;

    // undo what open allocated

    freePrivateData_av(data);
    free(file->private_data);

    return 0;
}

