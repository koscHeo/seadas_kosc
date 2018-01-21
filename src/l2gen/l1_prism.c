/*
 * l1_prism.c
 *
 *  Portable Remote Imaging SpectroMeter (PRISM)
 *
 *  Created on: June 11, 2015
 *      Author: Rick Healy SAIC
 *              NASA-GSFC OBPG
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "filehandle.h"
#include "l12_proto.h"
#include <proj_api.h>
#include "jplaeriallib.h"
#include "prism.h"

#define SKIP -9999

static const int maxBands = 285;
static double *lat, *lon;

prism_t* createPrivateData_pr (int numBands, int32_t nscan, int32_t npix) {

    int i;

    prism_t* data = (prism_t*)calloc(1,sizeof(prism_t));
    if(data == NULL) {
        fprintf(stderr,"-E- %s line %d: unable to allocate private data for prism\n",
                __FILE__,__LINE__);
        exit(1);
    }

    data->wave = (double *) malloc(numBands*sizeof(double) );
    data->fwhm = (double *) malloc(numBands*sizeof(double) );
    data->gain = (double *) malloc(numBands*sizeof(double) );
    if(data->wave==NULL || data->fwhm==NULL || data->gain==NULL) {
        fprintf(stderr,"-E- %s line %d: unable to allocate scale/offset data for prism\n",
                __FILE__,__LINE__);
        exit(1);
    }

    return data;
}


int  openl1_prism( filehandle *file ) {

    int16_t buffer[10];
    FILE *ptr;
    char tag[itemSize];
    char *val0;
    char val[itemSize];
    char *inbasename, *infile;
    static char *infile2;
    int i, j, k, status,pos;
    double *indata;
    float  *indataf;
    float rotation;
    int numBands, num;
    char* result;
    char line[itemSize],val1[itemSize];
    int count;
    int year, month, day, hour, minute, second, doy;
    float knts2mps = 0.51444444444; // knots to meters per second
    float ft2m     = 0.3048;        // feet to meters
    int isec = 0;
    double sec;

    char projStr[1024];
    static char *dupline;
    int cnt,linelength;

    printf("file=%s\n",file->name);
    inbasename = getinbasename(file->name);
    printf("Basename=%s\n",inbasename);
    pos = strlen(inbasename);
    if (pos<=0) {
        fprintf(stderr,"-E- %s line %d: Not a avalid prism file %s\n",
                __FILE__,__LINE__,file->name);
        return 1;
    }

    sscanf(inbasename, "prm%4d%2d%2dt%2d%2d%2d", &year, &month, &day,&hour,&minute,&second);

    prism_t* data = file->private_data = createPrivateData_pr(maxBands,file->nscan,file->npix);

//    if (year >= 92) year = year + 1900; else year = year + 2000;

    data->month = month;
    data->day   = day;

    ymdhms2ydmsec(year, month, day, hour, minute, second,
                       &data->year, &data->doy, &data->msec);

    printf("Date of prism flight: Y-%d M-%d D-%d %d %d %d %d\n",year,month,day, hour, minute,data->doy,data->msec);

    if ( (ptr = fopen(file->name,"r")) == NULL) {
        fprintf(stderr,"-E- %s line %d: unable to open %s\n",
                __FILE__,__LINE__,file->name);
        return 1;
    }

//    if (!infile2) {
//        linelength = 128;
//        infile2 = (char *) malloc(linelength*sizeof(char));
//    }
//    infile2 = malloc((pos+strlen("_rdn_ort"))*sizeof(char));

//                if ((dupline = (char *) malloc(1024)) == NULL) {
//                    fprintf(stderr,
//                        "-E- %s line %d: Memory allocation failure.\n",
//                        __FILE__,__LINE__);
//                    return(0);
//                }

    int numLinesNeeded   = 1;
    int numSamplesNeeded = 1;
    int bandsNeeded      = 1;
    int utmZoneNeeded    = 1;
    int eastingNeeded    = 1;
    int northingNeeded   = 1;
    int pixelSizeNeeded  = 1;
    int interleaveNeeded = 1;
    int rotationAngNeeded= 1;
    int EndTimeNeeded= 1;
    int altNeeded= 1;

    // loop metadata
    while( numLinesNeeded ||
            numSamplesNeeded ||
            bandsNeeded      ||
            pixelSizeNeeded ||
            utmZoneNeeded ||
            eastingNeeded ||
            northingNeeded ||
            rotationAngNeeded ||
            EndTimeNeeded ||
            altNeeded ||
            interleaveNeeded) {

        result = fgets(line, itemSize, ptr);
         if(result == NULL) {
            fprintf(stderr,"-E- %s line %d: unable to read all of the required metadata from prism file\n",
                    __FILE__,__LINE__);
            exit(1);
        }
        trimBlanks(line);

        if((val0 = checkTagLine(line,"lines"))) {
            numLinesNeeded = 0;
            file->nscan = atoi(val0);
            data->nscan = file->nscan;
            printf("lines=%d\n",data->nscan);
        }
        if((val0=checkTagLine(line,"samples"))) {
            numSamplesNeeded = 0;
            file->npix = atoi(val0);
            data->npix = file->npix;
            printf("samples=%d\n",data->npix);
       }
        if((val0=checkTagLine(line,"Alt"))) {
            altNeeded = 0;
            data->alt = atof(val0)*ft2m;
            printf("Altitude=%lf\n",data->alt);
       }
        if((val0=checkTagLine(line,"EndTime"))) {
            EndTimeNeeded = 0;
            sscanf(val0, "%2d%2d", &hour,&minute);
            printf("End hour=%d minute=%d\n",hour,minute);
            second = 0;
            ymdhms2ydmsec(year, month, day, hour, minute, second,
                               &data->eyear, &data->edoy, &data->emsec);

            printf("End Time=%d/%d\n",data->edoy,data->emsec);
       }
        if((val0=checkTagLine(line, "bands"))) {
            bandsNeeded = 0;
            numBands = atoi(val0);
            data->numBands = numBands;
            if (numBands > maxBands) {
                fprintf(stderr,"-E- %s line %d: number of bands (%d) from prism file > maxBands (%d)\n",
                        __FILE__,__LINE__,numBands,maxBands);
                exit(1);
            }
            printf("numBands=%d\n",data->numBands);
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
            printf("Interleave=%d\n",data->interleave);
       }
        if((val0=checkTagLine(line,"map info"))) {
            cnt = 0;
            linelength = strlen(line);
            if (dupline) free(dupline);
            if ((dupline = (char *) malloc(linelength*sizeof(char))) == NULL) {
                fprintf(stderr,
                    "-E- %s line %d: Memory allocation failure.\n",
                    __FILE__,__LINE__);
                return(0);
            }
            strcpy(dupline,line);
            result = strtok(dupline,",");
            while (result) {
                switch(cnt) {
                case 3:
                    data->easting   = atof(result);
                    eastingNeeded = 0;
                    break;
                case 4:
                    data->northing  = atof(result);
                    northingNeeded = 0;
                    break;
                case 5:
                    data->pixelSize = atof(result);
                    pixelSizeNeeded = 0;
                    break;
                case 7:
                    data->utmZone   = atoi(result);
                    utmZoneNeeded = 0;
                    break;
                default:
                    break;
                }
                printf(">>%d) %s\n",cnt,result);
                cnt ++;
                result = strtok(NULL,",");
            }
            if((val0=checknspTagLine(line,"rotation"))) {
                rotationAngNeeded = 0;
                rotation = atof(val0);
                if (rotation > 45)
                    data->eastbyscan = -1;
                 else if (rotation < -45)
                     data->eastbyscan = 1;
                 else
                     data->eastbyscan = 0;
            } else {
                printf("Rotation angle expected in line: %s\n",val0);
                exit(-1);
            }
            printf("Rotation/easting/northing/pixsize/utmzone=%f/%lf/%lf/%lf/%d\n",rotation,data->easting,data->northing,data->pixelSize,data->utmZone);
        }
//        if((val0=checkTagLine(line,"map info"))) {
//            sscanf(val0,"%*[^,], %*f, %*f, %lf, %lf, %lf, %lf, %d, %*[^,], %*[^,], %*[^,], %s}",&data->easting,&data->northing,&data->pixelSize,&data->pixelSize,&data->utmZone,val1);
//            if((val0=checknspTagLine(line,"rotation"))) {
//                rotationAngNeeded = 0;
//                rotation = atof(val0);
//                if (rotation > 45)
//                    data->eastbyscan = -1;
//                 else if (rotation < -45)
//                     data->eastbyscan = 1;
//                 else
//                     data->eastbyscan = 0;
//            } else {
//                printf("Rotation angle expected in line: %s\n",val0);
//                exit(-1);
//            }
//            printf("Rotation/easting/northing/pixsize/utmzone=%f/%lf/%lf/%lf/%d\n",rotation,data->easting,data->northing,data->pixelSize,data->utmZone);
//            pixelSizeNeeded = 0;
//            northingNeeded = 0;
//            eastingNeeded = 0;
//            utmZoneNeeded = 0;
//
//        }


    }

    fclose(ptr);


    // Get the sensor and solar data
    data->sena = (double **) malloc(data->nscan*sizeof(double*));
    data->senz = (double **) malloc(data->nscan*sizeof(double*));
    data->solz = (double **) malloc(data->nscan*sizeof(double*));
    data->sola = (double **) malloc(data->nscan*sizeof(double*));
    data->utc  = (double **) malloc(data->nscan*sizeof(double*));
    if(data->sena==NULL || data->senz==NULL || data->sola==NULL || data->solz==NULL ) {
        fprintf(stderr,"-E- %s line %d: unable to allocate sensor and solar angle data for prism\n",
                __FILE__,__LINE__);
        exit(1);
    }

    for (i=0;i<data->nscan; i++) {
        data->sena[i] = (double *) malloc(data->npix*sizeof(double));
        data->senz[i] = (double *) malloc(data->npix*sizeof(double));
        data->sola[i] = (double *) malloc(data->npix*sizeof(double));
        data->solz[i] = (double *) malloc(data->npix*sizeof(double));
        data->utc[i]  = (double *) malloc(data->npix*sizeof(double));
        if(data->sena[i]==NULL || data->senz[i]==NULL || data->sola[i]==NULL || data->solz[i]==NULL || data->utc[i]==NULL) {
            fprintf(stderr,"-E- %s line %d: unable to allocate sensor and solar angle data for prism\n",
                    __FILE__,__LINE__);
            exit(1);
        }
    }

    for (i=0;i<numBands;i++) {
        *(data->gain+i) = 1.0;
    }

   // free(infile);
    infile = malloc((pos+strlen("_rdn_ort"))*sizeof(char));
    strcpy(infile, inbasename);
    strcat(infile, "_rdn_ort");
    printf("Opening prism image file %s\n", infile);

    if ( (data->av_fp = fopen(infile,"rb")) == NULL) {
        fprintf(stderr,"-E- %s line %d: unable to open %s\n",
                __FILE__,__LINE__,infile);
        return 1;

    }

    sprintf(projStr, "+proj=utm +zone=%d +ellps=WGS84 +datum=WGS84 +units=m +no_defs ",
            data->utmZone);
//    printf("RJH: projStr=%s\n",projStr);
    if (!(data->pj_ortho = pj_init_plus(projStr))) {
        printf("Error - prism first projection failed to init\n");
        exit(1);
    }
    if (!(data->pj_latlong = pj_latlong_from_proj(data->pj_ortho))) {
        fprintf(stderr,"-E- %s line %d: prism latlon projection failed to init\n",
            __FILE__,__LINE__);
        exit(1);
    }

    lat = (double *)   malloc ( file->npix*sizeof(double) );
    lon = (double *)   malloc ( file->npix*sizeof(double) );

    return(0);

}

int readl1_prism( filehandle *file, int recnum, l1str *l1rec, int lonlat)
/*
 *  fill standard record with L1B line of data
 */
{
    int   status, min_msec = 86401 * 1000, max_msec = -1;

    static double last_good_hour=18;
    static int firstCall=1;
    float rel_sec;
    double sec,pos[3];
    float  epos[3], sunpos[3];
    int year, month, day,hour, minute ;
    double dist;
    double secondOfDay;
    float longitude, latitude, sunDist;
    int isec, i;
    int    npix=file->npix, ip, ib, ipb, msec;
    prism_t* data = (prism_t*)file->private_data;

    if (firstCall) {
        if (want_verbose)
            printf("file->nbands = %d, l1rec->nbands = %d npix=%d\n",
                    (int) file->nbands, (int) l1rec->nbands,npix);
        firstCall = 0;

        for (ip = 0; ip < npix; ip++) {
            l1rec->pixnum[ip] = ip;
            l1rec->flags[ip] = 0;
        }
    }

    //  set information about data
    l1rec->npix     = file->npix;
    l1rec->sensorID = file->sensorID;
    l1rec->alt      = data->alt;

    *(l1rec->year) = data->year;
    *(l1rec->day)  = data->doy;
    *(l1rec->msec) = data->msec + (data->emsec - data->msec)*recnum/(data->nscan-1);

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

    prism_proj4_convert(data, npix, lon, lat);

    if (lat[npix/2] > SKIP)
        latitude = lat[npix/2];
    else {
        fprintf(stderr,"-E- %s line %d: Don't have sensor latitude for geometry calculation\n",
            __FILE__,__LINE__);
        exit(1);
    }

    if (lon[npix/2] > SKIP)
        longitude = lon[npix/2];
    else {
        fprintf(stderr,"-E- %s line %d: Don't have sensor longitude for geometry calculation\n",
            __FILE__,__LINE__);
        exit(1);
    }

    getPosVec(latitude,longitude, data->alt, pos); // get position vector of sensor

    secondOfDay = *(l1rec->msec)/1000;
    l_sun_((l1rec->year),(l1rec->day),&secondOfDay,sunpos,&sunDist); // get position vector for the sun

    for (i=0; i < 3; i++) {
        sunpos[i] *= 1.496e8; //convert to km for call to get_zenaz
        epos[i]    = pos[i];
    }

    for (ip=0; ip<npix; ip++) {

        l1rec->pixnum[ip] = ip;

        if ( isnan(lat[ip]) ) lat[ip] = SKIP;
        if ( isnan(lon[ip]) ) lon[ip] = SKIP;
        l1rec->lat[ip]    = lat[ip];
        l1rec->lon[ip]    = lon[ip];
        //printf("ip=%d scan=%d lat=%f lon=%f\n",ip,recnum,lat[ip],lon[ip]);

        if (l1rec->lon[ip] < -181.0 || l1rec->lon[ip] > 181.0 ||
            l1rec->lat[ip] <  -91.0 || l1rec->lat[ip] >  91.0 )
          l1rec->navfail[ip] = 1;

        get_zenaz  (epos, lon[ip], lat[ip], &l1rec->senz[ip], &l1rec->sena[ip]);
        get_zenaz(sunpos, lon[ip], lat[ip], &l1rec->solz[ip], &l1rec->sola[ip]);

        //printf("RJH: %d %d senz=%f sena=%f solz=%f sola=%f\n",recnum, ip, l1rec->senz[ip],l1rec->sena[ip], l1rec->solz[ip],l1rec->sola[ip]);


    }

    readBinScanLine_float(l1rec->Lt, recnum, l1rec->npix, data->gain, l1rec->nbands, data->numBands, data->interleave, 0, data->av_fp);

    return(LIFE_IS_GOOD);
}


void prism_proj4_convert(prism_t * data, int numPoints, double *x, double *y)
{
    int i;
    if(pj_transform(data->pj_ortho, data->pj_latlong, numPoints, 1, x, y, NULL)) {
        fprintf(stderr,"-E- %s line %d: prism proj4 transformation blew up\n",
            __FILE__,__LINE__);
        exit(1);
    }
    for(i=0; i<numPoints; i++) {
        x[i] *= RAD_TO_DEG;
        y[i] *= RAD_TO_DEG;
    }
}

int closel1_prism(filehandle *file) {
    int ib;
    prism_t* data = (prism_t*) file->private_data;

    // undo what open allocated

    freePrivateData_pr(data);
    free(file->private_data);

    return 0;
}

void freePrivateData_pr(prism_t* data) {
    int k;
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

