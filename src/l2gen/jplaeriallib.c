/*
 * jplaeriallib.c
 *
 *  Created on: Jun 12, 2015
 *      Author: r. healy
 *      Various functions for BIP and BIL header files (AVIRIS and PRISM presently)
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "jplaeriallib.h"

char* getinbasename_av(char *file) {

    int pos,k;
    char *inbasename;

    pos = 0;
    k=strlen(file);
    while (k>0 && file[k] != '/') {
        if (file[k] == '_') pos = k;
        k--;
    }

    if (pos>0) {
        pos+=2; //file should have an _x
        inbasename = (char *) malloc(pos*sizeof(char));
        strncpy(inbasename, file,pos);
    } else
        inbasename = NULL;

    return(inbasename);

}
char* getinbasename(char *file) {

    int pos,k;
    char *inbasename;

    pos = 0;
    k=strlen(file);
    while (k>0 && file[k] != '/') {
        if (file[k] == '_') pos = k;
        k--;
    }

    printf("getinbasename=%s\n",file);
    if (pos>0) {
        printf("Inbase: pos=%d\n",pos);
        inbasename = (char *) malloc((pos+1)*sizeof(char));
        strncpy(inbasename, file,pos);
        inbasename[pos] = '\0';
    } else
        inbasename = NULL;

    printf("2)getinbasename=%s\n",inbasename);
    return(inbasename);

}
void readWavInfo_jpl(FILE* fp, char* tag, char* val) {
    char* result;
    char line[itemSize];
    int count;

    result = fgets(line, itemSize, fp);
    if(result == NULL) {
        fprintf(stderr,"-E- %s line %d: unable to read all of the required metadata from prism file\n",
                __FILE__,__LINE__);
        exit(1);
    }
    trimBlanks(line);

    count = sscanf(line, "%s %s", val,tag);
    if(count != 2) {
            // not found so return blank line
            tag[0] = 0;
            val[0] = 0;
            return;
    }

    trimBlanks(val);

}

void readNextLine_jpl(FILE* fp, char* tag, char* val) {
    char* result;
    char line[itemSize];
    int count;

    result = fgets(line, itemSize, fp);
    if(result == NULL) {
        fprintf(stderr,"-E- %s line %d: unable to read all of the required metadata from prism file\n",
            __FILE__,__LINE__);
        exit(1);
    }
    trimBlanks(line);

    count = sscanf(line, "%s = %s", tag, val);
    if(count != 2) {
        // not found so return blank line
        tag[0] = 0;
        val[0] = 0;
        return;
    }

    trimBlanks(val);
}


/*
 * A generic reader for 2 byte integer binary BIL/BIP files
 * @param out Lt array
 * @param in recnum - scan line number
 * @param in npix   - number of pixels in scan line
 * @param in gain array - factor to divide Lt by band number
 * @param in nbands     - number of bands in L2 sensor
 * @param in numBands   - number of bands in input file (numBands >= nbands)
 * @param in interleave - binary interleave (either BIP or BIL)
 * @param in ptr        - input file pointer
*/
int readBinScanLine_int2(float *Lt, int32_t recnum, int32_t npix, double *gain, int nbands, int numBands, int interleave, int swap, FILE *ptr) {

    static int32_t last_recnum=0;
    int i,j,ib,ibp, ipb_av,ip,ipb, status;
    static int firstCall=1;
    int16_t *ibuf;

    ibuf = (int16_t *) malloc(npix*numBands*sizeof(int16_t) );

    if (!ibuf) {
        fprintf(stderr,"-E- %s line %d: Out of Memory allocating buffer array.\n",
                __FILE__,__LINE__);
        exit(1);

    }

    fseek(ptr,sizeof(int16_t)*numBands*npix*(recnum-1),0);
    status = fread(ibuf,sizeof(int16_t),numBands*npix,ptr);
    if ( status != numBands*npix) {
        fprintf(stderr,"-E- %s line %d: Wrong Lt data read size: want %d got %d\n",
            __FILE__,__LINE__,npix*numBands,status);
        exit(1);
    }
    if ( swap == 1 )
        swapc_bytes((char *)ibuf    ,sizeof(int16_t),npix*numBands);

    for(ip=0;ip<npix;ip++) {
        for(ib = 0; ib < nbands; ib++) {
            ipb    = ip*nbands+ib;
            switch (interleave) {
                case BIP:
                    ipb_av = ip*numBands+ib;
                    break;
                case BIL:
                    ipb_av = ib*npix+ip;
                    break;
            }

            //ibp    = ib*npix+ip;

            if(ibuf[ipb_av] <= 0) {
                Lt[ipb] = BAD_FLT;                            // I assume this is outside the projected tile
                //                l1rec->navfail[ip] = 1;     // so set navigation failure  - Commented out for hyperspectal - RJH
            } else {
                Lt[ipb] = ibuf[ipb_av] / gain[ib];  //uWatts/cm^2/sr
            }
        }
    }

    free(ibuf);
    return(status);
}



/*
 * A generic reader for 4 byte float binary BIL/BIP files
 * @param out Lt array
 * @param in recnum - scan line number
 * @param in npix   - number of pixels in scan line
 * @param in gain array - factor to divide Lt by band number
 * @param in nbands     - number of bands in L2 sensor
 * @param in numBands   - number of bands in input file (numBands >= nbands)
 * @param in interleave - binary interleave (either BIP or BIL)
 * @param in ptr        - input file pointer
*/
int readBinScanLine_float(float *Lt, int32_t recnum, int32_t npix, double *gain, int nbands, int numBands, int interleave, int swap, FILE *ptr) {

    static int32_t last_recnum=0;
    int i,j,ib,ibp, ipb_av,ip,ipb, status;
    static int firstCall=1;
    float *ibuf;

    ibuf = (float *) malloc(npix*numBands*sizeof(float) );

    if (!ibuf) {
        fprintf(stderr,"-E- %s line %d: Out of Memory allocating buffer array.\n",
                __FILE__,__LINE__);
        exit(1);

    }

    fseek(ptr,4*numBands*npix*(recnum-1),0);
    status = fread(ibuf,4,numBands*npix,ptr);
    if ( status != numBands*npix) {
        fprintf(stderr,"-E- %s line %d: Wrong Lt data read size: want %d got %d\n",
            __FILE__,__LINE__,npix*numBands,status);
        exit(1);
    }

    if (swap)
    	swapc_bytes((char *)ibuf    ,4,npix*numBands);

    for(ip=0;ip<npix;ip++) {
        for(ib = 0; ib < nbands; ib++) {
            ipb    = ip*nbands+ib;
            switch (interleave) {
                case BIP:
                    ipb_av = ip*numBands+ib;
                    break;
                case BIL:
                    ipb_av = ib*npix+ip;
                    break;
            }

            if(ibuf[ipb_av] <= 0) {
                Lt[ipb] = BAD_FLT;                            // I assume this is outside the projected tile
                //                l1rec->navfail[ip] = 1;     // so set navigation failure  - Commented out for hyperspectal - RJH
            } else {
                Lt[ipb] = ibuf[ipb_av] / gain[ib];  //uWatts/cm^2/sr
            }
        }
    }

    free(ibuf);
    return(status);
}

double getValidAngle(double *ang, int32_t npix, int32_t skip) {
    int32_t i,k;

    for (i=0;i<npix && ang[i] <= skip; i++)
        k=i;

    return(ang[k]);

}

char* checkTagLine(char *linein, char* tag) {
    char* result;
    int count;
    char line[1024];

    strcpy(line,linein);
    count = 0;
    if (strstr(line,tag) ) {
        result = strtok(line,"=");
        if (result) {
            result = strtok(NULL,"=");
            count = 1;
        }
    }
    if (count < 1) {
        return(NULL);
    }

    trimBlanks(result);

    return (result);
}
/* Find the tag within a line that has multiple tokens
 * Return the value to the right of the token
 */
char* checknspTagLine(char *linein, char* tag) {

    char* result;
    int count;
    char line[1024];

    strcpy(line,linein);

    count = 0;
    if (strstr(line,tag) ) {
        result = strtok(line,"=");
        while (!strstr(result,tag)) {
             result = strtok(NULL,"=");
       }
        if (result) {
            result = strtok(NULL,"=");
            count = 1;
        }
    }
    if (count < 1) {
        return(NULL);
    }

    return (result);
}

void getPosVec(float lat,float lon, float alt, double *pos) {

    int i;
    double gv[3];
    double hv[3];
    double re = 6378.137;
    double f = 1 / 298.257;
    double omf2 = (1.0 - f) * (1.0 - f);
    double slat = lat * DEG_TO_RAD;
    double slon = lon * DEG_TO_RAD;

    // Compute geocentric radius
    double slatg = atan(tan(slat) * omf2);
    double r = re * (1.0 - f) / sqrt(1.0 - (2.0 - f) * f * pow(cos(slatg), 2));

    // Compute geocentric subsensor vector
    gv[0] = cos(slatg) * cos(slon);
    gv[1] = cos(slatg) * sin(slon);
    gv[2] = sin(slatg);

    // Compute geodetic height vector
    hv[0] = cos(slat) * cos(slon);
    hv[1] = cos(slat) * sin(slon);
    hv[2] = sin(slat);

    // Compute sensor position vector
    for (i = 0; i < 3; i++) {
        pos[i] = gv[i] * r + hv[i] * alt;
    }

}
