/*
 * leapsecond.c
 *
 *  Created on: Feb 6, 2015
 *      Author: swbaile1
 */
#define LEAPSEC_DAT "leapsec.dat"
#define JD_1993 2448988.5
#define SECONDS_IN_DAY 86400
#define MAX_LEAP_SECONDS 16
static double leap_seconds_js[MAX_LEAP_SECONDS] = {0};

#include <timeutils.h>
#include <stdlib.h>
#include <stdio.h>

int leapseconds_since_1993(double tai93time) {
    char *varRoot;
    char leapsecdat[FILENAME_MAX];

    if (!leap_seconds_js[0]) {
        if ((varRoot = getenv("OCVARROOT")) == NULL) {
            printf("-E- OCVARROOT environment variable is not defined.\n");
            return (-1);
        }
        sprintf(leapsecdat, "%s/modis/leapsec.dat", varRoot);
        FILE *file_h = fopen(leapsecdat, "rb");
        int i = 0, year;
        float jd;
        while (fgetc(file_h) != '\n')
            ;
        while (fscanf(file_h, " %d %*s %*d =JD %f%*[^\n]", &year, &jd) == 2) {
            if (year >= 1993) {
                leap_seconds_js[i] = ((jd - JD_1993 + 1) * SECONDS_IN_DAY);
                i++;
            }
        }
        fclose(file_h);
    }

    int i;
    for (i = 0; i < MAX_LEAP_SECONDS - 1 && leap_seconds_js[i]; i++) {
        if (tai93time < leap_seconds_js[i]) {
            return i;
        }
    }
    return i;
}
