/*
 * dataday.h
 *
 *  Created on: Jan 22, 2015
 *      Author: rhealy
 */

#ifndef SRC_L2BIN_DATADAY_H_
#define SRC_L2BIN_DATADAY_H_
#include <time.h>

#define DATELINE_NOT_CROSSED    0
#define DATELINE_CROSSED        1
#define DATELINE_NORTH_POLE     2
#define DATELINE_SOUTH_POLE     3
#define DATELINE_IRREGULAR      4
#define INVALID_COORD   99999.999

#define EPOCH_YEAR              1970
#define SECONDS_PER_YEAR        31536000
#define SECONDS_PER_LEAP_YEAR   31622400
#define MALLOC(ptr,typ,num) {                                           \
  (ptr) = (typ *)malloc((num) * sizeof(typ));                           \
  if((ptr) == NULL){                                                    \
    fprintf(stderr,"-E- %s line %d: Memory allocation failure.\n",      \
    __FILE__,__LINE__);                                                 \
    exit(EXIT_FAILURE);                                                 \
  }                                                                     \
}

enum hort {
    head, tail
};

struct coord {
    float lat;
    float lon;
    struct coord *next;
};
void get_coord_extrema(int isnight, int n,
        float *olat, float *olon, float *north, float *south,
        float *west, float *east, int *dateline);
void get_datadays(time_t starttime, float eqxhour, int isnight,
        int dateline, float west, float east,
        int32_t *dataday0, int32_t *dataday1);
extern void cdata_();
extern void l_sun_(int *iyr, int *iday, double *sec, float sunr[3], float *rs);
void push(struct coord **vrtx, enum hort hort, float lat, float lon);

int daynight_outlines(int32_t *year, int32_t *dayOfYear,
        int32_t *msecondOfDay, int32_t wid, int32_t hgt,
        float **lat, float **lon, int32_t n[2],
        float *olat[2], float *olon[2], int8_t **dorn);



#endif /* SRC_L2BIN_DATADAY_H_ */
