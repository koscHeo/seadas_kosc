#ifndef _XCAL_H
#define _XCAL_H

#define XTNTIME 300
#define XTNDET   40
#define XTNMSIDE  2
#define XTNORDER  6
#define XRVS      0
#define XM12      1
#define XM13      2

double *xcal_modis(l1str *l1rec, int type, int bandnum);

#endif
