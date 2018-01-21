/*---------------------------------------------------------------------------------  
 |        Copyright (C) 1999  Emergent IT Inc.  and Raytheon Systems Company      |    
 |  Permission to use, modify, and distribute this software and its documentation |
 |  for any purpose without fee is hereby granted, provided that the above        | 
 |  copyright notice appear in all copies and that both that copyright notice and |
 |  this permission notice appear in supporting documentation.                    |
 |                                                                                |
 |  BEGIN_FILE_PROLOG:                                                            |
 |  FILENAME:                                                                     |
 |  HE5_GCTPFUNC.H                                                                |
 |  DESCRIPTION:                                                                  |
 |	This file contains function prototypes that are specific to the GCTP          |
 |  AUTHOR:                                                                       |
 |	Alex Muslimov / Emergent Information Tecnologies, Inc.                        |
 |  HISTORY:                                                                      |
 |	14-May-01 AM Initial version                                                  |
 | END_FILE_PROLOG:                                                               |
 ----------------------------------------------------------------------------------
 */

#ifndef HE5_GctpFunc_h
#define HE5_GctpFunc_h

/*****************************************************************
    Function prototypes.
*****************************************************************/

int stplnfor(double lon, double lat, double *x, double *y);
int stplninv(double x, double y, double *lon, double *lat);
int stplnforint(long zone, long sphere, char *fn27, char *fn83);
int stplninvint(long zone, long sphere, char *fn27, char *fn83);
int alberfor(double lon, double lat, double *x, double *y);
int alberinv(double x, double y, double *lon, double *lat);
int alberforint(double r_maj, double r_min, double lat1, double lat2, 
                double lon0, double lat0, double false_east, 
                double false_north);
int alberinvint(double r_maj, double r_min, double lat1, double lat2, 
                double lon0, double lat0, double false_east, 
                double false_north);
int lamccfor(double lon, double lat, double *x, double *y);
int lamccinv(double x, double y, double *lon, double *lat);
int lamccforint(double r_maj, double r_min, double lat1, double lat2, 
                double c_lon, double c_lat, double false_east, 
                double false_north);
int lamccinvint(double r_maj, double r_min, double lat1, double lat2, 
                double c_lon, double c_lat, double false_east, 
                double false_north);
int merfor(double lon, double lat, double *x, double *y);
int merinv(double x, double y, double *lon, double *lat);
int merforint(double r_maj, double r_min, double center_lon, 
                double center_lat, double false_east, double false_north);
int merinvint(double r_maj, double r_min, double center_lon, 
                double center_lat, double false_east, double false_north);
int psfor(double lon, double lat, double *x, double *y);
int psinv(double x, double y, double *lon, double *lat);
int psforint(double r_maj, double r_min, double c_lon, double c_lat, 
                double false_east, double false_north);
int psinvint(double r_maj, double r_min, double c_lon, double c_lat, 
                double false_east, double false_north);
int polyfor(double lon, double lat, double *x, double *y);
int polyinv(double x, double y, double *lon, double *lat);
int polyforint(double r_maj, double r_min, double center_lon, 
                double center_lat, double false_east, double false_north);
int polyinvint(double r_maj, double r_min, double center_lon, 
                double center_lat, double false_east, double false_north);
int eqconfor(double lon, double lat, double *x, double *y);
int eqconinv(double x, double y, double *lon, double *lat);
int eqconforint(double r_maj, double r_min, double lat1, double lat2, 
                double center_lon, double center_lat, double false_east, 
                double false_north, long mode);
int eqconinvint(double r_maj, double r_min, double lat1, double lat2, 
                double center_lon, double center_lat, double false_east, 
                double false_north, long mode);
int tmfor(double lon, double lat, double *x, double *y);
int tminv(double x, double y, double *lon, double *lat);
int tmforint(double r_maj, double r_min, double scale_fact, 
                double center_lon, double center_lat, double false_east, 
                double false_north);
int tminvint(double r_maj, double r_min, double scale_fact, 
                double center_lon, double center_lat, double false_east, 
                double false_north);
int sterfor(double lon, double lat, double *x, double *y);
int sterinv(double x, double y, double *lon, double *lat);
int sterforint(double r_maj, double center_lon, double center_lat, 
                double false_east, double false_north);
int sterinvint(double r_maj, double center_lon, double center_lat, 
                double false_east, double false_north);
int lamazfor(double lon, double lat, double *x, double *y);
int lamazinv(double x, double y, double *lon, double *lat);
int lamazforint(double r, double center_long, double center_lat, 
                double false_east, double false_north);
int lamazinvint(double r, double center_long, double center_lat, 
                double false_east, double false_north);
int azimfor(double lon, double lat, double *x, double *y);
int aziminv(double x, double y, double *lon, double *lat);
int azimforint(double r_maj, double center_lon, double center_lat, 
                double false_east, double false_north);
int aziminvint(double r_maj, double center_lon, double center_lat, 
                double false_east, double false_north);
int gnomfor(double lon, double lat, double *x, double *y);
int gnominv(double x, double y, double *lon, double *lat);
int gnomforint(double r, double center_long, double center_lat, 
                double false_east, double false_north);
int gnominvint(double r, double center_long, double center_lat, 
                double false_east, double false_north);
int orthfor(double lon, double lat, double *x, double *y);
int orthinv(double x, double y, double *lon, double *lat);
int orthforint(double r_maj, double center_lon, double center_lat, 
                double false_east, double false_north);
int orthinvint(double r_maj, double center_lon, double center_lat, 
                double false_east, double false_north);
int gvnspfor(double lon, double lat, double *x, double *y);
int gvnspinv(double x, double y, double *lon, double *lat);
int gvnspforint(double r, double h, double center_long, double center_lat, 
                double false_east, double false_north);
int gvnspinvint(double r, double h, double center_long, double center_lat, 
                double false_east, double false_north);
int sinfor(double lon, double lat, double *x, double *y);
int sininv(double x, double y, double *lon, double *lat);
int sinforint(double r, double center_long, double false_east, 
                double false_north);
int sininvint(double r, double center_long, double false_east, 
                double false_north);
int equifor(double lon, double lat, double *x, double *y);
int equiinv(double x, double y, double *lon, double *lat);
int equiforint(double r_maj, double center_lon, double lat1, 
                double false_east, double false_north);
int equiinvint(double r_maj, double center_lon, double lat1, 
                double false_east, double false_north);
int millfor(double lon, double lat, double *x, double *y);
int millinv(double x, double y, double *lon, double *lat);
int millforint(double r, double center_long, double false_east, 
                double false_north);
int millinvint(double r, double center_long, double false_east, 
                double false_north);
int vandgfor(double lon, double lat, double *x, double *y);
int vandginv(double x, double y, double *lon, double *lat);
int vandgforint(double r, double center_long, double false_east, 
                double false_north);
int vandginvint(double r, double center_long, double false_east, 
                double false_north);
int omerfor(double lon, double lat, double *x, double *y);
int omerinv(double x, double y, double *lon, double *lat);
int omerforint(double r_maj, double r_min, double scale_fact, 
                double azimuth, double lon_orig, double lat_orig, 
                double false_east, double false_north, double lon1, 
                double lat1, double lon2, double lat2, long mode);
int omerinvint(double r_maj, double r_min, double scale_fact, 
                double azimuth, double lon_orig, double lat_orig, 
                double false_east, double false_north, double lon1, 
                double lat1, double lon2, double lat2, long mode);
int somfor(double lon, double lat, double *x, double *y);
int sominv(double x, double y, double *lon, double *lat);
int somforint(double r_major, double r_minor, long satnum, long path, 
                double alf_in, double lon, double false_east, 
                double false_north, double time, long start1, long flag, 
                double sat_ratio);
int sominvint(double r_major, double r_minor, long satnum, long path, 
                double alf_in, double lon, double false_east, 
                double false_north, double time, long start1, long flag, 
                double sat_ratio);
int hamfor(double lon, double lat, double *x, double *y);
int haminv(double x, double y, double *lon, double *lat);
int hamforint(double r, double center_long, double false_east, 
                double false_north);
int haminvint(double r, double center_long, double false_east, 
                double false_north);
int robfor(double lon, double lat, double *x, double *y);
int robinv(double x, double y, double *lon, double *lat);
int robforint(double r, double center_long, double false_east, 
                double false_north);
int robinvint(double r, double center_long, double false_east, 
                double false_north);
int goodfor(double lon, double lat, double *x, double *y);
int goodinv(double x, double y, double *lon, double *lat);
int goodforint(double r);
int goodinvint(double r);
int molwfor(double lon, double lat, double *x, double *y);
int molwinv(double x, double y, double *lon, double *lat);
int molwforint(double r, double center_long, double false_east, 
                double false_north);
int molwinvint(double r, double center_long, double false_east, 
                double false_north);
int imolwfor(double lon, double lat, double *x, double *y);
int imolwinv(double x, double y, double *lon, double *lat);
int imolwforint(double r);
int imolwinvint(double r);
int alconfor(double lon, double lat, double *x, double *y);
int alconinv(double x, double y, double *lon, double *lat);
int alconforint(double r_maj, double r_min, double false_east, 
                double false_north);
int alconinvint(double r_maj, double r_min, double false_east, 
                double false_north);
int wivfor(double lon, double lat, double *x, double *y);
int wivinv(double x, double y, double *lon, double *lat);
int wivforint(double r, double center_long, double false_east, 
                double false_north);
int wivinvint(double r, double center_long, double false_east, 
                double false_north);
int wviifor(double lon, double lat, double *x, double *y);
int wviiinv(double x, double y, double *lon, double *lat);
int wviiforint(double r, double center_long, double false_east, 
                double false_north);
int wviiinvint(double r, double center_long, double false_east, 
                double false_north);
int obleqfor(double lon, double lat, double *x, double *y);
int obleqinv(double x, double y, double *lon, double *lat);
int obleqforint(double r, double center_long, double center_lat, 
                double shape_m, double shape_n, double angle, 
                double false_east, double false_north);
int obleqinvint(double r, double center_long, double center_lat, 
                double shape_m, double shape_n, double angle, 
                double false_east, double false_north);
int isinusfor(double lon, double lat, double *x, double *y);
int isinusinv(double x, double y, double *lon, double *lat);
int isinusforinit(double sphere, double lon_cen_mer, double false_east,
                   double false_north, double dzone, double djustify);
int isinusinvinit(double sphere, double lon_cen_mer, double false_east,
                   double false_north, double dzone, double djustify);
int utmfor(double lon, double lat, double *x, double *y);
int utminv(double x, double y, double *lon, double *lat);
int utmforint(double r_maj, double r_min, double scale_fact, long zone);
int utminvint(double r_maj, double r_min, double scale_fact, long zone);
long calc_utm_zone(double lon);
int bceaforint(double r_maj, double r_min, double center_lon, 
               double center_lat, double false_east, double false_north);
int bceafor(double lon, double lat, double *x, double *y);
int bceainvint(double r_maj, double r_min, double center_lon, 
               double center_lat, double false_east, double false_north);
int bceainv(double x, double y, double *lon, double *lat);
void p_error(char *what, char *where);
void ptitle(char *A);
void tsincos(double val, double *sin_val, double *cos_val);
double msfnz(double eccent, double sinphi, double cosphi);
double qsfnz(double eccent, double sinphi, double cosphi);
double tsfnz(double eccent, double phi, double sinphi);
void radius2(double A, double B);
void radius(double A);
void stanparl(double A, double B);
void cenlonmer(double A);
void cenlon(double A);
void cenlat(double A);
void origin(double A);
void offsetp(double A, double B);
double adjust_lon(double x);
double phi1z(double eccent, double qs, long  *flag);
double phi2z(double eccent, double ts, long *flag);
double phi3z(double ml, double e0, double e1, double e2, 
             double e3, long *flag);
double phi4z(double eccent, double e0, double e1, double e2, 
             double e3, double a, double b, double *c, double *phi);
double asinz(double con);
int sign(double x);
double e0fn(double x);
double e1fn(double x);
double e2fn(double x);
double e3fn(double x);
double e4fn(double x);
double mlfn(double e0, double e1, double e2, double e3, double phi);
double paksz(double ang, long *iflg);
double pakcz(double pak);
void stparl1(double A);
void genrpt(double A, char *S);
void genrpt_long(long A, char *S);
void pblank();

#if defined(SGI64)
void for_init(int a, int b, double c[], int d, char *str1, char *str2, int *e, int (**func)(double, double, double *, double *)); 

void inv_init(int a,int b,double *c,int d,char *str1,char *str2,int *e,int (**func)(double, double, double*, double*));
#else
void for_init(long a, long b, double c[], long d, char *str1, char *str2, long *e, long (**func)(double, double, double *, double *)); 

void inv_init(long a,long b,double *c,long d,char *str1,char *str2,long *e,long (**func)(double, double, double*, double*));
#endif


#endif


