/*
 *  W. Robinson, SAIC, 10 Dec 2004  new for CZCS
 */
#ifndef L1_CZCSW_H
#define L1_CZCS_H

int openl1_czcs(filehandle *file);
int readl1_czcs(filehandle *file, int32_t recnum, l1str *l1rec);
int closel1_czcs(filehandle *file);
int cz_posll_2_satang( float *, int, float *, float *, float *, float * );
void matrix_mult( double[3], double[3][3], double[3] );
void cross_prod( double *, double *, double * );

#endif
