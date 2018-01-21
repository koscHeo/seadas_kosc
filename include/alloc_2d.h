#ifndef _ALLOC_2D_H
#define _ALLOC_2D_H

void free2d_float(float **p);
float **alloc2d_float(int w, int h);

void free2d_short(short **p);
float **alloc2d_short(int w, int h);

void free2d_char(char **p);
char **alloc2d_char(int w, int h);

#endif



