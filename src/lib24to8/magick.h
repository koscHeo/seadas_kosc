/*
  Include declarations
*/
#include <stdio.h>
#if defined(__STDC__) || defined(sgi) || defined(_AIX)
#include <stdlib.h>
#include <unistd.h>
#else
#ifndef vms
#include <malloc.h>
#include <memory.h>
#endif
#endif
#include <ctype.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/stat.h>
#undef index

/*
  Define declarations for the Display program.
*/
#if __STDC__ || defined(sgi) || defined(_AIX)
#define _Declare(formal_parameters) formal_parameters
#else
#define const 
#define _Declare(formal_parameters) ()
#endif
#define DownShift(x) ((int) ((x)+(1L << 15)) >> 16)
#define False  0
#define Max(x,y)  (((x) > (y)) ? (x) : (y))
#define Min(x,y)  (((x) < (y)) ? (x) : (y))
#define True  1
#define UpShift(x) ((x) << 16)
#define UpShifted(x) ((int) ((x)*(1L << 16)+0.5))
#ifdef vms
#define pclose(file)  exit(1)
#define popen(command,mode)  exit(1)
#define unlink(file)  remove(file)
#endif

/*
  Typedef declarations.
*/
typedef void
  (*ErrorHandler) _Declare((char *,char *));

/*
  Variable declarations.
*/
static char 
  Version[]="@(#)ImageMagick 3.0 94/05/01 cristy@dupont.com";
