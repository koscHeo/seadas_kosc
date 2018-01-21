/****************************************************************************
 * NCSA HDF                                                                 *
 * Software Development Group                                               *
 * National Center for Supercomputing Applications                          *
 * University of Illinois at Urbana-Champaign                               *
 * 605 E. Springfield, Champaign IL 61820                                   *
 *                                                                          *
 * For conditions of distribution and use, see the accompanying             *
 * hdf/COPYING file.                                                        *
 *                                                                          *
 ****************************************************************************/

/* tutils.h,v 1.12 1995/10/05 16:17:28 koziol Exp */

#ifndef _TUTILS_H
#define _TUTILS_H

/* Define these for use in all the tests */
#ifndef TESTMASTER
extern
#endif
int         num_errs
#ifdef TESTMASTER
= 0
#endif
,           Verbosity
#ifdef TESTMASTER
= 0
#endif
           ;

#ifdef TEST_PC
#define FAR far
#else
#ifndef FAR
#define FAR     /* */
#endif /* FAR */
#endif /* TEST_PC */

/* Use %ld to print the value because long could cover most cases. */
/* Used to make certain a return value _is_not_ a value */
#define CHECK(ret, val, where) \
do {if (Verbosity>9) printf("   Call to HDF routine: %15s at line %4d in %s returned %ld \n",where,(int)__LINE__,__FILE__,(long)ret);\
if(ret == val) {printf("*** UNEXPECTED RETURN from %s is %ld at line %4d in %s\n", where, (long)ret, (int)__LINE__,__FILE__); num_errs++;} \
} while(0)

/* Used to make certain a return value _is_ a value */
#define VERIFY(x, val, where) \
do {if (Verbosity>9) printf("   Call to HDF routine: %15s at line %4d in %s had value %ld \n",where,(int)__LINE__,__FILE__,(long)x);\
if(x != val) {printf("*** UNEXPECTED VALUE from %s is %ld at line %4d in %s\n", where, (long)x,(int)__LINE__,__FILE__); num_errs++;} \
} while(0)

#define RESULT(a) \
do { \
if (Verbosity>8) printf("   Call to HDF routine: %15s at line %4d in %s returned %ld \n",a,(int)__LINE__,__FILE__,(long)ret); \
if (Verbosity>9) HEprint(stdout,0); \
if(ret == FAIL) {printf("*** UNEXPECTED RETURN from %s is %ld at line %4d in %s\n", a, (long)ret,(int)__LINE__,__FILE__); num_errs++;} \
} while(0)

#define MESSAGE(v,a) {if (Verbosity>=v) {a}}

#endif /* _TUTILS_H */
