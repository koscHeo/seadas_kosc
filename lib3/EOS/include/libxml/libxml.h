/*
 * libxml.h: internal header only used during the compilation of libxml
 *
 * See COPYRIGHT for the status of this software
 *
 * Author: breese@users.sourceforge.net
 */

#ifndef __XML_LIBXML_H__
#define __XML_LIBXML_H__



/*
## ----------- ##
## confdefs.h. ##
## ----------- ##
*/

#define PACKAGE_NAME ""
#define PACKAGE_TARNAME ""
#define PACKAGE_VERSION ""
#define PACKAGE_STRING ""
#define PACKAGE_BUGREPORT ""
#define PACKAGE "libxml2"
#define VERSION "2.7.3"
#define PROTOTYPES 1
#define __PROTOTYPES 1
#define STDC_HEADERS 1
#define HAVE_SYS_TYPES_H 1
#define HAVE_SYS_STAT_H 1
#define HAVE_STDLIB_H 1
#define HAVE_STRING_H 1
#define HAVE_MEMORY_H 1
#define HAVE_STRINGS_H 1
#define HAVE_INTTYPES_H 1
#define HAVE_STDINT_H 1
#define HAVE_UNISTD_H 1
#define HAVE_STRING_H 1
#define HAVE_DLFCN_H 1
#define HAVE_ZLIB_H 1
#define HAVE_LIBZ 1
#define HAVE_DIRENT_H 1
#define STDC_HEADERS 1
#define HAVE_FCNTL_H 1
#define HAVE_UNISTD_H 1
#define HAVE_CTYPE_H 1
#define HAVE_DIRENT_H 1
#define HAVE_ERRNO_H 1
#define HAVE_MALLOC_H 1
#define HAVE_STDARG_H 1
#define HAVE_SYS_STAT_H 1
#define HAVE_SYS_TYPES_H 1
#define HAVE_STDINT_H 1
#define HAVE_TIME_H 1
#define HAVE_MATH_H 1
#define HAVE_LIMITS_H 1
#define HAVE_FLOAT_H 1
#define HAVE_STDLIB_H 1
#define HAVE_SYS_SOCKET_H 1
#define HAVE_NETINET_IN_H 1
#define HAVE_ARPA_INET_H 1
#define HAVE_NETDB_H 1
#define HAVE_SYS_TIME_H 1
#define HAVE_SYS_SELECT_H 1
#define HAVE_SYS_MMAN_H 1
#define HAVE_SYS_TIMEB_H 1
#define HAVE_SIGNAL_H 1
#define HAVE_ARPA_NAMESER_H 1
#define HAVE_RESOLV_H 1
#define HAVE_DLFCN_H 1
#define HAVE_STRFTIME 1
#define HAVE_STRDUP 1
#define HAVE_STRNDUP 1
#define HAVE_STRERROR 1
#define HAVE_FINITE 1
#define HAVE_STRFTIME 1
#define HAVE_LOCALTIME 1
#define HAVE_GETTIMEOFDAY 1
#define HAVE_FTIME 1
#define HAVE_STAT 1
#define HAVE_SIGNAL 1
#define HAVE_PRINTF 1
#define HAVE_SPRINTF 1
#define HAVE_FPRINTF 1
#define HAVE_SNPRINTF 1
#define HAVE_VFPRINTF 1
#define HAVE_VSPRINTF 1
#define HAVE_VSNPRINTF 1
#define HAVE_SSCANF 1
#define HAVE_VA_COPY 1
#define XML_SOCKLEN_T socklen_t
#define SUPPORT_IP6
#define HAVE_GETADDRINFO
#define HAVE_ISNAN
#define HAVE_ISINF
#define HAVE_DLOPEN
#define HAVE_LIBPTHREAD
#define HAVE_PTHREAD_H
#define ICONV_CONST 
#define LIBXML_THREAD_ENABLED




#ifndef NO_LARGEFILE_SOURCE
#ifndef _LARGEFILE_SOURCE
#define _LARGEFILE_SOURCE
#endif
#ifndef _FILE_OFFSET_BITS
#define _FILE_OFFSET_BITS 64
#endif
#endif

#if defined(macintosh)
#include "config-mac.h"
#elif defined(_WIN32_WCE)
/*
 * Windows CE compatibility definitions and functions
 * This is needed to compile libxml2 for Windows CE.
 * At least I tested it with WinCE 5.0 for Emulator and WinCE 4.2/SH4 target
 */
#include <win32config.h>
#include <libxml/xmlversion.h>
#else
#include "config.h"
#include <libxml/xmlversion.h>
#endif

#if defined(__Lynx__)
#include <stdio.h> /* pull definition of size_t */
#include <varargs.h>
int snprintf(char *, size_t, const char *, ...);
int vfprintf(FILE *, const char *, va_list);
#endif

#ifndef WITH_TRIO
#include <stdio.h>
#else
/**
 * TRIO_REPLACE_STDIO:
 *
 * This macro is defined if teh trio string formatting functions are to
 * be used instead of the default stdio ones.
 */
#define TRIO_REPLACE_STDIO
#include "trio.h"
#endif

/*
 * Internal variable indicating if a callback has been registered for
 * node creation/destruction. It avoids spending a lot of time in locking
 * function while checking if the callback exists.
 */
extern int __xmlRegisterCallbacks;
/*
 * internal error reporting routines, shared but not partof the API.
 */
void __xmlIOErr(int domain, int code, const char *extra);
void __xmlLoaderErr(void *ctx, const char *msg, const char *filename);
#ifdef LIBXML_HTML_ENABLED
/*
 * internal function of HTML parser needed for xmlParseInNodeContext
 * but not part of the API
 */
void __htmlParseContent(void *ctx);
#endif

/*
 * internal global initialization critical section routines.
 */
void __xmlGlobalInitMutexLock(void);
void __xmlGlobalInitMutexUnlock(void);
void __xmlGlobalInitMutexDestroy(void);

#ifdef IN_LIBXML
#ifdef __GNUC__
#ifdef PIC
#ifdef linux
#if (__GNUC__ == 3 && __GNUC_MINOR__ >= 3) || (__GNUC__ > 3)
#include "elfgcchack.h"
#endif
#endif
#endif
#endif
#endif
#endif /* ! __XML_LIBXML_H__ */
