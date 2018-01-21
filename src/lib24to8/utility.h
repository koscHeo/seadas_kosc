#ifndef vms
#ifndef NoDIRENT
#include <dirent.h>
#else
#include <sys/dir.h>
#ifndef dirent
#define dirent direct
#endif
#endif
#endif
#ifndef S_ISDIR
#define S_ISDIR(mode) (((mode) & S_IFMT) == S_IFDIR)
#endif

#include <stdint.h>

/*
  Utility define declarations.
*/
#define DirectorySeparator  "/"

/*
  Utilities routines.
*/
extern char
  **StringToList _Declare((char *)),
  **ListFiles _Declare((char *,char *,int *));

extern int
  GlobExpression _Declare((char *,char *)),
  ReadDataBlock _Declare((char *,FILE *));

extern unsigned int
  ReadData _Declare((char *,int,int,FILE *));

extern uint32_t 
  LSBFirstReadLong _Declare((FILE *)),
  MSBFirstReadLong _Declare((FILE *));

extern unsigned short
  LSBFirstReadShort _Declare((FILE *)),
  MSBFirstReadShort _Declare((FILE *));

extern void
  LSBFirstWriteLong _Declare((uint32_t,FILE *)),
  LSBFirstWriteShort _Declare((unsigned int,FILE *)),
  MSBFirstOrderLong _Declare((char *,unsigned int)),
  MSBFirstOrderShort _Declare((char *,unsigned int)),
  MSBFirstWriteLong _Declare((uint32_t,FILE *)),
  MSBFirstWriteShort _Declare((unsigned int,FILE *));
