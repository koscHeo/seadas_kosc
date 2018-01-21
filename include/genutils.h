#ifndef _GENUTILS_H
#define _GENUTILS_H

#include <timeutils.h>
#include <stdio.h>


#ifdef __cplusplus
extern "C" {
#endif

/**\file 
   This is a bunch of useful general utilities.
*/

/** Define a few useful BAD values
*/
#define BAD_FLT  -32767.0
#define BAD_INT    -32767
#define BAD_UINT    65535
#define BAD_BYTE     -128
#define BAD_UBYTE     255

/** Do we want to print out extra info to stdout */
extern int want_verbose;

/** Is this machine little endian
    \return 0 for BIG_ENDIAN machines, 1 for LITTLE_ENDIAN mahcines 
*/
int endianess(void);

void parse_file_name(const char *inpath, char *outpath);

/** Swap bytes in place 
    \param in pointer to memory to byte swap
    \param nbyte size of object to reverse
    \param ntime number of objects to swap
    \return 0 always
*/
int swapc_bytes(char *in, int nbyte, int ntime);

/** Swap bytes from in to out.  in and out should not be overlapping memory
    \param in pointer to source memory to byte swap
    \param out pointer to destination memory
    \param nbyte size of object to reverse
    \param ntime number of objects to swap
    \return 0 always
*/
int swapc_bytes2(const char *in, char *out, int nbyte, int ntime);

int isleap(int year);

/** read a binary file and swap bytes if necessary
    \param little_endian set to 0 for big endian files, 1 for little endian files
    \param ptr memory to read into
    \param size size of object to read in bytes
    \param nmemb number of objects to read
    \param stream file to read from
    \return number of objects read
*/
size_t fread_swap(int little_endian, void *ptr, size_t size, size_t nmemb, FILE *stream);

/** write a binary file and swap bytes if necessary
    \param little_endian set to 0 for big endian files, 1 for little endian files
    \param ptr memory to write
    \param size size of object to write in bytes
    \param nmemb number of objects to write
    \param stream file to write to
    \return number of objects written
*/
size_t fwrite_swap(int little_endian, const  void  *ptr,  size_t size, size_t nmemb, FILE *stream);

void spline(float [],float [],int,float,float,float []);
void splint(float [],float [],float [],int,float,float *);
void lspline(float xin [],float yin [],int32_t nin,
             float xout[],float yout[],int32_t nout);
float linterp(float xin [],float yin [],int32_t nin,float xout);
float bioBandShift(float xin [],float yin [],int32_t nin,float xout);

int32_t filesize_(const char *filename);
int32_t filesize(const char *filename);

int getlut_file(char *lut_file, short *rlut, short *glut, short * blut);

char *lowcase(char *instr);
char *upcase(char *instr);

int isValidInt(const char* str);

void* allocateMemory(size_t numBytes, const char* name);

void trimBlanks(char* str);
char* trimBlanksDup(const char* str);

int getFileFormatIndex(const char* str);
const char* getFileFormatName(const char* str);
const char* getFileFormatExtension(const char* str);


#ifdef __cplusplus
}
#endif

#endif
