#ifndef  NETCDFUTILS_H
#define  NETCDFUTILS_H

#ifdef __cplusplus
extern "C" {
#endif

void report_err(const int stat, const int line, const char *file);
void check_err(const int stat, const int line, const char *file);

int writeBinList_nc( int32_t deflate, int32_t grpid, int32_t nbins_to_write, 
		     const void *data);
int writeBinData_nc( int32_t deflate, int32_t grpid, int32_t nbins_to_write,
		     int32_t iprod, const char *prodname, const void *data);
int writeBinIndex_nc( int32_t grpid, int32_t nbins_to_write, const void *data);
int writeQuality_nc( int32_t deflate, int32_t grpid, int32_t nbins_to_write, 
        const void *data);



#ifdef __cplusplus
}
#endif

#endif
