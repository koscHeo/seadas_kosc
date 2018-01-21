/* 
 * File:   ancnrt_proto.h
 * Author: dshea
 *
 * Created on February 16, 2011, 11:49 AM
 */

#ifndef ANCNRT_PROTO_H
#define	ANCNRT_PROTO_H

#ifdef	__cplusplus
extern "C" {
#endif

/* fortran  prototypes */
void rdgrid_();
void gregor_(int *julday, int *year, int *month, int *day);

int readgrib2( char *file, int npix, int nlin, int rgmode,
    float *data );
int grib2_t( char *grib2_t_str, int *year, int *doy, int *hour, int *npix,
    int *nlin, int *h_fcst );
int readgrib2_3d( char *file, char *grib2_t_str, int npix, int nlin,
    float *data, int nprfl, int *year, int *month, int *day, int *hour );
int count_annot(char *filename);
void shift_180( float *in, int npix, int nlin, float fact, float *out );


/* prototypes from ANCroutines.c */
int startHDF(char *outfile, int32 *sdfid, int32 *fid, int32 mode);
int32 setupGrid(int32 fid, char *grpname);
int32 gridToGrid(int32 outergridid, int32 innergridid);
int32 writeGeom(int32 fid, int32 gridid, char *geomname, int32 bin_meth,
                int32 registration, float32 vsize, float32 hsize,
                float32 max_north, float32 max_south, float32 max_west,
                float32 max_east);
int32 findGeomId(int32 fid, char *geomname);
int32 linkGeom(int32 gridid, int32 geomid);
int32 detachGeom(int32 geomid);
int addAttr(int32 sdsid, char *dataattr, int32 datatype, char *dataunit);
int setSDSref(int32 sdsid, int32 gridid);
int deattachHDFgrid(int32 gridid);
int closeHDFstructs(int32 sdfid, int32 fid);
int32 wrtsds(int32 sdfid, int rank, int32 *shape, int32 datatype,
             char *datalabel, void *data);
int32 rewrtsds(int32 sdsid, int32 *shape, void *data);
int rdsds(char *filename, char *vgname, char *sdsname, int32 *dimsizes,
          void *inData);
int wrtattr(int32 dfile, struct annotation *annot, int numannarr);
int32 SDSinFile (char *sdsname, char *longname, char *units, char *datafmt,
    int32 datatype, int32 sdfid, int32 rank, int32 *shape, void *data, int32 gridid);


void pexit (char *string);
int pwarning (char *string);


#ifdef	__cplusplus
}
#endif

#endif	/* ANCNRT_PROTO_H */

