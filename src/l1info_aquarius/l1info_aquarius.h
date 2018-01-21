#include "FileStager.h"
#include "hdf5_Aquarius.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_linalg.h>

typedef struct input_struct {
  char ifile[FILENAME_MAX];
  char pfile[FILENAME_MAX];

} instr;


// subroutine prototypes
int getGeoNavSun( Hdf::hdf5_Aquarius *l1afile, 
		  float *cellon, float *cellat, float *celtht, float *celphi,
		  float *suntht, float *sunphi, float *gkxkib, float *glxlat,
		  double *zang);

// Update for rpy_adj  07/05/12  JMG
extern "C" void geolocation_(double *, double *, double *, double *, double *, 
			     float *, float *, float *, float *, float *, 
			     float *, float *, float *, float *, float *,
			     double *, 
			     double *, double *, double *, double *,
			     float *, float *,
			     double *, double *, double *, double *, double *);
