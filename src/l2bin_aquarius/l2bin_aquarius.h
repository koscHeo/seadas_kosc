#include "hdf5_Aquarius.h"
//#include <gsl/gsl_interp.h>

#define PRODSTRLEN     2048  /* String length limit for product spec */

#define MAXNFILES 525

//#define TRUE 1

#define RFI            0
#define RAIN           2
#define LAND           3
#define ICE            4
#define WIND           5
#define TEMP           6
#define FLUXD          7
#define FLUXR          8
#define SUNGLINT       9
#define MOON          10
#define GALACTIC      11
#define NAV           12
#define SAOVERFLOW    13
#define ROUGH         14
#define FLARE         15
#define POINTING      16
#define TBCONS        17
#define COLDWATER     18
#define TFTADIFF      19
#define SPARE1        20
#define REFL_1STOKES  21
#define SPARE2        22
#define RFI_REGION    23

typedef struct input_struct {
  char infile[FILENAME_MAX];
  char ofile[FILENAME_MAX];
  char pfile[FILENAME_MAX];

  char    flaguse[1024];
  char    l3bprod[PRODSTRLEN];
  //  int64_t l3bprodword;

  char    resolve[16];
  int32_t    beam;

  char    fit[16];
  float   filter_width;
  int32_t    fill;

  char pversion[64];

  char orbit_type[4];

  bool require;
} instr;

int parseInput( int argc, char* argv[], instr *l2binInput, string *procControl);

inline
int extParmWordValue( string sLine, string *sParmWord, string *sValue) {

  string::size_type posBeginIdx, posEndIdx;
  string::size_type ipos=0;
  const string      sDelim( "=" );

  // Extract parameter word
  posEndIdx = sLine.find_first_of( sDelim );
  *sParmWord = sLine.substr( ipos, posEndIdx );
  posBeginIdx = posEndIdx + 1;  // Beginning of next word (after '=')

  // Convert to uppercase
  for (size_t j=0; j<(*sParmWord).length(); j++)
    (*sParmWord)[j] = toupper((*sParmWord)[j]);

  // Extract parameter value
  *sValue  = sLine.substr( posBeginIdx);

  return 0;
}


inline
int expandEnvVar( string *sValue) {
  if ( (*sValue).find_first_of( "$" ) == string::npos) return 0;
  string::size_type posEndIdx = (*sValue).find_first_of( "/" );
  if ( posEndIdx == string::npos) return 0;
  char *envVar_str = getenv((*sValue).substr( 1, posEndIdx-1 ).c_str());
  if (envVar_str == 0x0) {
    printf("Environment variable: %s not defined.\n", envVar_str);
    exit(1);
  }
  *sValue = envVar_str + (*sValue).substr( posEndIdx);

  return 0;
}


