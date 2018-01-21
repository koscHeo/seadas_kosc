#include <stdio.h>
#include <math.h>

#include <timeutils.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <list>
#include <stack>

#include "hdf_bin.h"
#include "hdf5utils.h"
#include "passthebuck.h"
#include "l2bin_aquarius.h"

#include <gsl/gsl_multifit.h>

#define MISSIONCHARACTERISTICS "Nominal orbit: inclination=98.0 (Sun-synchronous); node=6PM (ascending); eccentricity=<0.002; altitude=657 km; ground speed=6.825 km/sec"

#define SENSORCHARACTERISTICS "Number of beams=3; channels per receiver=4; frequency 1.413 GHz; bits per sample=16; instatntaneous field of view=6.5 degrees; science data block period=1.44 sec"

#define NUMBER_OF_VALUES 10000

#define PI      3.141592653589793
#define VERSION "2.06"

//////////////////////////////////////////////////////////////////////////////

//  Revision 2.06  05/31/15
//  Add Start Orbit, End Orbit metadata
//  J. Gales

//  Revision 2.05  05/27/15
//  Add Sensor Characteristics, L2 Flag Names metadata
//  J. Gales

//  Revision 2.04  08/07/14
//  Revert back to original metadata
//  J. Gales

//  Revision 2.03  04/04/14
//  Add missing H5Tclose() statements to fix memory leaks
//  J. Gales

//  Revision 2.02  03/31/14
//  Exit if missing values are binned
//  J. Gales

//  Revision 2.01  01/31/14
//  Modify galactic REFL_1STOKES test
//  J. Gales

//  Revision 2.00  01/10/14
//  Add support for POINTING,TBCONS,COLDWATER,TFTADIFF,
//                  REFL_1STOKES,RFI_REGION flags
//  J. Gales

//  Revision 1.12  10/21/13
//  Add "Start Time" and "End Time" to metadata
//  J. Gales

//  Revision 1.11  09/20/13
//  Change SM Dense Vegetation flag check to 512
//  J. Gales

//  Revision 1.10  08/20/13
//  Initialize l2binInput.require to false
//  J. Gales

//  Revision 1.09  07/01/13
//  Change attribute name to Easternmost/Westernmost Longitude
//  J. Gales

//  Revision 1.08  04/30/13
//  Modifiy radiometer flags read for SM
//  J. Gales

//  Revision 1.07  04/08/13
//  Add support for additional flags (RAIN not yet supported)
//  Add support for flag binfiles (required parameter)
//  J. Gales

//  Revision 1.06  03/13/13
//  Bin soil moisture pixels with dense vegatation flag set
//  J. Gales

//  Revision 1.05  01/30/13
//  Add support for orbit type (A/D)
//  J. Gales

//  Revision 1.04  07/19/12
//  Expand input parameters in process control when using parfile
//  J. Gales

//  Revision 1.03  05/15/12
//  Add support for soil moisture
//  J. Gales

//  Revision 1.02  04/03/12
//  Add support for flag repackaging
//  J. Gales

//  Revision 1.01  03/03/12
//  Add support for NAV,ROUGH flagging
//  Add support for default products
//  Add support for PAR, SUITE command line parameters
//  Add check for cellon/cellat = -999
//  J. Gales

//  Revision 1.00  09/16/11
//  Add support for TEMP flagging
//  J. Gales

//  Revision 0.56  07/18/11
//  Use hdf_bin.h rather than hdf5_bin.h
//  J. Gales

//  Revision 0.55  05/16/11
//  Change H5T_NATIVE_B32 to H5T_STD_U32LE
//  J. Gales
//
//  Revision 0.52  04/11/11
//  Exit (early) with 1 if no filled bins
//  J. Gales
//
//  Revision 0.51  01/06/11
//  Support for pversion
//  J. Gales
//
//  Revision 0.50  06/03/10
//  Support for V3 flags
//  J. Gales
//
//  Revision 0.42  06/03/10
//  Fix metadata
//  J. Gales
//
//  Revision 0.40  05/25/10
//  Add support for beam_number=-1
//  J. Gales
//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace Hdf;


static int32 *numbin;
static int32 *basebin;
static uint32 nrows;
static uint16 *nobs;
static float32 **data_values;
static float32 **data_longitude;
static float32 **data_latitude;
static int16 **file_index;

int main (int argc, char* argv[])
{
  double blkSec;

  int32_t ibeam, iblk;
  uint32_t ifile;
  uint32_t nfiles;

  instr l2binInput;
  int32_t numBlks[MAXNFILES];
  double granstart[MAXNFILES], granstop[MAXNFILES];

  int32 i32;

  float32 latbin=0.0;

  int32 *bin_flag;
  int16 *nscenes, *lastfile;
  int16 *allocated_space;

  int32 n_allocperbin=10;
  int32 n_total_obs=0;
  int32 total_filled_bins=0;

  float32 data[NUMBER_OF_BEAMS];

  float32 *sum_data;
  float32 *sum2_data;
  float32 f32, wgt, sum, sum2;

  float32 northmost=-90.0, southmost=90.0, eastmost=-180.0, westmost=180.0;

  char **prodname;
  char **flagname;

  char *char_ptr1;

  char* fldname3[2];

  char buf[2048];

  time_t tnow;
  struct tm *tmnow;
  
  uint8 *a, *bin_indx;

  int32 getbinnum(float32, float32);

  FILE *fp;

  cout << "l2bin_aquarius " << VERSION << " (" 
       <<  __DATE__ << " " << __TIME__ << ")" << endl;

  if ( argc == 1) {
    cout << endl << "l2bin_aquarius infile=input_filelist ofile=L3B_output_file resolve=resolution l3bprod=product_list beam=beam_number flaguse=flaguse filter_width=filter_width fill=fill fit=fit pversion=pversion"
	 << endl;
    cout << "where input_filelist is the list of L2 files to bin" << endl;
    cout << "      L3B_output_filelist is the L3B output filename" << endl;
    cout << "      resolve is the bin resolution (typically '1D')" << endl;
    cout << "      product_list is the output products to bin" << endl;
    cout << "      beam_number is the radiometer beam (1-3)" << endl;
    cout << "      pversion is processing version" << endl;

    cout << " " << endl;
    cout << "The following parameters are for smoothing:" << endl;
    cout << " " << endl;
    cout << "      filter_width is the angular radius in degrees of the cutoff filter" << endl;
    cout << "      fill = {0,1} {0 - no fill of undefined bins, 1 - fill undefined bins}" << endl;
    cout << "      fit = {\"C\", \"L\", \"B\", \"Q\"}" << endl;
    cout << "      \"C\" - constant fit  (1 parameter)" << endl;
    cout << "      \"L\" - linear fit    (2 parameters)" << endl;
    cout << "      \"B\" - bilinear fit  (3 parameters)" << endl;
    cout << "      \"Q\" - quadratic fit (4 parameters)" << endl;
    return 0;
  }

  string processControl;
  parseInput(argc, argv, &l2binInput, &processControl);
  static Hdf::hdf5_Aquarius l2file[MAXNFILES];

  if (strcmp(l2binInput.resolve, "9") == 0) nrows = 2160/1;
  if (strcmp(l2binInput.resolve, "18") == 0) nrows = 2160/2;
  if (strcmp(l2binInput.resolve, "36") == 0) nrows = 2160/4;
  if (strcmp(l2binInput.resolve, "46") == 0) nrows = 2160/5;
  if (strcmp(l2binInput.resolve, "72") == 0) nrows = 2160/8;
  if (strcmp(l2binInput.resolve, "1D") == 0) nrows = 2160/12;
  if (strcmp(l2binInput.resolve, "144") == 0) nrows = 2160/16;

  int32_t nOrder = -1;
  bool smooth = true;
  float32 cs_filter;
  if (strcmp(l2binInput.fit, "C") == 0) nOrder = 1;
  if (strcmp(l2binInput.fit, "L") == 0) nOrder = 3;
  if (strcmp(l2binInput.fit, "B") == 0) nOrder = 4;
  if (strcmp(l2binInput.fit, "Q") == 0) nOrder = 6;
  if (nOrder != -1) {
    cs_filter = cos(l2binInput.filter_width*PI/180.);
    printf("nOrder: %d\n", nOrder);
  } else {
    smooth = false;
  }

  int32_t s_orb, e_orb;

  // Single HDF input
  // ----------------
  if (H5Fis_hdf5(l2binInput.infile) == TRUE) {
    nfiles = 1;
    l2file[0].openl2(l2binInput.infile, H5F_ACC_RDONLY, &numBlks[0],
		     &granstart[0], &granstop[0]);
    cout << granstop[0] - granstart[0] << endl;
    cout << "Single HDF input" << endl;

    hid_t gid[6];
    l2file[ifile].getH5gid( l2file[ifile].getH5fid(), gid);
    hid_t attr = H5Aopen_name(gid[0], "Orbit Number");
    H5Aread(attr, H5T_NATIVE_LONG, &s_orb);
    e_orb = s_orb;
    H5Aclose(attr);

  } else {

    // Filelist input - Determine number of input files
    // ------------------------------------------------
    nfiles = 0;
    fp = fopen(l2binInput.infile, "r");
    if (fp == NULL) {
      printf("Input listing file: \"%s\" not found.\n", l2binInput.infile);
      return -1;
    }
    while(fgets(buf, 256, fp) != NULL) nfiles++;
    fclose(fp);
    cout << nfiles << " input files" << endl;

    // Open L2 input files
    // -------------------
    fp = fopen(l2binInput.infile, "r");
    for (ifile=0; ifile<nfiles; ifile++) {

      fgets(buf, 256, fp);
      buf[strlen(buf)-1] = 0;
      l2file[ifile].openl2(buf, H5F_ACC_RDONLY, &numBlks[ifile],
			   &granstart[ifile], &granstop[ifile]);

      if ( ifile == 0) {
        hid_t gid[6];
        l2file[ifile].getH5gid( l2file[ifile].getH5fid(), gid);
        hid_t attr = H5Aopen_name(gid[0], "Orbit Number");
        H5Aread(attr, H5T_NATIVE_LONG, &s_orb);
        H5Aclose(attr);
      }

      if ( ifile == (nfiles-1)) {
        hid_t gid[6];
        l2file[ifile].getH5gid( l2file[ifile].getH5fid(), gid);
        hid_t attr = H5Aopen_name(gid[0], "Orbit Number");
        H5Aread(attr, H5T_NATIVE_LONG, &e_orb);
        H5Aclose(attr);
      }

    } // ifile loop
    fclose(fp);
  }


  // Check Data Type
  hid_t gid[6];
  char datatype[8];
  hid_t atype = H5Tcopy(H5T_C_S1);
  H5Tset_size(atype, 8);
  l2file[0].getH5gid( l2file[0].getH5fid(), gid);
  hid_t attr = H5Aopen_name(gid[0], "Data Type");
  H5Aread(attr, atype, datatype);
  H5Aclose(attr);
  H5Tclose(atype);


  /* Parse L3 Product list */
  /* --------------------- */
  uint32_t len = strlen(l2binInput.l3bprod);
  int32_t nProd = 0;
  for (size_t i=1; i<len; i++) 
    if (l2binInput.l3bprod[i] == ':') nProd++;

  prodname = (char **) calloc(nProd+1, sizeof(char *));

  size_t j = 0;
  for (size_t i=0; i<len; i++) {
    if (l2binInput.l3bprod[i] == ':') {
      prodname[j] = l2binInput.l3bprod + i + 1; 
      l2binInput.l3bprod[i] = 0;
      j++;
    }
  }

  for (size_t j=0; j<nProd; j++) cout << prodname[j] << endl;

  /* Parse FLAGUSE list */
  /* ------------------ */
  int32 len_flag = strlen(l2binInput.flaguse);
  int32 nFlag = 0;
  for (size_t i=1; i<len_flag; i++) 
    if (l2binInput.flaguse[i] == ',') nFlag++;

  flagname = (char **) calloc(nFlag+1, sizeof(char *));

  j = 0;
  for (size_t i=0; i<len_flag; i++) {
    if (l2binInput.flaguse[i] == ',') {
      flagname[j] = l2binInput.flaguse + i + 1; 
      l2binInput.flaguse[i] = 0;
      j++;
    }
  }

  int32 flagtest[4] = {0,0,0,0};
  for (int32_t i=0; i<nFlag; i++) {
    if ( strcmp( flagname[i], "LANDYELLOW") == 0) 
      flagtest[0] |= (uint32_t) pow( 2, LAND); 
    if ( strcmp( flagname[i], "LANDRED") == 0) 
      flagtest[1] |= (uint32_t) pow( 2, LAND); 
    if ( strcmp( flagname[i], "LAND") == 0) { 
      flagtest[0] |= (uint32_t) pow( 2, LAND); 
      flagtest[1] |= (uint32_t) pow( 2, LAND); 
    }

    if ( strcmp( flagname[i], "ICEYELLOW") == 0) 
      flagtest[0] |= (uint32_t) pow( 2, ICE); 
    if ( strcmp( flagname[i], "ICERED") == 0) 
      flagtest[1] |= (uint32_t) pow( 2, ICE); 
    if ( strcmp( flagname[i], "ICE") == 0) {
      flagtest[0] |= (uint32_t) pow( 2, ICE); 
      flagtest[1] |= (uint32_t) pow( 2, ICE); 
    }

    if ( strcmp( flagname[i], "WINDYELLOW") == 0) 
      flagtest[0] |= (uint32_t) pow( 2, WIND); 
    if ( strcmp( flagname[i], "WINDRED") == 0) 
      flagtest[1] |= (uint32_t) pow( 2, WIND); 
    if ( strcmp( flagname[i], "WIND") == 0) {
      flagtest[0] |= (uint32_t) pow( 2, WIND); 
      flagtest[1] |= (uint32_t) pow( 2, WIND); 
    }
    if ( strcmp( flagname[i], "WINDCONVERGE") == 0) 
      flagtest[2] |= (uint32_t) pow( 2, WIND); 
    if ( strcmp( flagname[i], "SCATRFI") == 0) 
      flagtest[3] |= (uint32_t) pow( 2, WIND); 


    if ( strcmp( flagname[i], "FLAREYELLOW") == 0) 
      flagtest[0] |= (uint32_t) pow( 2, FLARE); 
    if ( strcmp( flagname[i], "FLARERED") == 0) 
      flagtest[1] |= (uint32_t) pow( 2, FLARE); 
    if ( strcmp( flagname[i], "FLARE") == 0) {
      flagtest[0] |= (uint32_t) pow( 2, FLARE); 
      flagtest[1] |= (uint32_t) pow( 2, FLARE); 
    }

    if ( strcmp( flagname[i], "RFIYELLOW") == 0) 
      flagtest[0] |= (uint32_t) pow( 2, RFI); 
    if ( strcmp( flagname[i], "RFIRED") == 0) 
      flagtest[1] |= (uint32_t) pow( 2, RFI); 
    if ( strcmp( flagname[i], "RFI") == 0) {
      flagtest[0] |= (uint32_t) pow( 2, RFI); 
      flagtest[1] |= (uint32_t) pow( 2, RFI); 
    }


    if ( strcmp( flagname[i], "TEMPYELLOW") == 0) {
      flagtest[0] |= (uint32_t) pow( 2, TEMP); 
      flagtest[2] |= (uint32_t) pow( 2, TEMP); 
    }

    if ( strcmp( flagname[i], "TEMPRED") == 0) {
      flagtest[1] |= (uint32_t) pow( 2, TEMP); 
      flagtest[3] |= (uint32_t) pow( 2, TEMP); 
    }
    if ( strcmp( flagname[i], "TEMPVYELLOW") == 0) {
      flagtest[0] |= (uint32_t) pow( 2, TEMP); 
    }
    if ( strcmp( flagname[i], "TEMPHYELLOW") == 0) {
      flagtest[2] |= (uint32_t) pow( 2, TEMP); 
    }
    if ( strcmp( flagname[i], "TEMPVRED") == 0) {
      flagtest[1] |= (uint32_t) pow( 2, TEMP); 
    }
    if ( strcmp( flagname[i], "TEMPHRED") == 0) {
      flagtest[3] |= (uint32_t) pow( 2, TEMP); 
    }
    if ( strcmp( flagname[i], "TEMP") == 0) {
      flagtest[0] |= (uint32_t) pow( 2, TEMP); 
      flagtest[1] |= (uint32_t) pow( 2, TEMP); 
      flagtest[2] |= (uint32_t) pow( 2, TEMP); 
      flagtest[3] |= (uint32_t) pow( 2, TEMP); 
    }


    if ( strcmp( flagname[i], "MOONYELLOW") == 0) {
      flagtest[0] |= (uint32_t) pow( 2, MOON); 
      flagtest[2] |= (uint32_t) pow( 2, MOON); 
    }

    if ( strcmp( flagname[i], "MOONRED") == 0) {
      flagtest[1] |= (uint32_t) pow( 2, MOON); 
      flagtest[3] |= (uint32_t) pow( 2, MOON); 
    }
    if ( strcmp( flagname[i], "MOONVYELLOW") == 0) {
      flagtest[0] |= (uint32_t) pow( 2, MOON); 
    }
    if ( strcmp( flagname[i], "MOONHYELLOW") == 0) {
      flagtest[2] |= (uint32_t) pow( 2, MOON); 
    }
    if ( strcmp( flagname[i], "MOONVRED") == 0) {
      flagtest[1] |= (uint32_t) pow( 2, MOON); 
    }
    if ( strcmp( flagname[i], "MOONHRED") == 0) {
      flagtest[3] |= (uint32_t) pow( 2, MOON); 
    }
    if ( strcmp( flagname[i], "MOON") == 0) {
      flagtest[0] |= (uint32_t) pow( 2, MOON); 
      flagtest[1] |= (uint32_t) pow( 2, MOON); 
      flagtest[2] |= (uint32_t) pow( 2, MOON); 
      flagtest[3] |= (uint32_t) pow( 2, MOON); 
    }

    if ( strcmp( flagname[i], "SUNGLINTYELLOW") == 0) {
      flagtest[0] |= (uint32_t) pow( 2, SUNGLINT); 
      flagtest[2] |= (uint32_t) pow( 2, SUNGLINT); 
    }

    if ( strcmp( flagname[i], "SUNGLINTRED") == 0) {
      flagtest[1] |= (uint32_t) pow( 2, SUNGLINT); 
      flagtest[3] |= (uint32_t) pow( 2, SUNGLINT); 
    }
    if ( strcmp( flagname[i], "SUNGLINTVYELLOW") == 0) {
      flagtest[0] |= (uint32_t) pow( 2, SUNGLINT); 
    }
    if ( strcmp( flagname[i], "SUNGLINTHYELLOW") == 0) {
      flagtest[2] |= (uint32_t) pow( 2, SUNGLINT); 
    }
    if ( strcmp( flagname[i], "SUNGLINTVRED") == 0) {
      flagtest[1] |= (uint32_t) pow( 2, SUNGLINT); 
    }
    if ( strcmp( flagname[i], "SUNGLINTHRED") == 0) {
      flagtest[3] |= (uint32_t) pow( 2, SUNGLINT); 
    }
    if ( strcmp( flagname[i], "SUNGLINT") == 0) {
      flagtest[0] |= (uint32_t) pow( 2, SUNGLINT); 
      flagtest[1] |= (uint32_t) pow( 2, SUNGLINT); 
      flagtest[2] |= (uint32_t) pow( 2, SUNGLINT); 
      flagtest[3] |= (uint32_t) pow( 2, SUNGLINT); 
    }

    if ( strcmp( flagname[i], "FLUXDYELLOW") == 0) {
      flagtest[0] |= (uint32_t) pow( 2, FLUXD); 
      flagtest[2] |= (uint32_t) pow( 2, FLUXD); 
    }

    if ( strcmp( flagname[i], "FLUXDRED") == 0) {
      flagtest[1] |= (uint32_t) pow( 2, FLUXD); 
      flagtest[3] |= (uint32_t) pow( 2, FLUXD); 
    }
    if ( strcmp( flagname[i], "FLUXDVYELLOW") == 0) {
      flagtest[0] |= (uint32_t) pow( 2, FLUXD); 
    }
    if ( strcmp( flagname[i], "FLUXDHYELLOW") == 0) {
      flagtest[2] |= (uint32_t) pow( 2, FLUXD); 
    }
    if ( strcmp( flagname[i], "FLUXDVRED") == 0) {
      flagtest[1] |= (uint32_t) pow( 2, FLUXD); 
    }
    if ( strcmp( flagname[i], "FLUXDHRED") == 0) {
      flagtest[3] |= (uint32_t) pow( 2, FLUXD); 
    }
    if ( strcmp( flagname[i], "FLUXD") == 0) {
      flagtest[0] |= (uint32_t) pow( 2, FLUXD); 
      flagtest[1] |= (uint32_t) pow( 2, FLUXD); 
      flagtest[2] |= (uint32_t) pow( 2, FLUXD); 
      flagtest[3] |= (uint32_t) pow( 2, FLUXD); 
    }

    if ( strcmp( flagname[i], "FLUXRYELLOW") == 0) {
      flagtest[0] |= (uint32_t) pow( 2, FLUXR); 
      flagtest[2] |= (uint32_t) pow( 2, FLUXR); 
    }

    if ( strcmp( flagname[i], "FLUXRRED") == 0) {
      flagtest[1] |= (uint32_t) pow( 2, FLUXR); 
      flagtest[3] |= (uint32_t) pow( 2, FLUXR); 
    }
    if ( strcmp( flagname[i], "FLUXRVYELLOW") == 0) {
      flagtest[0] |= (uint32_t) pow( 2, FLUXR); 
    }
    if ( strcmp( flagname[i], "FLUXRHYELLOW") == 0) {
      flagtest[2] |= (uint32_t) pow( 2, FLUXR); 
    }
    if ( strcmp( flagname[i], "FLUXRVRED") == 0) {
      flagtest[1] |= (uint32_t) pow( 2, FLUXR); 
    }
    if ( strcmp( flagname[i], "FLUXRHRED") == 0) {
      flagtest[3] |= (uint32_t) pow( 2, FLUXR); 
    }
    if ( strcmp( flagname[i], "FLUXR") == 0) {
      flagtest[0] |= (uint32_t) pow( 2, FLUXR); 
      flagtest[1] |= (uint32_t) pow( 2, FLUXR); 
      flagtest[2] |= (uint32_t) pow( 2, FLUXR); 
      flagtest[3] |= (uint32_t) pow( 2, FLUXR); 
    }

    if ( strcmp( flagname[i], "GALACTICYELLOW") == 0) {
      flagtest[0] |= (uint32_t) pow( 2, GALACTIC); 
      flagtest[2] |= (uint32_t) pow( 2, GALACTIC); 
    }

    if ( strcmp( flagname[i], "GALACTICRED") == 0) {
      flagtest[1] |= (uint32_t) pow( 2, GALACTIC); 
      flagtest[3] |= (uint32_t) pow( 2, GALACTIC); 
    }
    if ( strcmp( flagname[i], "GALACTICVYELLOW") == 0) {
      flagtest[0] |= (uint32_t) pow( 2, GALACTIC); 
    }
    if ( strcmp( flagname[i], "GALACTICHYELLOW") == 0) {
      flagtest[2] |= (uint32_t) pow( 2, GALACTIC); 
    }
    if ( strcmp( flagname[i], "GALACTICVRED") == 0) {
      flagtest[1] |= (uint32_t) pow( 2, GALACTIC); 
    }
    if ( strcmp( flagname[i], "GALACTICHRED") == 0) {
      flagtest[3] |= (uint32_t) pow( 2, GALACTIC); 
    }
    if ( strcmp( flagname[i], "GALACTIC") == 0) {
      flagtest[0] |= (uint32_t) pow( 2, GALACTIC); 
      flagtest[1] |= (uint32_t) pow( 2, GALACTIC); 
      flagtest[2] |= (uint32_t) pow( 2, GALACTIC); 
      flagtest[3] |= (uint32_t) pow( 2, GALACTIC); 
    }


    if ( strcmp( flagname[i], "NAVROLL") == 0) 
      flagtest[0] |= (uint32_t) pow( 2, NAV); 
    if ( strcmp( flagname[i], "NAVPITCH") == 0) 
      flagtest[1] |= (uint32_t) pow( 2, NAV); 
    if ( strcmp( flagname[i], "NAVYAW") == 0) 
      flagtest[2] |= (uint32_t) pow( 2, NAV); 
    if ( strcmp( flagname[i], "NAVOFF") == 0) 
      flagtest[3] |= (uint32_t) pow( 2, NAV); 
    if ( strcmp( flagname[i], "NAV") == 0) {
      flagtest[0] |= (uint32_t) pow( 2, NAV); 
      flagtest[1] |= (uint32_t) pow( 2, NAV); 
      flagtest[2] |= (uint32_t) pow( 2, NAV); 
      flagtest[3] |= (uint32_t) pow( 2, NAV);
    }


    if ( strcmp( flagname[i], "SAOVERFLOW") == 0) 
      flagtest[0] |= (uint32_t) pow( 2, SAOVERFLOW);  

    if ( strcmp( flagname[i], "ROUGHOOB") == 0) 
      flagtest[0] |= (uint32_t) pow( 2, ROUGH);  
    if ( strcmp( flagname[i], "SWH") == 0) 
      flagtest[1] |= (uint32_t) pow( 2, ROUGH);  

    if ( strcmp( flagname[i], "POINTING") == 0) {
      flagtest[0] |= (uint32_t) pow( 2, POINTING); 
      flagtest[1] |= (uint32_t) pow( 2, POINTING);  
    }

    if ( strcmp( flagname[i], "TBCONS") == 0)
      flagtest[0] |= (uint32_t) pow( 2, TBCONS); 
    if ( strcmp( flagname[i], "TBCONS") == 0) 
      flagtest[1] |= (uint32_t) pow( 2, TBCONS); 

    if ( strcmp( flagname[i], "COLDWATERYELLOW") == 0) 
      flagtest[0] |= (uint32_t) pow( 2, COLDWATER); 
    if ( strcmp( flagname[i], "COLDWATERRED") == 0) 
      flagtest[1] |= (uint32_t) pow( 2, COLDWATER); 
    if ( strcmp( flagname[i], "COLDWATER") == 0) {
      flagtest[0] |= (uint32_t) pow( 2, COLDWATER); 
      flagtest[1] |= (uint32_t) pow( 2, COLDWATER); 
    }

    if ( strcmp( flagname[i], "TFTADIFFYELLOW") == 0) 
      flagtest[0] |= (uint32_t) pow( 2, TFTADIFF); 
    if ( strcmp( flagname[i], "TFTADIFFRED") == 0) 
      flagtest[1] |= (uint32_t) pow( 2, TFTADIFF); 
    if ( strcmp( flagname[i], "TFTADIFF") == 0) {
      flagtest[0] |= (uint32_t) pow( 2, TFTADIFF); 
      flagtest[1] |= (uint32_t) pow( 2, TFTADIFF); 
    }

    if ( strcmp( flagname[i], "REFL_1STOKESMOONYELLOW") == 0) {
      flagtest[0] |= (uint32_t) pow( 2, REFL_1STOKES);
    }
    if ( strcmp( flagname[i], "REFL_1STOKESMOONRED") == 0) {
      flagtest[1] |= (uint32_t) pow( 2, REFL_1STOKES);
    }
    if ( strcmp( flagname[i], "REFL_1STOKESGAL") == 0) {
      flagtest[2] |= (uint32_t) pow( 2, REFL_1STOKES);
    }
    if ( strcmp( flagname[i], "REFL_1STOKES") == 0) {
      flagtest[0] |= (uint32_t) pow( 2, REFL_1STOKES);
      flagtest[1] |= (uint32_t) pow( 2, REFL_1STOKES);
      flagtest[2] |= (uint32_t) pow( 2, REFL_1STOKES);
    }

    if ( strcmp( flagname[i], "RFI_REGION") == 0)
      flagtest[0] |= (uint32_t) pow( 2, RFI_REGION); 

  }


  // Exclude all flagged data but Dense Vegetation
  if ( strcmp(datatype, "SM") == 0) {
    // Dense Vegetation flag switched to 512
    //    flagtest[0] = 0xffff - 1024;
    flagtest[0] = 0xffff - 512; 
  }

  /* Compute numbin array (Number of bins in each row) */
  /* ------------------------------------------------- */
  numbin = (int32 *) calloc(nrows, sizeof(int32));
  for (size_t i=0; i<nrows; i++) {
    latbin = (i + 0.5) * (180.0 / nrows) - 90.0;
    numbin[i] = (int32) (cos(latbin * PI/180.0) * (2.0*nrows) + 0.5);
  }

  /* Compute basebin array (Starting bin of each row [1-based]) */
  /* ---------------------------------------------------------- */
  basebin = (int32 *) calloc(nrows+1, sizeof(int32));
  basebin[0] = 1;
  for (size_t i=1; i<=nrows; i++) {
    basebin[i] = basebin[i-1] + numbin[i-1];
  }
  printf("total number of bins: %d\n", basebin[nrows]-1);

  uint32_t n_bins_in_group = basebin[nrows]-1;


  // Open near-bins files
  int32 *near_bins_index;
  float32 *bins_ll;
  FILE *fp_near_bins, *fp_ll;
  char *tmp_str = getenv("OCDATAROOT");
  if (tmp_str == 0x0) {
    printf("Environment variable OCDATAROOT not defined.\n");
    exit(1);
  }

  if (smooth) {
    strcpy(buf, tmp_str);
    strcat(buf, "/aquarius/radiometer/near_1D_4D.bin");
    fp_near_bins = fopen(buf, "rb");
    if ( fp_near_bins == NULL) {
      printf("%s not found.\n", buf);
      exit(1);
    }
    near_bins_index = (int32 *) calloc(n_bins_in_group+1, sizeof(int32));
    fread(near_bins_index, n_bins_in_group+1, sizeof(int32), fp_near_bins);

    strcpy(buf, tmp_str);
    strcat(buf, "/aquarius/radiometer/1D_ll.bin");
    fp_ll = fopen(buf, "rb");
    if ( fp_ll == NULL) {
      printf("%s not found.\n", buf);
      exit(1);
    }
    bins_ll = (float32 *) calloc(n_bins_in_group*2, sizeof(float32));
    fread(bins_ll, n_bins_in_group*2, sizeof(float32), fp_ll);
    fclose(fp_ll);
  }


  /* Create output file */
  /* ------------------ */
  strcpy(buf, l2binInput.ofile);
  static Hdf::hdf5_bin output_binfile;
  output_binfile.create( buf, nrows);


  /* Allocate Arrays for Bin Index */
  /* ----------------------------- */
  bin_indx = (uint8 *) calloc(36 * nrows, 1);

  float *cellon, *cellat;
  double *zang;

  int32 n_filled_bins = 0;
  bin_flag = (int32 *) calloc(n_bins_in_group, sizeof(int32));
  nscenes  = (int16 *) calloc(n_bins_in_group, sizeof(int16));
  lastfile = (int16 *) calloc(n_bins_in_group, sizeof(int16));


  /* Allocate bin accumulator & data value arrays */
  /* -------------------------------------------- */
  nobs = (uint16 *) calloc(n_bins_in_group, sizeof(int16));
  allocated_space = (int16 *) calloc(n_bins_in_group, sizeof(int16));
  data_values = (float32 **) calloc(n_bins_in_group, sizeof(float32 *));

  if (smooth) {
    data_longitude = (float32 **) calloc(n_bins_in_group, sizeof(float32 *));
    data_latitude = (float32 **) calloc(n_bins_in_group, sizeof(float32 *));
  }

  file_index = (int16 **) calloc(n_bins_in_group, sizeof(int16 *));

  for (size_t i=0; i<n_bins_in_group; i++) {
      nobs[i] = 0;
      allocated_space[i] = 0;
      lastfile[i] = -1;
  }


  double minTime=1e30, maxTime=-1e30;

  // ============ Main loop =========
  for (ifile=0; ifile<nfiles; ifile++) {

    cellon = new float[numBlks[ifile]*NUMBER_OF_BEAMS];
    cellat = new float[numBlks[ifile]*NUMBER_OF_BEAMS];
    zang = new double[numBlks[ifile]];

    hid_t gid[6];
    l2file[ifile].getH5gid( l2file[ifile].getH5fid(), gid);

    l2file[ifile].readl2_nav(cellat, cellon, zang);

    cout << "Binning file: " << ifile << endl;

    // Loop over blocks
    for (iblk=0; iblk<numBlks[ifile]; iblk++) {

      // Skip if in orbit overlap region
      l2file[ifile].readl2sec( iblk, gid[1], &blkSec);
      if (blkSec <= granstart[ifile] || blkSec > granstop[ifile]) {
      	continue;
      }

      // Skip if wrong type of orbit (A/D)
      if ( (strcmp(l2binInput.orbit_type, "A") == 0 && (zang[iblk] > 180)) ||
	   (strcmp(l2binInput.orbit_type, "D") == 0 && (zang[iblk] < 180))) {
	continue;
      }

      int32_t rad_flags[NUMBER_OF_BEAMS][4];
      int32_t rad_flags_sm[NUMBER_OF_BEAMS];
      if ( strcmp(datatype, "SM") == 0) {
	l2file[ifile].readl2radflagsm( iblk, gid[2], &rad_flags_sm);
      } else {
	l2file[ifile].readl2radflag( iblk, gid[2], &rad_flags[0]);
      }

      for (ibeam=0; ibeam<NUMBER_OF_BEAMS; ibeam++) {

	// Skip if wrong beam
	if ( l2binInput.beam >= 0)
	  if ( l2binInput.beam != ibeam) continue;

	if ( strcmp(datatype, "SM") == 0) {
	    if ((rad_flags_sm[ibeam] & flagtest[0]) != 0) continue;
	} else {
	  // Skip if flagged
	  if (!l2binInput.require) {
	    if ((rad_flags[ibeam][0] & flagtest[0]) != 0) continue;
	    if ((rad_flags[ibeam][1] & flagtest[1]) != 0) continue;
	    if ((rad_flags[ibeam][2] & flagtest[2]) != 0) continue;
	    if ((rad_flags[ibeam][3] & flagtest[3]) != 0) continue;
	  } else {
	    if ((rad_flags[ibeam][0] & flagtest[0]) == 0 &&
		(rad_flags[ibeam][1] & flagtest[1]) == 0 &&
		(rad_flags[ibeam][2] & flagtest[2]) == 0 &&
		(rad_flags[ibeam][3] & flagtest[3]) == 0) continue;
	  }
	}

        if (cellon[NUMBER_OF_BEAMS*iblk+ibeam] == -999) continue;
        if (cellat[NUMBER_OF_BEAMS*iblk+ibeam] == -999) continue;

	// Set min/max times is applicable
	if ( blkSec < minTime) minTime = blkSec;
	if ( blkSec > maxTime) maxTime = blkSec;

	if ( cellat[NUMBER_OF_BEAMS*iblk+ibeam] > northmost) 
	northmost = cellat[NUMBER_OF_BEAMS*iblk+ibeam];
	if ( cellat[NUMBER_OF_BEAMS*iblk+ibeam] < southmost) 
	southmost = cellat[NUMBER_OF_BEAMS*iblk+ibeam];
	if ( cellon[NUMBER_OF_BEAMS*iblk+ibeam] > eastmost) 
	eastmost = cellon[NUMBER_OF_BEAMS*iblk+ibeam];
	if ( cellon[NUMBER_OF_BEAMS*iblk+ibeam] < westmost) 
	westmost = cellon[NUMBER_OF_BEAMS*iblk+ibeam];

	// Get bin number (1-based)
	int32_t bin = getbinnum(cellon[NUMBER_OF_BEAMS*iblk+ibeam], 
				cellat[NUMBER_OF_BEAMS*iblk+ibeam]);

	int32_t ibin = bin - basebin[0]; // 0-based bin number

	/* increment nscenes */
	/* ----------------- */
	if (ifile != lastfile[ibin]) {
	  nscenes[ibin]++;
	  lastfile[ibin] = ifile;
	}

	/* Allocate space for file index & bin data values */
	/* ----------------------------------------------- */
	if (file_index[ibin] == NULL) {

          //          cout << "Alloc: " << iblk << " " << ibeam << " " << 
          // ibin << " " << nobs[ibin] << endl;

	  file_index[ibin] = (int16 *) calloc(n_allocperbin, sizeof(int16));

	  if ( nProd != 0) {
	    data_values[ibin] = 
	      (float32 *) calloc(n_allocperbin*nProd, sizeof(float32));

	    if (data_values[ibin] == 0x0) {
	      perror(buf);
	      printf("Allocation failed for data_values[bin]: %d %s\n",
		     ibin,buf);
	      exit(-1);
	    }

	    if (smooth) {
	      data_longitude[ibin] = 
		(float32 *) calloc(n_allocperbin, sizeof(float32));

	      if (data_longitude[ibin] == 0x0) {
		perror(buf);
		printf("Allocation failed for data_longitude[bin]: %d %s\n",
		       ibin,buf);
		exit(-1);
	      }
	      
	      data_latitude[ibin] = 
		(float32 *) calloc(n_allocperbin, sizeof(float32));

	      if (data_latitude[ibin] == 0x0) {
		perror(buf);
		printf("Allocation failed for data_latitude[bin]: %d %s\n",
		       ibin,buf);
		exit(-1);
	      }
	    }
	  }

	  allocated_space[ibin] = n_allocperbin;
	}


	/* Set file_index for each observation */
	/* ----------------------------------- */
	file_index[ibin][nobs[ibin]] = ifile;


	/* Get data values for all L3 products */
	/* ----------------------------------- */
	for (int32_t jprod=0; jprod<nProd; jprod++) {

          l2file[ifile].readl2( iblk, gid[3], prodname[jprod], data);

	  data_values[ibin][nProd*nobs[ibin]+jprod] = data[ibeam];

	  if ( strcmp( prodname[jprod], "SSS") == 0 &&
	       data[ibeam] == -9999 && !l2binInput.require) {
	    cout << "Missing value binned for ifile: " << ifile
		 << " iblk: " << iblk << " ibeam: " << ibeam
		 << " product: " << prodname[jprod] << endl;
	    exit(1);
	  }

	} /* jprod loop */


	/* Store data longitude & latitude */
	/* ------------------------------- */ 
	if ( smooth) {
	  data_longitude[ibin][nobs[ibin]] = (PI / 180.0) * 
	    cellon[NUMBER_OF_BEAMS*iblk+ibeam];  
	  data_latitude[ibin][nobs[ibin]] =  (PI / 180.0) * 
	    cellat[NUMBER_OF_BEAMS*iblk+ibeam];  
	}

	/* Increment number of observations in bin */
	/* --------------------------------------- */
	nobs[ibin]++;


	/* Reallocate if necessary */
	/* ----------------------- */
	if (nobs[ibin] == allocated_space[ibin]) {

          //          cout << "Realloc: " << iblk << " " << ibeam << " " << 
          //ibin << " " << nobs[ibin] << endl;

	  file_index[ibin] = 
	    (int16 *) realloc(file_index[ibin], 
			      (nobs[ibin]+n_allocperbin) * sizeof(int16));

	  if ( nProd != 0) {
	    data_values[ibin] = 
	      (float32 *) realloc(data_values[ibin], 
				  (nobs[ibin]+n_allocperbin) * nProd * 
				  sizeof(float32));
	    if (data_values[ibin] == 0x0) {
	      perror(buf);
	      printf("Reallocation failed for data_values[bin]: %d\n",
		     ibin);
	      exit(-1);
	    }

	    if (smooth) {
	      data_longitude[ibin] = 
		(float32 *) realloc(data_longitude[ibin], 
				    (nobs[ibin]+n_allocperbin) * 
				    sizeof(float32));
	      if (data_longitude[ibin] == 0x0) {
		perror(buf);
		printf("Reallocation failed for data_longitude[bin]: %d\n",
		       ibin);
		exit(-1);
	      }
	      
	      data_latitude[ibin] = 
		(float32 *) realloc(data_latitude[ibin], 
				    (nobs[ibin]+n_allocperbin) * 
				    sizeof(float32));
	      if (data_latitude[ibin] == 0x0) {
		perror(buf);
		printf("Reallocation failed for data_latitude[bin]: %d\n",
		       ibin);
		exit(-1);
	      }
	    }
	  }

	  allocated_space[ibin] += n_allocperbin;
	} /* end reallocate */
      } // ibeam

    } // iblk

    delete[] cellon;
    delete[] cellat;
    delete[] zang;
  } // ifile



  /* Compute Total # of filled bins */
  /* ------------------------------ */
  for (size_t i=0; i<n_bins_in_group; i++) {
    if (nobs[i] != 0) {
      n_filled_bins++;
      n_total_obs += nobs[i];
    }
  } /* bin loop */



  /* ********** If filled bins ********** */
  /* ------------------------------------ */
  if (n_filled_bins > 0) {

    /* Fill "Bin List" vdata array */
    /* --------------------------- */
    Hdf::binListStruct_hdf5 *outBinList;
    outBinList =(Hdf::binListStruct_hdf5 *) 
      calloc(n_filled_bins, sizeof(Hdf::binListStruct_hdf5));

    if ( not smooth) {
      int32_t i = 0;
      for (size_t ibin=0; ibin<n_bins_in_group; ibin++) {

	if (nobs[ibin] != 0) {

	  // weights {=sqrt(# of L2 files in given bin)} */
	  wgt = 0.0;
	  for (ifile=0; ifile<=nfiles; ifile++) {
	    i32 = 0;
	    for (size_t j=0; j<nobs[ibin]; j++) {
	      if (file_index[ibin][j] == ifile) i32++;
	    }
	    wgt += sqrt(i32);
	  }  
	  
	  int32_t bin = ibin + basebin[0];

	  outBinList[i].bin_num = bin;
	  outBinList[i].nobs = nobs[ibin];
	  outBinList[i].nscenes = nscenes[ibin];
	  outBinList[i].weights = wgt;
	  outBinList[i].flags_set = bin_flag[ibin];

	  i++;
	} /* nobs[ibin] != 0 */
      } /* ibin loop */

      // If no output product write binlist here
      if ( nProd == 0)
	output_binfile.write( NULL, n_filled_bins, NULL, outBinList);

      //	output_binfile.writeBinList( n_filled_bins);

    } // if not smooth

    /* Print info on filled row group */
    /* ------------------------------ */
    printf("%-20s:%8d\n", "# bins", n_bins_in_group);
    printf("%-20s:%8d\n", "# filled bins", n_filled_bins);
    printf("%-20s:%8d\n", "# total obs", n_total_obs);
    printf("\n");


    /* Allocate sum & sum-squared arrays */
    /* --------------------------------- */
    sum_data  = (float32 *) calloc(n_bins_in_group, sizeof(float32));
    sum2_data = (float32 *) calloc(n_bins_in_group, sizeof(float32));

    /* Loop over all L3 products to fill sum arrays */
    /* -------------------------------------------- */
    for (size_t iprod=0; iprod < nProd; iprod++) {

      memset(sum_data,  0, n_bins_in_group * sizeof(float32));
      memset(sum2_data, 0, n_bins_in_group * sizeof(float32));

      fldname3[0] = (char *) calloc(strlen(prodname[iprod])+5,sizeof(char));
      fldname3[1] = (char *) calloc(strlen(prodname[iprod])+8,sizeof(char));

      char_ptr1 = strchr(prodname[iprod], '/');
      if (char_ptr1 != NULL) *char_ptr1 = '_';

      strcpy(fldname3[0], prodname[iprod]);
      strcpy(fldname3[1], prodname[iprod]);
      strcat(fldname3[0], "_sum");
      strcat(fldname3[1], "_sum_sq");

      if (char_ptr1 != NULL) *char_ptr1 = '/';

      /* Process bins */
      /* ------------ */
      for (size_t bin=0; bin<n_bins_in_group; bin++) {

	if ((bin % 10000) == 0) {
	  time(&tnow);
	  tmnow = localtime(&tnow);
	}

	if ( smooth) {
	  sum_data[bin] = -999;
	}

	if (nobs[bin] == (0 & (l2binInput.fill == 0))) continue;

	// smooth
	if (smooth) {

	  float64 x_val[NUMBER_OF_VALUES];
	  float64 y_val[NUMBER_OF_VALUES];
	  float64 d_val[NUMBER_OF_VALUES];
	  float64 w_val[NUMBER_OF_VALUES];

	  // bin center coordinates
	  float32 x0 = cos(bins_ll[2*bin+1]) * cos(bins_ll[2*bin]);
	  float32 y0 = cos(bins_ll[2*bin+1]) * sin(bins_ll[2*bin]);
	  float32 z0 = sin(bins_ll[2*bin+1]);

	  // rotation matrix to transform bin center to north poll
	  float32 rotmat[3][3];
	  rotmat[0][0] =  cos(bins_ll[2*bin]) * sin(bins_ll[2*bin+1]);
	  rotmat[0][1] =  sin(bins_ll[2*bin]) * sin(bins_ll[2*bin+1]);
	  rotmat[0][2] = -cos(bins_ll[2*bin+1]);
	  rotmat[1][0] = -sin(bins_ll[2*bin]);
	  rotmat[1][1] =  cos(bins_ll[2*bin]);
	  rotmat[1][2] = 0;
	  rotmat[2][0] =  cos(bins_ll[2*bin]) * cos(bins_ll[2*bin+1]);
	  rotmat[2][1] =  sin(bins_ll[2*bin]) * cos(bins_ll[2*bin+1]);
	  rotmat[2][2] =  sin(bins_ll[2*bin+1]);

	  // If rot_mat is correct then rot_x0 & rot_y0 = 0 & rot_z0 = 1
	  //float32 rot_x0 = 
	  //rotmat[0][0] * x0 + rotmat[0][1] * y0 + rotmat[0][2] * z0;
	  //float32 rot_y0 = 
	  //rotmat[1][0] * x0 + rotmat[1][1] * y0 + rotmat[1][2] * z0;
	  //float32 rot_z0 = 
	  //rotmat[2][0] * x0 + rotmat[2][1] * y0 + rotmat[2][2] * z0;
	  
	  fseek(fp_near_bins, 
		(n_bins_in_group+1+near_bins_index[bin])*sizeof(int32),
		SEEK_SET);

	  i32 = 0;
	  int32 jbin;
	  for (size_t k=near_bins_index[bin]; k<near_bins_index[bin+1]; k++) {

	    fread(&jbin, sizeof(int32), 1, fp_near_bins);

	    for (j=0; j<nobs[jbin]; j++) {

	      float32 lon = data_longitude[jbin][j];
	      float32 lat = data_latitude[jbin][j];

	      float32 x1 = cos(lat) * cos(lon);
	      float32 y1 = cos(lat) * sin(lon);
	      float32 z1 = sin(lat);

	      float32 rot_x1 = 
		rotmat[0][0] * x1 + rotmat[0][1] * y1 + rotmat[0][2] * z1;
	      float32 rot_y1 = 
		rotmat[1][0] * x1 + rotmat[1][1] * y1 + rotmat[1][2] * z1;
	      //float32 rot_z1 = 
	      //rotmat[2][0] * x1 + rotmat[2][1] * y1 + rotmat[2][2] * z1;

	      float32 cs = x0*x1 + y0*y1 + z0*z1;
	      
	      float data = data_values[jbin][j*nProd+iprod];

	      if (cs >= cs_filter && data > 0) {
		// J. Gales (projection of unit vector onto tangent plane)
		//x_val[i32] = rot_x1 / rot_z1;
		//y_val[i32] = rot_y1 / rot_z1;

		// J. Lilly 2/27/09 (projection upon x-y plane)
		x_val[i32] = rot_x1;
		y_val[i32] = rot_y1;

		d_val[i32] = (double) data;

		double arclen = acos((double) cs) * 180 / PI; // in degrees
		double ratio = arclen / l2binInput.filter_width;
		w_val[i32] = 1 - ratio*ratio;
		if (w_val[i32] < 0) w_val[i32] = 0;

		i32++;
		
		if ( i32 == NUMBER_OF_VALUES) {
		  printf("NUMBER_OF_VALUES parameter must be increased.\n");
		  printf("bin: %d  i32: %d\n", (int)bin, i32);
		  exit(110);
		}
	      }
	    } // for j
	  } // for k

	  // If enough data points perform lsf
	  if (i32 >= nOrder) {
	    //	  std::cout << bin << " " << i32 << std::endl;

	    double chisq;
	    gsl_matrix *X, *cov;
	    gsl_vector *y, *w, *c;

	    X = gsl_matrix_alloc (i32, nOrder);
	    y = gsl_vector_alloc (i32);
	    w = gsl_vector_alloc (i32);

	    c = gsl_vector_alloc (nOrder);
	    cov = gsl_matrix_alloc (nOrder, nOrder);

	    gsl_multifit_linear_workspace * work 
	      = gsl_multifit_linear_alloc (i32, nOrder);

	    for (size_t i=0; i<i32; i++) {
	      gsl_matrix_set (X, i, 0, 1.0);

	      switch( nOrder) {
	      case 3: 
		gsl_matrix_set (X, i, 1, x_val[i]);
		gsl_matrix_set (X, i, 2, y_val[i]);
		break;

	      case 4: 
		gsl_matrix_set (X, i, 1, x_val[i]);
		gsl_matrix_set (X, i, 2, y_val[i]);
		gsl_matrix_set (X, i, 3, x_val[i]*y_val[i]);
		break;

	      case 6: 
		gsl_matrix_set (X, i, 1, x_val[i]);
		gsl_matrix_set (X, i, 2, y_val[i]);
		gsl_matrix_set (X, i, 3, x_val[i]*y_val[i]);
		gsl_matrix_set (X, i, 4, 0.5*x_val[i]*x_val[i]);
		gsl_matrix_set (X, i, 5, 0.5*y_val[i]*y_val[i]);
		break;
	      }

	      gsl_vector_set (y, i, d_val[i]);
	      gsl_vector_set (w, i, w_val[i]);
	    }
	    gsl_multifit_wlinear (X, w, y, c, cov, &chisq, work);

	    //   if ( fabs(gsl_matrix_get(cov,0,0)/gsl_vector_get(c,0)) > 0.01)
	    //continue;

	    //	  std::cout << gsl_vector_get(c,0) << std::endl;
	    //std::cout << gsl_vector_get(c,1) << std::endl;
	    //std::cout << gsl_vector_get(c,2) << std::endl;

	    sum_data[bin] = gsl_vector_get(c,0);
	    sum2_data[bin] = sum_data[bin] * sum_data[bin];

	    gsl_multifit_linear_free (work);
	    gsl_vector_free(y);
	    gsl_vector_free(w);
	    gsl_vector_free(c);
	    gsl_matrix_free(X);
	    gsl_matrix_free(cov);
	  } // i32 >= nOrder)

	} else {

	  // Each observation has a weight given by the reciprocal of the
	  // square root of the total number of observations from a particular
	  // input swath (L2) file.  In this way multiple observations 
	  // of a given bin from a single L2 file are reduced in weight.

	  i32 = 1;
	  sum  = data_values[bin][0*nProd+iprod];
	  sum2 = sum * sum;
	  for (j=1; j<nobs[bin]; j++) {
	    if (file_index[bin][j] == file_index[bin][j-1]) {
	      i32++;
	      sum += data_values[bin][j*nProd+iprod];
	      sum2 += data_values[bin][j*nProd+iprod] *
		data_values[bin][j*nProd+iprod];
	    } else {
	      sum_data[bin]  += (sum  / sqrt(i32));
	      sum2_data[bin] += (sum2 / sqrt(i32));
	      
	      i32 = 1;
	      sum  = data_values[bin][j*nProd+iprod];
	      sum2 = sum * sum;
	    }
	  } /* observation loop */
	  sum_data[bin]  += (sum  / sqrt(i32));
	  sum2_data[bin] += (sum2 / sqrt(i32));
	}

      } /* bin loop */

      /* Write Product Vdatas */
      /* -------------------- */

      /* Determine number of smoothed bins */
      /* --------------------------------- */
      if ( smooth) {
	int n_smoothed_bins = 0;
	for (size_t bin=0; bin<n_bins_in_group; bin++) {
	  if (sum_data[bin] != -999) n_smoothed_bins++;
	}
	n_filled_bins = n_smoothed_bins;
      }

      a = (uint8 *) calloc(8 * n_filled_bins, 1);

      /* Fill bin data array */
      /* ------------------- */
      int i = 0;
      for (size_t ibin=0; ibin<n_bins_in_group; ibin++) {
	if ( smooth) {
	  if (sum_data[ibin] != -999) {
	    memcpy(&a[i*8  ], &sum_data[ibin],  4);
	    memcpy(&a[i*8+4], &sum2_data[ibin], 4);

	    int32 bin = ibin + basebin[0];

	    if ( iprod == 0) {
	      outBinList[i].bin_num = bin;
	      outBinList[i].nobs = nobs[ibin];
	      outBinList[i].nscenes = nscenes[ibin];
	      outBinList[i].weights = 1.0;
	      outBinList[i].flags_set = bin_flag[ibin];	      
	    }
	    i++;
	  }
	} else {
	  if (nobs[ibin] != 0) {
	    memcpy(&a[i*8  ], &sum_data[ibin],  4);
	    memcpy(&a[i*8+4], &sum2_data[ibin], 4);
	    i++;
	  }
	}
      }
      char_ptr1 = strchr(prodname[iprod], '/');
      if (char_ptr1 != NULL) *char_ptr1 = '_';

      if ( iprod == 0)
	output_binfile.write( prodname[iprod], n_filled_bins, 
			      (float *) a, outBinList);
      else
	output_binfile.write( prodname[iprod], n_filled_bins, 
			      (float *) a, NULL);

      if (char_ptr1 != NULL) *char_ptr1 = '/';


      free(a);

      free(fldname3[0]);
      free(fldname3[1]);
      
    } /* iprod loop */

    output_binfile.n_data_records += n_filled_bins;


    /* Free dynamic memory */
    /* ------------------- */
    if (sum_data  != NULL) free(sum_data);
    if (sum2_data != NULL) free(sum2_data);

    total_filled_bins += n_filled_bins;


  } /* ********** n_filled_bin > 0 ********** */
  time(&tnow);
  tmnow = localtime(&tnow);


  /* Free Dynamic Memory */
  /* ------------------- */
  time(&tnow);
  tmnow = localtime(&tnow);

  for (size_t i=0; i<n_bins_in_group; i++) {
    //    cout << i << endl;
    if (file_index[i] != NULL) free(file_index[i]);
    if (data_values[i] != NULL) free(data_values[i]);
  }

  time(&tnow);
  tmnow = localtime(&tnow);

  free(data_values);
  if (smooth) {
    free(data_longitude);
    free(data_latitude);
  }
  free(nobs);
  free(allocated_space);
  free(file_index);

  time(&tnow);
  tmnow = localtime(&tnow);

  printf("total_filled_bins: %d\n", total_filled_bins); 

  if (total_filled_bins == 0) {
    for (size_t i=0; i<nfiles; i++) l2file[i].closel2();
   strcpy(buf, "rm -f ");
   strcat(buf, l2binInput.ofile);
   strcat(buf, "*");
   printf("%s\n", buf);
   system(buf);
   return 1;
  }

  cout << "in close" << endl;
  for (size_t i=0; i<nfiles; i++) l2file[i].closel2();


  /* Read and write global attributes */
  /* -------------------------------- */
  printf("Writing Global Attributes\n");

  // Mark metadata as already written (hdf5_bin::close())
  strcpy( output_binfile.meta_l3b.product_name, "___");

  hid_t grp0 = output_binfile.get_grp0();

  PTB( SetScalarH5A (grp0, "Product Name", H5T_STRING, l2binInput.ofile));

  PTB( SetScalarH5A (grp0, "Title" , H5T_STRING, 
		     (VOIDP) "Aquarius Level-3 Binned Data"));

  PTB( SetScalarH5A (grp0, "Sensor" , H5T_STRING, (VOIDP) "Aquarius"));

  PTB( SetScalarH5A (grp0, "Sensor Characteristics", H5T_STRING, 
		       (VOIDP) SENSORCHARACTERISTICS));

  PTB( SetScalarH5A (grp0, "Mission", H5T_STRING, (VOIDP) "SAC-D Aquarius"));

  PTB( SetScalarH5A (grp0, "Mission Characteristics", H5T_STRING, 
		       (VOIDP) MISSIONCHARACTERISTICS));

  //  strcpy(buf, "ORIGINAL");
  //PTB( SetScalarH5A (grp0, "Replacement Flag", H5T_STRING, buf));

  strcpy(buf, "l2bin_aquarius");
  PTB( SetScalarH5A (grp0, "Software Name", H5T_STRING, buf));

  strcpy(buf, VERSION);
  PTB( SetScalarH5A (grp0, "Software ID", H5T_STRING, buf));

  PTB( SetScalarH5A (grp0, "Processing Version", H5T_STRING, 
		     (VOIDP) l2binInput.pversion));

  PTB( SetScalarH5A (grp0, "Start Orbit", H5T_STD_U32LE, (VOIDP) &s_orb));
  PTB( SetScalarH5A (grp0, "End Orbit",   H5T_STD_U32LE, (VOIDP) &e_orb));

  PTB( SetScalarH5A (grp0, "Processing Time", H5T_STRING, 
		     ydhmsf(time(NULL),'G')));

  PTB( SetScalarH5A (grp0, "Processing Control", H5T_STRING, 
		     processControl.c_str()));

  // Replace 0 with ','
  for (size_t i=0; i<len_flag; i++) {
    if (l2binInput.flaguse[i] == 0) l2binInput.flaguse[i] = ',';
  }
  l2binInput.flaguse[len_flag-1] = 0;

  PTB( SetScalarH5A (grp0, "L2 Flag Names", H5T_STRING, 
		     &l2binInput.flaguse[1]));

  /*
  status = SDsetattr(sd_id_w, "Input Parameters", DFNT_CHAR, strlen(l2binInput.parms)+1, 
		     l2binInput.parms);
  */

  /*
  strcpy(buf, l2_str[0].filename);
  if (nfiles > 1) strcat(buf, ",");
  for (ifile=1; ifile<nfiles-1; ifile++) {
    strcat(buf, l2_str[ifile].filename);
    strcat(buf, ",");
  }
  if (nfiles > 1) strcat(buf, l2_str[nfiles-1].filename);
  PTB( SetScalarH5A (grp0, "Input Files", H5T_STRING, buf));
  */

  /*

  */

  istringstream istr;
  
  int32_t itemp;
  string ydhmsf_str;
  int32_t millisec;

  ydhmsf_str = ydhmsf(minTime + DIFFJAN0680_1970,'G');
  //  cout << ydhmsf_str.c_str() << endl;
  PTB( SetScalarH5A (grp0, "Start Time", H5T_STRING, ydhmsf_str.c_str()));
  istr.str( ydhmsf_str.substr(0,4)); istr >> itemp;
  //cout << itemp << endl; 
  PTB( SetScalarH5A (grp0, "Start Year", H5T_NATIVE_SHORT, (VOIDP) &itemp));

  istr.clear(); istr.str( ydhmsf_str.substr(4,3)); istr >> itemp;
  PTB( SetScalarH5A (grp0, "Start Day", H5T_NATIVE_SHORT, (VOIDP) &itemp));

  millisec = get_millisec( &ydhmsf_str);
  PTB( SetScalarH5A (grp0, "Start Millisec", H5T_STD_U32LE, 
		     (VOIDP) &millisec));

  ydhmsf_str = ydhmsf(maxTime + DIFFJAN0680_1970,'G');
  PTB( SetScalarH5A (grp0, "End Time", H5T_STRING, ydhmsf_str.c_str()));
  istr.clear(); istr.str( ydhmsf_str.substr(0,4)); istr >> itemp;
  PTB( SetScalarH5A (grp0, "End Year", H5T_NATIVE_SHORT, (VOIDP) &itemp));
  istr.clear(); istr.str( ydhmsf_str.substr(4,3)); istr >> itemp;
  PTB( SetScalarH5A (grp0, "End Day", H5T_NATIVE_SHORT, (VOIDP) &itemp));
  millisec = get_millisec( &ydhmsf_str);
  PTB( SetScalarH5A (grp0, "End Millisec", H5T_STD_U32LE, 
		     (VOIDP) &millisec));

  strcpy(buf, "degrees North");
  PTB( SetScalarH5A (grp0, "Latitude Units", H5T_STRING, buf));

  strcpy(buf, "degrees East");
  PTB( SetScalarH5A (grp0, "Longitude Units", H5T_STRING, buf));

  PTB( SetScalarH5A (grp0, "Northernmost Latitude", H5T_NATIVE_FLOAT,
		     (VOIDP) &northmost));
  PTB( SetScalarH5A (grp0, "Southernmost Latitude", H5T_NATIVE_FLOAT,
		     (VOIDP) &southmost));
  PTB( SetScalarH5A (grp0, "Easternmost Longitude", H5T_NATIVE_FLOAT,
		     (VOIDP) &eastmost));
  PTB( SetScalarH5A (grp0, "Westernmost Longitude", H5T_NATIVE_FLOAT,
		     (VOIDP) &westmost));

  PTB( SetScalarH5A (grp0, "Data Bins", H5T_STD_U32LE, 
		     (VOIDP) &total_filled_bins));

  f32 = 100 * ((float32) total_filled_bins) / 
    (basebin[nrows-1]+numbin[nrows-1]-1);
  PTB( SetScalarH5A (grp0, "Percent Data Bins", H5T_NATIVE_FLOAT,
		     (VOIDP) &f32));


  /* Free Dynamic Memory */
  /* ------------------- */
  if (bin_flag != NULL) free(bin_flag);
  if (nscenes != NULL) free(nscenes);
  if (lastfile != NULL) free(lastfile);

  free(numbin);
  free(basebin);
  free(bin_indx);

  if (smooth) fclose(fp_near_bins);

  output_binfile.close();

  return 0;
}


int parseInput( int argc, char* argv[], instr *l2binInput, string *procControl)
{
  using namespace std;

  string::size_type ipos=0;
  string            sLine, sValue;
  string            sParmWord;
  istringstream istr;
  stack<string> parmStack, procCntlStack;

  l2binInput->beam = -1;
  strcpy(l2binInput->pversion, "");
  strcpy(l2binInput->orbit_type, "");
  l2binInput->require = false;

  // Push regular command line parameters
  for (int i=1; i<argc; i++) {
    extParmWordValue( sLine.assign(argv[i]), &sParmWord, &sValue);

    if (sParmWord.compare("PAR") == 0) continue;
    if (sParmWord.compare("SUITE") == 0) continue;

    parmStack.push( sLine);
  }

  // Push SUITE parameters
  for (int i=1; i<argc; i++) {
    extParmWordValue( sLine.assign(argv[i]), &sParmWord, &sValue);
    if (sParmWord.compare("SUITE") == 0) parmStack.push( sLine);
  }

  // Push PAR parameters
  for (int i=1; i<argc; i++) {
    extParmWordValue( sLine.assign(argv[i]), &sParmWord, &sValue);
    if (sParmWord.compare("PAR") == 0) parmStack.push( sLine);
  }

  // Push default parameters
  sValue = "$OCDATAROOT/aquarius/l2bin_aquarius_defaults.par";
  expandEnvVar( &sValue);
  sLine = "par=" + sValue;
  parmStack.push( sLine);

  // Pop stack loop
  while (!parmStack.empty()) {

    sLine = parmStack.top();
    parmStack.pop();
    procCntlStack.push( sLine);

    extParmWordValue( sLine, &sParmWord, &sValue);

    if (sParmWord.compare("INFILE") == 0) {
      int len = sValue.copy(l2binInput->infile, sValue.length());
      l2binInput->infile[len] = 0;
      continue;
    }

    if (sParmWord.compare("OFILE") == 0) {
      int len = sValue.copy(l2binInput->ofile, sValue.length());
      l2binInput->ofile[len] = 0;
      continue;
    }

    if (sParmWord.compare("L3BPROD") == 0) {
      size_t found = sValue.find("default:");
      if ( found == 0) {
	strcat( l2binInput->l3bprod, sValue.substr(8).c_str());
	strcat( l2binInput->l3bprod, ":");
	continue;
      } else {
	int len;
	l2binInput->l3bprod[0] = ':';
	len = sValue.copy(&l2binInput->l3bprod[1], sValue.length());
	l2binInput->l3bprod[len+1] = ':';
	l2binInput->l3bprod[len+2] = 0;
	if ( strcmp(l2binInput->l3bprod, "::") == 0) l2binInput->l3bprod[0] = 0;
	continue;
      }
    }

    if (sParmWord.compare("FLAGUSE") == 0) {
      int len;
      l2binInput->flaguse[0] = ',';
      len = sValue.copy(&l2binInput->flaguse[1], sValue.length());
      l2binInput->flaguse[len+1] = ',';
      l2binInput->flaguse[len+2] = 0;
      continue;
    }

    if (sParmWord.compare("REQUIRE") == 0) {
      // Note: "boolalpha" allows "true/false" to be used
      istringstream iStream(sValue);
      iStream >> boolalpha >> l2binInput->require;
      continue;
    }

    if (sParmWord.compare("RESOLVE") == 0) {
      int len = sValue.copy(l2binInput->resolve, sValue.length());
      l2binInput->resolve[len] = 0;
      continue;
    }

    if (sParmWord.compare("BEAM") == 0) {
      istringstream iStream(sValue);
      if (!(iStream >> l2binInput->beam)) {
	cout << "Improper beam number" << endl;
	exit(110);
      }
      l2binInput->beam--;
      continue;
    }

    if (sParmWord.compare("FIT") == 0) {
      int len = sValue.copy(l2binInput->fit, sValue.length());
      l2binInput->fit[len] = 0;
      continue;
    }

    if (sParmWord.compare("FILTER_WIDTH") == 0) {
      istringstream iStream(sValue);
      if (!(iStream >> l2binInput->filter_width)) {
	cout << "Improper filter width value" << endl;
	exit(110);
      }
      continue;
    }

    if (sParmWord.compare("FILL") == 0) {
      istringstream iStream(sValue);
      if (!(iStream >> l2binInput->fill)) {
	cout << "Improper fill value" << endl;
	exit(110);
      }
      continue;
    }

    if (sParmWord.compare("PVERSION") == 0) {
      int len = sValue.copy(l2binInput->pversion, sValue.length());
      l2binInput->pversion[len] = 0;
      continue;
    }

    if (sParmWord.compare("ORBIT_TYPE") == 0) {
      int len = sValue.copy(l2binInput->orbit_type, sValue.length());
      l2binInput->orbit_type[len] = 0;
      continue;
    }

    if (sParmWord.compare("PAR") == 0 ||
	sParmWord.compare("SUITE") == 0) {
      ifstream parfile;
      if ( sParmWord.compare("SUITE") == 0)
	sValue = 
	  "$OCDATAROOT/aquarius/l2bin_aquarius_defaults_" + sValue + ".par";
      if ( sValue.substr( 0, 1).compare("$") == 0) expandEnvVar( &sValue);
      parfile.open ( sValue.c_str(), ifstream::in);
      if ( parfile.fail() == true) {
	cout << "Cannot open: "  <<  sValue.c_str() << endl;
	exit(1);
      }

      while( !parfile.eof()) {
	getline( parfile, sLine);
	if ( sLine.size() != 0 && sLine.substr( 0, 1).compare("#") != 0)
	  parmStack.push( sLine);
      }
      parfile.close();
      continue;
    }

    cout << endl << 
      "Parameter '" << sParmWord.c_str() << "' not defined." << endl;
    exit(1);
  }

  // Write to process control string
  while (!procCntlStack.empty()) {
    sLine = procCntlStack.top();
    procCntlStack.pop();

    extParmWordValue( sLine, &sParmWord, &sValue);
    if (sParmWord.compare("PAR") != 0 && sParmWord.compare("SUITE") != 0) {
      //      cout << sLine.c_str() << endl;
      for (size_t j=0; j<sParmWord.length(); j++)
	sParmWord[j] = tolower(sParmWord[j]);
      sParmWord = " " + sParmWord + "=";
      if ( procControl->find( sParmWord) == string::npos)
	procControl->append( " " + sLine);
    }
  }
  *procControl = procControl->substr(1);

  return 0;
}


int32 getbinnum(float32 lon, float32 lat) {

  int32 bin_row; // 0-based
  int32 bin_col; // 0-based

  bin_row = (int32_t) ((90.0 + lat) * ((float64) nrows / (float64) 180.0));
  if (bin_row < 0 || bin_row >= nrows) {
    //    printf("bin_row: %d out of bounds. (%f)\n", bin_row, scan_lat);
    return -1;
  }

  bin_col = (long) ((float64) numbin[bin_row] * 
		    (lon + 180.0) / (float64) 360.0);
 
  // bin is 1-based
  int32 bin = basebin[bin_row] + bin_col;

  return (bin);
}

