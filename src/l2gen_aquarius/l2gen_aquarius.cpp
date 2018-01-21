#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <iomanip>
#include <math.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stack>
#include <list>
#include <libgen.h>
#include "l2gen_aquarius.h"
#include <clo.h>

#ifndef MACINTOSH
#include <malloc.h>
#endif


#define DEG2RAD 0.017453293
#define TB_COS  3.0  
#define TEFF    290.0

#define VERSION "4.00"

//    Modification history:
// Programmer  Organization  Date      Version  Description of change
// ----------  ------------  ----      -------  ---------------------
// Joel Gales  FutureTech    04/02/08           Original development
// Joel Gales  FutureTech    08/02/09           Use gpstai2unix
// Joel Gales  FutureTech    08/10/09           GranStart/Stop reflect
//                                              actual first/last record times.
// Joel Gales  FutureTech    01/08/10           Add initial support for V2
// Joel Gales  FutureTech    03/17/10           Add new navigation fields
// Joel Gales  FutureTech    04/01/10           Switch to rad landmask,
//                                              Add ice flagging
// Joel Gales  FutureTech    04/02/10           Add scatterometer navigation
//                                              fields
// Joel Gales  FutureTech    05/05/10           Add support for browse flag
// Joel Gales  FutureTech    05/13/10           Trap no active blocks error
// Joel Gales  FutureTech    07/17/10           Fix readl1_caltemps rbe temps
// Joel Gales  FutureTech    07/17/10           Add "Input Files" &
//                                              "Processing Control" metadata
// Joel Gales  FutureTech    08/13/10           Add matchup delta lat/lon &
//                                              minimum distance
// Joel Gales  FutureTech    11/23/10           Initial update to V3
// Joel Gales  FutureTech    12/07/10           Implement new scatter products
// Joel Gales  FutureTech    12/08/10           Fix time bug in orb interp
// Joel Gales  FutureTech    12/09/10           Reduce size of scat prod arrays
// Joel Gales  FutureTech    12/13/10           Add support for scatterometer 
//                                              flags
// Joel Gales  FutureTech    12/15/10           Add support for WIND,TEMP,MOON
//                                              rad flags and scat_samples
// Joel Gales  FutureTech    12/21/10           Add Kpc_ANT fields
//                                              Fix rad_rfi_flags dim (2->4)
//                                              Add support for rad_samples
// Joel Gales  FutureTech    01/12/11           Check for improper products
// Joel Gales  FutureTech    01/28/11           Add support for gice
// Joel Gales  FutureTech    01/31/11           Fix RED WIND flag bug
// Joel Gales  FutureTech    02/16/11           Add sclon, sclat, scalt
// Joel Gales  FutureTech    02/16/11           Convert tasun, tagal, tamon
//                                              from C.S. to VH3 before write
// Joel Gales  FutureTech    02/18/11           Add support for SUNGLINT, 
//                                              GALACTIC flags
// Joel Gales  FutureTech    02/21/11           Add support for Solar Flux flag
// Joel Gales  FutureTech    03/11/11           Fix ta expected bugs
// Joel Gales  FutureTech    03/17/11           Make smoothed ice product
//                                              standard rad ice product
// Joel Gales  FutureTech    03/28/11           Use gland instead of fland
//                                              Use gice instead of fice
// Joel Gales  FutureTech    03/30/11           Process all blocks for
//                                              orbStart/orbStop = -999
// Joel Gales  FutureTech    04/13/11  0.965    Adjust lower yellow LAND
//                                              flag for gland mask,
//                                              Add ancillary files to metadata
// Joel Gales  FutureTech    05/07/11  0.970    Add parameter/suite file 
//                                              support
// Joel Gales  FutureTech    05/16/11  0.980    Add TOA/SSS land correction
// Joel Gales  FutureTech    06/02/11  0.985    Add scatterometer ancillary
//                                              files attribute
// Joel Gales  FutureTech    06/21/23  0.990    Add scatterometer process cntl
//                                              attribute, improved radiometer
//                                              process control, support for 
//                                              scat wind use in SSS
//
// Joel Gales  FutureTech    07/15/11  0.991    Use tbsur rather than tb_sur
//                                              for ta_expected calc
// Joel Gales  FutureTech    08/10/11  1.000    Update to V 1.00
// Joel Gales  FutureTech    08/10/11  1.010    If yancfile 2 & 3 not defined
//                                              then fill with yancfile1
// Joel Gales  FutureTech    08/22/11  1.020    Trap missing RAD frames
// Joel Gales  FutureTech    08/24/11  1.030    Add rfi-filtered TF data
//                                              products
// Joel Gales  FutureTech    08/25/11  1.040    Change index limit from 6 to
//                                              5 in rad_samples calc
// Joel Gales  FutureTech    08/25/11  1.050    Fix rad_rfi_flag write bug
//                                              (writel2radrfi() in l2_io.cpp)
// Joel Gales  FutureTech    09/02/11  1.060    Change cellonfoot & sclon
//                                              to -180 -> +180
// Joel Gales  FutureTech    09/02/11  1.060    Fix iopt_rfi to read 
//                                              "true/false" values correctly.
// Joel Gales  FutureTech    09/02/11  1.060    Fix spurious IHOUR problem in
//                                              get_ancillary_data() causing
//                                              navigation step to fail.
// Joel Gales  FutureTech    09/05/11  1.060    Report and change to -9999
//                                              NaN calibration temperatures
//                                              (readl1_caltemps() in 
//                                              l1_io.cpp)
// Joel Gales  FutureTech    09/09/11  1.010    Update to V1.1 SSS algorithm
//                                              Fix "Input Files" attribute
//                                              Remove "Product Center Time"
//                                              Report static files used 
// Joel Gales  FutureTech    09/09/11  1.11     Fix completely bonehead coding
//                                              mistake for cellonfoot range
//                                              fix.  (See 09/02/11) 
// Joel Gales  FutureTech    09/09/11  1.12     For TEMP flagging change
//                                              comparison from tb_toa to
//                                              ta_beam, correct index from 0
//                                              to 1 in taexpected[] for H-pol
// Joel Gales  FutureTech    09/27/11  1.13     Incorporate fixes in fd_tausq
//                                              and find_refl_tot (stop to
//                                              write) to allow navigation to
//                                              be generated.
// Joel Gales  FutureTech    10/12/11  1.14     Add gain/offset fields to
//                                              Converted Telemetry Group
//                                              Add NAV flag
// Joel Gales  FutureTech    10/20/11  1.15     Add support for rad ta drift
//                                              correction
//                                              Add support for updated 
//                                              code
//                                              JPL_scat_sw_delivery_10.17.2011
//                                              Fix delta_Tb write order bug
// Joel Gales  FutureTech    12/05/11  1.16     Fix indexing error in iflag_rfi
//                                              icnt -> iblk for
//                                              RFI flags/rad_samples
// Joel Gales  FutureTech    02/07/12  1.17     Initial implementation of DR cal
// Joel Gales  FutureTech    02/29/12  1.18     Fix cTND bug in DR cal
//                                              Roughness correction SSS V1.3
//                                              _ -> _nolc, _lc -> _
//                                              Add ROUGH flag
// Joel Gales  FutureTech    03/07/12  1.19     Remove unneeded gain/offset
//                                              fields
// Joel Gales  FutureTech    03/15/12  1.20     Add DR Kcoefficient attributes
//                                              rad_Tbx for rc,lc; rc,no_lc;
//                                              norc,lc; no_rc,no_lc
//                                              Comment out scatterometer flag 
//                                              check for roughness correction
// Joel Gales  FutureTech    03/22/12  1.21     Skip ancillary and SSS
//                                              processing for blocks with
//                                              off-earth navigation
//                                              Fix off-earth NAV bug
//                                              (Not set for beams 2 & 3)
// Joel Gales  FutureTech    03/31/12  1.22     Use NCEP if scat flag 29/31 set
//                                              or roughness corr table oob
//                                              Repackage the radiometer flag
//                                              Fix galactic flag limits
// Joel Gales  FutureTech    04/19/12  1.23     Add iopt_nosa1
//                                              Fix oob boresight (range) bug
//                                              in fd_cell_parameters()
// Joel Gales  FutureTech    04/19/12  1.24     Add support for SWE
//                                              Fix V/H reflect solar flag
// Joel Gales  FutureTech    05/04/12  1.25     Add support for scatterometer
//                                              V 1.3.1 code
//                                              Fix mis-dimensioned
//                                              tasun_bak_tab
//                                              (in fd_sun_backscatter)
// Joel Gales  FutureTech    05/27/12  1.26     Change FLUXD, FLUXR flagging
// Joel Gales  FutureTech    06/05/12  1.26     Add iopt_zero option
// Joel Gales  FutureTech    06/27/12  1.26     Add rpy_adj parameter
// Joel Gales  FutureTech    07/13/12  1.27     Implement rad flag limits table
//                                              and updated flag limits
// Joel Gales  FutureTech    07/20/12  1.28     Add and flag ACS mode in NAV
// Joel Gales  FutureTech    07/25/12  1.29     Add support for c_deltaTND
//                                              input, fix c_deltaTND[][]
//                                              ordering
// Joel Gales  Futuretech    08/08/12  1.30     Add support for solar flare flag
// Joel Gales  Futuretech    08/10/12  1.31     Modify orbit_phase_dif for V2.0
//                                              Adjust rad_samples, 
//                                              minSA for nosa1
//                           08/29/12           Add V2.0 part 2
// Joel Gales  Futuretech    09/06/12  1.32     Rearrange code to compute dtbc
//                                              for L2_CAL processing
// Joel Gales  Futuretech    09/11/12  1.33     Fix Ta/RFI bug in count_to_ta
//                                              Add new scatterometer model
//                                              functions
//                                              Add gpp,gmm,op,om fields to
//                                              Converted Telemetry
//                                              Remove taexpected goto for CAL
//                                              Add Tb expected
// Joel Gales  Futuretech    09/14/12  1.34     Add TA0/TF0 to SCI processing
//                                              Fix # of rad_samples bug
//                                              84/60 flipped for sa1/nosa1
//                                              Fix percentRFI
// Joel Gales  Futuretech    11/05/12  1.35     Add clo routines
//                                              Add support for subsurface
//                                              temp ancillary fields
// Joel Gales  Futuretech    11/05/12  1.35     Add support for radiometer
//                                              offset correction (day)
//                                              Apply sigma0 patch for
//                                              rad roughness corrrection
// Joel Gales  Futuretech    11/30/12  1.36     Add support for radiometer
//                                              offset correction (orbit)
// Joel Gales  Futuretech    12/12/12  1.37     Add gsl_interp_init()
//                                              Turn off gsl error handler
//                                              for attitude interpolation
// Joel Gales  Futuretech    01/25/13  2.00     Add anomaly_status param
//                                              Add CF metadata
// Joel Gales  Futuretech    02/11/13  2.01     Add support for solar flux
//                                              attributes and solar xray
//                                              block attribute
// Joel Gales  Futuretech    02/28/13  2.02     Add xrayfile parameters to xml
// Joel Gales  Futuretech    04/05/13  2.10     Add support for V2.1
//                                              HH/HHH winds, etc
// Joel Gales  Futuretech    04/08/13  2.11     Make climate salinity code
//                                              independent of WOA
// Joel Gales  Futuretech    04/18/13  2.12     Fix gice/gland test in HH/HHH
//                                              so it won't break out of ibeam
//                                              loop
// Joel Gales  Futuretech    05/05/13  2.13     Use HH rather than HHH winds
//                                              in TA_EXP to remove dependence
//                                              on TF.
// Joel Gales  Futuretech    08/02/13  2.20     Incorporate V2.3 RSS code
// Joel Gales  Futuretech    10/03/13  2.21     Add missing surtep/TEFF factor
//                                              for dtbc_hh, Add ta_expected
//                                              using HHH winds, Add symmetric
//                                              galactic dTb correction
// Joel Gales  Futuretech    12/02/13  2.30     Support for new flags:
//                                              POINTING, TBCONS, COLDWATER,
//                                              TFTADIFF, REFL_1STOKES,
//                                              RFI_REGION
// Joel Gales  Futuretech    12/13/13  2.31     Move POINTING flag code
//                                              before OOB continue
// Joel Gales  Futuretech    01/09/14  2.32     Fix WIND,TEMP,SA_OVERFLOW issues
//                                              Add anc_swh product
// Joel Gales  Futuretech    01/09/14  2.40     Modify galactic REFL_1STOKES
//                                              flag to use both HH winds and
//                                              1st stokes reflected galaxy
//                                              temperature
// Joel Gales  Futuretech    01/29/14  2.41     Fix bug in TFTADIFF and
//                                              galactic REFL_1STOKES flags
// Joel Gales  Futuretech    02/04/14  2.42     Set LAND/ICE severe flags
//                                              for mask values (0.5)
// Joel Gales  Futuretech    02/12/14  2.43     Add check for bad emissivity
//                                              of TBCONS flag, Add Tf-Ta>0.3
//                                              to TFTADIFF red flag
// Joel Gales  Futuretech    02/18/14  2.44     Compute HH/HHH winds for
//                                              gland/ice < 0.1
//                                              Add HH/HHH non-covergence to
//                                              wind flag, Add Tf-Ta>0.3 to
//                                              TFTADIFF yellow flag
// Joel Gales  Futuretech    02/19/14  2.45     Reverse HHH & HH R.C. to fix
//                                              NaN em0 values
// Joel Gales  Futuretech    02/20/14  2.46     Both HHH & HH R.C. must be good
//                                              for ioob=0
// Joel Gales  Futuretech    03/06/14  2.47     Add tagal_ref_GO products
// Joel Gales  Futuretech    03/11/14  2.48     Add dtagal_ref products
// Joel Gales  Futuretech    03/30/14  3.00     Refactor WIND and ROUGH flags
//                                              to better reflect error
//                                              conditions
// Joel Gales  Futuretech    03/31/14  3.01     Put back flagging on NaN em0
//                                              values in TBCONS
// Joel Gales  Futuretech    04/02/14  3.02     Fix LAND mask (> 0.5) flagcheck
//                                              Test for SWH = -999 with
//                                              floating point tolerance
// Joel Gales  Futuretech    04/14/14  3.03     Add support for geographic
//                                              rfi ta_nominal values
// Joel Gales  Futuretech    04/28/14  3.04     Add SSS[SST] adjustment
// Joel Gales  Futuretech    07/29/14  3.10     Add support for L1B l_acc
//                                              instrument-based correction
// Joel Gales  Futuretech    10/03/14  3.20     Add support for
//                                              instrument-based drift 
//                                              correction (TA_hat)
//                                              Change iopt_cal to iopt_l1b
//                                              Add support for generation of
//                                              instrument-based calibration 
//                                              L2 granules. 
// Joel Gales  Futuretech    10/25/14  3.21     Add support for
//                                              instrument-based drift 
//                                              correction (TA_hat) with offset
//                                              term
// Joel Gales  Futuretech    11/12/14  3.22     Add support for
//                                              instrument-based drift 
//                                              correction w/o gain
// Joel Gales  Futuretech    11/21/14  3.23     Copy dtbc to dtbc_hh if HH
//                                              wind retrieval succeeds but
//                                              HHH wind retrieval fails
// Joel Gales  Futuretech    12/18/14  3.24     Change scat_toi to scat_ant
// Joel Gales  Futuretech    02/10/15  3.30     Fix SWH bug:
//                                              if (abs(swh-(-999.0))<0.001) 
//                                              Add support for V3.4
//                                              demiss_sst_wspd correction
// Joel Gales  Futuretech    02/25/15  3.40     Rename to V3.40
// Joel Gales  Futuretech    02/26/15  3.41     Update igrf11syn.f for 2015
// Joel Gales  Futuretech    02/27/15  3.42     Add support for
//                                              rad_dtb_sst_wspd_V &
//                                              rad_dtb_sst_wspd_H fields
// Joel Gales  Futuretech    03/19/15  3.43     Don't extrapolate in l_acc_corr
// Joel Gales  Futuretech    04/14/15  3.44     Implement density and
//                                              spiciness code, Non-linear
//                                              UI coupling
// Joel Gales  Futuretech    04/22/15  3.45     Add support for uncertainties,
//                                              UI coupling code for 
//                                              ta_expected,
//                                              Add RFI TD (A & D) files
// Joel Gales  Futuretech    05/15/15  3.46     Add UI coupling code for
//                                              ta_expected HHH winds,
//                                              Add support for static
//                                              bias adjustment to SSS,
//                                              Fix asc/desc determination
//                                              for SSS unc
// Joel Gales  Futuretech    05/18/15  4.00     Tag Version 4.0

extern "C" double zulu2unix( char *zulu);
extern "C" double gsw_sa_from_sp( double sp, double p, double lon, double lat);
extern "C" double gsw_ct_from_t(double sa, double t, double p);
extern "C" double gsw_rho_t_exact(double sa, double t, double p);
extern "C" double gsw_spiciness0(double sa, double ct);

using namespace std;

// Comparison function for qsort
int compare (const void * a, const void * b)
{
  float v = *(float*)a - *(float*)b;
  if ( v < 0.0) return -1;
  if ( v > 0.0) return +1;
  return 0;
}

int main (int argc, char* argv[])
{
  char tmpStr[256];
  sprintf(tmpStr, "%s %s (%s %s)", "l2gen_aquarius", 
	  VERSION, __DATE__, __TIME__);
  clo_setVersion(tmpStr);

  if ( argc == 1) {
    cout << tmpStr << endl;
    return 0;
  }

  float zerof = 0.0;

  // blkSec is seconds from 01/06/1980
  double blkSec;

  int32_t iblk, ibeam, first_blk, last_blk;

  // Coefficients need to do a very small adjustment to incidence angle to get 
  // effective incidence angle
  static float tht_coef[3] = {1.0017714, 1.0018588, 1.0014806};
  float freq_aq = 1.413;

  static float dir_coef[4*2] = {0,0,0,0,0,0,0,0};

#include "l2gen_aquarius_prod.inc"

  instr l2genInput;
  uint32_t nProd;
  string processControl;
  parseInput(argc, argv, &l2genInput, &nProd, &processControl);

  if ( l2genInput.iopt_rfi == 0) {
    cout << endl <<
      "Non-RFI filtered antenna temperatures used for SSS retrieval" << endl;
  } else {
    cout << endl <<
      "RFI filtered antenna temperatures used for SSS retrival" << endl;
  }

  // If yancfile2 and yancfile3 not define then fill with yancfile1
  if ( strcmp(l2genInput.yancfile2, "") == 0) {
    strcpy( l2genInput.yancfile2, l2genInput.yancfile1);
  }

  if ( strcmp(l2genInput.yancfile3, "") == 0) {
    strcpy( l2genInput.yancfile3, l2genInput.yancfile1);
  }


  // Open L1A input file
  static Hdf::hdf5_Aquarius l1afile;
  l1afile.openl1(l2genInput.ifile, H5F_ACC_RDONLY);

  int32_t nBlks = l1afile.nBlks();
  if ( nBlks > MAXCYC) {
    cout << "nBlks: " << nBlks << " greater than MAXCYC: " << MAXCYC << endl;
    exit(7);
  }

  string inputFiles;
  // Change to l2genInput.ifile from "ifile="  JMG  09/09/11 
  inputFiles = basename( l2genInput.ifile);


  //  double granuleStart, granuleStop;
  //l1afile.getGranuleTimes( &granuleStart, &granuleStop);

  double orbitStart, orbitStop;
  int32_t orbitNumber, cycleNumber, passNumber;
  l1afile.getOrbitTimes( &orbitStart, &orbitStop, 
			 &orbitNumber, &cycleNumber, &passNumber);

  time_t rawtime;
  tm * ptm;
  char buffer [80];

  rawtime = (uint32_t) gpstai2unix( orbitStart);
  ptm = gmtime ( &rawtime );
  strftime (buffer, 80, "%Y-%m-%d %H:%M:%S", ptm);
  cout << endl;
  cout << "OrbitStart: " << buffer << endl;

  rawtime = (uint32_t) gpstai2unix (orbitStop);
  ptm = gmtime ( &rawtime );
  strftime (buffer, 80, "%Y-%m-%d %H:%M:%S", ptm);
  cout << "OrbitStop:  " << buffer << endl;
  cout << endl;

  double nodeCrossingTime;
  float nodeLongitude;
  l1afile.getNodeInfo( &nodeCrossingTime, &nodeLongitude);

 
  // Generate ACS mode
  getHKT( &l1afile, acsmode);
 
  // Generate geolocation
  double *Pos, *Vel, *rpy;
  Pos = new double[MAXCYC*3];
  Vel = new double[MAXCYC*3];
  rpy = new double[MAXCYC*3];
  getGeoNavSun( &l1afile, l2genInput.rpy_adj, 
		cellon, cellat, celtht, celphi,
		suntht, sunphi, sunglt, moonglt,
		glxlon, glxlat,
		zang, sun_zenith, sclon, sclat, scalt,
		cellonfoot, cellatfoot,
		sund, sunr, moond, moonr, bore_sight,
		Pos, Vel, rpy);
  //  int32_t landmask = 0;

  //  char *ocdataroot_str = getenv("OCDATAROOT");
  //if (ocdataroot_str == 0x0) {
  //  printf("Environment variable OCDATAROOT not defined.\n");
  //  exit(1);
  //}

  bool inout[MAXSCAN];

  // Get extract limits if matchup
  if ( l2genInput.matchup_lat != -999 &&
       l2genInput.matchup_lon != -999) {

    float clat2 = l2genInput.matchup_lat * DEG2RAD;
    float clon2 = l2genInput.matchup_lon * DEG2RAD;

    for (iblk=0; iblk<nBlks; iblk++) {
      //      printf("iblk: %d\n", iblk);
      // Skip empty blks
      l1afile.readl1_radiometer(iblk, 0, 0, &blkSec, NULL,
				NULL, NULL, NULL);
      inout[iblk] = false;

      // Skip if outside orbit
      if ( blkSec == -1 || zang[iblk] == -999 || 
	   blkSec < orbitStart || blkSec >= orbitStop) continue;

      // Check blk time is within 1000 sec of matchup time;
      //      if ( l2genInput.matchup_time != -1) {
      //	float diff = blkSec - l2genInput.matchup_time;
      //	if ( diff < -1000 || diff > 1000) 
      //	  continue;
      //}

      float north=-90;
      float south=+90;
      float east=-180;
      float west=+180;
      for (ibeam=0; ibeam<NUMBER_OF_BEAMS; ibeam++) {
	for (size_t i=0; i<4; i++) {
	  float lat = cellatfoot[4*(ibeam+NUMBER_OF_BEAMS*iblk)+i];
	  float lon = cellonfoot[4*(ibeam+NUMBER_OF_BEAMS*iblk)+i];
	  if ( lon > 180) lon -= 360;

	  if ( lat > north) north = lat;
	  if ( lat < south) south = lat;
	  if ( lon > east)  east  = lon;
	  if ( lon < west)  west  = lon;
	}
      } // ibeam loop 


      north += l2genInput.matchup_delta_lat;
      south -= l2genInput.matchup_delta_lat;
      east  += l2genInput.matchup_delta_lon;
      west  -= l2genInput.matchup_delta_lon;

      if ( (east - west) > 90) {
	if ( l2genInput.matchup_lat > +80 || l2genInput.matchup_lat < -80) {
	  east = -180;
	  west = +180;
	}
	if ( l2genInput.matchup_lat <= north &&
	     l2genInput.matchup_lat >= south &&
	     (l2genInput.matchup_lon >= east ||
	     l2genInput.matchup_lon <= west)) {
	  cout << right << fixed;
	  cout.precision( 3);
	  cout << "lat: " << l2genInput.matchup_lat << " lon: " << 
	    l2genInput.matchup_lon;
	  cout <<  " -box: " << north << " " << south << " "
	       << east << " " << west << " " << endl;

	  // Check for minimum distance
	  for (ibeam=0; ibeam<NUMBER_OF_BEAMS; ibeam++) {
	    int ibeamblk = ibeam+NUMBER_OF_BEAMS*iblk;

	    float clat1 = cellat[ibeamblk] * DEG2RAD;
	    float clon1 = cellon[ibeamblk] * DEG2RAD;

	    float dis = acos(sin(clat1)*sin(clat2) + 
			     cos(clat1)*cos(clat2)*cos(clon2-clon1))*6371;

	    if ( dis < l2genInput.matchup_min_dist) 
	      inout[iblk] = true;
	  }
	}
     } else {
	if ( l2genInput.matchup_lat > +80 || l2genInput.matchup_lat < -80) {
	  east = +180;
	  west = -180;
	}
	if ( l2genInput.matchup_lat <= north &&
	     l2genInput.matchup_lat >= south &&
	     l2genInput.matchup_lon <= east &&
	     l2genInput.matchup_lon >= west) {
	  cout << right << fixed;
	  cout.precision( 3);
	  cout << "lat: " << l2genInput.matchup_lat << " lon: " << 
	    l2genInput.matchup_lon;
	  cout << " iblk: " << iblk;
	  cout << " +box: " << north << " " << south << " "
	       << east << " " << west << " " << endl;

	  // Check for minimum distance
	  for (ibeam=0; ibeam<NUMBER_OF_BEAMS; ibeam++) {
	    int ibeamblk = ibeam+NUMBER_OF_BEAMS*iblk;

	    float clat1 = cellat[ibeamblk] * DEG2RAD;
	    float clon1 = cellon[ibeamblk] * DEG2RAD;

	    float dis = acos(sin(clat1)*sin(clat2) + 
			     cos(clat1)*cos(clat2)*cos(clon2-clon1))*6371;

	    if ( dis < l2genInput.matchup_min_dist) 
	      inout[iblk] = true;
	  }
	}
      }
    } // iblk loop
    bool first = true;
    for (iblk=0; iblk<nBlks; iblk++) {
      if ( inout[iblk] == true) {
	if ( first == true) first_blk = iblk;
	last_blk = iblk;
	first = false;
      }
    }
    if ( first == true) {
      cout << endl << "Matchup not found." << endl;
      exit(110);
    }

    // Pad ends by 10 blks
    first_blk -= 10;
    if ( first_blk < 0) first_blk = 0;
    last_blk += 10;
    if ( last_blk >= nBlks) last_blk = nBlks-1;
    for (iblk=first_blk; iblk<=last_blk; iblk++) {
      inout[iblk] = true;
    }


    // Create output limits (1-based) text file
    if ( l2genInput.matchup_lim_file[0] != 0) {
      string str = l2genInput.matchup_lim_file;
      ofstream fout;
      fout.open(str.c_str());
      fout << first_blk+1 << endl;
      fout << last_blk+1 << endl;
      fout.close();

      l1afile.closel1();
      return 0;
    }

    // End matchup section
  } else {
    // Not Matchup run
    // Clip off "wings" at ends of orbit
    for (iblk=0; iblk<nBlks; iblk++) {
      l1afile.readl1_radiometer(iblk, 0, 0, &blkSec, NULL,
				NULL, NULL, NULL);
      if ( orbitStart == -999 && orbitStop == -999)
	inout[iblk] = true; 
      else if ( blkSec == -1 || zang[iblk] == -999 || 
		blkSec < orbitStart || blkSec >= orbitStop)
	inout[iblk] = false; 
      else
	inout[iblk] = true; 
    }
  }


  // Get number of active blocks
  // Find first record time of L2 granule
  int32_t nactBlks = 0;
  bool notFound = true;
  for (iblk=0; iblk<nBlks; iblk++) {
    if ( inout[iblk] == true) {
      if ( notFound) {
	l1afile.readl1_radiometer(iblk, 0, 0, &blkSec, NULL, NULL, NULL, NULL);
	notFound = false;
      }
      nactBlks++;
    }
  } 

  // If no active blocks then exit with error
  if ( nactBlks == 0) {
    cout << "*** No active blocks found." << endl;
    exit(1);
  }


  // Get scatter_basename
  string str;
  string scatter_basename;
  string scatter_textfile;
  if ( l2genInput.scatter_basename[0] == 0) {
    str = basename( l2genInput.ifile);
    size_t pos = str.find("."); 
    scatter_basename = str.substr( 0, pos);  
  } else {
    scatter_basename = l2genInput.scatter_basename;
  }

  // Get Scatterometer process control string
  ifstream in_SCAT_PROCCTRL;
  string sScatProcCtrl="";
  scatter_textfile = scatter_basename + "_L2_process_control.txt";
  in_SCAT_PROCCTRL.open ( scatter_textfile.c_str(), ifstream::in);
  if ( in_SCAT_PROCCTRL.fail() == false) {
    getline( in_SCAT_PROCCTRL, sScatProcCtrl);
    in_SCAT_PROCCTRL.close();
  }

  // Create L2 output file
  static Hdf::hdf5_Aquarius l2file;
  double granStart = blkSec;
  l2file.createl2(l2genInput.ofile, nactBlks, nProd, 
		  l2genInput.l2activeprodflag, granStart,
		  nodeCrossingTime, nodeLongitude, orbitNumber,
		  cycleNumber, passNumber,
		  (char *) inputFiles.c_str(), 
		  (char *) processControl.c_str(),
		  (char *) sScatProcCtrl.c_str(),
		  (char *) VERSION, l2genInput.pversion, l2genInput.iopt_l1b);

  hid_t gid[7];
  l2file.getH5gid( l2file.getH5fid(), gid);


  ///////////////////////////////////////////////////////////////////////
  ////////////////////////// Read Rad Flag Limits ///////////////////////
  ///////////////////////////////////////////////////////////////////////
  istringstream istr;
  ifstream radFlag;

  radFlag.open ( l2genInput.radflaglimitsfile, ifstream::in);
  if ( radFlag.fail() != 0) {
    cout << endl << "Radiometer Flag Limits File: " <<
      l2genInput.radflaglimitsfile << " not found." << endl;
    exit(1);
  }
  char txtbuf[512];
  ostringstream radLimits;
  radLimits.str("");

  long limitRFI[2];
  radFlag.getline( txtbuf, 512); istr.str( string(txtbuf) + " ");
  istr >> limitRFI[0];
  istr >> limitRFI[1];
  radLimits << "RFI: " << limitRFI[0] << " " << limitRFI[1] << " ";

  long limitRain[2];
  radFlag.getline( txtbuf, 512); istr.str( string(txtbuf) + " ");
  //  istr >> limitRain[0];
  //istr >> limitRain[1];

  float limitLand[2];
  radFlag.getline( txtbuf, 512); istr.str( string(txtbuf) + " ");
  istr >> limitLand[0];
  istr >> limitLand[1];
  radLimits << "Land: " << limitLand[0] << " " << limitLand[1] << " ";

  float limitIce[2];
  radFlag.getline( txtbuf, 512); istr.str( string(txtbuf) + " ");
  istr >> limitIce[0];
  istr >> limitIce[1];
  radLimits << "Ice: " << limitIce[0] << " " << limitIce[1] << " ";

  float limitWind[2];
  radFlag.getline( txtbuf, 512); istr.str( string(txtbuf) + " ");
  istr >> limitWind[0];
  istr >> limitWind[1];
  radLimits << "Wind: " << limitWind[0] << " " << limitWind[1] << " ";

  float limitTemp[2];
  radFlag.getline( txtbuf, 512); istr.str( string(txtbuf) + " ");
  istr >> limitTemp[0];
  istr >> limitTemp[1];
  radLimits << "Temp: " << limitTemp[0] << " " << limitTemp[1] << " ";

  float limitFluxD[2];
  radFlag.getline( txtbuf, 512); istr.str( string(txtbuf) + " ");
  istr >> limitFluxD[0];
  istr >> limitFluxD[1];
  radLimits << "FluxD: " << limitFluxD[0] << " " << limitFluxD[1] << " ";

  float limitFluxR[2];
  radFlag.getline( txtbuf, 512); istr.str( string(txtbuf) + " ");
  istr >> limitFluxR[0];
  istr >> limitFluxR[1];
  radLimits << "FluxR: " << limitFluxR[0] << " " << limitFluxR[1] << " ";

  float limitGlint[2];
  radFlag.getline( txtbuf, 512); istr.str( string(txtbuf) + " ");
  istr >> limitGlint[0];
  istr >> limitGlint[1];
  radLimits << "Glint: " << limitGlint[0] << " " << limitGlint[1] << " ";

  float limitMoon[2];
  radFlag.getline( txtbuf, 512); istr.str( string(txtbuf) + " ");
  istr >> limitMoon[0];
  istr >> limitMoon[1];
  radLimits << "Moon: " << limitMoon[0] << " " << limitMoon[1] << " ";

  float limitGal[2];
  radFlag.getline( txtbuf, 512); istr.str( string(txtbuf) + " ");
  istr >> limitGal[0];
  istr >> limitGal[1];
  radLimits << "Gal: " << limitGal[0] << " " << limitGal[1] << " ";

  float limitRPY[3];
  radFlag.getline( txtbuf, 512); istr.str( string(txtbuf) + " ");
  istr >> limitRPY[0];
  istr >> limitRPY[1];
  istr >> limitRPY[2];
  radLimits << "RPY: " << limitRPY[0] << " " << limitRPY[1] << " "
	    << limitRPY[2] << " ";

  float limitFlare[2];
  radFlag.getline( txtbuf, 512); istr.str( string(txtbuf) + " ");
  istr >> limitFlare[0];
  istr >> limitFlare[1];
  radLimits << "Flare: " << limitFlare[0] << " " << limitFlare[1] << " ";

  float limitTbcons;
  radFlag.getline( txtbuf, 512); istr.str( string(txtbuf) + " ");
  istr >> limitTbcons;
  radLimits << "Tb cons: " << limitTbcons << " ";

  float limitColdWater[2];
  radFlag.getline( txtbuf, 512); istr.str( string(txtbuf) + " ");
  istr >> limitColdWater[0];
  istr >> limitColdWater[1];
  radLimits << "Cold Water: " << limitColdWater[0] << " " << limitColdWater[1] << " ";

  float limitTFTA[4];
  radFlag.getline( txtbuf, 512); istr.str( string(txtbuf) + " ");
  istr >> limitTFTA[0];
  istr >> limitTFTA[1];
  radFlag.getline( txtbuf, 512); istr.str( string(txtbuf) + " ");
  istr >> limitTFTA[2];
  istr >> limitTFTA[3];
  radLimits << "TF-TA: " 
	    << limitTFTA[0] << " " << limitTFTA[1] << " "
	    << limitTFTA[2] << " " << limitTFTA[3] << " ";

  float limit1stStokes[5];
  radFlag.getline( txtbuf, 512); istr.str( string(txtbuf) + " ");
  istr >> limit1stStokes[0];
  istr >> limit1stStokes[1];
  radFlag.getline( txtbuf, 512); istr.str( string(txtbuf) + " ");
  istr >> limit1stStokes[2];
  istr >> limit1stStokes[3];
  radFlag.getline( txtbuf, 512); istr.str( string(txtbuf) + " ");
  istr >> limit1stStokes[4];
  radLimits << "Refl 1st Stokes (moon): " 
	    << limit1stStokes[0] << " " << limit1stStokes[1] << " "
	    << "Refl 1st Stokes (galaxy): " 
	    << limit1stStokes[2] << " " << limit1stStokes[3] << " "
	    << "(wind): " 
	    << limit1stStokes[4];

  radFlag.close();


  ///////////////////////////////////////////////////////////////////////
  ///////////////////////// Initialize Static Data //////////////////////
  ///////////////////////////////////////////////////////////////////////
  cout << "Initializing static data" << endl;
  //static static_data staticData;
  static_data *staticData = new static_data;

  initialize_static_data( &l2genInput, staticData);


  ///////////////////////////////////////////////////////////////////////
  /////////////////////// Load SSS static bias adj //////////////////////
  ///////////////////////////////////////////////////////////////////////
  string bias_adj_file = l2genInput.static_bias_adj_file;

  double sss_static_bias_adj[2][3][180][360];

  if ( bias_adj_file != "") {

    hid_t h5fid, grp;
    string filename = bias_adj_file;
    expandEnvVar( &filename);
    h5fid =  H5Fopen( filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if ( h5fid == -1) {
      cout << "Unable to open SSS bias adjustment file: " << 
        filename.c_str() << endl;
      exit(1);
    }

    grp = H5Gopen1(h5fid,"/");

    hsize_t start[4]={0,0,0,0}, count[4]={2,3,180,360};

    Hdf::h5d_read(grp, (char*) "bias_adjustment", 
                  &sss_static_bias_adj[0][0][0][0],
                  4, start, count);

    H5Gclose(grp);
    H5Fclose(h5fid);
  }

  ///////////////////////////////////////////////////////////////////////
  /////////////////////////// Get Ancillary Data ////////////////////////
  ///////////////////////////////////////////////////////////////////////
  cout << endl << "Getting ancillary data" << endl;

  //     solar_flux = solar flux (sfu)
  //     fland  = fraction of land within 3db footprint (0 to 1.0)
  //     gland  = fraction of land over entire gain pattern (0 to 1.0)
  //     fice   = fraction of sea ice at boresight (0 to 1.0)
  //     gice   = fraction of sea ice over entire gain pattern (0 to 1.0)
  //     sss_reference = salinity HYCOM model (psu)
  //     surtep = ncep surface temp (kelvin)
  //     subtep = ncep sub-surface temp (kelvin)
  //     cwat = ncep cloud water (kg m-2)
  //     sm = ncep soil moisture
  //     winspd_ncep = ncep wind speed (m/s) value has been increased by 3% to 
  //     match buoys and satellite wind retreivals 
  //     windir = ncep wind direction in meteorlogical sense (out of) (deg)
  //     tran   = one-way atmospheric transmission (0 to 1)
  //     tbup   = upwelling   atmospheric brightness temperature (k)
  //     tbdw   = downwelling atmospheric brightness temperature (k)

  // Open first solar xray file
  string filename = l2genInput.xrayfile1;
  hid_t h5fid, grp;
  hsize_t zero=0, xray_cnt=1440;
  float xray[1440];
  float solar_attr[4]={-9999,0,-9999,0};
  int32_t n_solar=0, n_xray=0;

  if ( filename != "") {
    expandEnvVar( &filename);
    cout << setw(35) << "Reading Solar Xray file: " << filename.c_str() << endl;

    h5fid =  H5Fopen( filename.c_str(),  H5F_ACC_RDONLY, H5P_DEFAULT);
    if ( h5fid == -1) {
      cout << "Solar Xray file: " << filename.c_str() << " not found." << endl;
      exit(1);
    }
    grp = H5Gopen1(h5fid,"/");
    Hdf::h5d_read(grp, (char*) "1min_xray_irradiance", xray, 
		  1, &zero, &xray_cnt);
    H5Gclose(grp);
    H5Fclose(h5fid);
  } else {
    cout << "No solar xray file1 provided" << endl;
    for (size_t i=0; i<xray_cnt; i++) xray[i] = -9999.0;
  }


  // Get seconds from 01/01/2000 & ancillary data
  double *sec2000;
  int lyear, idayjl, isecdy, idum;
  float lsecdy=-1;
  double secdy, dumd;
  sec2000 = new double[nBlks];
  float percentWater = 0.0;
  for (iblk=0; iblk<nBlks; iblk++) {

    if ( l2genInput.iopt_l1b) continue;

    l1afile.readl1_radiometer(iblk, 0, 0, &blkSec, NULL,
			      NULL, NULL, NULL);
    if ( blkSec == -1) continue;

    sec2000[iblk] = gpstai2utc2000( blkSec);
    fd_date_2000_(&sec2000[iblk], &dumd, &lyear, &idayjl, 
		  &idum, &idum, &secdy);
    isecdy = (int) (secdy + 0.5);

    // Load new solar xray file if day change
    if ( secdy < lsecdy) {
      filename = l2genInput.xrayfile2;
      if ( filename != "") {
	expandEnvVar( &filename);
	cout << setw(35) << "Reading Solar Xray file: " << 
	  filename.c_str() << endl;
	h5fid =  H5Fopen( filename.c_str(),  H5F_ACC_RDONLY, H5P_DEFAULT);
	if ( h5fid == -1) {
	  cout << "Solar Xray file: " << filename.c_str() << 
	    " not found." << endl;
	  exit(1);
	}
	grp = H5Gopen1(h5fid,"/");
	Hdf::h5d_read(grp, (char*) "1min_xray_irradiance", xray, 
		      1, &zero, &xray_cnt);
	H5Gclose(grp);
	H5Fclose(h5fid);
      } else {
	cout << "No solar xray file2 provided" << endl;
	for (size_t i=0; i<xray_cnt; i++) xray[i] = -9999.0;
      }
    }
    lsecdy = secdy;

    for (ibeam=0; ibeam<NUMBER_OF_BEAMS; ibeam++) {

      int ibeamblk = ibeam+NUMBER_OF_BEAMS*iblk;

      float clat = cellat[ibeamblk];
      float clon = cellon[ibeamblk];
      float ctht = celtht[ibeamblk];

      float clon360 = clon;
      if ( clon != -999)
	if (clon360 < 0) clon360 += 360;

      if ( clon == -999 || clat == -999) continue;

      int32_t irad = ibeam + 1;

      get_ancillary_data_(&lyear, &idayjl, &isecdy, 
			  &clat, &clon360, &ctht,
			  &surtep[ibeamblk], &surp[ibeamblk], &cwat[ibeamblk], 
			  &swe[ibeamblk], &sm[ibeamblk], &subtep[ibeamblk], 
			  &winspd_ncep[ibeamblk], &windir[ibeamblk], 
			  &tran[ibeamblk], &tbup[ibeamblk], &tbdw[ibeamblk], 
			  &solar_flux[iblk],
			  &irad, &zang[iblk], &sclon[iblk],    // in
			  &staticData->fpt_lnd[0][0][0],        // in
			  &staticData->frc_lnd[0][0][0],        // in
			  &staticData->landcorr[0][0][0][0][0], // in
			  &fland[ibeamblk], &gland[ibeamblk],  // out
			  &fice[ibeamblk], &gice[ibeamblk],    // out
			  &swh[ibeamblk],                      // out
			  &sss_reference[ibeamblk],            // out
			  &tb_landcorr_vpol[ibeamblk],         // out
			  &tb_landcorr_hpol[ibeamblk],         // out
			  l2genInput.yancfile1,
			  l2genInput.yancfile2,
			  l2genInput.yancfile3);

      percentWater += (1 - gland[ibeamblk]);

      //      cout << "fland/fice: " << 
      //  setw(12) << setprecision(3) << fixed << secdy+0.72 <<
      //  setw(4) << right << ibeam << 
      //  setw(10) << setprecision(5) << fland[ibeamblk] << 
      //  setw(10) << setprecision(5) << fice[ibeamblk] << endl;

    } // ibeam loop

    // Assign block solar_xray value
    solar_xray[iblk] = xray[isecdy/60];

  } // iblk loop
  percentWater /= nBlks * NUMBER_OF_BEAMS;



  //////////////////////////////////////////////////////////////////////////
  ////////////////////////// Calibration ///////////////////////////////////
  //////////////////////////////////////////////////////////////////////////
  cout << endl << "Calibrating Instrument" << endl;

  uint16_t radiomEar[RADIOMETER_SUBCYCLES][RADIOMETER_SIGNALS_PER_SUBCYCLE]
    [NUMBER_OF_BEAMS][RADIOMETER_POLARIZATIONS];

  uint16_t radiomCnd[RADIOMETER_SUBCYCLES][NUMBER_OF_BEAMS]
    [RADIOMETER_POLARIZATIONS];

  uint16_t radiomCal[RADIOMETER_LONG_ACCUM][NUMBER_OF_BEAMS]
    [RADIOMETER_POLARIZATIONS];

  float *S_acc_raw;
  float *L_acc_raw;
  float *S_acc;
  float *L_acc;
  S_acc_raw = new float[MAXCYC*(NUMRADIOEAR+NUMRADIOCND)];
  S_acc = new float[MAXCYC*(NUMRADIOEAR+NUMRADIOCND)];

  L_acc_raw = new float[MAXCYC*NUMRADIOCAL];
  L_acc = new float[MAXCYC*NUMRADIOCAL];
  
  float *calTemps;
  calTemps = new float[MAXCYC*NUMCALTEMPS];

#ifndef MACINTOSH
  struct mallinfo minfo;
  minfo = mallinfo();
  //printf("Used space: %10d\n", minfo.arena+minfo.hblkhd);
#endif

  int32_t atc_frmnum, rad_frmnum;
  uint8_t atc_subframe, rad_subframe;

  //  int32_t last_subframe;
  //int32_t last_frmnum;

  if ( l2genInput.iopt_nosa1 == true) {
    cout << "No SA1 counts used." << endl;
  }


  int k = 0, l = 0;
  for (iblk=0; iblk<nBlks; iblk++) {
    //    cout << "iblk: " << iblk << endl;

    l1afile.readl1_frame_info( iblk,
			       &atc_frmnum, &atc_subframe,
			       &rad_frmnum, &rad_subframe);

    if ( rad_frmnum == -1) continue;

    l1afile.readl1_radiometer( iblk, rad_frmnum, rad_subframe, NULL, NULL,
			       &radiomEar[0][0][0][0],
			       &radiomCnd[0][0][0], 
			       &radiomCal[0][0][0]); 

    for (size_t irad=0; irad<NUMBER_OF_BEAMS; irad++) {
      for (size_t ipol=0; ipol<RADIOMETER_POLARIZATIONS; ipol++) {

	// Re-arrange short accum
	for (size_t isub=0; isub<RADIOMETER_SUBCYCLES; isub++) {
	  for (size_t i=0; i<5; i++) {
                if ( l2genInput.iopt_nosa1 == 1 && i == 0)
              S_acc_raw[k++] = -1.0;
            else if (i < 2)
              S_acc_raw[k++] = radiomEar[isub][i][irad][ipol] / 2.0;
            else
              S_acc_raw[k++] = radiomEar[isub][i][irad][ipol];
    	  }
    	  S_acc_raw[k++] = radiomCnd[isub][irad][ipol];
	}

	// Re-arrange long accum
	for (size_t ilng=0; ilng<RADIOMETER_LONG_ACCUM; ilng++) {
          if (ilng < 4)
            L_acc_raw[l++] = radiomCal[ilng][irad][ipol] / 10.0;
          else
            L_acc_raw[l++] = radiomCal[ilng][irad][ipol] / 2.0;
	}
      } // ipol loop
    }  // irad loop

    l1afile.readl1_caltemps( iblk, rad_frmnum, rad_subframe, 
			     &calTemps[iblk*NUMCALTEMPS]);

    //    last_subframe = rad_subframe;
    //last_frmnum = rad_frmnum;
  } // iblk loop


  // V/P/M/H
  rfi_data rfi;
  float receiver_noise=74.6; // Kelvin
  float bw_tau=474.34; // sqrt(bandwidth*integration time)=sqrt( 25 MHz * 0.009 sec)

  if ( l2genInput.ta_nominal_files[0] != 0) {
    ifstream in_ta_nominal_filelist;
    ifstream in_ta_nominal;
    string line;

    in_ta_nominal_filelist.open( l2genInput.ta_nominal_files, ifstream::in);
    if ( in_ta_nominal_filelist.fail() == false) {
      for (size_t i=0; i<NUMBER_OF_BEAMS*RADIOMETER_POLARIZATIONS; i++) {
        int irad = i/RADIOMETER_POLARIZATIONS;
        int ipol = i - irad*RADIOMETER_POLARIZATIONS;

        getline( in_ta_nominal_filelist, line);

        if ( line.substr( 0, 1).compare("$") == 0) expandEnvVar( &line);

        cout << "Reading " << line.c_str() << endl;
        in_ta_nominal.open( line.c_str(), ifstream::in);
        if ( in_ta_nominal.fail() == false) {
          // read the entire line into a string
          for (size_t ilon=0; ilon<360; ilon++) {
            getline( in_ta_nominal, line);

            // now we'll use a stringstream to separate the fields in the line
            stringstream ss( line);
            string field;
            int ilat=0;
            while (getline( ss, field, ',' )) {
              // for each field we wish to convert it to a float
              stringstream fs( field);
              float f;
              fs >> f;
              rfi.stdta_rad[irad][ipol][ilon][ilat] = 
                (f + receiver_noise) / bw_tau;
              //cout << j << " " << i << " " << f << endl;
              ilat++;
            }
          }
        } else {
          cout << "Cannot open TA nominal file: " << line.c_str() << endl;
          exit(1);
        }
        in_ta_nominal.close();
      } // file loop
    } else {
      cout << "Cannot open TA nominal filelist: " << 
        l2genInput.ta_nominal_files << endl;
      exit(1);
    }
    in_ta_nominal_filelist.close();
  } else {

    // Load RFI ocean nom values
    for (size_t i=0; i<NUMBER_OF_BEAMS*RADIOMETER_POLARIZATIONS; i++) {

      int irad = i/RADIOMETER_POLARIZATIONS;
      int ipol = i - irad*RADIOMETER_POLARIZATIONS;
      for (size_t ilat=0; ilat<181; ilat++) {
        for (size_t ilon=0; ilon<360; ilon++) {
          rfi.stdta_rad[irad][ipol][ilon][ilat] = 
            (l2genInput.ta_ocean_nom[3*ipol+irad] + receiver_noise) / bw_tau;
        }
      }
    }
  }

  itemp instrTemps;
  gainoff gainOff;
  float c_deltaTND[NUMBER_OF_BEAMS][2];

  // Generate L_acc correction if l_acc_corr_files specified
  // This is the instrument-based wiggle correction provied by S. Misra
  filename = l2genInput.l_acc_corr_files;
  if ( filename != "") {
    cout << endl << 
      "Applying instrument-based corrections to long accumulation counts" 
         << endl;
    //    l2file.writel2_s_acc( nBlks, "S_acc_raw", "Raw Antenna Counts",
    //                    S_acc_raw);
    if ( l2genInput.iopt_l1b)
      l2file.writel2_l_acc( nBlks, "L_acc_raw", 
                            "Long Accumulations Raw",
                            L_acc_raw);

    // Compute TA0/TF0 (without exponential drift correction)
    if ( l2genInput.c_deltaTND[0][0][0] == 0.0) {
      for (size_t irad=0; irad<NUMBER_OF_BEAMS; irad++) {
        for (size_t ipol=0; ipol<2; ipol++) { // 0=V,1=H
          c_deltaTND[irad][ipol] = 0.0;
        }
      }

      count_to_ta_( l2genInput.coeff_loss_file, l2genInput.coeff_nl_file, 
                    &nBlks, calTemps, &instrTemps, &gainOff, &rfi, 
                    S_acc_raw, L_acc_raw, S_acc, L_acc, 
                    cellat, cellon, zang, &c_deltaTND[0][0], TA_hat_drift_corr,
                    TA_hat, TA0, TF_hat, TF0);
    }

    // Apply l_acc_corr to L_acc_raw
    l_acc_corr( filename, L_acc_raw);

    if ( l2genInput.iopt_l1b)
      l2file.writel2_l_acc( nBlks, "L_acc_load_corr", 
                            "Long Accumulations Load_Corrected",
                            L_acc_raw);
  }  // if ( l2genInput.l_acc_corr_files != "")

#if 0
  // Compute instrument-based gain (drift) correction from S. Misra
  cout << "Applying instrument-based drift correction to TA_hat" << endl;
  for (size_t icyc=0; icyc<MAXCYC; icyc++) {
    for (size_t irad=0; irad<NUMBER_OF_BEAMS; irad++) {
      for (size_t ipol=0; ipol<RADIOMETER_POLARIZATIONS; ipol++) { // V/P/M/H
        if ( ipol == 1 || ipol == 2) continue;

        int indx = ipol / 3;

        if ( icyc == 0) {
          if ( ipol == 0) {
            cout << "V" << (irad+1) << " (non-lin): " << scientific
                 << l2genInput.instrument_gain_corr[0][irad][indx] << endl;
            cout << "V" << (irad+1) << " (offset):  " << scientific
                 << l2genInput.instrument_gain_corr[1][irad][indx] << endl;
          }
          if ( ipol == 3)  {
            cout << "H" << (irad+1) << " (non-lin): " << scientific
                 << l2genInput.instrument_gain_corr[0][irad][indx] << endl;
            cout << "H" << (irad+1) << " (offset) : " << scientific
                 << l2genInput.instrument_gain_corr[1][irad][indx] << endl;
          }
        }

        int off = RADIOMETER_LONG_ACCUM * 
          (RADIOMETER_POLARIZATIONS * 
           (NUMBER_OF_BEAMS*icyc + irad) + ipol);

        // Antenna minus reference
        float a_minus_r = L_acc_raw[off+4] - L_acc_raw[off+0];

        float rnd_minus_r;
        if ( ipol == 0) // V
          rnd_minus_r = 
            0.5*(L_acc_raw[off+1]+L_acc_raw[off+2]) - L_acc_raw[off+0];
        else if ( ipol == 3) // H
          rnd_minus_r = 
            0.5*(L_acc_raw[off+3]+L_acc_raw[off+2]) - L_acc_raw[off+0];

        //        TA_hat_drift_corr[2*(icyc*NUMBER_OF_BEAMS+irad)+indx] = 
        //((l2genInput.instrument_gain_corr[0][irad][indx] * a_minus_r) +
        // l2genInput.instrument_gain_corr[1][irad][indx]) / rnd_minus_r;

        // Run 5  JMG  11/12/14
        TA_hat_drift_corr[2*(icyc*NUMBER_OF_BEAMS+irad)+indx] = 
          ((l2genInput.instrument_gain_corr[0][irad][indx] * a_minus_r) +
           l2genInput.instrument_gain_corr[1][irad][indx]);
      }
    }
  }
#endif

  // Compute TA0/TF0 (without exponential drift correction)
  if ( l2genInput.c_deltaTND[0][0][0] != 0.0) {
    for (size_t irad=0; irad<NUMBER_OF_BEAMS; irad++) {
      for (size_t ipol=0; ipol<2; ipol++) { // 0=V,1=H
        c_deltaTND[irad][ipol] = 0.0;
      }
    }
    count_to_ta_( l2genInput.coeff_loss_file, l2genInput.coeff_nl_file, 
                  &nBlks, calTemps, &instrTemps, &gainOff, &rfi, 
                  S_acc_raw, L_acc_raw, S_acc, L_acc, 
                  cellat, cellon, zang, &c_deltaTND[0][0], TA_hat_drift_corr,
                  TA_hat, TA0, TF_hat, TF0);


    // Compute c_deltaTND from exponential fit
    if ( l2genInput.c_deltaTND[0][0][0] != 0.0) {
      // c_deltaTND(orbitNumber)
      cout << "Computing c_deltaTND from exponential fit" << endl;
      for (size_t irad=0; irad<NUMBER_OF_BEAMS; irad++) {
        for (size_t ipol=0; ipol<2; ipol++) { // 0=V,1=H
          c_deltaTND[irad][ipol] = 
            l2genInput.c_deltaTND[irad][ipol][0] -
            l2genInput.c_deltaTND[irad][ipol][1] *
            exp( -l2genInput.c_deltaTND[irad][ipol][2] * orbitNumber);
        }
      }
    }
  }

  // TA_hat is antenna temperatures before loss correction [kelvin] 
  // TA     is antenna temperatures after  loss correction [kelvin]
  // (1=V, 2=H, 3=3rd Stokes)

  // Note: TA_hat_drift_corr is multiplied by TND within count_to_ta  
  count_to_ta_( l2genInput.coeff_loss_file, l2genInput.coeff_nl_file, 
		&nBlks, calTemps, &instrTemps, &gainOff, &rfi, 
		S_acc_raw, L_acc_raw, S_acc, L_acc, 
                cellat, cellon, zang, &c_deltaTND[0][0], TA_hat_drift_corr,
		TA_hat, TA, TF_hat, TF);


  // Open calibration offset correction file
  filename = l2genInput.rad_offset_corr_file;
  float rad_offset_corr[6]={0.0,0.0,0.0,0.0,0.0,0.0};
  if ( filename != "") {
    expandEnvVar( &filename);
    h5fid =  H5Fopen( filename.c_str(),  H5F_ACC_RDONLY, H5P_DEFAULT);
    if ( h5fid == -1) {
      cout << "Unable to open offset correction file: " << 
	filename.c_str() << endl;
      exit(1);
    }

    grp = H5Gopen1(h5fid,"/");

    hsize_t start[2]={0,0}, count[2]={1,6};
    start[0] = orbitNumber-1;
    Hdf::h5d_read(grp, (char*) "rad_offset_corr", rad_offset_corr, 
		  2, start, count);
    H5Gclose(grp);
    H5Fclose(h5fid);

    cout << "Applying radiometer offset correction" << endl;

    // Subtract calibration offset correction
    for (iblk=0; iblk<nBlks; iblk++) {
      for (size_t irad=0; irad<NUMBER_OF_BEAMS; irad++) {
	int ibeamblk = irad+NUMBER_OF_BEAMS*iblk;
	for (size_t ipol=0; ipol<2; ipol++) { // 0=V,1=H
	  if ( rad_offset_corr[irad*2+ipol] != -9999.0) {
	    TA[3*ibeamblk+ipol] -= rad_offset_corr[irad*2+ipol];
	    TF[3*ibeamblk+ipol] -= rad_offset_corr[irad*2+ipol];
	  }
	}
      }
    }
  } // if (filename != "")


  //  filename = l2genInput.l_acc_corr_files;
  if ( l2genInput.iopt_l1b) {
    l2file.writel2_l_acc( nBlks, "L_acc_nl_corr", 
                          "Long Accumulations NL_Corrected",
                          L_acc);

    l2file.writel2_sec_TA( nBlks, &l1afile, TA);
  }

  // Note: TA_hat_drift_corr has been multiplied by TND
  //       within count_to_ta  
  if ( l2genInput.iopt_l1b == 1) {
    l2file.writel2_drift_corr( nBlks, TA_hat_drift_corr);
  }

  delete[] S_acc_raw;
  delete[] L_acc_raw;
  delete[] S_acc;
  delete[] L_acc;
  delete[] TA_hat_drift_corr;

  // Open scatterometer output txt files
  bool scatfilesexist = true;
  ifstream in_ANT, in_TOA, in_EXP, in_KPC_ANT, in_KPC_TOA;
  ifstream in_sWind, in_sLIfrac, in_sCenLL, in_sEdgeLL, in_sPolRoll;
  ifstream in_tot_TOA, in_tot_KPC_TOA;
  ifstream in_sWindUnc, in_exSurT, in_exSurT_unc;
  ifstream in_sFlags, in_sRFIflags, in_sSamples;

  cout << "Reading scatterometer output files" << endl;


  // Antenna
  scatter_textfile = scatter_basename + "_L2_sigma_ANT.txt";
  in_ANT.open ( scatter_textfile.c_str(), ifstream::in);
  if ( in_ANT.fail() == true) {
    cout << "Cannot open: "  << scatter_textfile.c_str() << endl;
    scatfilesexist = false;
  }

  // TOA
  scatter_textfile = scatter_basename + "_L2_sigma_TOA.txt";
  in_TOA.open ( scatter_textfile.c_str(), ifstream::in);
  if ( in_TOA.fail() == true) {
    cout << "Cannot open: "  << scatter_textfile.c_str() << endl;
  }

  // EXP
  scatter_textfile = scatter_basename + "_L2_sigma_exp.txt";
  in_EXP.open ( scatter_textfile.c_str(), ifstream::in);
  if ( in_EXP.fail() == true) {
    cout << "Cannot open: "  << scatter_textfile.c_str() << endl;
  }

  // tot_TOA
  scatter_textfile = scatter_basename + "_L2_tot_sigma_TOA.txt";
  in_tot_TOA.open ( scatter_textfile.c_str(), ifstream::in);
  if ( in_tot_TOA.fail() == true) {
    cout << "Cannot open: "  << scatter_textfile.c_str() << endl;
  }

  // ANT (KPC)
  scatter_textfile = scatter_basename + "_L2_KPC_ANT.txt";
  in_KPC_ANT.open ( scatter_textfile.c_str(), ifstream::in);
  if ( in_KPC_ANT.fail() == true) {
    cout << "Cannot open: "  << scatter_textfile.c_str() << endl;
  }

  // TOA (KPC)
  scatter_textfile = scatter_basename + "_L2_KPC_TOA.txt";
  in_KPC_TOA.open ( scatter_textfile.c_str(), ifstream::in);
  if ( in_KPC_TOA.fail() == true) {
    cout << "Cannot open: "  << scatter_textfile.c_str() << endl;
  }

  // tot_TOA (KPC)
  scatter_textfile = scatter_basename + "_L2_tot_KPC_TOA.txt";
  in_tot_KPC_TOA.open ( scatter_textfile.c_str(), ifstream::in);
  if ( in_tot_KPC_TOA.fail() == true) {
    cout << "Cannot open: "  << scatter_textfile.c_str() << endl;
  }

  // Scat Wind
  scatter_textfile = scatter_basename + "_L2_wind_speed.txt";
  in_sWind.open ( scatter_textfile.c_str(), ifstream::in);
  if ( in_sWind.fail() == true) {
    cout << "Cannot open: "  << scatter_textfile.c_str() << endl;
  }

  // Scat Land/Ice fraction
  scatter_textfile = scatter_basename + "_L2_land_ice_fractions.txt";
  in_sLIfrac.open ( scatter_textfile.c_str(), ifstream::in);
  if ( in_sLIfrac.fail() == true) {
    cout << "Cannot open: "  << scatter_textfile.c_str() << endl;
  }

  // Scat Center Lon/Lat
  scatter_textfile = scatter_basename + "_L2_center_lon_lat.txt";
  in_sCenLL.open ( scatter_textfile.c_str(), ifstream::in);
  if ( in_sCenLL.fail() == true) {
    cout << "Cannot open: "  << scatter_textfile.c_str() << endl;
  }

  // Scat Edge Lon/Lat
  scatter_textfile = scatter_basename + "_L2_edges_lon_lat.txt";
  in_sEdgeLL.open ( scatter_textfile.c_str(), ifstream::in);
  if ( in_sEdgeLL.fail() == true) {
    cout << "Cannot open: "  << scatter_textfile.c_str() << endl;
  }

  // Scat Polarization Roll angle
  scatter_textfile = scatter_basename + "_L2_polarization_roll.txt";
  in_sPolRoll.open ( scatter_textfile.c_str(), ifstream::in);
  if ( in_sPolRoll.fail() == true) {
    cout << "Cannot open: "  << scatter_textfile.c_str() << endl;
  }

  // Scat Wind Uncertainty
  scatter_textfile = scatter_basename + "_L2_wind_speed_err.txt";
  in_sWindUnc.open ( scatter_textfile.c_str(), ifstream::in);
  if ( in_sWindUnc.fail() == true) {
    cout << "Cannot open: "  << scatter_textfile.c_str() << endl;
  }

  // Excess Surface Emissivity
  scatter_textfile = scatter_basename + "_L2_excess_surface_emissivity.txt";
  in_exSurT.open ( scatter_textfile.c_str(), ifstream::in);
  if ( in_exSurT.fail() == true) {
    cout << "Cannot open: "  << scatter_textfile.c_str() << endl;
  }

  // Excess Surface Emissivity Uncertainty
  scatter_textfile = scatter_basename + "_L2_excess_surface_emissivity_err.txt";
  in_exSurT_unc.open ( scatter_textfile.c_str(), ifstream::in);
  if ( in_exSurT_unc.fail() == true) {
    cout << "Cannot open: "  << scatter_textfile.c_str() << endl;
  }

  // Scatterometer Flags
  scatter_textfile = scatter_basename + "_L2_flags.txt";
  in_sFlags.open ( scatter_textfile.c_str(), ifstream::in);
  if ( in_sFlags.fail() == true) {
    cout << "Cannot open: "  << scatter_textfile.c_str() << endl;
  }

  // Scatterometer RFI Flags
  scatter_textfile = scatter_basename + "_L2_scat_RFI_flags.txt";
  in_sRFIflags.open ( scatter_textfile.c_str(), ifstream::in);
  if ( in_sRFIflags.fail() == true) {
    cout << "Cannot open: "  << scatter_textfile.c_str() << endl;
  }

  // Scatterometer Samples
  scatter_textfile = scatter_basename + "_L2_samples.txt";
  in_sSamples.open ( scatter_textfile.c_str(), ifstream::in);
  if ( in_sSamples.fail() == true) {
    cout << "Cannot open: "  << scatter_textfile.c_str() << endl;
  }

 skipscat:

  ///////////////
  // Main Loop //
  ///////////////
  int32_t icnt = 0;
  double granStop;
  string nominalNav="TRUE";
  uint32_t tot_rad_samples[2]={0,0};
  uint32_t minSA=1;
  if ( l2genInput.iopt_nosa1 == true) {
    minSA = 2;
  }

  // If SSS algorithm "SIGMA0_HHH" then use SST/SSS in emiss function 
  staticData->emiss_sst_sss = 0;
  if ( strcmp( l2genInput.sss_algorithm, "SIGMA0_HHH") == 0) {
    staticData->emiss_sst_sss = 1;
  }

  // Compute GPS pointing anomoly event times
  double pointingTime[32];
  int32_t nPointing = 0;
  if (strcmp(l2genInput.pointing_anomaly_file, "") != 0) {
    cout << "Reading pointing anomaly file: " << l2genInput.pointing_anomaly_file << endl;
    ifstream inPointingAnomoly;
    char txtbuf[512];
    char zulu[16];
    inPointingAnomoly.open( l2genInput.pointing_anomaly_file, ifstream::in);
    while(1) {
      inPointingAnomoly.getline( txtbuf, 512);
      if( inPointingAnomoly.eof()) break;
      memset(zulu, 0, 16);
      strncpy( zulu, &txtbuf[9], 7);
      strncat( zulu, &txtbuf[17], 2);
      strncat( zulu, &txtbuf[20], 2);
      strncat( zulu, &txtbuf[23], 2);
      strncat( zulu, &txtbuf[26], 2);
      pointingTime[nPointing++] = unix2gpstai(zulu2unix( zulu));

      memset(zulu, 0, 16);
      strncpy( zulu, &txtbuf[39], 7);
      strncat( zulu, &txtbuf[47], 2);
      strncat( zulu, &txtbuf[50], 2);
      strncat( zulu, &txtbuf[53], 2);
      strncat( zulu, &txtbuf[57], 2);
      pointingTime[nPointing++] = unix2gpstai(zulu2unix( zulu));
    }
  }
  int32_t currentPointing = 0;

  for (iblk=0; iblk<nBlks; iblk++) {

    // If not matchup then monitor iblk index
    if ( l2genInput.matchup_lat == -999 && l2genInput.matchup_lon == -999) {
      if ((iblk % 500) == 0) cout << "iblk: " << iblk << endl;
    }

    // Set SSS, etc to missing value
    for (ibeam=0; ibeam<NUMBER_OF_BEAMS; ibeam++) {
      sss_aquarius_nolc[NUMBER_OF_BEAMS*iblk+ibeam] = -9999;
      sss_aquarius[NUMBER_OF_BEAMS*iblk+ibeam] = -9999;
      sss_aquarius_err_ran[NUMBER_OF_BEAMS*iblk+ibeam] = -9999;
      sss_aquarius_err_sys[NUMBER_OF_BEAMS*iblk+ibeam] = -9999;
      sss_aquarius_adj[NUMBER_OF_BEAMS*iblk+ibeam] = -9999;

      density[NUMBER_OF_BEAMS*iblk+ibeam] = -9999;
      spiciness[NUMBER_OF_BEAMS*iblk+ibeam] = -9999;
      
      tb_sur_nolc[2*(NUMBER_OF_BEAMS*iblk+ibeam)+0] = -9999;
      tb_sur_nolc[2*(NUMBER_OF_BEAMS*iblk+ibeam)+1] = -9999;
      tb_sur_rc_nolc[2*(NUMBER_OF_BEAMS*iblk+ibeam)+0] = -9999;
      tb_sur_rc_nolc[2*(NUMBER_OF_BEAMS*iblk+ibeam)+1] = -9999;
      tb_sur_rc[2*(NUMBER_OF_BEAMS*iblk+ibeam)+0] = -9999;
      tb_sur_rc[2*(NUMBER_OF_BEAMS*iblk+ibeam)+1] = -9999;
      tb_sur[2*(NUMBER_OF_BEAMS*iblk+ibeam)+0] = -9999;
      tb_sur[2*(NUMBER_OF_BEAMS*iblk+ibeam)+1] = -9999;
      ta_expected[3*(NUMBER_OF_BEAMS*iblk+ibeam)+0] = -9999;
      ta_expected[3*(NUMBER_OF_BEAMS*iblk+ibeam)+1] = -9999;
      ta_expected[3*(NUMBER_OF_BEAMS*iblk+ibeam)+2] = -9999;

      winspd_hh[NUMBER_OF_BEAMS*iblk+ibeam] = -9999; 
      winspd_hhh[NUMBER_OF_BEAMS*iblk+ibeam] = -9999; 
    }

    if ( inout[iblk] == false) {
      continue;
    }

    if ( l2genInput.iopt_l1b) continue;

    // blkSec - GPS seconds from 01/06/1980
    l1afile.readl1_radiometer(iblk, 0, 0, &blkSec, NULL,
			      NULL, NULL, NULL);

    /*
      inout[] (set above) should take care of bad entries!?
    // Skip processing for blocks with negative time
    */
    if ( blkSec == -1 || zang[iblk] == -999) {
      cout << "Empty frame at block " << iblk << endl;
      //      cout << iblk << " " << blkSec << " " << zang[iblk] << endl;
      //      icnt++; Removed JMG 05/01/09
      continue;
    }

    l1afile.readl1_frame_info( iblk,
			       &atc_frmnum, &atc_subframe,
			       &rad_frmnum, &rad_subframe);
    if ( rad_frmnum == -1) continue;

    // If pointing anomaly file exist then go to next bracketing
    // period if necessary
    if (strcmp(l2genInput.pointing_anomaly_file, "") != 0 &&
	currentPointing < nPointing-2) {
      while ( blkSec > pointingTime[currentPointing+1]) currentPointing += 2;
    }

    // Writes blkSec as both secs from start of day and blkSec
    // Note: secs from start is mid-point time, blkSec is start-block time
    l2file.writel2sec(icnt, gid[1], (VOIDP) &blkSec);
    l2file.writeSolarXrayFlux(icnt, gid[1], (VOIDP) &solar_xray[iblk]);

    double secdy, secyr;
    int lyear, idayjl, isecdy, idum;
    fd_date_2000_(&sec2000[iblk], &secyr, &lyear, &idayjl, 
		  &idum, &idum, &secdy);
    isecdy = (int) secdy; // MUST BE A TRUNCATION

    float orbit_angle = zang[iblk];
    if (orbit_angle < 0)   orbit_angle += 360;
    if (orbit_angle > 360) orbit_angle -= 360;


    //////////////////////////////////////////////////////////////////////////
    //////////////////////// Read Scatterometer Data /////////////////////////
    //////////////////////////////////////////////////////////////////////////

    // HV,VV,VH,HH 
    char scat_txtbuf[20][512];
    string scat_data[20];
    int32_t icds_blknum, scat_icds_blknum;
    //    istringstream istr;
    l1afile.readl1_icds_blknum( iblk, &icds_blknum);

    enum { inANT, inTOA, inEXP, inTOTTOA, inSWIND, inSLIFRAC, inSCENLL, 
	   inSEDGELL, inSPOLROLL, inKPCANT, inKPCTOA, inTOTKPCTOA, inSWINDUNC, 
	   inEXSURT, inEXSURTUNC, inSFLAGS, inSRFIFLAGS, inSAMPLES, nSCATFILES};

    while (1) {
      if ( !scatfilesexist) break;

      in_ANT.getline( scat_txtbuf[inANT], 512);
      in_TOA.getline( scat_txtbuf[inTOA], 512);
      in_EXP.getline( scat_txtbuf[inEXP], 512);
      in_tot_TOA.getline( scat_txtbuf[inTOTTOA], 512);
      in_sWind.getline( scat_txtbuf[inSWIND], 512);
      in_sLIfrac.getline( scat_txtbuf[inSLIFRAC], 512);
      in_sCenLL.getline( scat_txtbuf[inSCENLL], 512);
      in_sEdgeLL.getline( scat_txtbuf[inSEDGELL], 512);
      in_sPolRoll.getline( scat_txtbuf[inSPOLROLL], 512);

      in_KPC_ANT.getline( scat_txtbuf[inKPCANT], 512);
      in_KPC_TOA.getline( scat_txtbuf[inKPCTOA], 512);
      in_tot_KPC_TOA.getline( scat_txtbuf[inTOTKPCTOA], 512);

      in_sWindUnc.getline( scat_txtbuf[inSWINDUNC], 512);
      in_exSurT.getline( scat_txtbuf[inEXSURT], 512);
      in_exSurT_unc.getline( scat_txtbuf[inEXSURTUNC], 512);

      in_sFlags.getline( scat_txtbuf[inSFLAGS], 512);
      in_sRFIflags.getline( scat_txtbuf[inSRFIFLAGS], 512);

      in_sSamples.getline( scat_txtbuf[inSAMPLES], 512);

      for (size_t i=0; i<nSCATFILES; i++) {    
	scat_data[i] = scat_txtbuf[i];
      }

      istr.clear(); istr.str( scat_data[0].substr( 0, 12));
      istr >> scat_icds_blknum;  

      if ( scat_icds_blknum == icds_blknum) {
	for (size_t i=0; i<NUMBER_OF_BEAMS; i++) {    

	  // Read Scatterometer ANT & TOA & EXP & KPC_ANT & KPC_TOA
	  for (size_t j=0; j<4; j++) {    
	    istr.clear(); istr.str( scat_data[inANT].substr( 32+i*4*14+j*14, 
							     14));
	    istr >> scat_ant[i*4+j];

	    istr.clear(); istr.str( scat_data[inTOA].substr( 32+i*4*14+j*14, 
							     14));
	    istr >> scat_toa[i*4+j];

	    istr.clear(); istr.str( scat_data[inEXP].substr( 32+i*4*14+j*14, 
							     14));
	    istr >> scat_exp[i*4+j];

	    istr.clear(); istr.str( scat_data[inKPCANT].substr( 32+i*4*14+j*14,
								14));
	    istr >> scat_kpc_ant[i*4+j];

	    istr.clear(); istr.str( scat_data[inKPCTOA].substr( 32+i*4*14+j*14,
								14));
	    istr >> scat_kpc_toa[i*4+j];
	  }

	  // Read Scatterometer tot_TOA
	  istr.clear(); istr.str( scat_data[inTOTTOA].substr( 32+i*14, 14));
	  istr >> scat_tot_toa[i];

	  // Read Scatterometer tot_KPC_TOA
	  istr.clear(); istr.str( scat_data[inTOTKPCTOA].substr( 32+i*14, 14));
	  istr >> scat_tot_kpc_toa[i];

	  // Read Scatterometer Wind
	  istr.clear(); istr.str( scat_data[inSWIND].substr( 32+i*14, 14));
	  istr >> scat_swind[i];

	  // Read Scatterometer Wind Uncertainty
	  istr.clear(); istr.str( scat_data[inSWINDUNC].substr( 32+i*14, 14));
	  istr >> scat_swindunc[i];

	  // Read Scatterometer excess emissivity and uncertainty
	  // Fix ordering bug  JMG  10/20/11
	  for (size_t j=0; j<2; j++) {
	    istr.clear();
	    istr.str( scat_data[inEXSURT].substr( 32+i*2*14+j*14, 14));
	    istr >> scat_esurf[i*2+j];

	    istr.clear(); 
	    istr.str( scat_data[inEXSURTUNC].substr( 32+i*2*14+j*14, 14));
	    istr >> scat_esurf_unc[i*2+j];
	  }

	  // Read Scatterometer Land/Ice Fraction
	  for (size_t j=0; j<2; j++) { 
	    istr.clear(); 
	    istr.str( scat_data[inSLIFRAC].substr( 32+j*3*14+i*14, 14));
	    istr >> scat_sLIfrac[i*2+j];
	  }

	  // Read Scatterometer Flags
	  istr.clear(); 
	  istr.str( scat_data[inSFLAGS].substr( 32+i*20, 20));
	  istr >> scat_sFlags[i];

	  // Read Scatterometer RFI Flags
	  uint32_t tmp;
	  for (size_t j=0; j<2; j++) {
	    istr.clear(); 
	    istr.str( scat_data[inSRFIFLAGS].substr( 32+j*3*10+i*10, 10));
	    //	    cout << scat_data[inSRFIFLAGS].substr( 32+j*3*10+i*10, 10).c_str() << endl;
	    // Need to read into 32-bit int tmp, then assign to byte array
	    istr >> tmp;
	    scat_sRFI_Flags[i*2+j] = tmp;
	  }

	  // Read Scatterometer Samples
	  istr.clear(); istr.str( scat_data[inSAMPLES].substr( 32+i*10, 10));
	  istr >> scat_samples[i];

	  // Read Scatterometer Center Lon/Lat
	  // Beam 1 (lon/lat): *  * 
	  // Beam 2 (lon/lat): *  * 
	  // Beam 3 (lon/lat): *  *
	  istr.clear(); istr.str( scat_data[inSCENLL].substr( 32+i*14*2, 14));
	  istr >> scat_sCenLon[iblk*NUMBER_OF_BEAMS+i];
	  if ( scat_sCenLon[iblk*NUMBER_OF_BEAMS+i] > 180)
	    scat_sCenLon[iblk*NUMBER_OF_BEAMS+i] -= 360.0;

	  istr.clear(); 
	  istr.str( scat_data[inSCENLL].substr( 32+i*14*2+14, 14));
	  istr >> scat_sCenLat[iblk*NUMBER_OF_BEAMS+i];

	  // Read Scatterometer Edge Lon/Lat
	  // Beam 1 (lon/lat): * * * *   * * * *
	  // Beam 2 (lon/lat): * * * *   * * * *
	  // Beam 3 (lon/lat): * * * *   * * * *
	  for (size_t j=0; j<4; j++) {    
	    istr.clear(); 
	    istr.str( scat_data[inSEDGELL].substr( 32+i*4*14*2+j*14, 14));
	    istr >> scat_sEdgeLon[iblk*NUMBER_OF_BEAMS*4+i*4+j];
	    if ( scat_sEdgeLon[iblk*NUMBER_OF_BEAMS*4+i*4+j] > 180)
	      scat_sEdgeLon[iblk*NUMBER_OF_BEAMS*4+i*4+j] -= 360.0;

	    istr.clear(); 
	    istr.str( scat_data[inSEDGELL].substr( 32+i*4*14*2+4*14+j*14, 14));
	    istr >> scat_sEdgeLat[iblk*NUMBER_OF_BEAMS*4+i*4+j];
	  }

	  // Read Scatterometer Polarization Roll angle
	  istr.clear(); istr.str( scat_data[inSPOLROLL].substr( 32+i*14, 14));
	  istr >> scat_polarization_roll[iblk*NUMBER_OF_BEAMS+i];
	  if ( scat_polarization_roll[iblk*NUMBER_OF_BEAMS+i] > 180)
	    scat_polarization_roll[iblk*NUMBER_OF_BEAMS+i] -= 360.0;
	}  
	break;
      } // if scat_icds_blknum >= icds_blknum
    } // if scatfilesexist
    //////////////////////////////////////////////////////////////////////////
    //////////////////////// End Scatterometer Data //////////////////////////
    //////////////////////////////////////////////////////////////////////////


    if ( ! l2genInput.iopt_l1b) { 
      l2file.writel2samples( icnt, gid[1], (char *) "scat_samples", 
			     (VOIDP) scat_samples);
    }

    // tagal_dir = direct    galaxy radiation in terms of classical stokes (K)
    // tagal_ref = reflected galaxy radiation in terms of classical stokes (K)
    // tasun_dir = direct     solar radiation in terms of classical stokes (K)
    // tasun_ref = reflected  solar radiation in terms of classical stokes (K)
    fd_ta_sun_( &sec2000[iblk], &zang[iblk], &staticData->time_sun[0], 
		&staticData->tasun_dir_tab[0][0][0][0], 
		&staticData->tasun_ref_tab[0][0][0][0], 
		&tasun_dir[iblk*NUMBER_OF_BEAMS*3],
		&tasun_ref[iblk*NUMBER_OF_BEAMS*3]);
    for (size_t i=0; i<NUMBER_OF_BEAMS*3; i++) {
      tasun_dir[iblk*NUMBER_OF_BEAMS*3+i] *= solar_flux[iblk];
      tasun_ref[iblk*NUMBER_OF_BEAMS*3+i] *= solar_flux[iblk];
      //     fprintf(fp, "%15.2f, %12.8f %12.8f\n", sec2000[iblk], 
      //      tasun_dir[iblk*NUMBER_OF_BEAMS*3+i],
      //      tasun_ref[iblk*NUMBER_OF_BEAMS*3+i]);
    }

    // Update 1415 MHz flux peak and mean
    if ( solar_flux[iblk] > solar_attr[0]) solar_attr[0] = solar_flux[iblk];
    if ( solar_flux[iblk] > 0) {
      solar_attr[1] += solar_flux[iblk];
      n_solar++;
    }

    // Update solar xray peak and mean
    if ( solar_xray[iblk] > solar_attr[2]) solar_attr[2] = solar_xray[iblk];
    if ( solar_xray[iblk] > 0) {
      solar_attr[3] += solar_xray[iblk];
      n_xray++;
    }

    float ta_beam[3];
    uint32_t rad_flag[NUMBER_OF_BEAMS*4];

    uint8_t rd_sts_oplut;
    l1afile.readl1_dpu_status_tlm( rad_frmnum, rad_subframe, 10, 
				   &rd_sts_oplut);

    // Beam (Horn) Loop
    for (ibeam=0; ibeam<NUMBER_OF_BEAMS; ibeam++) {

      tb_toi[3*ibeam+0] = -9999;
      tb_toi[3*ibeam+1] = -9999;
      tb_toi[3*ibeam+2] = -9999;
      tb_toa[2*ibeam+0] = -9999;
      tb_toa[2*ibeam+1] = -9999;
      tb_toa_nolc[2*ibeam+0] = -9999;
      tb_toa_nolc[2*ibeam+1] = -9999;

      int32_t irad = ibeam + 1;

      int ibeamblk = ibeam+NUMBER_OF_BEAMS*iblk;

      float clon = cellon[ibeamblk];
      float clat = cellat[ibeamblk];

      float clon360 = clon;
      if ( clon != -999)
	if (clon360 < 0) clon360 += 360;

      rad_flag[4*ibeam+0] = 0;
      rad_flag[4*ibeam+1] = 0;
      rad_flag[4*ibeam+2] = 0;
      rad_flag[4*ibeam+3] = 0;

      // POINTING flag
      if ((blkSec > pointingTime[currentPointing]) && 
	  (blkSec < pointingTime[currentPointing+1])) {
	rad_flag[4*ibeam+0] |= (uint32_t) pow( 2, POINTING);
	nominalNav = "FALSE";
      }

      // POINTING flag (bad ACS mode)
      if ( (acsmode[iblk] != 5) && (acsmode[iblk] != 6)) {
	rad_flag[4*ibeam+1] |= (uint32_t) pow( 2, POINTING);
	nominalNav = "FALSE";
      }

      // NAV flag
      // Bit 1 (attitude)
      if ( fabs(rpy[3*iblk]) > limitRPY[0]) {
	rad_flag[4*ibeam+0] |= (uint32_t) pow( 2, NAV);
	nominalNav = "FALSE";
      }
      if ( fabs(rpy[3*iblk+1]) > limitRPY[1]) {
	rad_flag[4*ibeam+1] |= (uint32_t) pow( 2, NAV);
	nominalNav = "FALSE";
      }
      if ( fabs(rpy[3*iblk+2]) > limitRPY[2]) {
	rad_flag[4*ibeam+2] |= (uint32_t) pow( 2, NAV);
	nominalNav = "FALSE";
      }

      // NAV flag (off-earth beam)
      if ( clon == -999 || clat == -999) {
	rad_flag[4*ibeam+3] |= (uint32_t) pow( 2, NAV);
	nominalNav = "FALSE";
	continue;
      }

      // LAND flag
      if ( gland[ibeamblk] > limitLand[0] && gland[ibeamblk] <= limitLand[1])
	rad_flag[4*ibeam] |= (uint32_t) pow( 2, LAND);
      else if (gland[ibeamblk] > limitLand[1])
	rad_flag[4*ibeam+1] |= (uint32_t) pow( 2, LAND);

      if (gland[ibeamblk] > 0.5)
	rad_flag[4*ibeam+2] |= (uint32_t) pow( 2, LAND);


      // ICE flag
      if ( gice[ibeamblk] > limitIce[0] && gice[ibeamblk] <= limitIce[1])
	rad_flag[4*ibeam] |= (uint32_t) pow( 2, ICE);
      else if (gice[ibeamblk] > limitIce[1])
	rad_flag[4*ibeam+1] |= (uint32_t) pow( 2, ICE);
       else if (gice[ibeamblk] > 0.5)
	rad_flag[4*ibeam+2] |= (uint32_t) pow( 2, ICE);

      // FLARE flag
      if ( solar_xray[iblk] == -9999) {
	rad_flag[4*ibeam+2] |= (uint32_t) pow( 2, FLARE);
      } else {
	if ( solar_xray[iblk] > limitFlare[0] && 
	     solar_xray[iblk] <= limitFlare[1])
	  rad_flag[4*ibeam] |= (uint32_t) pow( 2, FLARE);
	else if (solar_xray[iblk] > limitFlare[1]) 
	  rad_flag[4*ibeam+1] |= (uint32_t) pow( 2, FLARE);
      }

      // Short accumulation counts overflow
      int32_t overflow_chan[3]={8,9,10};
      uint8_t overflow;
      l1afile.readl1_radiom_nrt_tlm( rad_frmnum, rad_subframe, 
				     overflow_chan[ibeam],  &overflow);
      if ( overflow == 1) {
	rad_flag[4*ibeam] |= (uint32_t) pow( 2, SAOVERFLOW);
	for (size_t i=0; i<2; i++) {
	  TA[3*ibeamblk+i] = -9999; 
	  TF[3*ibeamblk+i] = -9999; 
	}
	continue;
      }

      // Trap non-nominal command conditions
      if ( rd_sts_oplut != 0) {
	rad_flag[4*ibeam] |= (uint32_t) pow( 2, NONNOMCMD);
	continue;
      }


      // {HV,VV,VH,HH}
      float xsigma0_vv = scat_toa[ibeam*4+1];
      float xsigma0_hh = scat_toa[ibeam*4+3];

      // Adjust for sigma0 adjustment 11/2012
      if ( strcmp( l2genInput.sss_algorithm, "SIGMA0_HHH") != 0) {
	if ( ibeam >= 1) xsigma0_vv -= 0.2;
	if ( ibeam == 2) xsigma0_hh -= 0.3;
      }

      if ( (scat_sFlags[ibeam] & (1 << 24)) != 0) xsigma0_hh=0.0;
      if ( (scat_sFlags[ibeam] & (1 << 22)) != 0) xsigma0_vv=0.0;
      //if ( (scat_sFlags[ibeam] & (1 << 23)) != 0) xsigma0_vh=0.0;
      //if ( (scat_sFlags[ibeam] & (1 << 22)) != 0) xsigma0_vv=0.0;

      // convert from dB to real units
      float xsigma0_vv_ru=powf(10, 0.1*xsigma0_vv); 
      float xsigma0_hh_ru=powf(10, 0.1*xsigma0_hh); 

      // relative azimuth angle, 0=upwind
      float phir = celphi[ibeamblk] - windir[ibeamblk];

      float wspd_scat, chisq_scat;
      int32_t iflag_wspd;

      // It is used in DEW, whcih is poportional to E0
      // In the derivation of the DEW model function I have averaged over 
      // a large data set so effectively SSS=35.0 seems appropriate
      float ysss = 35.0; 
      float ysst =  surtep[ibeamblk] - 273.15;

      ///////////////////////////////////////////////////////////////////////
      ///////////////////////// HH winds ////////////////////////////////////
      ///////////////////////////////////////////////////////////////////////

      // If scat flags 29/31 set then skip
      if ( strcmp( l2genInput.sss_algorithm, "SIGMA0_HHH") == 0 &&
	   (scat_sFlags[ibeam] & (1 << 29)) == 0 &&
	   (scat_sFlags[ibeam] & (1 << 31)) == 0) {

	if ( (gice[ibeamblk] < 0.1) && (gland[ibeamblk] < 0.1)) {
	  FD_WSPD( &irad, &phir, 
		   &winspd_ncep[ibeamblk], 
		   NULL,
		   &xsigma0_vv_ru, 
		   &xsigma0_hh_ru,
		   &ysst,
		   &ysss,
		   staticData,
		   &winspd_hh[ibeamblk], 
		   &chisq_scat, 
		   &iflag_wspd);

	  // wind non-convergence
	  if ( iflag_wspd == 1) {
	    rad_flag[4*ibeam+2] |= (uint32_t) pow( 2, WIND);
	  }
	}
      } else {
	// WIND flag (scat RFI)
	rad_flag[4*ibeam+3] |= (uint32_t) pow( 2, WIND);
      }

      // WIND flag
      float winspd_chk;
      if ( winspd_hh[ibeamblk] != -9999)
	winspd_chk = winspd_hh[ibeamblk];
      else
	winspd_chk = winspd_ncep[ibeamblk];
      if ( (winspd_chk > limitWind[0]) && (winspd_chk <= limitWind[1]))
	rad_flag[4*ibeam] |= (uint32_t) pow( 2, WIND);
      else if (winspd_chk > limitWind[1])
	rad_flag[4*ibeam+1] |= (uint32_t) pow( 2, WIND);

      // ta_beam in V/H polarization at this point, changed to Stokes in ta_earth()
      if ( l2genInput.iopt_rfi == 0) {
	for (size_t i=0; i<3; i++) {
	  ta_beam[i] = TA[3*ibeamblk+i];
	}
      } else {
	for (size_t i=0; i<3; i++) {
	  ta_beam[i] = TF[3*ibeamblk+i];
	}
      }


      // TF-TA flag
      for (size_t i=0; i<2; i++) {
	float tfta = TF[3*ibeamblk+i]-TA[3*ibeamblk+i];

	// red
	if (tfta < limitTFTA[2] || tfta > limitTFTA[3]) 
	  rad_flag[4*ibeam+2*i+1] |= (uint32_t) pow( 2, TFTADIFF);

	// yellow
	if ( tfta >= limitTFTA[0] && tfta < limitTFTA[1])
	  rad_flag[4*ibeam+2*i] |= (uint32_t) pow( 2, TFTADIFF);

	if ( tfta > limitTFTA[3])
	  rad_flag[4*ibeam+2*i] |= (uint32_t) pow( 2, TFTADIFF);
      }


      // RFI mask flag
      uint8_t iasc, mask_value;
      if ( zang[iblk] < 180) iasc = 0; else iasc = 1;
      mask_value = 
	staticData->rfi_mask[(int) (clat+90)/2][(int) (clon360/2)][iasc];
      if ( mask_value == 1)
	rad_flag[4*ibeam] |= (uint32_t) pow( 2, RFI_REGION);


      uint8_t iopt[4]={1,1,1,1};

      float winspd = winspd_hh[ibeamblk];
      if ( winspd == -9999) winspd = winspd_ncep[ibeamblk];

      get_taearth( idayjl, ibeamblk, staticData,
		   cellat, clon360, celtht, phir, 
		   fland, gland, gice, surtep, tran, sss_reference, &winspd, 
                   sm, solar_flux, 
		   sec2000, zang, sun_zenith, bore_sight, moonr,
		   ta_beam, ta_earth,
		   tagal_dir, tagal_ref, tasun_dir,
		   tasun_ref, tasun_bak, tamon_ref, tagal_ref_GO,
		   dtagal_ref);

      //     ===============================================================
      //     ===== apply apc table to convert from ta_earth to toi tb =====
      //     ===============================================================

      // tb_toi - top of ionosphere brightness temperature in terms of 
      // vpol, hpol, 3rd stokes (K)
      for (size_t istokes=0; istokes<3; istokes++) {
	tb_toi[3*ibeam+istokes] = 0;
	for (size_t j=0; j<3; j++) {
	  tb_toi[3*ibeam+istokes] += 
	    staticData->apc_matrix[ibeam][istokes][j] * ta_earth[3*ibeam+j];
	}
      }


      // TM 04/08/2015: added array coeffs_dI_U
      if ( strcmp(l2genInput.dI_U_coeff_file, "") != 0) {
        float uval = ta_earth[3*ibeam+2];
        float dI;
        find_di_u_( &irad, &uval, &staticData->coeffs_dI_U[0][0], &dI);
        dI = 2.0*dI; // DI/2 to DI 
        tb_toi[3*ibeam] -= dI;
      }

      if ( abs(tb_toi[3*ibeam+2]) > 1.0e-8 || 
           abs(tb_toi[3*ibeam+1]) > 1.0e-8) {
        faraday_deg[ibeamblk] = 
          0.5*atan2f(tb_toi[3*ibeam+2], tb_toi[3*ibeam+1]) / DEG2RAD;
      } else {
        faraday_deg[ibeamblk] = -9999.0;
      }

      // tb_toa(_nolc) - top of atmosphere brightness temperature in terms of 
      // vpol, hpol
      tb_toa_nolc[2*ibeam+0] = tb_toi[3*ibeam+0];
      // faraday rotation correction
      tb_toa_nolc[2*ibeam+1] = sqrt(tb_toi[3*ibeam+1] * tb_toi[3*ibeam+1] + 
			       tb_toi[3*ibeam+2] * tb_toi[3*ibeam+2]);

      float tb_toa_tmp[3];
      tb_toa_tmp[0] = tb_toa_nolc[2*ibeam];
      tb_toa_tmp[1] = tb_toa_nolc[2*ibeam+1];
      
      stokes_to_vh_( &ta_earth[3*ibeam]);
      stokes_to_vh_( &tb_toi[3*ibeam]);
      stokes_to_vh_( tb_toa_tmp);

      tb_toa_nolc[2*ibeam+0] = tb_toa_tmp[0];
      tb_toa_nolc[2*ibeam+1] = tb_toa_tmp[1];

      //  3/5/2013: land corr only applied for gland>0.0005 (ATBD section 3.8)
      if ( gland[ibeamblk] > 0.0005) {
	tb_toa[2*ibeam+0] = tb_toa_nolc[2*ibeam+0] - tb_landcorr_vpol[ibeamblk];
	tb_toa[2*ibeam+1] = tb_toa_nolc[2*ibeam+1] - tb_landcorr_hpol[ibeamblk];
      } else {
	tb_toa[2*ibeam+0] = tb_toa_nolc[2*ibeam+0];
	tb_toa[2*ibeam+1] = tb_toa_nolc[2*ibeam+1];
      }
      
      // apply a very small adjustment to incidence angle to get 
      // effective incidence angle
      float thtadj = tht_coef[ibeam] * celtht[ibeamblk];

      float tbdown = tbdw[ibeamblk] + tran[ibeamblk] * TB_COS;

      float surtb_nolc[2];
      float surtb[2];
      surtb_nolc[0] = 
	(tb_toa_nolc[2*ibeam+0] - tbup[ibeamblk]) / tran[ibeamblk];
      surtb_nolc[1] = 
	(tb_toa_nolc[2*ibeam+1] - tbup[ibeamblk]) / tran[ibeamblk];
      surtb[0] = (tb_toa[2*ibeam+0] - tbup[ibeamblk]) / tran[ibeamblk];
      surtb[1] = (tb_toa[2*ibeam+1] - tbup[ibeamblk]) / tran[ibeamblk];

      float emiss_nolc[2];
      float emiss[2];
      emiss_nolc[0] = (surtb_nolc[0] - tbdown) / (surtep[ibeamblk] - tbdown);
      emiss_nolc[1] = (surtb_nolc[1] - tbdown) / (surtep[ibeamblk] - tbdown);
      emiss[0] = (surtb[0] - tbdown) / (surtep[ibeamblk] - tbdown);
      emiss[1] = (surtb[1] - tbdown) / (surtep[ibeamblk] - tbdown);

      float surf_emission_nolc[2], surf_emission0_nolc[2];
      float surf_emission[2], surf_emission0[2];
      surf_emission_nolc[0] = emiss_nolc[0] * surtep[ibeamblk];
      surf_emission_nolc[1] = emiss_nolc[1] * surtep[ibeamblk];
      surf_emission[0] = emiss[0] * surtep[ibeamblk];
      surf_emission[1] = emiss[1] * surtep[ibeamblk];

      // tb_sur - surface brightness temperature, 
      // ie emissivity times temperature, in terms of vpol, hpol
      // Stored as rad_TbV, rad_TbH in L2 file
      tb_sur_nolc[2*ibeamblk+0] = surf_emission_nolc[0];
      tb_sur_nolc[2*ibeamblk+1] = surf_emission_nolc[1];

      // No roughness/directional correction
      tb_sur[2*ibeamblk+0] = surf_emission[0];
      tb_sur[2*ibeamblk+1] = surf_emission[1];

      float em0[2], em_meas[2], dew_meas[2];
      float sss_clim;

      ///////////////////////////////////////////////////////////////////////
      //////////////////////// HHH winds ////////////////////////////////////
      ///////////////////////////////////////////////////////////////////////

      // If scat flags 29/31 set then skip
      if ( strcmp( l2genInput.sss_algorithm, "SIGMA0_HHH") == 0 &&
	   (scat_sFlags[ibeam] & (1 << 29)) == 0 &&
	   (scat_sFlags[ibeam] & (1 << 31)) == 0) {

	if ( (gice[ibeamblk] < 0.1) && (gland[ibeamblk] < 0.1)) {

	  float xsurtep = surtep[ibeamblk];
	  float sst = xsurtep-273.15;
	  int32_t isecyr = (int32_t) secyr;
	  fd_sss_clm_( &isecyr, &clat, &clon, staticData, &sss_clim);

	  MEISSNER_WENTZ_SALINITY( &freq_aq, &thtadj, 
				   &sst, &sss_clim, 
				   em0);

	  em_meas[0] = surf_emission[0] / xsurtep;
	  em_meas[1] = surf_emission[1] / xsurtep;
	  dew_meas[0] = (em_meas[0] - em0[0])*TEFF;
	  dew_meas[1] = (em_meas[1] - em0[1])*TEFF;

	  FD_WSPD( &irad, &phir, 
		   &winspd_ncep[ibeamblk], 
		   dew_meas,
		   &xsigma0_vv_ru, 
		   &xsigma0_hh_ru,
		   &ysst,
		   &ysss,
		   staticData,
		   &winspd_hhh[ibeamblk], 
		   &chisq_scat, 
		   &iflag_wspd);

	  // wind non-convergence
	  if ( iflag_wspd == 1) {
	    rad_flag[4*ibeam+2] |= (uint32_t) pow( 2, WIND);
	  }
	}
      }

      ///////////////////////////////////////////////////////////////////////
      //////////////////////// Roughness correction /////////////////////////
      ///////////////////////////////////////////////////////////////////////
      float dtbc[2];
      float dtbc_hh[2];
      float dtb_sst_wspd[2];

      int32_t itype;
      int32_t iwdir=1;
      int32_t ioob=0;

      // If scat flags 29/31 set or no HHH wind use NCEP
      if ( (scat_sFlags[ibeam] & (1 << 29)) != 0 ||
	   (scat_sFlags[ibeam] & (1 << 31)) != 0 ||
	   ( winspd_hhh[ibeamblk] == -9999)) {

	itype = 0;

	// Roughness correction with NCEP
	fd_dtb_roughness_( &itype, &iwdir, &irad, &winspd_ncep[ibeamblk], 
			   &phir, &xsigma0_vv_ru, &swh[ibeamblk],
			   &ysst, &ysss,
			   staticData, &ioob, &dtbc[0], &dtb_sst_wspd[0]);

        // Set dtbc_hh to be used in ta_expected
        dtbc_hh[0] = dtbc[0];
        dtbc_hh[1] = dtbc[1];

      } else {
	// Roughness correction with sigma 0
	itype = 1;
        // JMG 02/10/15 change '<' to '>'
        if ( abs(swh[ibeamblk] - (-999.0)) > 0.001) itype = 2;

	int32_t ioob1, ioob2;
	// Compute R.C. using HH wind for ta_exp calc
	fd_dtb_roughness_( &itype, &iwdir, &irad, &winspd_hh[ibeamblk], 
			   &phir, &xsigma0_vv_ru, &swh[ibeamblk],
			   &ysst, &ysss,
			   staticData, &ioob1, &dtbc_hh[0], NULL);
	
	// Compute R.C. using HHH wind
	fd_dtb_roughness_( &itype, &iwdir, &irad, &winspd_hhh[ibeamblk], 
			   &phir, &xsigma0_vv_ru, &swh[ibeamblk],
			   &ysst, &ysss,
			   staticData, &ioob2, &dtbc[0], &dtb_sst_wspd[0]);

	// Both HH and HHH wind R.C. corrections must be OK
	ioob = ioob1 | ioob2;
	
	if ( ioob == 1) {
	  rad_flag[4*ibeam] |= (uint32_t) pow( 2, ROUGH);

	  itype = 0;
	  fd_dtb_roughness_( &itype, &iwdir, &irad, &winspd_ncep[ibeamblk], 
			     &phir, &xsigma0_vv_ru, &swh[ibeamblk],
			     &ysst, &ysss,
			     staticData, &ioob, &dtbc[0], &dtb_sst_wspd[0]);
	}
      }

      // Set SWH flag if not value in ancillary map
      if ( abs(swh[ibeamblk] - (-999.0)) < 0.001) {
	rad_flag[4*ibeam+1] |= (uint32_t) pow( 2, ROUGH);
      }

      // the dtbc is actually a demiss in V7 and 
      // needs to be converted to a TBsurf
      // V2.0  08/29/12
      dtbc[0] = dtbc[0] * surtep[ibeamblk] / TEFF;
      dtbc[1] = dtbc[1] * surtep[ibeamblk] / TEFF;

      // Fixed 10/03/13  JMG
      dtbc_hh[0] = dtbc_hh[0] * surtep[ibeamblk] / TEFF;
      dtbc_hh[1] = dtbc_hh[1] * surtep[ibeamblk] / TEFF;

      surf_emission0_nolc[0] = surf_emission_nolc[0] - dtbc[0]; 
      surf_emission0_nolc[1] = surf_emission_nolc[1] - dtbc[1]; 
      
      tb_sur_rc_nolc[2*ibeamblk+0] = surf_emission0_nolc[0];
      tb_sur_rc_nolc[2*ibeamblk+1] = surf_emission0_nolc[1];

      // Land corrected
      surf_emission0[0] = surf_emission[0] - dtbc[0]; 
      surf_emission0[1] = surf_emission[1] - dtbc[1]; 
      
      tb_sur_rc[2*ibeamblk+0] = surf_emission0[0];
      tb_sur_rc[2*ibeamblk+1] = surf_emission0[1];

      // estimate sea-surface salinity
      float surtep_celsius = surtep[ibeamblk] - 273.15;
      if ( surtep_celsius < -5.000) surtep_celsius = -5.000;
      if ( surtep_celsius > 39.999) surtep_celsius = 39.999;

      // COLDWATER flag
      if (surtep_celsius < limitColdWater[0] && 
	  surtep_celsius >= limitColdWater[1])
	rad_flag[4*ibeam+0] |= (uint32_t) pow( 2, COLDWATER);
      else if (surtep_celsius < limitColdWater[1])
	rad_flag[4*ibeam+1] |= (uint32_t) pow( 2, COLDWATER);


      // If land,ice < 0.5 then compute SSS
      if ( (gland[ibeamblk] < 0.5) && (gice[ibeamblk] < 0.5)) {// &&
	   //	   (scat_sFlags[ibeam] & (1 << 29)) == 0 &&
	   //(scat_sFlags[ibeam] & (1 << 31)) == 0) {

	float chisq_sss;
	int32_t iflag_sss;

        //        double sss_sst_adj = 
        // l2genInput.sss_sst_adj_parms[0]*surtep[ibeamblk]*surtep[ibeamblk] +
        //l2genInput.sss_sst_adj_parms[1]*surtep[ibeamblk] +
        //l2genInput.sss_sst_adj_parms[2];

	// No land-corrected SSS 
	if ( strcmp( l2genInput.sss_algorithm, "SIGMA0_HHH") == 0) {
	  FD_SSS( &irad, &thtadj, &surtep_celsius, 
		  &surf_emission0_nolc, 
		  &sss_aquarius_nolc[ibeamblk], 
		  &chisq_sss, &iflag_sss);

	  MEISSNER_WENTZ_SALINITY( &freq_aq, &thtadj, 
				   &surtep_celsius, 
				   &sss_aquarius_nolc[ibeamblk], 
				   em0);
	} else {
	  fd_sss_sigma0_( &thtadj, &surtep_celsius, &surf_emission0_nolc, em0,
			  &irad, &sss_aquarius_nolc[ibeamblk]);
	}

	tb_con_nolc[ibeamblk] = 
	  sqrt( pow(surf_emission0_nolc[0]-surtep[ibeamblk]*em0[0], 2) + 
		pow(surf_emission0_nolc[1]-surtep[ibeamblk]*em0[1], 2));

	// Land-corrected SSS
	if ( strcmp( l2genInput.sss_algorithm, "SIGMA0_HHH") == 0) {
	  FD_SSS( &irad, &thtadj, &surtep_celsius, 
		  &surf_emission0, 
		  &sss_aquarius[ibeamblk], 
		  &chisq_sss, &iflag_sss);

	  MEISSNER_WENTZ_SALINITY( &freq_aq, &thtadj, 
				   &surtep_celsius, 
				   &sss_aquarius[ibeamblk], 
				   em0);
	} else {
	  fd_sss_sigma0_( &thtadj, &surtep_celsius, &surf_emission0, em0,
			  &irad, &sss_aquarius[ibeamblk]);
	}

        // Apply SSS_static_bias_adj
        if ( abs(sss_aquarius[ibeamblk] - (-9999.0)) > 0.001 &&
             bias_adj_file != "") {
          int ilat = (int) (cellat[ibeamblk] + 90);
          int ilon = (int) (cellon[ibeamblk] + 180);
          sss_aquarius_adj[ibeamblk] = sss_aquarius[ibeamblk] - 
            sss_static_bias_adj[iasc][ibeam][ilat][ilon];
        }

	// Compute tb_con and set TBCONS flag is necessary
	tb_con[ibeamblk] = 
	  sqrt( pow(surf_emission0[0]-surtep[ibeamblk]*em0[0], 2) + 
		pow(surf_emission0[1]-surtep[ibeamblk]*em0[1], 2));

	if (tb_con[ibeamblk] > limitTbcons) {
	  rad_flag[4*ibeam+0] |= (uint32_t) pow( 2, TBCONS);
	}
	if (isnan(em0[0])) {
	  rad_flag[4*ibeam+1] |= (uint32_t) pow( 2, TBCONS);
	}


        // Compute random and systematic SSS error
        // Thomas Meissner, Remote Sensing Systems, April 20, 2015.
        if ( strcmp(l2genInput.l2_unc_maps_file, "") != 0) {
          int32_t iasc;
          // iasc = 1 ascending, iasc = 2 descending
          if (0 <= zang[iblk] && zang[iblk] <= 180) iasc = 1; else iasc = 2;
          sec2000[iblk] = gpstai2utc2000( blkSec);
          fd_date_2000_(&sec2000[iblk], &dumd, &lyear, &idayjl, 
                        &idum, &idum, &secdy);
          interpolate_l2_error_(l2genInput.l2_unc_maps_file, &idayjl, 
                                &clat, &clon, &iasc,
                                &sss_aquarius_err_ran[ibeamblk], 
                                &sss_aquarius_err_sys[ibeamblk]);
        }


        // Compute Density and Spiciness
        if ( sss_aquarius[ibeamblk] != -9999) {
          // Density and Spiciness
          // Julian J Schanze <jschanze@esr.org>

          double aqS, aqT, aqSA, aqCT, aqlat, aqlon, aqrho, aqtau;

          aqlon = (double) clon;
          aqlat = (double) clat;
          aqS = (double) sss_aquarius[ibeamblk];
          aqT = (double) surtep_celsius;

          // In this first step, absolute salinity is computed from Aquarius 
          // salinity (PSS-78), pressure (always 0 dbar in the case of Aquarius)
          // and the longitude and latitude for each data point 
          // (which is needed to determine absolute salinity).
          aqSA = gsw_sa_from_sp( aqS, 0, aqlon, aqlat);

          // This computes conservative temperature (CT) from in-situ 
          // temperature. 
          // CT is a required variable in the spiciness computation.
          // aqrho=gsw_rho(aqSA,aqCT,0); // inputs: SA, CT, pressure
          aqCT = gsw_ct_from_t( aqSA, aqT, 0);

          // inputs: SA, in-situ temp, pressure
          // Note that you have a choice here of computing density (rho) from 
          // using gsw_rho_t_exact, which uses the full equation of state or 
          // using gsw_rho, which uses a computationally much more efficient 
          // (albeit less accurate) code. 
          // I am not sure how significant the speed impact would be on the 
          // total processing system would be. If you can use gsw_rho_ct_exact 
          // without too much penalty that would be the preferred option, but 
          // gsw_rho would be good enough for the T and S accuracies we have.
          aqrho = gsw_rho_t_exact( aqSA, aqT, 0); 

          // This computes spiciness from absolute salinity (SA) and
          // conservative temperature (CT) at the surface (0dB pressure). 
          aqtau = gsw_spiciness0( aqSA, aqCT);

          density[ibeamblk] = (float) aqrho;
          spiciness[ibeamblk] = (float) aqtau;
        }

      } // If land,ice < 0.5 then compute SSS

    ta_expected:
      // expected ta
      thtadj = tht_coef[ibeam] * celtht[ibeamblk];
      phir = celphi[ibeamblk] - windir[ibeamblk];
      tbdown = tbdw[ibeamblk]+ tran[ibeamblk] * TB_COS;

      // Repeat for L2_CAL generation
      if ( surtep_celsius < -5.000) surtep_celsius = -5.000;
      if ( surtep_celsius > 39.999) surtep_celsius = 39.999;
      
      // Set windspeed to 0 to get tbsur0
      float tbsur0[2], tbsur[2];
      float refl0[2], refl[2];
      find_refl_tot_( &idayjl, 
		      &cellat[ibeamblk], &clon360, &thtadj,
		      &phir, &fland[ibeamblk], &fice[ibeamblk], 
		      &sss_reference[ibeamblk],
		      &surtep_celsius, &zerof, &sm[ibeamblk],
		      staticData, refl0, tbsur0);

      winspd = winspd_hh[ibeamblk];
      if ( winspd == -9999) {
	winspd = winspd_ncep[ibeamblk];
	dtbc_hh[0] = dtbc[0];
	dtbc_hh[1] = dtbc[1];
      }

      if ( (rad_flag[4*ibeam] & (uint32_t) pow( 2, ROUGH)) != 0) {
	find_refl_tot_( &idayjl, 
			&cellat[ibeamblk], &clon360, &thtadj,
			&phir, &fland[ibeamblk], &fice[ibeamblk], 
			&sss_reference[ibeamblk],
			&surtep_celsius, &winspd, &sm[ibeamblk],
			staticData, refl, tbsur);
      } else {
	tbsur[0] = tbsur0[0] + dtbc_hh[0];
	tbsur[1] = tbsur0[1] + dtbc_hh[1];

	refl[0] = 1.0 - (tbsur[0]/surtep[ibeamblk]);
	refl[1] = 1.0 - (tbsur[1]/surtep[ibeamblk]);
      }


      // Compute tb_expected, tb_expected0
      tb_expected[2*ibeamblk+0] = tbsur[0];
      tb_expected[2*ibeamblk+1] = tbsur[1];
      tb_expected0[2*ibeamblk+0] = tbsur0[0];
      tb_expected0[2*ibeamblk+1] = tbsur0[1];

      float tbtoa_expected[3];
      // Use tbsur rather than tb_sur JMG 07/15/11
      for (size_t i=0; i<2; i++) {
	tbtoa_expected[i] = tbup[ibeamblk] + 
	  tran[ibeamblk] * (tbsur[i] + refl[i] * tbdown);
      }
      vh_to_stokes_( tbtoa_expected);
      
      float tbtoi_expected[3];
	
      tbtoi_expected[0] = tbtoa_expected[0];
      tbtoi_expected[1] = tbtoa_expected[1] * 
	cos(2*faraday_deg[ibeamblk]*DEG2RAD);
      tbtoi_expected[2] = tbtoa_expected[1] * 
	sin(2*faraday_deg[ibeamblk]*DEG2RAD);

      // TM 04/08/2015: added array coeffs_dI_U
      if ( strcmp(l2genInput.dI_U_coeff_file, "") != 0) {
        float uval = 
          staticData->apc_inverse[ibeam][2][0] * tbtoi_expected[0] +
	  staticData->apc_inverse[ibeam][2][1] * tbtoi_expected[1] +
	  staticData->apc_inverse[ibeam][2][2] * tbtoi_expected[2];
        float dI;
        find_di_u_( &irad, &uval, &staticData->coeffs_dI_U[0][0], &dI);
        dI = 2.0*dI; // DI/2 to DI 
        tbtoi_expected[0] += dI;
      }

      float taexpected[3];
      for (size_t i=0; i<3; i++) {
	taexpected[i] = 
	  staticData->apc_inverse[ibeam][i][0] * tbtoi_expected[0] +
	  staticData->apc_inverse[ibeam][i][1] * tbtoi_expected[1] +
	  staticData->apc_inverse[ibeam][i][2] * tbtoi_expected[2];

	taexpected[i] += tagal_dir[3*ibeamblk+i] + tagal_ref[3*ibeamblk+i] + 
	  iopt[0]*tasun_dir[3*ibeamblk+i] + iopt[1]*tasun_ref[3*ibeamblk+i] +
	  iopt[2]*tamon_ref[3*ibeamblk+i] + iopt[3]*tasun_bak[3*ibeamblk+i];
      }
      stokes_to_vh_( taexpected);

      for (size_t i=0; i<3; i++) {
	ta_expected[3*ibeamblk+i] = taexpected[i];
      }


      // TEMP flag
      // Change from tb_toa to ta_beam for comparison JMG  09/14/2011
      // Convert ta_beam back to VH
      stokes_to_vh_( ta_beam);

      // V Pol
      if ( fabs(ta_beam[0] - taexpected[0]) >  limitTemp[0] && 
	   fabs(ta_beam[0] - taexpected[0]) <= limitTemp[1])
	rad_flag[4*ibeam+0] |= (uint32_t) pow( 2, TEMP);
      else if (fabs(ta_beam[0] - taexpected[0]) > limitTemp[1])
	rad_flag[4*ibeam+1] |= (uint32_t) pow( 2, TEMP);

      // H Pol
      if ( fabs(ta_beam[1] - taexpected[1]) >  limitTemp[0] && 
	   fabs(ta_beam[1] - taexpected[1]) <= limitTemp[1])
	rad_flag[4*ibeam+2] |= (uint32_t) pow( 2, TEMP);
      else if (fabs(ta_beam[1] - taexpected[1]) > limitTemp[1])
	rad_flag[4*ibeam+3] |= (uint32_t) pow( 2, TEMP);


      // Repeat with HHH winds
      winspd = winspd_hhh[ibeamblk];
      if ( winspd == -9999) {
	winspd = winspd_ncep[ibeamblk];
      }

      if ( (rad_flag[4*ibeam] & (uint32_t) pow( 2, ROUGH)) != 0) {
	find_refl_tot_( &idayjl, 
			&cellat[ibeamblk], &clon360, &thtadj,
			&phir, &fland[ibeamblk], &fice[ibeamblk], 
			&sss_reference[ibeamblk],
			&surtep_celsius, &winspd, &sm[ibeamblk],
			staticData, refl, tbsur);
      } else {
	tbsur[0] = tbsur0[0] + dtbc[0];
	tbsur[1] = tbsur0[1] + dtbc[1];

	refl[0] = 1.0 - (tbsur[0]/surtep[ibeamblk]);
	refl[1] = 1.0 - (tbsur[1]/surtep[ibeamblk]);
      }

      // Use tbsur rather than tb_sur JMG 07/15/11
      for (size_t i=0; i<2; i++) {
	tbtoa_expected[i] = tbup[ibeamblk] + 
	  tran[ibeamblk] * (tbsur[i] + refl[i] * tbdown);
      }
      vh_to_stokes_( tbtoa_expected);
      
      tbtoi_expected[0] = tbtoa_expected[0];
      tbtoi_expected[1] = tbtoa_expected[1] * 
	cos(2*faraday_deg[ibeamblk]*DEG2RAD);
      tbtoi_expected[2] = tbtoa_expected[1] * 
	sin(2*faraday_deg[ibeamblk]*DEG2RAD);

      // TM 04/08/2015: added array coeffs_dI_U
      if ( strcmp(l2genInput.dI_U_coeff_file, "") != 0) {
        float uval = 
          staticData->apc_inverse[ibeam][2][0] * tbtoi_expected[0] +
	  staticData->apc_inverse[ibeam][2][1] * tbtoi_expected[1] +
	  staticData->apc_inverse[ibeam][2][2] * tbtoi_expected[2];
        float dI;
        find_di_u_( &irad, &uval, &staticData->coeffs_dI_U[0][0], &dI);
        dI = 2.0*dI; // DI/2 to DI 
        tbtoi_expected[0] += dI;
      }

      for (size_t i=0; i<3; i++) {
	taexpected[i] = 
	  staticData->apc_inverse[ibeam][i][0] * tbtoi_expected[0] +
	  staticData->apc_inverse[ibeam][i][1] * tbtoi_expected[1] +
	  staticData->apc_inverse[ibeam][i][2] * tbtoi_expected[2];

	taexpected[i] += tagal_dir[3*ibeamblk+i] + tagal_ref[3*ibeamblk+i] + 
	  iopt[0]*tasun_dir[3*ibeamblk+i] + iopt[1]*tasun_ref[3*ibeamblk+i] +
	  iopt[2]*tamon_ref[3*ibeamblk+i] + iopt[3]*tasun_bak[3*ibeamblk+i];
      }
      stokes_to_vh_( taexpected);


      for (size_t i=0; i<3; i++) {
	ta_expected_hhh[3*ibeamblk+i] = taexpected[i];
      }
      // End repeat with HHH winds

      // MOON flag
      float ta_moon_tmp[3];
      ta_moon_tmp[0] = tamon_ref[3*ibeamblk+0];
      ta_moon_tmp[1] = tamon_ref[3*ibeamblk+1];
      ta_moon_tmp[2] = tamon_ref[3*ibeamblk+2];

      // MOON 1st Stokes flag
      if (ta_moon_tmp[0] > limit1stStokes[0] && 
	  ta_moon_tmp[0] <= limit1stStokes[1])
	rad_flag[4*ibeam+0] |= (uint32_t) pow( 2, REFL_1STOKES);
      else if (ta_moon_tmp[0] > limit1stStokes[1])
	rad_flag[4*ibeam+1] |= (uint32_t) pow( 2, REFL_1STOKES);

      stokes_to_vh_( ta_moon_tmp);

      // V Pol
      if ( (ta_moon_tmp[0] > limitMoon[0]) && (ta_moon_tmp[0] <= limitMoon[1]))
	rad_flag[4*ibeam+0] |= (uint32_t) pow( 2, MOON);
      else if (ta_moon_tmp[0] > limitMoon[1])
	rad_flag[4*ibeam+1] |= (uint32_t) pow( 2, MOON);

      // H Pol
      if ( (ta_moon_tmp[1] > limitMoon[0]) && (ta_moon_tmp[1] <= limitMoon[1]))
	rad_flag[4*ibeam+2] |= (uint32_t) pow( 2, MOON);
      else if (ta_moon_tmp[1] > limitMoon[1])
	rad_flag[4*ibeam+3] |= (uint32_t) pow( 2, MOON);


      // SUNGLINT flag
      float ta_sun_tmp[3];
      ta_sun_tmp[0] = tasun_bak[3*ibeamblk+0];
      ta_sun_tmp[1] = tasun_bak[3*ibeamblk+1];
      ta_sun_tmp[2] = tasun_bak[3*ibeamblk+2];
      stokes_to_vh_( ta_sun_tmp);

      // V Pol
      if ( (ta_sun_tmp[0] > limitGlint[0]) && (ta_sun_tmp[0] <= limitGlint[1]))
	rad_flag[4*ibeam+0] |= (uint32_t) pow( 2, SUNGLINT);
      else if (ta_sun_tmp[0] > limitGlint[1])
	rad_flag[4*ibeam+1] |= (uint32_t) pow( 2, SUNGLINT);

      // H Pol
      if ( (ta_sun_tmp[1] > limitGlint[0]) && (ta_sun_tmp[1] <= limitGlint[1]))
	rad_flag[4*ibeam+2] |= (uint32_t) pow( 2, SUNGLINT);
      else if (ta_sun_tmp[1] > limitGlint[1])
	rad_flag[4*ibeam+3] |= (uint32_t) pow( 2, SUNGLINT);


      // FLUXD flag
      ta_sun_tmp[0] = tasun_dir[3*ibeamblk+0];
      ta_sun_tmp[1] = tasun_dir[3*ibeamblk+1];
      ta_sun_tmp[2] = tasun_dir[3*ibeamblk+2];
      stokes_to_vh_( ta_sun_tmp);

      // V Pol
      if ( (ta_sun_tmp[0] > limitFluxD[0]) && (ta_sun_tmp[0] <= limitFluxD[1]))
	rad_flag[4*ibeam+0] |= (uint32_t) pow( 2, FLUXD);
      else if (ta_sun_tmp[0] > limitFluxD[1])
	rad_flag[4*ibeam+1] |= (uint32_t) pow( 2, FLUXD);

      // H Pol
      if ( (ta_sun_tmp[1] > limitFluxD[0]) && (ta_sun_tmp[1] <= limitFluxD[1]))
	rad_flag[4*ibeam+2] |= (uint32_t) pow( 2, FLUXD);
      else if (ta_sun_tmp[1] > limitFluxD[1])
	rad_flag[4*ibeam+3] |= (uint32_t) pow( 2, FLUXD);


      // FLUXR flag
      ta_sun_tmp[0] = tasun_ref[3*ibeamblk+0];
      ta_sun_tmp[1] = tasun_ref[3*ibeamblk+1];
      ta_sun_tmp[2] = tasun_ref[3*ibeamblk+2];
      stokes_to_vh_( ta_sun_tmp);

      // V Pol
      if ( (ta_sun_tmp[0] > limitFluxR[0]) && (ta_sun_tmp[0] <= limitFluxD[1]))
	rad_flag[4*ibeam+0] |= (uint32_t) pow( 2, FLUXR);
      else if (ta_sun_tmp[0] > limitFluxR[1])
	rad_flag[4*ibeam+1] |= (uint32_t) pow( 2, FLUXR);

      // H Pol
      if ( (ta_sun_tmp[1] > limitFluxR[0]) && (ta_sun_tmp[1] <=limitFluxR[1]))
	rad_flag[4*ibeam+2] |= (uint32_t) pow( 2, FLUXR);
      else if (ta_sun_tmp[1] > limitFluxR[1])
	rad_flag[4*ibeam+3] |= (uint32_t) pow( 2, FLUXR);



      // GALACTIC flag
      float ta_gal_tmp[3];
      ta_gal_tmp[0] = tagal_ref[3*ibeamblk+0];
      ta_gal_tmp[1] = tagal_ref[3*ibeamblk+1];
      ta_gal_tmp[2] = tagal_ref[3*ibeamblk+2];

      // GALACTIC 1st Stokes flag
      if (ta_gal_tmp[0] > limit1stStokes[2] || 
	  (ta_gal_tmp[0] > limit1stStokes[3] && 
	   winspd_hh[ibeamblk] < limit1stStokes[4])) 
	rad_flag[4*ibeam+2] |= (uint32_t) pow( 2, REFL_1STOKES);

      stokes_to_vh_( ta_gal_tmp);

      // V Pol
      if ( (ta_gal_tmp[0] > limitGal[0]) && (ta_gal_tmp[0] <= limitGal[1]))
	rad_flag[4*ibeam+0] |= (uint32_t) pow( 2, GALACTIC);
      else if (ta_gal_tmp[0] > limitGal[1])
	rad_flag[4*ibeam+1] |= (uint32_t) pow( 2, GALACTIC);

      // H Pol
      if ( (ta_gal_tmp[1] > limitGal[0]) && (ta_gal_tmp[1] <= limitGal[1]))
	rad_flag[4*ibeam+2] |= (uint32_t) pow( 2, GALACTIC);
      else if (ta_gal_tmp[1] > limitGal[1])
	rad_flag[4*ibeam+3] |= (uint32_t) pow( 2, GALACTIC);


      // Convert tagal, tasun, tamon from Classical Stokes to VH3
      stokes_to_vh_( &tagal_dir[3*ibeamblk]);
      stokes_to_vh_( &tagal_ref[3*ibeamblk]);
      stokes_to_vh_( &tagal_ref_GO[3*ibeam]);
      stokes_to_vh_( &dtagal_ref[3*ibeam]);
      stokes_to_vh_( &tasun_dir[3*ibeamblk]);
      stokes_to_vh_( &tasun_ref[3*ibeamblk]);
      stokes_to_vh_( &tasun_bak[3*ibeamblk]);
      stokes_to_vh_( &tamon_ref[3*ibeamblk]);

    } // ibeam loop


    // Write radiometer RFI flags/rad_samples
    // Note: polarization order: V/P/M/H
    uint8_t rad_rfi_flag[NUMBER_OF_BEAMS][RADIOMETER_POLARIZATIONS]
      [RADIOMETER_SUBCYCLES];
    for (ibeam=0; ibeam<NUMBER_OF_BEAMS; ibeam++) {
      for (size_t ipol=0; ipol<RADIOMETER_POLARIZATIONS; ipol++) {

	// Initialize radiometer samples
	if ( l2genInput.iopt_nosa1 == true) {
	  rad_samples[ibeam][ipol] = 60;
	  tot_rad_samples[1] += 60;
	} else {
	  rad_samples[ibeam][ipol] = 84;
	  tot_rad_samples[1] += 84;
	}

	for (size_t isub=0; isub<RADIOMETER_SUBCYCLES; isub++) {
	  uint8_t byt = 0;
	  if ( rfi.iflag_rfi_CND[iblk][ibeam][ipol][isub] == 1) {
	    byt = 2;
	  }
	  // Change from 6 to 5 JMG 08/25/11
	  for (size_t i=minSA-1; i<5; i++) {
	    if ( rfi.iflag_rfi[iblk][ibeam][ipol][isub][i] == 1) {
	      byt += (2 << (i+1));

	      if (i <= 1) 
		rad_samples[ibeam][ipol] -= 2;
	      else
		rad_samples[ibeam][ipol] -= 1;
	    }
	  }
	  rad_rfi_flag[ibeam][ipol][isub] = byt;
	} // isub

	// Set RFI flag
	if ( rad_samples[ibeam][ipol] <  limitRFI[0] && 
	     rad_samples[ibeam][ipol] >= limitRFI[1])
	  rad_flag[4*ibeam+ipol] |= (uint32_t) pow( 2, RFI);
	if ( rad_samples[ibeam][ipol] < limitRFI[1])
	  rad_flag[4*ibeam+ipol] |= (uint32_t) pow( 2, RFI+1);

	tot_rad_samples[0] += rad_samples[ibeam][ipol];
      } // ipol
    } // ibeam

    // If not CAL then write samples and flags
    if ( ! l2genInput.iopt_l1b) { 
      l2file.writel2radrfi(icnt, gid[2], (VOIDP) rad_rfi_flag);
      l2file.writel2samples( icnt, gid[1], (char *) "rad_samples", 
			     (VOIDP) rad_samples);

      // Write scatterometer RFI flags
      l2file.writel2scatrfi(icnt, gid[2], (VOIDP) scat_sRFI_Flags);


      // Write radiometer flags
      l2file.writel2radflag( icnt, gid[2], (VOIDP) rad_flag);

      // Write scatterometer flags
      l2file.writel2(icnt, gid[2], (char *) "scatterometer_flags", 
		     (VOIDP) scat_sFlags);
    }


    // Write L2 products
    float datbuf[NUMBER_OF_BEAMS];
    float *dptr;
    for (size_t i=0; i<MAXNAQPROD; i++) {
      if ( l2genInput.l2activeprodflag[i] != 0) {

	// Product does not exist
	if ( prodPtr[i] == 0x0) {
	  cout << "Bad product pointer for " << 
	    shortProdName[i].c_str() << endl;
	  exit(1);
	}

	for (ibeam=0; ibeam<NUMBER_OF_BEAMS; ibeam++) {
	  if ( prodMultiplicity[i] > 0) {
	    dptr = prodPtr[i] + iblk*NUMBER_OF_BEAMS*prodMultiplicity[i] +
	      ibeam*prodMultiplicity[i] + prodOffset[i];
	  } else {
	    dptr = prodPtr[i] - ibeam*prodMultiplicity[i] + prodOffset[i];
	  }
	  memcpy( &datbuf[ibeam], dptr, sizeof(float));
	}
	l2file.writel2(icnt, gid[3], 
		       (char *) shortProdName[i].c_str(), 
		       (VOIDP) datbuf);
      } // if active prod
    } // product loop

    // Write calTemps and gainOff
    l2file.writel2_caltemps( icnt, &calTemps[iblk*NUMCALTEMPS]);
    l2file.writel2_gainoff( icnt, &gainOff, l2genInput.iopt_l1b);


    granStop = blkSec;

    icnt++;
  } // iblk (End of Main Loop)

  // Writing Navigation Data
  cout << "Writing Navigation Data" << endl;

  // Convert cellonfoot, sclon to -180 - +180
  // Fix bonehead mistake (upper limit for i set to argc!?!)  JMG  09/09/11
  for (iblk=0; iblk<nBlks; iblk++) {
    for (ibeam=0; ibeam<NUMBER_OF_BEAMS; ibeam++) {
      for (size_t i=0; i<4; i++) {
	if ( cellonfoot[4*(ibeam+NUMBER_OF_BEAMS*iblk)+i] > 180)
	  cellonfoot[4*(ibeam+NUMBER_OF_BEAMS*iblk)+i] -= 360;
      }
    }
    if ( sclon[iblk] > 180) sclon[iblk] -= 360;
  }

  if ( ! l2genInput.iopt_l1b) {
    l2file.writel2_nav( nBlks, inout, Pos, Vel, rpy, 
			cellat, cellon,
			celtht, celphi, suntht, sunphi, sunglt, moonglt,
			glxlat, glxlon, cellatfoot, cellonfoot, 
			sund, sunr, moond, zang, 
			sclon, sclat, scalt,
			scat_sCenLon, scat_sCenLat,
			scat_sEdgeLon, scat_sEdgeLat, scat_polarization_roll,
			acsmode, l2genInput.browse, false);
  } else {
    l2file.writel2_nav( nBlks, inout, Pos, Vel, rpy, 
			cellat, cellon,
			celtht, celphi, suntht, sunphi, sunglt, moonglt,
			glxlat, glxlon, cellatfoot, cellonfoot, 
			sund, sunr, moond, zang, 
			sclon, sclat, scalt,
			scat_sCenLon, scat_sCenLat,
			scat_sEdgeLon, scat_sEdgeLat, scat_polarization_roll,
			acsmode, true, l2genInput.iopt_l1b);
  }

  delete[] Pos;
  delete[] Vel;
  delete[] rpy;


  l1afile.closel1();

  float percentRFI = 100 *
    (1 - ((float) tot_rad_samples[0])/tot_rad_samples[1]);

  // Determine Ancillary Files Used
  string rad_anc_files[3];
  rad_anc_files[0].assign( l2genInput.yancfile1);
  rad_anc_files[1].assign( l2genInput.yancfile2);
  rad_anc_files[2].assign( l2genInput.yancfile3);

  // Generate Radiometer calibration files string
  string rad_calfiles;
  rad_calfiles.assign( basename(l2genInput.coeff_loss_file));
  rad_calfiles.append( ",");
  rad_calfiles.append( basename(l2genInput.coeff_nl_file)); 

  // Generate Radiometer data tables string
  string rad_tables;
  rad_tables.assign( basename(l2genInput.rad_landtables_file));
  rad_tables.append( ",");
  rad_tables.append( basename(l2genInput.rad_gainice_file)); 
  rad_tables.append( ",");
  rad_tables.append( basename(l2genInput.rad_tausq_file)); 
  rad_tables.append( ",");
  rad_tables.append( basename(l2genInput.rad_oceanrefl_file)); 
  rad_tables.append( ",");
  rad_tables.append( basename(l2genInput.rad_sun_file)); 
  rad_tables.append( ",");
  rad_tables.append( basename(l2genInput.rad_sunbak_file)); 
  rad_tables.append( ",");
  rad_tables.append( basename(l2genInput.rad_galwind_file)); 
  rad_tables.append( ",");
  rad_tables.append( basename(l2genInput.rad_apc_file));
  rad_tables.append( ",");
  rad_tables.append( basename(l2genInput.rad_landcorr_file));
  rad_tables.append( ",");
  rad_tables.append( basename(l2genInput.wind_errortab_file));
  rad_tables.append( ",");
  rad_tables.append( basename(l2genInput.emiss_coeff_harm_file));
  rad_tables.append( ",");
  rad_tables.append( basename(l2genInput.scat_coeff_harm_file));
  rad_tables.append( ",");
  rad_tables.append( basename(l2genInput.dtbw_win_sigma_file));
  rad_tables.append( ",");
  rad_tables.append( basename(l2genInput.dtbw_win_wav_file));
  rad_tables.append( ",");
  rad_tables.append( basename(l2genInput.climate_sal_file));
  rad_tables.append( ",");
  rad_tables.append( basename(l2genInput.rad_dta_gal_file));

  ifstream in_SCAT_ANCFILES;
  string sScatAncFiles="";
  scatter_textfile = scatter_basename + "_L2_ancillary_files.txt";
  in_SCAT_ANCFILES.open ( scatter_textfile.c_str(), ifstream::in);
  if ( in_SCAT_ANCFILES.fail() == false) {
    string sLine;
    while(1) {
      getline( in_SCAT_ANCFILES, sLine);
      if( in_SCAT_ANCFILES.eof()) break;
      sScatAncFiles.append( "," + sLine); 
    }
    sScatAncFiles = sScatAncFiles.substr( 1, string::npos);
    in_SCAT_ANCFILES.close();
  }

  solar_attr[1] /= n_solar;
  solar_attr[3] /= n_xray;
  l2file.closel2( granStop, percentWater, percentRFI, solar_attr,
		  rad_anc_files, &sScatAncFiles, &rad_calfiles, &rad_tables,
		  l2genInput.anomaly_status, &c_deltaTND[0][0],
		  radLimits.rdbuf()->str().c_str(), nominalNav.c_str(),
		  rad_offset_corr);

  // Close scatterometer output text files
  in_TOA.close();
  in_ANT.close();
  in_tot_TOA.close();
  in_sWind.close();
  in_sLIfrac.close();
  in_sCenLL.close();
  in_sEdgeLL.close();
  in_sPolRoll.close();
  in_KPC_ANT.close();
  in_KPC_TOA.close();
  in_tot_KPC_TOA.close();
  in_sWindUnc.close();
  in_exSurT.close();
  in_exSurT_unc.close();
  in_sFlags.close();
  in_sRFIflags.close();
  in_sSamples.close();

  delete[] cellon;
  delete[] cellat;
  delete[] celtht;
  delete[] celphi;
  delete[] suntht;
  delete[] sunphi;
  delete[] glxlon;
  delete[] glxlat;
  delete[] zang;
  delete[] cellonfoot;
  delete[] cellatfoot;
  delete[] sund;
  delete[] sunr;
  delete[] moond;

  delete[] solar_flux;
  delete[] solar_xray;
  delete[] fland;
  delete[] gland;
  delete[] fice;
  delete[] gice;
  delete[] swh;
  delete[] sss_reference;
  delete[] surtep;
  delete[] subtep;
  delete[] surp;
  delete[] cwat;
  delete[] sm;
  delete[] winspd_ncep;
  delete[] winspd_hh;
  delete[] winspd_hhh;
  delete[] windir;
  delete[] tran;
  delete[] tbup;
  delete[] tbdw;
  delete[] calTemps;
  delete staticData;

  cout << "Normal Completion" << endl << endl;

  return 0;
}


int getGeoNavSun( Hdf::hdf5_Aquarius *l1afile, double *rpy_adj,
		  float *cellon, float *cellat, float *celtht, float *celphi,
		  float *suntht, float *sunphi, float *sunglt, float *moonglt,
		  float *glxlon, float *glxlat,
		  double *zang, double *sun_zenith,
		  double *sclon, double *sclat, double *scalt, 
		  float *cellonfoot, float *cellatfoot,
		  double *sund, double *sunr, double *moond, double *moonr, 
		  double *bore_sight,
		  double *Pos, double *Vel, double *rpy)
{

  // Read Orbital Data
  uint32_t numOrbVec;
  l1afile->readl1_eph( &numOrbVec);
  double *orbTime, *orbPos, *orbVel, *dum;
  orbTime = new double[numOrbVec];
  orbPos = new double[numOrbVec*3];
  orbVel = new double[numOrbVec*3];
  dum = new double[numOrbVec];
  l1afile->readl1_eph( orbTime, orbPos, orbVel);

  double *Time;
  Time = new double[l1afile->nBlks()];

  double blkSec;

  // Cubic Interpolation Using orbPos and orbVel
  gsl_matrix *time_mat = gsl_matrix_alloc(4,4);
  gsl_vector *b = gsl_vector_alloc(4);

  gsl_permutation *p = gsl_permutation_alloc(4);
  int s;
  double C, D;

  for (size_t i=0; i<l1afile->nBlks(); i++) {
    l1afile->readl1_radiometer(i, 0, 0, &blkSec, NULL, NULL, NULL, NULL);

    // Add 0.72 sec to get mid-block time
    Time[i] = blkSec + 0.7200;
    if ( blkSec == -1) continue;
    for (size_t j=0; j<numOrbVec-1; j++) {
      if ( Time[i] < orbTime[0]) {
	Time[i] = -1;
	break;
      }
      if ( orbTime[j] <= Time[i] && orbTime[j+1] > Time[i]) {
	double time1_1 = orbTime[j] - Time[i];
	double time2_1 = orbTime[j+1] - Time[i];

	double time1_2 = time1_1 * time1_1;
	double time2_2 = time2_1 * time2_1;
	double time1_3 = time1_2 * time1_1;
	double time2_3 = time2_2 * time2_1;
	
	gsl_matrix_set( time_mat, 0, 0, time1_3);
	gsl_matrix_set( time_mat, 0, 1, time1_2);
	gsl_matrix_set( time_mat, 0, 2, time1_1);
	gsl_matrix_set( time_mat, 0, 3, 1.0);

	gsl_matrix_set( time_mat, 1, 0, time2_3);
	gsl_matrix_set( time_mat, 1, 1, time2_2);
	gsl_matrix_set( time_mat, 1, 2, time2_1);
	gsl_matrix_set( time_mat, 1, 3, 1.0);

	gsl_matrix_set( time_mat, 2, 0, 3*time1_2);
	gsl_matrix_set( time_mat, 2, 1, 2*time1_1);
	gsl_matrix_set( time_mat, 2, 2, 1.0);
	gsl_matrix_set( time_mat, 2, 3, 0.0);

	gsl_matrix_set( time_mat, 3, 0, 3*time2_2);
	gsl_matrix_set( time_mat, 3, 1, 2*time2_1);
	gsl_matrix_set( time_mat, 3, 2, 1.0);
	gsl_matrix_set( time_mat, 3, 3, 0.0);

	gsl_linalg_LU_decomp( time_mat, p, &s);


	// X, Vx
	gsl_vector_set( b, 0, orbPos[3*j]);
	gsl_vector_set( b, 1, orbPos[3*(j+1)]);
	gsl_vector_set( b, 2, orbVel[3*j]);
	gsl_vector_set( b, 3, orbVel[3*(j+1)]);
	gsl_linalg_LU_svx( time_mat, p, b);

	C = gsl_vector_get( b, 2);
	D = gsl_vector_get( b, 3);

	Pos[3*i] = D;
	Vel[3*i] = C;

	// Y, Vy
	gsl_vector_set( b, 0, orbPos[3*j+1]);
	gsl_vector_set( b, 1, orbPos[3*(j+1)+1]);
	gsl_vector_set( b, 2, orbVel[3*j+1]);
	gsl_vector_set( b, 3, orbVel[3*(j+1)+1]);
	gsl_linalg_LU_svx( time_mat, p, b);

	C = gsl_vector_get( b, 2);
	D = gsl_vector_get( b, 3);

	Pos[3*i+1] = D;
	Vel[3*i+1] = C;

	// Z, Vz
	gsl_vector_set( b, 0, orbPos[3*j+2]);
	gsl_vector_set( b, 1, orbPos[3*(j+1)+2]);
	gsl_vector_set( b, 2, orbVel[3*j+2]);
	gsl_vector_set( b, 3, orbVel[3*(j+1)+2]);
	gsl_linalg_LU_svx( time_mat, p, b);

	C = gsl_vector_get( b, 2);
	D = gsl_vector_get( b, 3);

	Pos[3*i+2] = D;
	Vel[3*i+2] = C;
      }
    }
  }


  // Read Attitude Data
  uint32_t numAttSamp;
  l1afile->readl1_att( &numAttSamp);
  double *attTime, *attAng;
  attTime = new double[numAttSamp];
  attAng = new double[numAttSamp*3];
  dum = new double[numAttSamp];
  l1afile->readl1_att( attTime, attAng, NULL);

  // Set up interpolation
  gsl_interp_accel *acc;
  gsl_interp *interp;
  acc = gsl_interp_accel_alloc ();
  interp = gsl_interp_alloc (gsl_interp_linear, numAttSamp);

  // Turn off error handler to avoid potential aborts at last attTime
  gsl_error_handler_t *old_gsl_handler = gsl_set_error_handler_off();

  // Roll
  for (size_t i=0; i<numAttSamp; i++) dum[i] = attAng[3*i];
  gsl_interp_init (interp, attTime, dum, numAttSamp);
  for (size_t i=0; i<l1afile->nBlks(); i++)
    rpy[3*i] = gsl_interp_eval (interp, attTime, dum, Time[i], acc);

  // Pitch
  for (size_t i=0; i<numAttSamp; i++) dum[i] = attAng[3*i+1];
  gsl_interp_init (interp, attTime, dum, numAttSamp);
  for (size_t i=0; i<l1afile->nBlks(); i++)
    rpy[3*i+1] = gsl_interp_eval (interp, attTime, dum, Time[i], acc);

  // Yaw
  for (size_t i=0; i<numAttSamp; i++) dum[i] = attAng[3*i+2];
  gsl_interp_init (interp, attTime, dum, numAttSamp);
  for (size_t i=0; i<l1afile->nBlks(); i++)
    rpy[3*i+2] = gsl_interp_eval (interp, attTime, dum, Time[i], acc);

  gsl_set_error_handler( old_gsl_handler);

  gsl_interp_accel_free (acc);
  gsl_interp_free (interp);

  cout << "Compute lon/lat" << endl;
  double time_utc_2000;
  for (size_t i=0; i<l1afile->nBlks(); i++) {
    // Check for -1 (need to subtract half-block time)
    if ( (Time[i] - 0.7200) == -1) continue;

    //    if ((i % 500) == 0) cout << i << endl;

    time_utc_2000 = gpstai2utc2000 (Time[i]);

    geolocation_( rpy_adj, &time_utc_2000, &Pos[3*i], &Vel[3*i], &rpy[3*i], 
		  &cellat[NUMBER_OF_BEAMS*i], &cellon[NUMBER_OF_BEAMS*i], 
		  &celtht[NUMBER_OF_BEAMS*i], &celphi[NUMBER_OF_BEAMS*i], 
		  &suntht[NUMBER_OF_BEAMS*i], &sunphi[NUMBER_OF_BEAMS*i], 
		  &sunglt[NUMBER_OF_BEAMS*i], &moonglt[NUMBER_OF_BEAMS*i], 
		  &glxlat[NUMBER_OF_BEAMS*i], &glxlon[NUMBER_OF_BEAMS*i],
		  &zang[i], &sun_zenith[i], &sclon[i], &sclat[i], &scalt[i], 
		  &cellatfoot[NUMBER_OF_BEAMS*4*i], 
		  &cellonfoot[NUMBER_OF_BEAMS*4*i],
		  &sund[NUMBER_OF_BEAMS*i], &sunr[NUMBER_OF_BEAMS*i], 
		  &moond[NUMBER_OF_BEAMS*i], &moonr[NUMBER_OF_BEAMS*i],
		  &bore_sight[NUMBER_OF_BEAMS*3*i]);

  }

  for (size_t i=0; i<NUMBER_OF_BEAMS*l1afile->nBlks(); i++) {
    if ( cellon[i] != -999) 
      cellon[i] = 
	cellon[i]*(cellon[i] <= 180.0) + (cellon[i]-360.)*(cellon[i] > 180.0);
  }

  delete[] orbTime;
  delete[] orbPos;
  delete[] orbVel;
  delete[] attTime;
  delete[] attAng;
  delete[] dum;
  delete[] Time;

  return 0;
}

int getHKT( Hdf::hdf5_Aquarius *l1afile, uint8_t *acsMode)
{
  uint32_t numHktSamp;
  l1afile->readl1_hkt( &numHktSamp);

  uint32_t *hktTime;
  hktTime = new uint32_t[numHktSamp];
  l1afile->readl1_hkttme( numHktSamp, hktTime);

  uint8_t *acsmode;
  acsmode = new uint8_t[numHktSamp];
  l1afile->readl1_hktacsmode( numHktSamp, acsmode);

  double *Time;
  Time = new double[l1afile->nBlks()];

  double blkSec;

  for (size_t i=0; i<l1afile->nBlks(); i++) {
    l1afile->readl1_radiometer(i, 0, 0, &blkSec, NULL, NULL, NULL, NULL);

    // Add 0.72 sec to get mid-block time
    Time[i] = blkSec + 0.7200;
    if ( blkSec == -1) continue;
    for (size_t j=0; j<numHktSamp-1; j++) {
      if ( Time[i] < hktTime[0]) {
	Time[i] = -1;
	break;
      }
      if ( hktTime[j] <= Time[i] && hktTime[j+1] > Time[i]) {
	acsMode[i] = acsmode[j];
      }
    }
  }

  delete[] hktTime;
  delete[] acsmode;
  delete[] Time;

  return 0;
}


int parseInput( int argc, char* argv[], instr *l2genInput, uint32_t *nProd,
		string *procControl)
{
  using namespace std;

  string::size_type posEndIdx;
  string::size_type ipos=0;
  string            sLine, sValue;
  string            sParmWord;
  istringstream istr;
  stack<string> parmStack, procCntlStack;
  bool iHelp=false;
  bool iXml=false;

  char xmlfile[256];

  procControl->clear();

  // Push regular command line parameters
  for (int i=1; i<argc; i++) {
    extParmWordValue( sLine.assign(argv[i]), &sParmWord, &sValue);

    if (sParmWord.compare("PAR") == 0) continue;
    if (sParmWord.compare("SUITE") == 0) continue;

    parmStack.push( sLine);
  }

  // Push PAR parameters
  for (int i=1; i<argc; i++) {
    extParmWordValue( sLine.assign(argv[i]), &sParmWord, &sValue);
    if (sParmWord.compare("PAR") == 0) parmStack.push( sLine);
  }

  // Push SUITE parameters
  for (int i=1; i<argc; i++) {
    extParmWordValue( sLine.assign(argv[i]), &sParmWord, &sValue);
    if (sParmWord.compare("SUITE") == 0) parmStack.push( sLine);
  }

  // Push default parameters
  sValue = "$OCDATAROOT/aquarius/l2gen_aquarius_defaults.par";
  expandEnvVar( &sValue);
  sLine = "par=" + sValue;
  parmStack.push( sLine);


  // Pop stack loop
  while (!parmStack.empty()) {

    sLine = parmStack.top();
    parmStack.pop();
    procCntlStack.push( sLine);

    extParmWordValue( sLine, &sParmWord, &sValue);

    if (sParmWord.compare("IFILE") == 0) {
      int len = sValue.copy(l2genInput->ifile, sValue.length());
      l2genInput->ifile[len] = 0;
      continue;
    }

    if (sParmWord.compare("OFILE") == 0) {
      int len = sValue.copy(l2genInput->ofile, sValue.length());
      l2genInput->ofile[len] = 0;
      continue;
    }

    if (sParmWord.compare("INCALFILELIST") == 0) {
      int len = sValue.copy(l2genInput->incalfilelist, sValue.length());
      l2genInput->incalfilelist[len] = 0;
      continue;
    }

    if (sParmWord.compare("L2PROD") == 0) {
      size_t found = sValue.find("default:");
      if ( found == 0) {
	strcat( l2genInput->l2prod, ",");
	strcat( l2genInput->l2prod, sValue.substr(8).c_str());
	continue;
      } else {
	int len = sValue.copy(l2genInput->l2prod, sValue.length());
	l2genInput->l2prod[len] = 0;
	continue;
      }
    }

    if (sParmWord.compare("MATCHUP_TIME") == 0) {
      istr.clear(); istr.str( sValue); istr >> l2genInput->matchup_time;
      continue;
    }

    if (sParmWord.compare("MATCHUP_LAT") == 0) {
      istr.clear(); istr.str( sValue); istr >> l2genInput->matchup_lat;
      continue;
    }

    if (sParmWord.compare("MATCHUP_LON") == 0) {
      istr.clear(); istr.str( sValue); istr >> l2genInput->matchup_lon;
      if ( l2genInput->matchup_lon > 180)
	l2genInput->matchup_lon -= 360;
      continue;
    }

    if (sParmWord.compare("MATCHUP_DELTA_LAT") == 0) {
      istr.clear(); istr.str( sValue); istr >> l2genInput->matchup_delta_lat;
      continue;
    }

    if (sParmWord.compare("MATCHUP_DELTA_LON") == 0) {
      istr.clear(); istr.str( sValue); istr >> l2genInput->matchup_delta_lon;
      continue;
    }

    if (sParmWord.compare("MATCHUP_MIN_DISTANCE") == 0) {
      istr.clear(); istr.str( sValue); istr >> l2genInput->matchup_min_dist;
      continue;
    }

    if (sParmWord.compare("MATCHUP_LIMITS_FILE") == 0) {
      int len = sValue.copy(l2genInput->matchup_lim_file, sValue.length());
      l2genInput->matchup_lim_file[len] = 0;
      continue;
    }

    if (sParmWord.compare("SCATTER_BASENAME") == 0) {
      int len = sValue.copy(l2genInput->scatter_basename, sValue.length());
      l2genInput->scatter_basename[len] = 0;
      continue;
    }

    if (sParmWord.compare("BROWSE") == 0) {
      istr.clear(); istr.str( sValue); istr >> l2genInput->browse;
      continue;
    }

    if (sParmWord.compare("IOPT_RFI") == 0) {
      // Note: "boolalpha" allows "true/false" to be used
      istr.clear(); istr.str( sValue);
      istr >> boolalpha >> l2genInput->iopt_rfi;
      continue;
    }

    if (sParmWord.compare("IOPT_L1B") == 0) {
      // Note: "boolalpha" allows "true/false" to be used
      istr.clear(); istr.str( sValue);
      istr >> boolalpha >> l2genInput->iopt_l1b;
      continue;
    }

    if (sParmWord.compare("IOPT_NOSA1") == 0) {
      // Note: "boolalpha" allows "true/false" to be used
      istr.clear(); istr.str( sValue);
      istr >> boolalpha >> l2genInput->iopt_nosa1;
      continue;
    }

    if (sParmWord.compare("IOPT_ZERO") == 0) {
      // Note: "boolalpha" allows "true/false" to be used
      istr.clear(); istr.str( sValue);
      istr >> boolalpha >> l2genInput->iopt_zero;
      continue;
    }

    if (sParmWord.compare("PVERSION") == 0) {
      int len = sValue.copy(l2genInput->pversion, sValue.length());
      l2genInput->pversion[len] = 0;
      continue;
    }

    if (sParmWord.compare("SSS_ALGORITHM") == 0) {
      int len = sValue.copy(l2genInput->sss_algorithm, sValue.length());
      l2genInput->sss_algorithm[len] = 0;
      continue;
    }

    if (sParmWord.compare("RAD_LANDTABLES_FILE") == 0) {
      int len = sValue.copy(l2genInput->rad_landtables_file, sValue.length());
      l2genInput->rad_landtables_file[len] = 0;
      continue;
    }

    if (sParmWord.compare("RAD_LANDCORR_FILE") == 0) {
      int len = sValue.copy(l2genInput->rad_landcorr_file, sValue.length());
      l2genInput->rad_landcorr_file[len] = 0;
      continue;
    }

    if (sParmWord.compare("RAD_GAINICE_FILE") == 0) {
      int len = sValue.copy(l2genInput->rad_gainice_file, sValue.length());
      l2genInput->rad_gainice_file[len] = 0;
      continue;
    }

    if (sParmWord.compare("RAD_TAUSQ_FILE") == 0) {
      int len = sValue.copy(l2genInput->rad_tausq_file, sValue.length());
      l2genInput->rad_tausq_file[len] = 0;
      continue;
    }

    if (sParmWord.compare("RAD_APC_FILE") == 0) {
      int len = sValue.copy(l2genInput->rad_apc_file, sValue.length());
      l2genInput->rad_apc_file[len] = 0;
      continue;
    }

    if (sParmWord.compare("RAD_SSSALGO_FILE") == 0) {
      int len = sValue.copy(l2genInput->rad_sssalgo_file, sValue.length());
      l2genInput->rad_sssalgo_file[len] = 0;
      continue;
    }

    if (sParmWord.compare("RAD_GALWIND_FILE") == 0) {
      int len = sValue.copy(l2genInput->rad_galwind_file, sValue.length());
      l2genInput->rad_galwind_file[len] = 0;
      continue;
    }

    if (sParmWord.compare("RAD_SUN_FILE") == 0) {
      int len = sValue.copy(l2genInput->rad_sun_file, sValue.length());
      l2genInput->rad_sun_file[len] = 0;
      continue;
    }

    if (sParmWord.compare("RAD_SUNBAK_FILE") == 0) {
      int len = sValue.copy(l2genInput->rad_sunbak_file, sValue.length());
      l2genInput->rad_sunbak_file[len] = 0;
      continue;
    }

    if (sParmWord.compare("RAD_OCEANREFL_FILE") == 0) {
      int len = sValue.copy(l2genInput->rad_oceanrefl_file, sValue.length());
      l2genInput->rad_oceanrefl_file[len] = 0;
      continue;
    }

    if (sParmWord.compare("WIND_ERRORTAB_FILE") == 0) {
      int len = sValue.copy(l2genInput->wind_errortab_file, sValue.length());
      l2genInput->wind_errortab_file[len] = 0;
      continue;
    }

    if (sParmWord.compare("EMISS_COEFF_HARM_FILE") == 0) {
      int len = sValue.copy(l2genInput->emiss_coeff_harm_file, 
			    sValue.length());
      l2genInput->emiss_coeff_harm_file[len] = 0;
      continue;
    }

    if (sParmWord.compare("SCAT_COEFF_HARM_FILE") == 0) {
      int len = sValue.copy(l2genInput->scat_coeff_harm_file, sValue.length());
      l2genInput->scat_coeff_harm_file[len] = 0;
      continue;
    }

    if (sParmWord.compare("DTBW_WIN_SIGMA_FILE") == 0) {
      int len = sValue.copy(l2genInput->dtbw_win_sigma_file, sValue.length());
      l2genInput->dtbw_win_sigma_file[len] = 0;
      continue;
    }

    if (sParmWord.compare("DTBW_WIN_WAV_FILE") == 0) {
      int len = sValue.copy(l2genInput->dtbw_win_wav_file, sValue.length());
      l2genInput->dtbw_win_wav_file[len] = 0;
      continue;
    }

    if (sParmWord.compare("DTB_EMISS_WSPD_FILE") == 0) {
      int len = sValue.copy(l2genInput->dtb_emiss_wspd_file, sValue.length());
      l2genInput->dtb_emiss_wspd_file[len] = 0;
      continue;
    }

    if (sParmWord.compare("RAD_DTA_GAL_FILE") == 0) {
      int len = sValue.copy(l2genInput->rad_dta_gal_file, sValue.length());
      l2genInput->rad_dta_gal_file[len] = 0;
      continue;
    }

    if (sParmWord.compare("COEFF_LOSS_FILE") == 0) {
      if ( sValue.substr( 0, 1).compare("$") == 0) expandEnvVar( &sValue);
      int len = sValue.copy(l2genInput->coeff_loss_file, sValue.length());
      l2genInput->coeff_loss_file[len] = 0;
      continue;
    }

    if (sParmWord.compare("COEFF_NL_FILE") == 0) {
      if ( sValue.substr( 0, 1).compare("$") == 0) expandEnvVar( &sValue);
      int len = sValue.copy(l2genInput->coeff_nl_file, sValue.length());
      l2genInput->coeff_nl_file[len] = 0;
      continue;
    }

    if (sParmWord.compare("YANCFILE1") == 0) {
      int len = sValue.copy(l2genInput->yancfile1, sValue.length());
      l2genInput->yancfile1[len] = 0;
      continue;
    }

    if (sParmWord.compare("YANCFILE2") == 0) {
      int len = sValue.copy(l2genInput->yancfile2, sValue.length());
      l2genInput->yancfile2[len] = 0;
      continue;
    }

    if (sParmWord.compare("YANCFILE3") == 0) {
      int len = sValue.copy(l2genInput->yancfile3, sValue.length());
      l2genInput->yancfile3[len] = 0;
      continue;
    }

    if (sParmWord.compare("RPY_ADJ") == 0) {
      istr.clear(); istr.str( sValue);
      istr >> l2genInput->rpy_adj[0];
      istr >> l2genInput->rpy_adj[1];
      istr >> l2genInput->rpy_adj[2];
      continue;
    }

    if (sParmWord.compare("C_DELTA_TND") == 0) {
      istr.clear(); istr.str( sValue);
      for (size_t irad=0; irad<NUMBER_OF_BEAMS; irad++) {
	for (size_t ipol=0; ipol<2; ipol++) {
	  istr >> l2genInput->c_deltaTND[irad][ipol][0];
	  istr >> l2genInput->c_deltaTND[irad][ipol][1];
	  istr >> l2genInput->c_deltaTND[irad][ipol][2];
	}
      }
      continue;
    }

    if (sParmWord.compare("TA_OCEAN_NOM") == 0) {
      istr.clear(); istr.str( sValue);
      size_t i=0;
      for (size_t irad=0; irad<NUMBER_OF_BEAMS; irad++) {
	for (size_t ipol=0; ipol<RADIOMETER_POLARIZATIONS; ipol++) {
	  istr >> l2genInput->ta_ocean_nom[i++];
	}
      }
      l2genInput->ta_nominal_files[0] = 0;
      continue;
    }

    if (sParmWord.compare("TA_NOMINAL_FILES") == 0) {
      if ( sValue.substr( 0, 1).compare("$") == 0) expandEnvVar( &sValue);
      int len = sValue.copy(l2genInput->ta_nominal_files, sValue.length());
      l2genInput->ta_nominal_files[len] = 0;
      l2genInput->ta_ocean_nom[0] = -1;
      continue;
    }

    if (sParmWord.compare("RADFLAGLIMITSFILE") == 0) {
      if ( sValue.substr( 0, 1).compare("$") == 0) expandEnvVar( &sValue);
      int len = sValue.copy(l2genInput->radflaglimitsfile, sValue.length());
      l2genInput->radflaglimitsfile[len] = 0;
      continue;
    }

    if (sParmWord.compare("XRAYFILE1") == 0) {
      int len = sValue.copy(l2genInput->xrayfile1, sValue.length());
      l2genInput->xrayfile1[len] = 0;
      continue;
    }

    if (sParmWord.compare("XRAYFILE2") == 0) {
      int len = sValue.copy(l2genInput->xrayfile2, sValue.length());
      l2genInput->xrayfile2[len] = 0;
      continue;
    }

    if (sParmWord.compare("RAD_OFFSET_CORR_FILE") == 0) {
      int len = sValue.copy(l2genInput->rad_offset_corr_file, sValue.length());
      l2genInput->rad_offset_corr_file[len] = 0;
      continue;
    }

    if (sParmWord.compare("CLIMATE_SAL_FILE") == 0) {
      int len = sValue.copy(l2genInput->climate_sal_file, sValue.length());
      l2genInput->climate_sal_file[len] = 0;
      continue;
    }

    if (sParmWord.compare("ANOMALY_STATUS") == 0) {
      int len = sValue.copy(l2genInput->anomaly_status, sValue.length());
      l2genInput->anomaly_status[len] = 0;
      continue;
    }

    if (sParmWord.compare("POINTING_ANOMALY_FILE") == 0) {
      int len = sValue.copy(l2genInput->pointing_anomaly_file, sValue.length());
      l2genInput->pointing_anomaly_file[len] = 0;
      continue;
    }

    if (sParmWord.compare("RFI_MASK_FILE") == 0) {
      int len = sValue.copy(l2genInput->rfi_mask_file, sValue.length());
      l2genInput->rfi_mask_file[len] = 0;
      continue;
    }

    if (sParmWord.compare("DI_U_COEFF_FILE") == 0) {
      int len = sValue.copy(l2genInput->dI_U_coeff_file, sValue.length());
      l2genInput->dI_U_coeff_file[len] = 0;
      continue;
    }

    if (sParmWord.compare("L2_UNCERTAINTY_MAPS_FILE") == 0) {
      if ( sValue.substr( 0, 1).compare("$") == 0) expandEnvVar( &sValue);
      int len = sValue.copy(l2genInput->l2_unc_maps_file, sValue.length());
      l2genInput->l2_unc_maps_file[len] = 0;
      continue;
    }

    if (sParmWord.compare("STATIC_BIAS_ADJ_FILE") == 0) {
      int len = sValue.copy(l2genInput->static_bias_adj_file, sValue.length());
      l2genInput->static_bias_adj_file[len] = 0;
      continue;
    }

    //    if (sParmWord.compare("SSS_SST_ADJ_PARMS") == 0) {
    //istr.clear(); istr.str( sValue);
    //istr >> l2genInput->sss_sst_adj_parms[0];
    //istr >> l2genInput->sss_sst_adj_parms[1];
    //istr >> l2genInput->sss_sst_adj_parms[2];
    //continue;
    //}

    if (sParmWord.compare("L_ACC_CORR_FILES") == 0) {
      int len = sValue.copy(l2genInput->l_acc_corr_files, sValue.length());
      l2genInput->l_acc_corr_files[len] = 0;
      continue;
    }

    if (sParmWord.compare("INSTRUMENT_GAIN_CORR") == 0) {
      istr.clear(); istr.str( sValue);
      // Order: V1 H1 V2 H2 V2 H2
      // itype=0: non-linear gain, itype=1:offset
      for (size_t itype=0; itype<2; itype++) {
        for (size_t irad=0; irad<NUMBER_OF_BEAMS; irad++) {
          for (size_t ipol=0; ipol<2; ipol++) {
            istr >> l2genInput->instrument_gain_corr[itype][irad][ipol];
          }
        }
      }
      continue;
    }

    if ((sParmWord.compare("-H") == 0) ||
	(sParmWord.compare("-HELP") == 0)) {
      iHelp = true;
      continue;
    }

    if ((sParmWord.compare( "-DUMP_OPTIONS_XMLFILE") == 0) ||
	(sParmWord.compare("--DUMP_OPTIONS_XMLFILE") == 0)) {
      int len = sValue.copy(xmlfile, sValue.length());
      xmlfile[len] = 0;
      iXml = true;
      continue;
    }

    if (sParmWord.compare("PAR") == 0 ||
	sParmWord.compare("SUITE") == 0) {
      ifstream parfile;
      if ( sParmWord.compare("SUITE") == 0)
	sValue = 
	  "$OCDATAROOT/aquarius/l2gen_aquarius_defaults_" + sValue + ".par";
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

  if ( iHelp || iXml) {
    // Setup clo list
    clo_optionList_t* list;
    list = (clo_optionList_t*) malloc(sizeof (clo_optionList_t));
    list->storageSize = 0;
    list->numOptions = 0;
    list->options = NULL;
    initOptions(list, l2genInput);
    cout << endl;
    if ( iHelp) clo_printUsage(list);
    if ( iXml) clo_writeXmlFile(list, xmlfile);
    exit(0);
  }


  // Determine number of products
  sValue.assign(l2genInput->l2prod);
  string::size_type found;
  *nProd = 1;
  found = sValue.find(",", 0);
  while (found < sValue.length()) {
    (*nProd)++;
    found = sValue.find(",", found+1);
  }

  ipos = 0;
  for (size_t i=0; i<MAXNAQPROD; i++) l2genInput->l2activeprodflag[i] = 0;
  for (size_t i=0; i<*nProd; i++) {
    posEndIdx = sValue.find(',', ipos);
    sParmWord  = sValue.substr( ipos, posEndIdx-ipos );
    //    cout << "Product " << i+1 << ": " << sParmWord.c_str() << endl;
    ipos = posEndIdx + 1;

    bool prod_found = false;
    for (size_t k=0; k<MAXNAQPROD; k++) {
      if (sParmWord.compare( shortProdName[k].c_str()) == 0) {
	l2genInput->l2activeprodflag[k] = 1;
	prod_found = true;
	break;
      }
    }
    if ( prod_found == false) {
      cout << "Improper product: " << sParmWord.c_str() << endl;
      exit(1);
    }
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


int initOptions( clo_optionList_t* list, instr *l2genInput) {

  char tmpStr[2048];

  clo_addOption(list, (char *) "par", CLO_TYPE_IFILE, NULL, 
		(char *) "input parameter file");
  strcpy(tmpStr, "processing version string");

  strcpy(tmpStr, "input L1 file name");
  clo_addOption(list, "ifile", CLO_TYPE_IFILE, NULL, (char*) tmpStr);

  strcpy(tmpStr, "output L2 file name,\n");
  clo_addOption(list, "ofile", CLO_TYPE_OFILE, NULL, tmpStr);

  strcpy(tmpStr, "L2 products to be included in ofile\n");
  strcat(tmpStr, "   l2prod = L2 products to be included in ofile");
  clo_addOption(list, "l2prod", CLO_TYPE_STRING, NULL, tmpStr);

  strcpy(tmpStr, "radiometer ancillary file 1");
  clo_addOption(list, "yancfile1", CLO_TYPE_IFILE, NULL, tmpStr);
  strcpy(tmpStr, "radiometer ancillary file 2");
  clo_addOption(list, "yancfile2", CLO_TYPE_IFILE, NULL, tmpStr);
  strcpy(tmpStr, "radiometer ancillary file 3");
  clo_addOption(list, "yancfile3", CLO_TYPE_IFILE, NULL, tmpStr);
  strcpy(tmpStr, "radiometer xray ancillary file 2");
  clo_addOption(list, "xrayfile2", CLO_TYPE_IFILE, NULL, tmpStr);

  return 0;
}


int initialize_static_data( instr *l2genInput, static_data *staticData)
{
  hid_t h5fid, grp;
  hsize_t start[6]={0,0,0,0,0,0}, count[6];
  string filename;

  H5E_auto_t old_func;
  void *old_client_data;
  H5Eget_auto(H5E_DEFAULT , &old_func, &old_client_data);

  hid_t fid;

  // Ocean Reflectance
  filename = l2genInput->rad_oceanrefl_file;
  expandEnvVar( &filename);
  cout << left << setw(35) << "Reading ocean reflectance file: " << 
    filename.c_str() << endl;
  h5fid =  H5Fopen( filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if ( h5fid == -1) {
    cout << "Ocean reflectance file: " << filename.c_str() << 
      " not found." << endl;
    exit(1);
  }

  // Sun file
  filename = l2genInput->rad_sun_file;
  expandEnvVar( &filename);
  cout << setw(35) << "Reading Sun file: " << filename.c_str() << endl;
  h5fid =  H5Fopen( filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if ( h5fid == -1) {
    cout << "Sun file: " << filename.c_str() << " not found." << endl;
    exit(1);
  }

  grp = H5Gopen1(h5fid,"/");

  count[0] = NZANG;
  Hdf::h5d_read(grp, (char*) "time_sun", &staticData->time_sun, 
		1, start, count);

  count[0] = NUMBER_OF_BEAMS;
  count[1] = NZANG;
  Hdf::h5d_read(grp, (char*) "eia_sun", &staticData->eia_sun, 2, start, count);

  count[0] = NUMBER_OF_BEAMS;
  count[1] = 3;
  count[2] = NZANG;
  count[3] = NOMEGA;
  Hdf::h5d_read(grp, (char*) "tasun_dir_tab", &staticData->tasun_dir_tab, 
		4, start, count);
  Hdf::h5d_read(grp, (char*) "tasun_ref_tab", &staticData->tasun_ref_tab, 
		4, start, count);

  H5Gclose(grp);
  H5Fclose(h5fid);


  // Sun bak file
  filename = l2genInput->rad_sunbak_file;
  expandEnvVar( &filename);
  cout << setw(35) << "Reading Sun Backscatter file: " << 
    filename.c_str() << endl;
  h5fid =  H5Fopen( filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if ( h5fid == -1) {
    cout << "Sun backscatter file: " << filename.c_str() << 
      " not found." << endl;
    exit(1);
  }

  grp = H5Gopen1(h5fid,"/");

  count[0] = NUMBER_OF_BEAMS;
  count[1] = 3;
  count[2] = 26;
  count[3] = 161;
  Hdf::h5d_read(grp, (char*) "tasun_bak_tab", &staticData->tasun_bak_tab, 
		4, start, count);

  H5Gclose(grp);
  H5Fclose(h5fid);


  // Galaxy file
  filename = l2genInput->rad_galwind_file;
  expandEnvVar( &filename);
  cout << setw(35) << "Reading Galaxy file: " << filename.c_str() << endl;
  h5fid =  H5Fopen( filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if ( h5fid == -1) {
    cout << "Galaxy file: " << filename.c_str() << " not found." << endl;
    exit(1);
  }

  grp = H5Gopen1(h5fid,"/");

  count[0] = NZANG;
  Hdf::h5d_read(grp, (char*) "time_galaxy", &staticData->time_galaxy, 
		1, start, count);
  double sidereal_year = 365.25636;
  //  double orbit_phase_dif= -0.84;
  double orbit_phase_dif= 3.00; // V2.0  08/10/12

  for (size_t i=0; i<NZANG; i++)
    staticData->time_galaxy[i] -= sidereal_year*(orbit_phase_dif+0.84);

  count[0] = NUMBER_OF_BEAMS;
  count[1] = NZANG;
  Hdf::h5d_read(grp, (char*) "eia_galaxy", &staticData->eia_galaxy, 
		2, start, count);

  count[0] = NUMBER_OF_BEAMS;
  count[1] = 3;
  count[2] = NZANG;
  count[3] = NOMEGA;
  Hdf::h5d_read(grp, (char*) "tagal_dir_tab", &staticData->tagal_dir_tab, 
		4, start, count);

  count[0] = 5;
  count[1] = NUMBER_OF_BEAMS;
  count[2] = 3;
  count[3] = NZANG;
  count[4] = NOMEGA;
  Hdf::h5d_read(grp, (char*) "tagal_ref_tab", &staticData->tagal_ref_tab, 
		5, start, count);

  H5Gclose(grp);
  H5Fclose(h5fid);


  // Land file
  filename = l2genInput->rad_landtables_file;
  expandEnvVar( &filename);
  cout << setw(35) << "Reading Land file: " << filename.c_str() << endl;
  h5fid =  H5Fopen( filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if ( h5fid == -1) {
    cout << "Land file: " << filename.c_str() << " not found." << endl;
    exit(1);
  }

  grp = H5Gopen1(h5fid,"/");

  count[0] = NUMBER_OF_BEAMS;
  count[1] = NZANG_LND;
  count[2] = NLON_LND;
  Hdf::h5d_read(grp, (char*) "fpt_lnd", &staticData->fpt_lnd, 3, start, count);

  count[0] = NUMBER_OF_BEAMS;
  count[1] = NZANG_LND;
  count[2] = NLON_LND;
  Hdf::h5d_read(grp, (char*) "frc_lnd", &staticData->frc_lnd, 3, start, count);

  H5Gclose(grp);
  H5Fclose(h5fid);


  // Land Correction file
  filename = l2genInput->rad_landcorr_file;
  expandEnvVar( &filename);
  cout << setw(35) << "Reading Land Correction file: " << 
    filename.c_str() << endl;
  h5fid =  H5Fopen( filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if ( h5fid == -1) {
    cout << "Land correction file: " << filename.c_str() << 
      " not found." << endl;
    exit(1);
  }

  grp = H5Gopen1(h5fid,"/");

  count[0] = NUMBER_OF_BEAMS;
  count[1] = 2;
  count[2] = 12;
  count[3] = 1440;
  count[4] = 1440;
  Hdf::h5d_read(grp, (char*) "land_dtb", &staticData->landcorr, 
		5, start, count);

  H5Gclose(grp);
  H5Fclose(h5fid);


  // Gain Ice file
  filename = l2genInput->rad_gainice_file;
  expandEnvVar( &filename);
  cout << setw(35) << "Reading Gain Ice file: " << filename.c_str() << endl;
  h5fid =  H5Fopen( filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if ( h5fid == -1) {
    cout << "Gain ice: " << filename.c_str() << " not found." << endl;
    exit(1);
  }

  grp = H5Gopen1(h5fid,"/");

  count[0] = NUMBER_OF_BEAMS;
  count[1] = 151;
  Hdf::h5d_read(grp, (char*) "gain_ice", &staticData->gain_ice, 
		2, start, count);

  H5Gclose(grp);
  H5Fclose(h5fid);


  // Tau sq file
  filename = l2genInput->rad_tausq_file;
  expandEnvVar( &filename);
  cout << setw(35) << "Reading Tau Sq file: " << filename.c_str() << endl;
  h5fid =  H5Fopen( filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if ( h5fid == -1) {
    cout << "Tau sq: " << filename.c_str() << " not found." << endl;
    exit(1);
  }

  grp = H5Gopen1(h5fid,"/");

  count[0] = 91;
  count[1] = 12;
  count[2] = 180;
  count[3] = 360;
  Hdf::h5d_read(grp, (char*) "itau", &staticData->itau, 4, 
		start, count);

  H5Gclose(grp);
  H5Fclose(h5fid);

  // coeffs_DI_U
  filename = l2genInput->dI_U_coeff_file;
  if ( filename != "") {
    expandEnvVar( &filename);
    cout << setw(35) << "Reading IU coupling file: " << 
      filename.c_str() << endl;
    h5fid =  H5Fopen( filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if ( h5fid == -1) {
      cout << "Non-Linear IU coupling file: " << filename.c_str() << 
        " not found." << endl;
      exit(1);
    }

    grp = H5Gopen1(h5fid,"/");

    count[0] = NUMBER_OF_BEAMS;
    count[1] = 4;
    Hdf::h5d_read(grp, (char*) "coeffs_DI_U", &staticData->coeffs_dI_U, 2, 
                  start, count);

    H5Gclose(grp);
    H5Fclose(h5fid);
  }

  // APC matrix
  filename = l2genInput->rad_apc_file;
  expandEnvVar( &filename);
  cout << setw(35) << "Reading APC file: " << filename.c_str() << endl;
  h5fid =  H5Fopen( filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if ( h5fid == -1) {
    cout << "APC matrix file: " << filename.c_str() << " not found." << endl;
    exit(1);
  }

  grp = H5Gopen1(h5fid,"/");

  count[0] = 3;
  count[1] = 3;
  count[2] = 3;
  Hdf::h5d_read(grp, (char*) "apc_matrix", &staticData->apc_matrix, 
		3, start, count);

  H5Gclose(grp);
  H5Fclose(h5fid);

  double a[3][3], a_inv[3][3], det;
  for (size_t i=0; i<3; i++) {

    for (size_t j=0; j<3; j++) {
      for (size_t k=0; k<3; k++) {
	a[j][k] = (double) staticData->apc_matrix[i][j][k];
      }
    }

    invert_3by3_( a, a_inv, &det);

    for (size_t j=0; j<3; j++) {
      for (size_t k=0; k<3; k++) {
	staticData->apc_inverse[i][j][k] = (float) a_inv[j][k];
      }
    }
  }


  // Wind estimated error array
  filename = l2genInput->wind_errortab_file;
  if ( filename != "") {
    expandEnvVar( &filename);
    cout << setw(35) << "Reading wind est error array file: " << 
      filename.c_str() << endl;
    h5fid =  H5Fopen( filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if ( h5fid == -1) {
      cout << "wind estimated error array file: " << filename.c_str() << 
	" not found." << endl;
      exit(1);
    }

    grp = H5Gopen1(h5fid,"/");

    count[0] = 40;
    count[1] = 7;
    Hdf::h5d_read(grp, (char*) "estimated_error_array", 
		  &staticData->estimated_error_array, 2, start, count);

    H5Gclose(grp);
    H5Fclose(h5fid);
  }

  // Harmonic coefficients, Roughness parameters for SIGMA0 SSS algorithm
  // Emiss coeff harmonics file (acoef)
  filename = l2genInput->emiss_coeff_harm_file;
  expandEnvVar( &filename);
  cout << setw(35) << "Reading emiss coeff harm file: " << 
    filename.c_str() << endl;
  h5fid =  H5Fopen( filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if ( h5fid == -1) {
    cout << "emiss coeff harm file: " << filename.c_str() << 
      " not found." << endl;
    exit(1);
  }

  grp = H5Gopen1(h5fid,"/");

  count[0] = 5;
  count[1] = NUMBER_OF_BEAMS;
  count[2] = 2;
  count[3] = 3;
  Hdf::h5d_read(grp, (char*) "acoef", &staticData->acoef, 
		4, start, count);

  // Turn off error handling
  H5Eset_auto( H5E_DEFAULT, NULL, NULL); 
  fid = H5Dopen1( grp, "wspd_max_a");
  if ( fid != -1) {
    H5Dclose( fid);
    // Restore previous error handler
    H5Eset_auto( H5E_DEFAULT, old_func, old_client_data); 
    
    count[0] = NUMBER_OF_BEAMS;
    count[1] = 2;
    count[2] = 3;
    Hdf::h5d_read(grp, (char*) "wspd_max_a", &staticData->wspd_max_a,
		  3, start, count);
  } else {
    double w0 = 28.5; // V2.0 linear extrapolation points 
    for (size_t irad=0; irad<NUMBER_OF_BEAMS; irad++) {
      for (size_t ipol=0; ipol<2; ipol++) { // 0=V,1=H
	for (size_t iharm=0; iharm<3; iharm++) {
	  staticData->wspd_max_a[irad][ipol][iharm] = w0;
	}
      }
    }
  }

  H5Gclose(grp);
  H5Fclose(h5fid);


  // Scat coeff harmonics file (bcoef)
  filename = l2genInput->scat_coeff_harm_file;
  expandEnvVar( &filename);
  cout << setw(35) << "Reading scat coeff harm file: " << 
    filename.c_str() << endl;
  h5fid =  H5Fopen( filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if ( h5fid == -1) {
    cout << "scat coeff harm file: " << filename.c_str() << 
      " not found." << endl;
    exit(1);
  }

  grp = H5Gopen1(h5fid,"/");

  count[0] = 5;
  count[1] = NUMBER_OF_BEAMS;
  count[2] = 4;
  count[3] = 3;
  Hdf::h5d_read(grp, (char*) "bcoef", &staticData->bcoef, 
		4, start, count);
    
  // Turn off error handling
  H5Eset_auto( H5E_DEFAULT, NULL, NULL); 
  fid = H5Dopen1( grp, "wspd_max_b");
  if ( fid != -1) {
    H5Dclose( fid);
    // Restore previous error handler
    H5Eset_auto( H5E_DEFAULT, old_func, old_client_data); 
    
    count[0] = NUMBER_OF_BEAMS;
    count[1] = 4;
    count[2] = 3;
    Hdf::h5d_read(grp, (char*) "wspd_max_b", &staticData->wspd_max_b,
		  3, start, count);
  } else {
    double w2[3] = {25.5, 22.5, 22.5};
    for (size_t irad=0; irad<NUMBER_OF_BEAMS; irad++) {
      for (size_t ipol=0; ipol<4; ipol++) {
	for (size_t iharm=0; iharm<3; iharm++) {
	  staticData->wspd_max_b[irad][ipol][iharm] = w2[iharm];
	}
      }
    }
  }

  H5Gclose(grp);
  H5Fclose(h5fid);


  // Roughness parameters
  filename = l2genInput->dtbw_win_sigma_file;
  expandEnvVar( &filename);
  cout << setw(35) << "Reading dtbw_win_sigma file: " << 
    filename.c_str() << endl;
  h5fid =  H5Fopen( filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if ( h5fid == -1) {
    cout << " dtbw_win_sigma file: " << filename.c_str() << 
      " not found." << endl;
    exit(1);
  }

  int32_t nbin_w=60;
  grp = H5Gopen1(h5fid,"/");
  count[0] = nbin_w;
  count[1] = nbin_w;
  count[2] = NUMBER_OF_BEAMS;
  count[3] = 2;
  count[4] = 3;
  Hdf::h5d_read(grp, (char*) "arr_dtbw", &staticData->arr_dtbw, 
		5, start, count);


  grp = H5Gopen1(h5fid,"/");
  count[0] = NUMBER_OF_BEAMS;
  count[1] = 2;
  Hdf::h5d_read(grp, (char*) "dsigma", &staticData->dsigma, 
		2, start, count);
    

  grp = H5Gopen1(h5fid,"/");
  count[0] = nbin_w;
  count[1] = nbin_w;
  count[2] = NUMBER_OF_BEAMS;
  count[3] = 2;
  count[4] = 3;
  Hdf::h5d_read(grp, (char*) "iflag_dtbw", &staticData->iflag_dtbw, 
		5, start, count);

  H5Gclose(grp);
  H5Fclose(h5fid);

  // WIN_WAV
  filename = l2genInput->dtbw_win_wav_file;
  if ( filename != "") {
    expandEnvVar( &filename);
    cout << setw(35) << "Reading dtbw_win_wav file: " << 
      filename.c_str() << endl;
    h5fid =  H5Fopen( filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if ( h5fid == -1) {
      cout << " dtbw_win_wav file: " << filename.c_str() << 
	" not found." << endl;
      exit(1);
    }

    int32_t nbin_winx=40;
    int32_t nbin_wavx=40;
    grp = H5Gopen1(h5fid,"/");
    count[0] = nbin_wavx;
    count[1] = nbin_winx;
    count[2] = NUMBER_OF_BEAMS;
    count[3] = 2;

    Hdf::h5d_read(grp, (char*) "arr_dtbwx", &staticData->arr_dtbwx, 
		  4, start, count);

    grp = H5Gopen1(h5fid,"/");
    count[0] = nbin_wavx;
    count[1] = nbin_winx;
    count[2] = NUMBER_OF_BEAMS;
    count[3] = 2;

    Hdf::h5d_read(grp, (char*) "iflag_winx_wavx", &staticData->iflag_winx_wavx, 
		  4, start, count);

    H5Gclose(grp);
    H5Fclose(h5fid);
  }


  // DTB_EMISS_WSPD
  filename = l2genInput->dtb_emiss_wspd_file;
  if ( filename != "") {
    expandEnvVar( &filename);
    cout << setw(35) << "Reading dtb_emiss_wspd file: " << 
      filename.c_str() << endl;
    h5fid =  H5Fopen( filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if ( h5fid == -1) {
      cout << " dtb_emiss_wspd file: " << filename.c_str() << 
	" not found." << endl;
      exit(1);
    }

    int32_t nsst_arr = 3001;
    grp = H5Gopen1(h5fid,"/");
    count[0] = nsst_arr;
    count[1] = NUMBER_OF_BEAMS;
    count[2] = 2;
    Hdf::h5d_read(grp, (char*) "yarr", &staticData->yarr, 
		  3, start, count);


    count[0] = NUMBER_OF_BEAMS;
    count[1] = 2;
    Hdf::h5d_read(grp, (char*) "T_STITCH", &staticData->T_STITCH, 
		  2, start, count);

    float misc[4];
    count[0] = 4;
    Hdf::h5d_read(grp, (char*) "misc", misc, 1, start, count);
    staticData->W_STITCH = misc[0]; 
    staticData->sst0 = misc[1]; 
    staticData->sst1 = misc[2]; 
    staticData->sst_step = misc[3]; 

    H5Gclose(grp);
    H5Fclose(h5fid);
  } else {
    staticData->W_STITCH = -1;
    cout << "No dtb_emiss_wspd file specified" << endl;
  }


  // Climatology salinity
  filename = l2genInput->climate_sal_file;
  if ( filename != "") {
    expandEnvVar( &filename);
    cout << setw(35) << "Reading Climatology salinity file: " << 
      filename.c_str() << endl;
    h5fid =  H5Fopen( filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if ( h5fid == -1) {
      cout << "Climatology salinity file: " << filename.c_str() << 
	" not found." << endl;
      exit(1);
    }

    grp = H5Gopen1(h5fid,"/");

    count[0] = 12;
    count[1] = 180;
    count[2] = 360;
    Hdf::h5d_read(grp, (char*) "climate salinity", &staticData->climate_sal, 
		  3, start, count);

    H5Gclose(grp);
    H5Fclose(h5fid);
  }


  // RFI mask file
  filename = l2genInput->rfi_mask_file;
  if ( filename != "") {
    expandEnvVar( &filename);
    cout << setw(35) << "Reading RFI mask file: " << filename.c_str() << endl;
    h5fid =  H5Fopen( filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if ( h5fid == -1) {
      cout << "RFI mask file: " << filename.c_str() << " not found." << endl;
      exit(1);
    }

    grp = H5Gopen1(h5fid,"/");

    count[0] = 90;
    count[1] = 180;
    count[2] = 2;

    Hdf::h5d_read(grp, (char*) "rfi_mask", 
    		  &staticData->rfi_mask, 3, start, count);

    H5Gclose(grp);
    H5Fclose(h5fid);
  }


  // dTa Galaxy Ref file                                                        
  filename = l2genInput->rad_dta_gal_file;
  if ( filename != "") {
    expandEnvVar( &filename);
    cout << setw(35) << "Reading dTa Galaxy file: " << filename.c_str() << endl;
    h5fid =  H5Fopen( filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if ( h5fid == -1) {
      cout << "dTa Galaxy file: " << filename.c_str() << " not found." << endl;
      exit(1);
    }

    grp = H5Gopen1(h5fid,"/");

    count[0] = NZANG;
    count[1] = NOMEGA;
    count[2] = NUMBER_OF_BEAMS;
    count[3] = 2;

    Hdf::h5d_read(grp, (char*) "dtagal_symm_ref_tab",
                  &staticData->dtagal_ref_tab, 4, start, count);

    H5Gclose(grp);
    H5Fclose(h5fid);
  } else {
    staticData->dtagal_ref_tab[0][0][0][0] = -9999.0;
  }

  return 0;
}


int comp_c_deltaTND( char* incalfilelist, float *c_delta)
{
  uint32_t nfiles = 0;
  FILE *fp = fopen( incalfilelist, "r");
  if (fp == NULL) {
    printf("Input listing file: \"%s\" not found.\n", incalfilelist);
    return -1;
  }
  char buf[2048];
  while(fgets(buf, 256, fp) != NULL) nfiles++;
  fclose(fp);

  typedef float c_delta_array[NUMBER_OF_BEAMS][2];
  c_delta_array *c_delta_file;
  c_delta_file = ( c_delta_array *) calloc( nfiles, sizeof(c_delta_array));


  // Open L2_CAL files
  // -----------------
  hid_t atype = H5Tcopy(H5T_C_S1);
  H5Tset_size(atype, 17);

  fp = fopen( incalfilelist, "r");
  for ( uint32_t ifile=0; ifile<nfiles; ifile++) {

    if ( fgets(buf, 256, fp) != NULL) {
      buf[strlen(buf)-1] = 0;
    } else {
      cout << "Error reading " << incalfilelist << endl;
      exit(1);
    }
      
    hid_t h5fid =  H5Fopen( buf, H5F_ACC_RDONLY, H5P_DEFAULT);
    if ( h5fid < 0) {
      cout << buf << " not found.";
      exit(1);
    }
    hid_t grp0 = H5Gopen1(h5fid, "/");
    hid_t attr;

    int32_t blknum;
    attr = H5Aopen_name(grp0, "Number of Blocks");
    H5Aread(attr, H5T_STD_I32LE, &blknum);
    H5Aclose(attr);

    float Texp_v[MAXCYC][NUMBER_OF_BEAMS];
    float Texp_h[MAXCYC][NUMBER_OF_BEAMS];
    float *Texp[2]={ &Texp_v[0][0], &Texp_h[0][0]};

    hid_t grp3 = H5Gopen1(h5fid, "/Aquarius Data");
    H5LTread_dataset_float( grp3, "rad_exp_TaV", &Texp_v[0][0]);
    H5LTread_dataset_float( grp3, "rad_exp_TaH", &Texp_h[0][0]);

    float Tf_v[MAXCYC][NUMBER_OF_BEAMS];
    float Tf_h[MAXCYC][NUMBER_OF_BEAMS];
    float *Tf[2]={ &Tf_v[0][0], &Tf_h[0][0]};
    H5LTread_dataset_float( grp3, "rad_TfV", &Tf_v[0][0]);
    H5LTread_dataset_float( grp3, "rad_TfH", &Tf_h[0][0]);

    float ice[MAXCYC][NUMBER_OF_BEAMS];
    float land[MAXCYC][NUMBER_OF_BEAMS];
    H5LTread_dataset_float( grp3, "rad_ice_frac", &ice[0][0]);
    H5LTread_dataset_float( grp3, "rad_land_frac", &land[0][0]);
    H5Gclose(grp3);

    float *calTemps;
    calTemps = ( float *) calloc(MAXCYC*NUMCALTEMPS, sizeof(float));

    hid_t grp5 = H5Gopen1(h5fid, "/Converted Telemetry");
    H5LTread_dataset_float( grp5, "rad_caltemps", calTemps);
    H5Gclose(grp5);

    int32_t idxTDL[3][2]={{23,21},{27,25},{31,29}};

    // <(Tf-Texp)/(Tf-Tdl)> 

    float *f;
    for (size_t irad=0; irad<NUMBER_OF_BEAMS; irad++) {
      for (size_t ipol=0; ipol<2; ipol++) { // 0=V,1=H
	f = (float *) calloc( blknum, sizeof(float));

	uint32_t k=0;
	for (int32_t i=0; i<blknum; i++) {
	  if ( land[i][irad] < 0.005 && ice[i][irad] == 0.0) {

	    float num_value = Tf[ipol][i*NUMBER_OF_BEAMS+irad] - 
	      Texp[ipol][i*NUMBER_OF_BEAMS+irad];

	    int32_t idx = i*NUMCALTEMPS+idxTDL[irad][ipol];
	    float Tdl = ((calTemps[idx] + calTemps[idx+1])/2.0) + 273.15;
	    float den_value = Tf[ipol][i*NUMBER_OF_BEAMS+irad] - Tdl;
	
	    f[k++] = num_value / den_value;
	  } // no land & no ice
	} // blk loop

	qsort( f, k, sizeof(float), compare);
	uint32_t m = k/2;
	if ( (k % 2) == 0) {
	  c_delta_file[ifile][irad][ipol] = (f[m] + f[m+1])/2;
	} else {
	  c_delta_file[ifile][irad][ipol] = f[m];
	}

	free (f);
      } // ipol loop
    } // irad loop

    H5Gclose(grp0);

    free( calTemps);

  } // ifile loop
  fclose(fp);

  double c0;
  uint32_t k=0;
  double *c_delta_d;
  for (size_t irad=0; irad<NUMBER_OF_BEAMS; irad++) {
    for (size_t ipol=0; ipol<2; ipol++) {

      c_delta_d = ( double *) calloc( nfiles, sizeof(double));
      for (size_t i=0; i<nfiles; i++)
	c_delta_d[i] = (double) c_delta_file[i][irad][ipol];
      c0 = gsl_stats_mean( c_delta_d, 1, nfiles);

      c_delta[k++] = (float) c0;

      free( c_delta_d);
    } // irad_loop
  } // ipol loop

  free( c_delta_file);

  return 0;
}


int get_taearth( int idayjl, int ibeamblk, static_data *staticData,
		 float *cellat, float clon360, float *celtht, float phir, 
		 float *gland, float *fland, float *fice, 
		 float *surtep, float *tran,
		 float *sss_reference, float *winspd, float *sm,
		 float *solar_flux, 
		 double *sec2000, double *zang, double *sun_zenith,
		 double *bore_sight, double *moonr,
		 float *ta_beam, float *ta_earth,
		 float *tagal_dir, float *tagal_ref, float *tasun_dir,
		 float *tasun_ref, float *tasun_bak, float *tamon_ref,
		 float *tagal_ref_GO, float *dtagal_ref)
{

  int32_t ibeam = ibeamblk % 3;
  int32_t iblk = ibeamblk / 3;
  int32_t irad = ibeam + 1;

  float surtep_celsius = surtep[ibeamblk] - 273.15;


  float transq, transq0=1, SST0, SSS0;
  uint8_t iopt[4]={1,1,1,1};

  // For the following call use celtht rather than thtadj because 
  // galaxy reflection tables are reference to boresight eia rather 
  // than effective eia
  float refl[2], tbsur[2], tbsur0[2];

  find_refl_tot_( &idayjl, 
		  &cellat[ibeamblk], &clon360, &celtht[ibeamblk],
		  &phir, &fland[ibeamblk], &fice[ibeamblk], 
		  &sss_reference[ibeamblk],
		  &surtep_celsius, winspd, &sm[ibeamblk],
		  staticData, refl, tbsur);

  transq = tran[ibeamblk] * tran[ibeamblk];

  float eia_gal;
  float nowind = 0.0;
  float tagal_spec[3];
  fd_ta_galaxy_( &irad, &sec2000[iblk], &zang[iblk], &nowind,
		 &staticData->eia_galaxy[0][0], 
		 &staticData->tagal_dir_tab[0][0][0][0],
		 &staticData->tagal_ref_tab[0][0][0][0][0], 
		 &staticData->time_galaxy[0],  
		 &eia_gal, &tagal_dir[3*ibeamblk], tagal_spec);

  fd_ta_galaxy_( &irad, &sec2000[iblk], &zang[iblk], winspd,
		 &staticData->eia_galaxy[0][0], 
		 &staticData->tagal_dir_tab[0][0][0][0],
		 &staticData->tagal_ref_tab[0][0][0][0][0], 
		 &staticData->time_galaxy[0],  
		 &eia_gal, &tagal_dir[3*ibeamblk], &tagal_ref[3*ibeamblk]);
  
  // empirical correction     
  tagal_ref_GO[3*ibeam+0] = tagal_ref[3*ibeamblk+0];
  tagal_ref_GO[3*ibeam+1] = tagal_ref[3*ibeamblk+1];
  tagal_ref_GO[3*ibeam+2] = tagal_ref[3*ibeamblk+2];

  if ( staticData->dtagal_ref_tab[0][0][0][0] != -9999.0) {
    fd_dta_symm_( &irad, &sec2000[iblk], &zang[iblk], 
    		  &staticData->time_galaxy[0], 
    		  &staticData->dtagal_ref_tab[0][0][0][0], 
		  &dtagal_ref[3*ibeam]);
    tagal_ref[3*ibeamblk+0] += dtagal_ref[3*ibeam+0];
    tagal_ref[3*ibeamblk+1] += dtagal_ref[3*ibeam+1];
  }

  // over land use tagal_spec 
  if ( gland[ibeamblk] > 0.5) {
    tagal_ref[3*ibeamblk+0] = tagal_spec[0];
    tagal_ref[3*ibeamblk+1] = tagal_spec[1];
    tagal_ref[3*ibeamblk+2] = tagal_spec[2];

    tagal_ref_GO[3*ibeam+0] = tagal_spec[0];
    tagal_ref_GO[3*ibeam+1] = tagal_spec[1];
    tagal_ref_GO[3*ibeam+2] = tagal_spec[2];
  }
         
  //cout << 0.5 * (tagal_spec[0]+tagal_spec[1]) << endl;
  //cout << 0.5 * (tagal_spec[0]-tagal_spec[1]) << endl;
  //cout << 0.5 * (tagal_ref[3*ibeamblk+0]+tagal_ref[3*ibeamblk+1]) << endl;
  //cout << 0.5 * (tagal_ref[3*ibeamblk+0]-tagal_ref[3*ibeamblk+1]) << endl;

  fd_ta_moon_( &irad, &bore_sight[3*ibeamblk], &moonr[ibeamblk], 
	       refl, &transq, 
	       &tamon_ref[3*ibeamblk]);

  fd_sun_backscatter_( &irad, &solar_flux[iblk], winspd, 
		       &sun_zenith[iblk], 
		       &staticData->tasun_bak_tab[0][0][0][0], 
		       &tasun_bak[3*ibeamblk]);

  // Convert ta_beam to Stokes
  vh_to_stokes_( ta_beam);

  // iopt[0]  set to 0 (1) to turn off (on) direct sun radiation
  // iopt[1]  set to 0 (1) to turn off (on) reflected sun radiation
  for (size_t i=0; i<3; i++) {
    ta_earth[3*ibeam+i] = ta_beam[i] - tagal_dir[3*ibeamblk+i] - 
      iopt[0]*tasun_dir[3*ibeamblk+i] - iopt[1]*tasun_ref[3*ibeamblk+i];
  }

  // refl coef assume for galaxy table, -999 means return istropic value
  // Note: sss & sst hard-coded values are used.
  float phir_iso = -999.0;
  SSS0 = 35.0;
  SST0 = 20.0;
  float refl0[2];
  fd_water_refl_exact_( &eia_gal, &SSS0, &SST0, winspd, &phir_iso,
			    staticData, refl0);

  // subtract reflected galaxy
  float ta_ref_adj[3];
  adjust_tagal_ref_( &irad, refl0, refl, 
		     &transq0, &transq, 
		     &ta_earth[3*ibeam], &tagal_ref[3*ibeamblk],
		     &staticData->apc_matrix[0][0][0], 
		     &staticData->apc_inverse[0][0][0], 
		     ta_ref_adj);

  for (size_t i=0; i<3; i++) {
    tagal_ref[3*ibeamblk+i] = ta_ref_adj[i];
  }

  // Do same for tagal_GO
  adjust_tagal_ref_( &irad, refl0, refl, 
		     &transq0, &transq, 
		     &ta_earth[3*ibeam], &tagal_ref_GO[3*ibeam],
		     &staticData->apc_matrix[0][0][0], 
		     &staticData->apc_inverse[0][0][0], 
		     ta_ref_adj);

  for (size_t i=0; i<3; i++) {
    tagal_ref_GO[3*ibeam+i] = ta_ref_adj[i];
  }


  adjust_tagal_ref_( &irad, refl0, refl, 
		     &transq0, &transq, 
		     &ta_earth[3*ibeam], &tasun_bak[3*ibeamblk],
		     &staticData->apc_matrix[0][0][0], 
		     &staticData->apc_inverse[0][0][0], 
		     ta_ref_adj);

  for (size_t i=0; i<3; i++) {
    tasun_bak[3*ibeamblk+i] = ta_ref_adj[i];
  }

  for (size_t i=0; i<3; i++) {
    ta_earth[3*ibeam+i] -= (tagal_ref[3*ibeamblk+i] + 
			    iopt[2]*tamon_ref[3*ibeamblk+i] + 
			    iopt[3]*tasun_bak[3*ibeamblk+i]);
  }
  return 0;
}


// Apply l_acc_corr to L_acc_raw
int l_acc_corr( string filename, float *L_acc_raw)
{
  cout << "No extrapolation applied" << endl;

  // Read L_acc correction file/data if specified
  typedef double l_acc_corr_array[2][10000];
  l_acc_corr_array *l_acc_corr_lut;

  l_acc_corr_lut = 
    ( l_acc_corr_array *) calloc( 6, sizeof(l_acc_corr_array));

  size_t found = filename.find(":");
  if (found == string::npos) {
    cout << "Delimiter not found in l_acc_corr_files parameter" << endl;
    exit(1);
  }

  const char *l_acc_corr_tag[6] = {"H1", "H2", "H3", "V1", "V2", "V3"};

  hid_t h5fid, grp;
  string filename1 = filename.substr(0, found);
  expandEnvVar( &filename1);
  h5fid =  H5Fopen( filename1.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if ( h5fid == -1) {
    cout << "Unable to open l_acc correction file: " << 
      filename1.c_str() << endl;
    exit(1);
  }

  grp = H5Gopen1(h5fid,"/");

  hsize_t start=0, count[6];
  for (size_t i=0; i<6; i++) {
    hid_t dset_id = H5Dopen1( grp, l_acc_corr_tag[i]);
    hid_t dspace_id = H5Dget_space(dset_id);
    herr_t status = H5Sget_simple_extent_dims( dspace_id, &count[i], NULL);
    
    Hdf::h5d_read(grp, l_acc_corr_tag[i], &l_acc_corr_lut[i][0][0], 
                  1, &start, &count[i]);

    H5Dclose(dset_id);
  }

  H5Gclose(grp);
  H5Fclose(h5fid);


  string filename2 = filename.substr(found+1);
  expandEnvVar( &filename2);
  h5fid =  H5Fopen( filename2.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if ( h5fid == -1) {
    cout << "Unable to open l_acc correction file: " << 
      filename2.c_str() << endl;
    exit(1);
  }

  grp = H5Gopen1(h5fid,"/");

  for (size_t i=0; i<6; i++) {
    hid_t dset_id = H5Dopen1( grp, l_acc_corr_tag[i]);

    Hdf::h5d_read(grp, l_acc_corr_tag[i], &l_acc_corr_lut[i][1][0], 
                  1, &start, &count[i]);

    H5Dclose(dset_id);
  }

  H5Gclose(grp);
  H5Fclose(h5fid);

  // Set up interpolation
  gsl_interp_accel *acc[6];
  gsl_interp *interp[6];

  for (size_t i=0; i<6; i++) {
    acc[i] = gsl_interp_accel_alloc ();
    interp[i] = gsl_interp_alloc (gsl_interp_linear, count[i]);
    gsl_interp_init (interp[i], &l_acc_corr_lut[i][0][0], 
                     &l_acc_corr_lut[i][1][0], count[i]);
  }

  // L_acc_raw(iacc,ipol,irad,icyc)

  // Turn off error handler to avoid potential aborts at last attTime
  gsl_error_handler_t *old_gsl_handler = gsl_set_error_handler_off();

  double l_acc_corr;
  int status;

  for (size_t icyc=0; icyc<MAXCYC; icyc++) {
    for (size_t irad=0; irad<NUMBER_OF_BEAMS; irad++) {
      for (size_t ipol=0; ipol<RADIOMETER_POLARIZATIONS; ipol++) { // V/P/M/H
        if ( ipol == 1 || ipol == 2) continue;

        int indx = (3 - ipol) + irad;

        int off = 
          RADIOMETER_LONG_ACCUM * 
          (RADIOMETER_POLARIZATIONS * 
           (NUMBER_OF_BEAMS*icyc + irad) + ipol);

        if ( L_acc_raw[off] != -9999) {
          double L_acc_double;

          double *xptr, *yptr, slope;

          if ( ipol == 0) {
            L_acc_double = (double) L_acc_raw[off+0];

            status = gsl_interp_eval_e (interp[indx], 
                                        &l_acc_corr_lut[indx][0][0],
                                        &l_acc_corr_lut[indx][1][0], 
                                        L_acc_double, acc[indx], 
                                        &l_acc_corr);

            if ( status == GSL_EDOM) {
              if ( L_acc_double < l_acc_corr_lut[indx][0][0]) {
                xptr = &l_acc_corr_lut[indx][0][0]; 
                yptr = &l_acc_corr_lut[indx][1][0]; 
              } else {
                xptr = &l_acc_corr_lut[indx][0][count[indx]-2]; 
                yptr = &l_acc_corr_lut[indx][1][count[indx]-2]; 
              }
              slope = (*(yptr+1) - *yptr) / (*(xptr+1) - *xptr);
              l_acc_corr = *yptr - slope * (*xptr - L_acc_double);
              l_acc_corr = (double) 0.0;  // no extrapolation for 031015
            }
            L_acc_raw[off+0] -= l_acc_corr;                             


            L_acc_double = (double) L_acc_raw[off+3];

            status = gsl_interp_eval_e (interp[indx], 
                                        &l_acc_corr_lut[indx][0][0],
                                        &l_acc_corr_lut[indx][1][0], 
                                        L_acc_double, acc[indx], 
                                        &l_acc_corr);

            if ( status == GSL_EDOM) {
              if ( L_acc_double < l_acc_corr_lut[indx][0][0]) {
                xptr = &l_acc_corr_lut[indx][0][0]; 
                yptr = &l_acc_corr_lut[indx][1][0]; 
              } else {
                xptr = &l_acc_corr_lut[indx][0][count[indx]-2]; 
                yptr = &l_acc_corr_lut[indx][1][count[indx]-2]; 
              }
              slope = (*(yptr+1) - *yptr) / (*(xptr+1) - *xptr);
              l_acc_corr = *yptr - slope * (*xptr - L_acc_double);
              l_acc_corr = (double) 0.0;  // no extrapolation for 031015
            }

            L_acc_raw[off+3] -= l_acc_corr;                             

          } else if ( ipol == 3) {
            L_acc_double = (double) L_acc_raw[off+0];

            status = gsl_interp_eval_e (interp[indx], 
                                        &l_acc_corr_lut[indx][0][0],
                                        &l_acc_corr_lut[indx][1][0], 
                                        L_acc_double, acc[indx], 
                                        &l_acc_corr);

            if ( status == GSL_EDOM) {
              if ( L_acc_double < l_acc_corr_lut[indx][0][0]) {
                xptr = &l_acc_corr_lut[indx][0][0]; 
                yptr = &l_acc_corr_lut[indx][1][0]; 
              } else {
                xptr = &l_acc_corr_lut[indx][0][count[indx]-2]; 
                yptr = &l_acc_corr_lut[indx][1][count[indx]-2]; 
              }
              slope = (*(yptr+1) - *yptr) / (*(xptr+1) - *xptr);
              l_acc_corr = *yptr - slope * (*xptr - L_acc_double);
              l_acc_corr = (double) 0.0;  // no extrapolation for 031015
            }

            L_acc_raw[off+0] -= l_acc_corr;         
            
            L_acc_double = (double) L_acc_raw[off+1];

            status = gsl_interp_eval_e (interp[indx], 
                                        &l_acc_corr_lut[indx][0][0],
                                        &l_acc_corr_lut[indx][1][0], 
                                        L_acc_double, acc[indx], 
                                        &l_acc_corr);

            if ( status == GSL_EDOM) {
              if ( L_acc_double < l_acc_corr_lut[indx][0][0]) {
                xptr = &l_acc_corr_lut[indx][0][0]; 
                yptr = &l_acc_corr_lut[indx][1][0]; 
              } else {
                xptr = &l_acc_corr_lut[indx][0][count[indx]-2]; 
                yptr = &l_acc_corr_lut[indx][1][count[indx]-2]; 
              }
              slope = (*(yptr+1) - *yptr) / (*(xptr+1) - *xptr);
              l_acc_corr = *yptr - slope * (*xptr - L_acc_double);
              l_acc_corr = (double) 0.0;  // no extrapolation for 031015
            }

            L_acc_raw[off+1] -= l_acc_corr;   
          }
        }
      }
    }
  }

  for (size_t i=0; i<6; i++) {
    gsl_interp_accel_free (acc[i]);
    gsl_interp_free (interp[i]);
  }

  gsl_set_error_handler( old_gsl_handler);

  return 0;
}
