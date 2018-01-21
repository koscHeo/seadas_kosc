#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stack>
#include <list>
#include <libgen.h>
#include "l2gen_aquarius.h"


#define VERSION "3.01"

//    Modification history:
// Programmer  Organization  Date      Version  Description of change
// ----------  ------------  ----      -------  ---------------------
// Joel Gales  FutureTech    05/20/12  0.10     Original development
// Joel Gales  FutureTech    05/20/12  0.90     Add additional metadata
// Joel Gales  FutureTech    11/07/12  0.91     Update for 121031 update
// Joel Gales  FutureTech    03/12/13  1.00     Add zang to navigation group
// Joel Gales  FutureTech    04/26/13  1.10     Add rpy to navigation group
//                                              Add anc_swe, anc_sm, gland,
//                                              gice, scat_toa
// Joel Gales  FutureTech    07/25/13  1.11     Add additional attributes for
//                                              NSICD compliance
// Joel Gales  FutureTech    09/30/14  3.00     Update version label to 3.0
// Joel Gales  FutureTech    10/02/14  3.01     Add "Input L2 File" attribute

extern "C" double zulu2unix( char *zulu);

using namespace std;

int main (int argc, char* argv[])
{
  cout << "l2gen_sm " 
       << VERSION << " (" 
       <<  __DATE__ << " " 
       << __TIME__ << ")" << endl;

  if ( argc == 1) {
    cout << endl << 
      "l2gen_sm ifile=L2_SM_textfile ofile=L2_SM_outputfile pversion=pversion"
	 << endl;
    return 0;
  }

  // blkSec is seconds from 01/06/1980
  double blkSec;

  int32_t iblk, ibeam, first_blk, last_blk;
  uint32_t rad_flag[NUMBER_OF_BEAMS];

#include "l2gen_aquarius_prod.inc"

  double *rpy;
  rpy = new double[MAXCYC*3];

  instr l2genInput;
  uint32_t nProd;
  string processControl;
  parseInput(argc, argv, &l2genInput, &nProd, &processControl);


  char starttime[14], endtime[14];
  char txtbuf[512];
  char input_L2_file[320];

  ifstream in_SM;
  in_SM.open ( l2genInput.ifile, ifstream::in);
  if ( in_SM.fail() == true) {
    cout << "Cannot open: "  << l2genInput.ifile << endl;
    exit(1);
  }

  istringstream istr;

  double nodeCrossingTime;
  float nodeLongitude;
  int32_t orbitNumber, cycleNumber, passNumber;
  string nominalNav, anomaly_status, ancillary_files;
  string strbuf;

  iblk = 0;
  in_SM.getline( input_L2_file, 320);
  strbuf.assign( basename(input_L2_file));
  size_t found = strbuf.find_first_of( " ", 0);
  strbuf.substr( 0, found);
  strcpy( input_L2_file, strbuf.c_str());
  input_L2_file[found] = 0;

  in_SM.getline( starttime, 14);
  in_SM.getline( endtime, 14);

  in_SM.getline( txtbuf, 512);
  txtbuf[15] = 0;
  nodeCrossingTime = unix2gpstai(zulu2unix( txtbuf));

  in_SM.getline( txtbuf, 512);
  istr.clear(); istr.str( txtbuf);
  istr >> nodeLongitude;

  in_SM.getline( txtbuf, 512);
  istr.clear(); istr.str( txtbuf);
  istr >> orbitNumber;

  in_SM.getline( txtbuf, 512);
  istr.clear(); istr.str( txtbuf);
  istr >> cycleNumber;

  in_SM.getline( txtbuf, 512);
  istr.clear(); istr.str( txtbuf);
  istr >> passNumber;

  unsigned int last_nonblk;
  in_SM.getline( txtbuf, 512);
  nominalNav = txtbuf;
  last_nonblk = nominalNav.find_last_not_of(" ");
  nominalNav = nominalNav.substr( 0, last_nonblk+1);

  in_SM.getline( txtbuf, 512);
  anomaly_status = txtbuf;
  last_nonblk = anomaly_status.find_last_not_of(" ");
  anomaly_status = anomaly_status.substr( 0, last_nonblk+1);

  in_SM.getline( txtbuf, 512);
  ancillary_files = txtbuf;
  last_nonblk = ancillary_files.find_last_not_of(" ");
  ancillary_files = ancillary_files.substr( 0, last_nonblk+1);


  while( !in_SM.eof()) {
    in_SM.getline( txtbuf, 512);
    iblk++;
  }
  in_SM.close();
  iblk--;
  iblk /= 2;

  // Create L2 output file
  static Hdf::hdf5_Aquarius l2file;
  double granStart = blkSec;
  int32_t nactBlks = iblk / 3;
  string inputFiles = basename( l2genInput.ifile);
  int32_t nBlks ;
  l2file.createl2sm(l2genInput.ofile, nactBlks, 
		    starttime, endtime,
		    nodeCrossingTime, nodeLongitude, orbitNumber,
		    cycleNumber, passNumber,
		    (char *) inputFiles.c_str(), 
		    (char *) processControl.c_str(),
		    (char *) VERSION, l2genInput.pversion,
		    (char *) ancillary_files.c_str(),
		    (char *) nominalNav.c_str(),
		    (char *) anomaly_status.c_str(),
                    (char *) input_L2_file);

  hid_t gid[6];
  l2file.getH5gid( l2file.getH5fid(), gid);

  in_SM.open ( l2genInput.ifile, ifstream::in);


  ///////////////
  // Main Loop //
  ///////////////
  int32_t icnt = 0;

  double granStop;

  // Skip header entries
  for (size_t i=0; i<11; i++) in_SM.getline( txtbuf, 512);

  while ( icnt < nactBlks) {

    if ((icnt % 500) == 0) cout << "icnt: " << icnt << endl;

    ibeam = -1;

    while ( ibeam < 2) {
      in_SM.getline( txtbuf, 512);

      strbuf = txtbuf;
      istr.clear(); istr.str( strbuf.substr( 0, 5));
      istr >> ibeam;
      ibeam--;

      istr.clear(); istr.str( strbuf.substr( 6, 15));
      istr >> blkSec;

      istr.clear(); istr.str( strbuf.substr( 21, 10));
      istr >> cellat[ibeam+NUMBER_OF_BEAMS*icnt];

      istr.clear(); istr.str( strbuf.substr( 31, 10));
      istr >> cellon[ibeam+NUMBER_OF_BEAMS*icnt];

      istr.clear(); istr.str( strbuf.substr( 41, 10));
      istr >> zang[icnt];
    
      istr.clear(); istr.str( strbuf.substr( 51, 10));
      istr >> celphi[ibeam+NUMBER_OF_BEAMS*icnt];

      istr.clear(); istr.str( strbuf.substr( 61, 12));
      istr >> tb_sur[1+2*(ibeam+NUMBER_OF_BEAMS*icnt)];

      istr.clear(); istr.str( strbuf.substr( 73, 12));
      istr >> tb_sur[0+2*(ibeam+NUMBER_OF_BEAMS*icnt)];

      istr.clear(); istr.str( strbuf.substr( 85, 12));
      istr >> surtep[ibeam+NUMBER_OF_BEAMS*icnt];

      istr.clear(); istr.str( strbuf.substr( 97, 12));
      istr >> subtep[ibeam+NUMBER_OF_BEAMS*icnt];

      istr.clear(); istr.str( strbuf.substr(109, 12));
      istr >> vol_sm[ibeam+NUMBER_OF_BEAMS*icnt];

      istr.clear(); istr.str( strbuf.substr(121, 8));
      istr >> rad_flag[ibeam];


      in_SM.getline( txtbuf, 512);
      strbuf = txtbuf;

      if ( ibeam == 0) {
	istr.clear(); istr.str( strbuf.substr( 0, 10));
	istr >> rpy[0+3*icnt];
	istr.clear(); istr.str( strbuf.substr(10, 10));
	istr >> rpy[1+3*icnt];
	istr.clear(); istr.str( strbuf.substr(20, 10));
	istr >> rpy[2+3*icnt];
      }

      istr.clear(); istr.str( strbuf.substr(30, 10));
      istr >> gland[ibeam+NUMBER_OF_BEAMS*icnt];

      istr.clear(); istr.str( strbuf.substr(40, 10));
      istr >> gice[ibeam+NUMBER_OF_BEAMS*icnt];

      istr.clear(); istr.str( strbuf.substr(50, 10));
      istr >> swe[ibeam+NUMBER_OF_BEAMS*icnt];

      istr.clear(); istr.str( strbuf.substr(60, 10));
      istr >> sm[ibeam+NUMBER_OF_BEAMS*icnt];

      istr.clear(); istr.str( strbuf.substr(70, 10));
      istr >> scat_toa[0+4*ibeam];

      istr.clear(); istr.str( strbuf.substr(80, 10));
      istr >> scat_toa[1+4*ibeam];

      istr.clear(); istr.str( strbuf.substr(90, 10));
      istr >> scat_toa[2+4*ibeam];

      istr.clear(); istr.str( strbuf.substr(100, 10));
      istr >> scat_toa[3+4*ibeam];

    /*
      write(3,'(i5,f15.3,4f10.2,5f12.4,i8)') kk,secGPS,lat,lon,zang,phi, &
      tbh,tbv,lkt2d(kk,k),lst2d(kk,k),mv,sm_flag
      write(3,'(5f10.5,6f10.2)') roll,pitch,yaw,lfr,ifr,swe,vsm2d(kk,k),hh,hv,vh,vv

      ! kk - Beam Number
      ! lat - Latitude
      ! lon - Longitude
      ! zang - intra-orbit angle
      ! phi - Azimuth angle (used to determine the direction of the orbit (Asc/Dsc) phi > 270 degrees is Asc
      ! tbh - Aquarius h pol radiometer observations (K)
      ! tbv - Aquarius v pol radiometer observations (K)
      ! lkt2d - NCEP surface temperature (K)
      ! lst2d - NCEP subsurface (0-10 cm) (temperature (K)
      ! mv - Volumetric Soil Moisture (m3/m3)
      ! sm_flag - Bit flag for soil moisture retrievals

      ! roll,pitch,yaw - Spacecraft Roll,Pitch,Yaw
      ! lfr - gland fraction
      ! ifr - gice fraction
      ! swe - snow water equ
      ! vsm - ancillary soil moisture
      ! hh,hv,vh,vv - TOA Scatterometer NRCS
    */
    } // while (ibeam)

    l2file.writel2sec(icnt, gid[1], (VOIDP) &blkSec);

    // Write radiometer flags
    l2file.writel2radflagsm( icnt, gid[2], (VOIDP) rad_flag);

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
	    dptr = prodPtr[i] + icnt*NUMBER_OF_BEAMS*prodMultiplicity[i] +
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

    icnt++;
  } // (End of Main Loop)


  l2file.writel2_navsm( nactBlks, cellat, cellon, zang, rpy);
  granStop = blkSec;

  l2file.closel2();

  // Close scatterometer output text files
  in_SM.close();

  free( rpy);

  cout << "Normal Completion" << endl << endl;

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

  procControl->clear();

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
  //  sValue = "$OCDATAROOT/aquarius/l2gen_aquarius_defaults.par";
  //expandEnvVar( &sValue);
  //sLine = "par=" + sValue;
  //parmStack.push( sLine);


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


    if (sParmWord.compare("PVERSION") == 0) {
      int len = sValue.copy(l2genInput->pversion, sValue.length());
      l2genInput->pversion[len] = 0;
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


  // Determine number of products
  strcpy( l2genInput->l2prod, 
	  "anc_surface_temp,anc_subsurf_temp,anc_swe,anc_sm,rad_land_frac,rad_ice_frac,");
  strcat( l2genInput->l2prod, 
	  "rad_sm,rad_TbH,rad_TbV,scat_VV_toa,scat_HH_toa,scat_HV_toa,scat_VH_toa");

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



