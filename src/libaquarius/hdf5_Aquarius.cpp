#include "hdf5_Aquarius.h"

namespace Hdf {
  hdf5_Aquarius::hdf5_Aquarius() {
    h5fid = -1;
  }


  hdf5_Aquarius::~hdf5_Aquarius() {
    
  }


  hid_t hdf5_Aquarius::getH5fid() {
    return h5fid;
  }

  int hdf5_Aquarius::getH5gid(hid_t h5fid, hid_t *gid) {
    gid[0] = grp0;
    gid[1] = grp1;
    gid[2] = grp2;
    gid[3] = grp3;
    gid[4] = grp4;
    gid[5] = grp5;
    gid[6] = grp6;
    return 0;
  }

  uint32_t hdf5_Aquarius::nBlks() {
    return blknum;
  }

  int hdf5_Aquarius::incBlks() {
    blknum++;
    return 0;
  }

  int hdf5_Aquarius::setFrameNums( int32_t out_rad_frmnum, 
				   int32_t out_atc_frmnum) {
    rad_frmnum = out_rad_frmnum;
    atc_frmnum = out_atc_frmnum;
    return 0;
  }


  int hdf5_Aquarius::getGranuleTimes(double *granStart, double *granStop) {
    if ( granStart != NULL) *granStart = granuleStart;
    if ( granStop != NULL)  *granStop = granuleStop;
    return 0;
  }

  int hdf5_Aquarius::getOrbitTimes(double *orbStart, double *orbStop,
				   int32_t *orbNumber, int32_t *cycNumber,
				   int32_t *psNumber) {
    *orbStart = orbitStart;
    *orbStop = orbitStop;
    *orbNumber = orbitNumber;
    *cycNumber = cycleNumber;
    *psNumber = passNumber;
    return 0;
  }

  int hdf5_Aquarius::getNodeInfo(double *nodeCrossTime, 
				 float *nodeLong) {
    *nodeCrossTime = nodeCrossingTime;
    *nodeLong = nodeLongitude;
    return 0;
  }
}


