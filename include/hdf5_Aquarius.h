//-*- Mode: C++; -*-

#ifndef hdf5_Aquarius_h
#define hdf5_Aquarius_h

#include "hdf5.h"
#define VOIDP void*

#define BLOCKS_PER_FRAME 4
#define RADIOMETER_SUBCYCLES 12
#define RADIOMETER_SIGNALS_PER_SUBCYCLE 5
#define NUMBER_OF_BEAMS 3
#define RADIOMETER_POLARIZATIONS 4
#define RADIOMETER_LONG_ACCUM 8
#define SCATTEROMETER_POLARIZATIONS 6
#define SCATTEROMETER_SUBCYCLES 8

#define SCATTEROMETER_TLM_OFFSET 164
#define SCATTEROMETER_SCI_OFFSET 201

#define RADIOMETER_SCI_OFFSET 537
#define RADIOMETER_ACCUM_OFFSET 599  // 537+62

#define DPLY_TLM_OFFSET  12
#define ICDS_STS_OFFSET  17  
#define ICDS_TLM_OFFSET  29
#define EXTT_TLM_OFFSET  53
#define APDU_TLM_OFFSET 123
#define ATC_TLM_OFFSET  128
#define ATC_TLM_SIZE     36

#define NUMCALTEMPS 85

// Pass/Cycle parameters
// 2011-08-25T00:25:39.8 corresponds to 998267154.8

#define NORBITCYCLE       103
#define FIRSTORBIT       1110
#define FIRSTTIME 998267154.8

#define SWAP_2(x) ( (((x) & 0xff) << 8) | ((unsigned short)(x) >> 8) )

#define SWAP_4(x) ( ((x) << 24) | \
         (((x) << 8) & 0x00ff0000) | \
         (((x) >> 8) & 0x0000ff00) | \
         ((x) >> 24) )

#define MAXSCAN 5200

#include <string>

#define MAXNAQPROD    200

enum { TAV, TAH, TA3, 
       TFV, TFH, TF3, 
       TAV0, TAH0, TA30, 
       TFV0, TFH0, TF30, 
       TOIV, TOIH, TOI3, 
       TOAV, TOAH, 
       FARTAH, 
       TAGALDV, TAGALDH, TAGALD3, 
       TAGALRV, TAGALRH, TAGALR3, 
       TAGALRGOV, TAGALRGOH,
       DTAGALV, DTAGALH, 
       TASUNDV, TASUNDH, TASUND3, 
       TASUNRV, TASUNRH, TASUNR3, 
       TASUNBV, TASUNBH, TASUNB3, 
       TAMONRV, TAMONRH, TAMONR3, 
       TBRAIN,
       TBV, TBH, TBVRC, TBHRC, 
       TBVRCLC, TBHRCLC, 
       TBVNORCLC, TBHNORCLC, 
       TBCON, TBCONLC,
       SSS, SSSLC, SSSERRRAN, SSSERRSYS, SSSADJ, 
       DENSITY, SPICE,
       TOAVLC, TOAHLC, HHWINDSPD, HHHWINDSPD,
       DTBSSTWSPDV, DTBSSTWSPDH,
       SRGHCORR,
       VVANT, HHANT, HVANT, VHANT,
       VVEXP, HHEXP, HVEXP, VHEXP,
       VVTOA, HHTOA, HVTOA, VHTOA, TOTTOA,
       VVKPCANT, HHKPCANT, HVKPCANT, VHKPCANT,
       VVKPCTOA, HHKPCTOA, HVKPCTOA, VHKPCTOA, TOTKPCTOA,
       SWINDUNC,
       SWINDSPD, ESURFV, ESURFH, ESURFVUNC, ESURFHUNC,
       AWINDSPD, AWINDDIR, ASURTEMP, ASURP, ACWAT, ASM, ASUBTEMP,
       ASWE, ASSS, AFLARE,
       ATRANS, ATBUP, ATBDW, ASWH,
       RLANDFRC, RICEFRC, RGICE, SLANDFRC, SICEFRC, 
       EXPTAV, EXPTAH, EXPTA3, EXPTAV_HHH, EXPTAH_HHH, EXPTA3_HHH, 
       EXPTBV, EXPTBH, EXPTBV0, EXPTBH0, VSM};

#ifdef AQUARIUS_L2PROD_DEFINE
std::string shortProdName[MAXNAQPROD] = 
  {"rad_TaV", "rad_TaH", "rad_Ta3", 
   "rad_TfV", "rad_TfH", "rad_Tf3", 
   "rad_TaV0", "rad_TaH0", "rad_Ta30", 
   "rad_TfV0", "rad_TfH0", "rad_Tf30", 
   "rad_toi_V", "rad_toi_H", "rad_toi_3",
   "rad_toa_V_nolc", "rad_toa_H_nolc",
   "rad_far_TaH", 
   "rad_galact_Ta_dir_V", "rad_galact_Ta_dir_H", "rad_galact_Ta_dir_3",
   "rad_galact_Ta_ref_V", "rad_galact_Ta_ref_H", "rad_galact_Ta_ref_3",
   "rad_galact_Ta_ref_GO_V", "rad_galact_Ta_ref_GO_H",
   "rad_galact_dTa_V", "rad_galact_dTa_H",
   "rad_solar_Ta_dir_V", "rad_solar_Ta_dir_H", "rad_solar_Ta_dir_3",
   "rad_solar_Ta_ref_V", "rad_solar_Ta_ref_H", "rad_solar_Ta_ref_3",
   "rad_solar_Ta_bak_V", "rad_solar_Ta_bak_H", "rad_solar_Ta_bak_3",
   "rad_moon_Ta_ref_V", "rad_moon_Ta_ref_H", "rad_moon_Ta_ref_3",
   "rad_rain_Tb", 
   "rad_TbV_nolc", "rad_TbH_nolc", "rad_TbV_rc_nolc", "rad_TbH_rc_nolc",
   "rad_TbV_rc", "rad_TbH_rc",
   "rad_TbV", "rad_TbH",
   "rad_Tb_consistency_nolc", "rad_Tb_consistency",
   "SSS_nolc", "SSS", "SSS_unc_ran", "SSS_unc_sys", "SSS_bias_adj",  
   "density", "spiciness",
   "rad_toa_V", "rad_toa_H",
   "rad_hh_wind_speed", "rad_hhh_wind_speed",
   "rad_dtb_sst_wspd_V", "rad_dtb_sst_wspd_H",
   "scat_rough_corr", 
   "scat_VV_ant", "scat_HH_ant", "scat_HV_ant", "scat_VH_ant", 
   "scat_VV_exp", "scat_HH_exp", "scat_HV_exp", "scat_VH_exp", 
   "scat_VV_toa", "scat_HH_toa", "scat_HV_toa", "scat_VH_toa",
   "scat_tot_toa",
   "Kpc_VV_ant", "Kpc_HH_ant", "Kpc_HV_ant", "Kpc_VH_ant",
   "Kpc_VV_toa", "Kpc_HH_toa", "Kpc_HV_toa", "Kpc_VH_toa",
   "Kpc_total",
   "wind_uncertainty",
   "scat_wind_speed", "scat_esurf_V", "scat_esurf_H", 
   "scat_esurf_V_uncertainty", "scat_esurf_H_uncertainty",
   "anc_wind_speed", "anc_wind_dir", "anc_surface_temp",
   "anc_surface_pressure", "anc_cwat", "anc_sm", "anc_subsurf_temp",
   "anc_swe", "anc_SSS", "anc_flare",
   "anc_trans", "anc_Tb_up", "anc_Tb_dw", "anc_swh", 
   "rad_land_frac", "rad_ice_frac", "rad_gice", 
   "scat_land_frac", "scat_ice_frac",
   "rad_exp_TaV", "rad_exp_TaH", "rad_exp_Ta3", 
   "rad_exp_TaV_hhh", "rad_exp_TaH_hhh", "rad_exp_Ta3_hhh", 
   "rad_exp_TbV", "rad_exp_TbH", "rad_exp_TbV0", "rad_exp_TbH0",
   "rad_sm"};
#endif

#define MAXCYC      5000 // Compatibility with Wentz fortran arrays

//      Gain/Offset structure
typedef struct gain_off_temp_struct {
  float gvv[MAXCYC][NUMBER_OF_BEAMS];
  float ghh[MAXCYC][NUMBER_OF_BEAMS];
  float ov [MAXCYC][NUMBER_OF_BEAMS];
  float oh [MAXCYC][NUMBER_OF_BEAMS];
  float op [MAXCYC][NUMBER_OF_BEAMS];
  float om [MAXCYC][NUMBER_OF_BEAMS];
  float gpv[MAXCYC][NUMBER_OF_BEAMS];
  float gph[MAXCYC][NUMBER_OF_BEAMS];
  float gmv[MAXCYC][NUMBER_OF_BEAMS];
  float gmh[MAXCYC][NUMBER_OF_BEAMS];
  float gpp[MAXCYC][NUMBER_OF_BEAMS];
  float gmm[MAXCYC][NUMBER_OF_BEAMS];
  float gpU[MAXCYC][NUMBER_OF_BEAMS];
  float gmU[MAXCYC][NUMBER_OF_BEAMS];
  float hpU[MAXCYC][NUMBER_OF_BEAMS];
  float hmU[MAXCYC][NUMBER_OF_BEAMS];
} gainoff;

float crs_dspl_temp_10( uint16_t tel_16);
float crs_dspl_temp_12( uint16_t tel_16);
float pls_dspl_temp_12( uint16_t tel_16);
float mns_dspl_temp_12( uint16_t tel_16);
float pls_dspl_ohms_12( uint16_t tel_16);
float mns_dspl_ohms_12( uint16_t tel_16);
float crs_dspl_ohms_12( uint16_t tel_16);

float fnep_ohms( float tel_16, float *extram);
float fnen_ohms( float tel_16, float *extram);
float finep_temp_12( uint16_t tel_16, float *extram,
		     float const_a, float const_b, float const_r0);
float finen_temp_12( uint16_t tel_16, float *extram,
		     float const_a, float const_b, float const_r0);
float coarse_temp_12( uint16_t tel_16, float *extram);
float coarse_temp_10( uint16_t tel_16, float *extram);

float ie_volt( uint8_t tel_8);

float resist_16bit_norm_disp( uint16_t tel_16);
float resist_16bit_ext_disp( uint16_t tel_16);

float temp_8bit_ext_fine( uint8_t tel_8,
			  float const_a, float const_b, float const_r0);
float temp_8bit_norm_fine( uint8_t tel_8,
			   float const_a, float const_b, float const_r0);

float temp_8bit_ext( uint8_t tel_8);
float temp_8bit_rad6k( uint8_t tel_8);
float temp_8bit_norm( uint8_t tel_8);
float temp_apdu_8bit( uint8_t tel_8, float const_ap);
float temp_8bit_norm_abr0( uint8_t tel_8, 
			   float const_a, float const_b, float const_r0);
float temp_8bit_ext_abr0( uint8_t tel_8, 
			  float const_a, float const_b, float const_r0);
float temp_8bit_norm_scat( uint8_t tel_8, 
			   float const_a, float const_b, float const_r0);
float temp_16bit_norm_disp( uint16_t tel_16, float gain, float offset,
			    float const_a, float const_b, float const_r0);
float temp_16bit_ext_disp( uint16_t tel_16, float gain, float offset,
			   float const_a, float const_b, float const_r0,
			   float const_h_r);
float temp_16bit_norm_fine( uint16_t tel_16, uint16_t *extram_et,
			    float const_a, float const_b, float const_r0);
float temp_16bit_ext_fine( uint16_t tel_16, uint16_t *extram_et,
			   float const_a, float const_b, float const_r0,
			   float const_h_r);
float temp_resist( float a, float b, float r0, float resist);
float dwnlk_fltg_pt( uint16_t tel_16);
uint16_t in32_out16( uint32_t i32, uint16_t nbits, uint16_t mask);


namespace Hdf {

  class hdf5_Aquarius {
    hid_t h5fid;
    hid_t grp0, grp1, grp2, grp3, grp4, grp5, grp6;
    uint32_t maxscan;
    int32_t blknum;
    uint8_t offset_corr;
    int32_t rad_frmnum;
    int32_t atc_frmnum;
    unsigned openFlags;
    int32_t nProd;

    bool first_write;
    uint32_t first_icds_blknum;
    uint8_t first_atc_subframe;

    double granuleStart;
    double granuleStop;
    double orbitStart;
    double orbitStop;
    int32_t orbitNumber;
    int32_t cycleNumber;
    int32_t passNumber;
    double nodeCrossingTime;
    float nodeLongitude;
  public:
    hdf5_Aquarius();
    ~hdf5_Aquarius();

    int createl1( char* l1_filename, int32_t numBlocks,
		  int32_t num_SACD_HKT, double granStart, 
		  char* inputFiles, char* processControl,
		  char *softwareId);
    int writel1( char* ptr);
    int readl1_radiometer( int32_t blknum, int32_t frmnum, uint8_t subframe, 
			   double *blkSec, 
		 	   unsigned short *radiomHdr, 
			   unsigned short *radiomEar, 
			   unsigned short *radiomCnd, 
			   unsigned short *radiomCal);
    int readl1_frame_info( int32_t blknum, 
			   int32_t *atc_frmnum, uint8_t *atc_subframe,
			   int32_t *rad_frmnum, uint8_t *rad_subframe);
    int writel1_frame_info( int32_t blknum, 
			    int32_t *atc_frmnum, uint8_t *atc_subframe,
			    int32_t *rad_frmnum, uint8_t *rad_subframe);
    int readl1_icds_blknum( int32_t blknum, int32_t *icds_blknum);
    int readl1_dpu_status_tlm( int32_t rad_frmnum, int32_t rad_subframe,
			       int32_t status_byte, uint8_t *dpu_status_tlm);
    int readl1_radiom_nrt_tlm( int32_t rad_frmnum, int32_t rad_subframe,
			       int32_t status_byte, uint8_t *radiom_nrt_tlm);
    int closel1();
    int closel1( double granuleStop);
    int closel1( double granuleStart, double granuleStop,
		 double nodeCrossingTime, float nodeLongitude,
		 int32_t *oplut_accum, float percent_non_default,
		 int32_t missingBlks);
    int openl1( char* l1_filename, unsigned flags);   
    int createl1_eph( int32_t numOrbVec);
    int writel1_eph( int32_t numOrbVec, double *time, 
		     double *pos, double *vel);
    int setOrbVec( int32_t numOrbVec);
    int createl1_att( int32_t numAttSamp);
    int writel1_att( int32_t numAttSamp, double *time, 
		     double *ang, double *quat, 
		     uint32_t *flags);
    int setAttSamp( int32_t numAttSamp);

    int readl1_caltemps( int32_t blknum, int32_t frmnum, uint8_t subframe, 
			 float *caltemps);
    int readl1_eph( uint32_t *numOrbVec);
    int readl1_eph( double *time, double *pos, double *vel);

    int readl1_att( uint32_t *numAttSamp);
    int readl1_att( double *time, double *ang, double *quat);

    int readl1_hkt( uint32_t *numHktRec);
    int readl1_hkttme( uint32_t numHktRec, uint32_t *time);
    int readl1_hktacsmode( uint32_t numHktRec, uint8_t *acsmode);
    int writel1_hkt( uint32_t recnum, uint8_t *hktrec);

    int createl2( char* l2_filename, int32_t numBlocks,
		  int32_t nprod, int8_t *l2activeprodflag,
		  double granStart,
		  double nodeCrossingTime, float nodeLongitude,
		  int32_t orbitNumber, int32_t cycleNumber, 
		  int32_t passNumber,
		  char* inputFiles, char* processControl,
		  char* processControlScat,
		  char *softver, char *pversion, bool iopt_l1b);

    int createl2sm( char* l2_filename, int32_t numBlocks,
		    char* starttime, char* endtime,   
		    double nodeCrossingTime, float nodeLongitude,
		    int32_t orbitNumber, int32_t cycleNumber, 
		    int32_t passNumber,
		    char* inputFiles, char* processControl,
		    char *softver, char *pversion, char *ancillary_files,
		    char *nominalNav, char *anomaly_status, 
                    char *input_L2_file);

    int openl2( char* l2_filename, unsigned flags, int32_t *numBlks,
		double *granstart, double *granstop);
    int writel2_nav( int32_t nBlks, bool *inout, 
		     double *pos, double *vel,
		     double *ang,
		     float *clat, float *clon,
		     float *celtht, float *celphi, 
		     float *suntht, float *sunphi, 
		     float *sunglt, float *moonglt, 
		     float *glxlat, float *glxlon,
		     float *cellatfoot, float *cellonfoot,
		     double *sund, double *sunr, double *moond,
		     double *zang,
		     double *sclon, double *sclat, double *scalt,
		     float *scat_clon, float *scat_clat,
		     float *scat_elon, float *scat_elat,
		     float *scat_polarization_roll,
		     uint8_t *acsmode,
		     bool browse, bool iopt_l1b);
    int writel2_navsm( int32_t nBlks, float *clat, float *clon,
		       double *zang, double *ang);
    int readl2_nav( float *beam_clat, float *beam_clon, double *zang);
    int readl2_nav( float *beam_clat, float *beam_clon,
		    float *cellatfoot, float *cellonfoot);
    int readl2radflag( int32_t iblk, hid_t grp, VOIDP data);
    int readl2radflagsm( int32_t iblk, hid_t grp, VOIDP data);

    int writel2( int32_t numBlocks, hid_t grp, char *prodName, VOIDP data);
    int writel2radflag( int32_t iblk, hid_t grp, VOIDP data);
    int writel2radflagsm( int32_t iblk, hid_t grp, VOIDP data);
    int writel2radrfi( int32_t iblk, hid_t grp, VOIDP data);
    int writel2scatrfi( int32_t iblk, hid_t grp, VOIDP data);
    int writel2samples( int32_t iblk, hid_t grp, 
			char *prodName, VOIDP data);
    int writel2sec( int32_t iblk, hid_t grp, VOIDP data);
    int writeSolarXrayFlux( int32_t iblk, hid_t grp, VOIDP data);
    int writel2_caltemps( int32_t iblk,  VOIDP data);
    int writel2_gainoff( int32_t iblk, gainoff *gainOff, bool iopt_l1b);
    int writel2_s_acc( int32_t nBlks, 
                       const char *name, const char *description,
                       float *s_acc);
    int writel2_l_acc( int32_t nBlks,
                       const char *name, const char *description,
                       float *l_acc);
    int writel2_sec_TA( int32_t nBlks, 
                        Hdf::hdf5_Aquarius *l1afile, 
                        float *TA);
    int writel2_drift_corr( int32_t nBlks, 
                            float *TA_hat_drift_corr);
    int readl2sec( int32_t iblk, hid_t grp, VOIDP data);

    int readl2( int32_t numBlocks, hid_t grp, char *prodName, VOIDP data);

    int copyl1( int32_t recIdx, int32_t frmIdx, int32_t wrtIdx,
			     hid_t grp[2], char *dsName);
    int closel2();
    int closel2( double granuleStop, float percentWater,
		 float percentRFI, float *solar_attr,
		 std::string *rad_anc_files, std::string *scat_anc_files,
		 std::string *rad_calfiles, std::string *rad_tables,
		 char *anomaly_status, float *c_deltaTND, 
		 const char *radLimitsStr, const char *nominalNav,
		 float *offset_corr);
    hid_t getH5fid();
    int getH5gid( hid_t, hid_t*);
    uint32_t nBlks();
    int incBlks();
    int setFrameNums( int32_t out_rad_frmnum, int32_t out_atc_frmnum);
    int getGranuleTimes(double *granStart, double *granStop);
    int getOrbitTimes(double *orbStart, double *orbStop, int32_t *orbNumber,
		      int32_t *cycNumber, int32_t *psNumber);
    int getNodeInfo(double *, float *);
    int writeOrbitTimes(hid_t grp, double orbStart, double orbStop);
  };

}

#endif

