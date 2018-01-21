#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_statistics.h>

#define  AQUARIUS_L2PROD_DEFINE 1

#include "hdf5_Aquarius.h"
#include "hdf5utils.h"
#include "H5LTpublic.h"
#include <clo.h>

#define NUMBEAMPOL  NUMBER_OF_BEAMS*RADIOMETER_POLARIZATIONS
#define NUMRADIOEAR NUMBEAMPOL*RADIOMETER_SIGNALS_PER_SUBCYCLE*RADIOMETER_SUBCYCLES
#define NUMRADIOCND NUMBEAMPOL*RADIOMETER_SUBCYCLES
#define NUMRADIOCAL NUMBEAMPOL*RADIOMETER_LONG_ACCUM

#define MAXCYC      5000 // Compatibility with Wentz fortran arrays

#define PRODSTRLEN     2048  /* String length limit for product spec */

#define NZANG       1441
#define NOMEGA      1441
#define NLON_LND    2160
#define NZANG_LND   2881

#define NBIN_W      60
#define NWIN_SC     40
#define NBIN_WINX   40
#define NBIN_WAVX   40

#define ISSTMAX     30
#define NSST_ARR    (ISSTMAX*100+1)

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
#define NONNOMCMD     20
#define REFL_1STOKES  21
#define SPARE2        22
#define RFI_REGION    23

typedef struct input_struct {
  char ifile[FILENAME_MAX];
  char ofile[FILENAME_MAX];
  //  char pfile[FILENAME_MAX];
  char incalfilelist[FILENAME_MAX];

  char rad_landtables_file[FILENAME_MAX];
  char rad_landcorr_file[FILENAME_MAX];
  char rad_gainice_file[FILENAME_MAX];
  char rad_tausq_file[FILENAME_MAX];

  char rad_apc_file[FILENAME_MAX];
  char rad_sssalgo_file[FILENAME_MAX];
  char rad_galwind_file[FILENAME_MAX];
  char rad_sun_file[FILENAME_MAX];
  char rad_sunbak_file[FILENAME_MAX];
  char rad_oceanrefl_file[FILENAME_MAX];
  char rad_dtb_wind_file[FILENAME_MAX];

  char coeff_loss_file[FILENAME_MAX];
  char coeff_nl_file[FILENAME_MAX];

  char emiss_coeff_harm_file[FILENAME_MAX];
  char scat_coeff_harm_file[FILENAME_MAX];
  char dtbw_win_sigma_file[FILENAME_MAX];
  char dtbw_win_wav_file[FILENAME_MAX];
  char dtb_emiss_wspd_file[FILENAME_MAX];

  char wind_errortab_file[FILENAME_MAX];

  char yancfile1[FILENAME_MAX];
  char yancfile2[FILENAME_MAX];
  char yancfile3[FILENAME_MAX];

  char l2prod[PRODSTRLEN];
  int8_t l2activeprodflag[MAXNAQPROD];
  double matchup_time;
  float matchup_lat;
  float matchup_lon;
  float matchup_delta_lat;
  float matchup_delta_lon;
  float matchup_min_dist;
  char matchup_lim_file[FILENAME_MAX];
  char scatter_basename[FILENAME_MAX];
  char pversion[64];
  char sss_algorithm[64];

  bool browse;
  bool iopt_rfi;
  //  bool iopt_drift;
  bool iopt_l1b;
  bool iopt_nosa1;
  bool iopt_zero;
  //  bool use_scatwind; // Remove JMG 02/24/12

  double rpy_adj[3];
  float c_deltaTND[3][2][3];
  float ta_ocean_nom[NUMBER_OF_BEAMS*RADIOMETER_POLARIZATIONS];
  char ta_nominal_files[FILENAME_MAX];

  char radflaglimitsfile[FILENAME_MAX];

  char xrayfile1[FILENAME_MAX];
  char xrayfile2[FILENAME_MAX];

  char rad_offset_corr_file[FILENAME_MAX];
  char climate_sal_file[FILENAME_MAX];

  char rad_dta_gal_file[FILENAME_MAX];

  char anomaly_status[FILENAME_MAX];
  char pointing_anomaly_file[FILENAME_MAX];

  char rfi_mask_file[FILENAME_MAX];

  char dI_U_coeff_file[FILENAME_MAX];

  char l2_unc_maps_file[FILENAME_MAX];

  double sss_sst_adj_parms[3];

  char l_acc_corr_files[FILENAME_MAX];

  char static_bias_adj_file[FILENAME_MAX];

  float instrument_gain_corr[2][NUMBER_OF_BEAMS][2];
} instr;


//      Instrument Temperature Structure
typedef struct instrument_temp_struct {
  float T1 [MAXCYC][NUMBER_OF_BEAMS];
  float T2A[MAXCYC][NUMBER_OF_BEAMS];
  float T2B[MAXCYC][NUMBER_OF_BEAMS];

  float T3[MAXCYC][NUMBER_OF_BEAMS][2];
  float T4[MAXCYC][NUMBER_OF_BEAMS][2];
  float T5[MAXCYC][NUMBER_OF_BEAMS][2];
  float T0[MAXCYC][NUMBER_OF_BEAMS][2];

  float TCND_P[MAXCYC][NUMBER_OF_BEAMS][2];
  float TND_P [MAXCYC][NUMBER_OF_BEAMS][2];
  float TCND  [MAXCYC][NUMBER_OF_BEAMS][2];

  float TND          [MAXCYC][NUMBER_OF_BEAMS][2];
  float TND_offset   [MAXCYC][NUMBER_OF_BEAMS][2];

  float TDL          [MAXCYC][NUMBER_OF_BEAMS][2];
  float TDL_offset   [MAXCYC][NUMBER_OF_BEAMS][2];

  float Tdet[MAXCYC][NUMBER_OF_BEAMS][4];

} itemp;


typedef struct static_data_struct {
  float apc_matrix[NUMBER_OF_BEAMS][3][3];
  float apc_inverse[NUMBER_OF_BEAMS][3][3];
  float sss_coef[251][451][4];
  double time_sun[NZANG];
  float eia_sun[NUMBER_OF_BEAMS][NZANG];
  float tasun_dir_tab[NUMBER_OF_BEAMS][3][NZANG][NOMEGA];
  float tasun_ref_tab[NUMBER_OF_BEAMS][3][NZANG][NOMEGA];
  float tasun_bak_tab[NUMBER_OF_BEAMS][3][26][161];
  double time_galaxy[NZANG];
  float eia_galaxy[NUMBER_OF_BEAMS][NZANG];
  float tagal_dir_tab[NUMBER_OF_BEAMS][3][NZANG][NOMEGA];
  float tagal_ref_tab[5][NUMBER_OF_BEAMS][3][NZANG][NOMEGA];
  double dtagal_ref_tab[NZANG][NOMEGA][NUMBER_OF_BEAMS][2];
  float fpt_lnd[NUMBER_OF_BEAMS][NZANG_LND][NLON_LND];
  float frc_lnd[NUMBER_OF_BEAMS][NZANG_LND][NLON_LND];
  uint16_t landcorr[NUMBER_OF_BEAMS][2][12][1440][1440];
  float gain_ice[NUMBER_OF_BEAMS][151];
  uint16_t itau[91][12][180][360];
  double acoef[5][3][2][3];
  double bcoef[5][3][4][3];
  double arr_dtbw[NBIN_W][NBIN_W][NUMBER_OF_BEAMS][2][3];
  double arr_dtbwx[NBIN_WAVX][NBIN_WINX][NUMBER_OF_BEAMS][2];
  float  dsigma[NUMBER_OF_BEAMS][2];
  uint8_t iflag_dtbw[NBIN_W][NBIN_W][NUMBER_OF_BEAMS][2][3];
  uint8_t iflag_winx_wavx[NBIN_WAVX][NBIN_WINX][NUMBER_OF_BEAMS][2];
  double wspd_max_a[NUMBER_OF_BEAMS][2][3];
  double wspd_max_b[NUMBER_OF_BEAMS][4][3];
  float coeffs_dI_U[NUMBER_OF_BEAMS][4];

  float sst_step, sst0, sst1;
  double yarr[NSST_ARR][NUMBER_OF_BEAMS][2];
  double T_STITCH[NUMBER_OF_BEAMS][2];
  float W_STITCH;

  float estimated_error_array[NWIN_SC][7];
  uint16_t climate_sal[12][180][360];
  uint8_t rfi_mask[90][180][2];
  uint8_t emiss_sst_sss;
} static_data;



//      RFI structure
typedef struct rfi_struct {
  uint8_t iflag_rfi[MAXCYC][NUMBER_OF_BEAMS][RADIOMETER_POLARIZATIONS]
  [RADIOMETER_SUBCYCLES][6];
  int32_t num_rfi[MAXCYC][NUMBER_OF_BEAMS][RADIOMETER_POLARIZATIONS];
  uint8_t iflag_glitch[MAXCYC][NUMBER_OF_BEAMS][RADIOMETER_POLARIZATIONS];
  uint8_t iflag_rfi_CND[MAXCYC][NUMBER_OF_BEAMS][RADIOMETER_POLARIZATIONS]
  [RADIOMETER_SUBCYCLES];

  int32_t wm[360][181];
  int32_t wd[360][181];
  float   tm[360][181];
  float   td[360][181];
  float   td_a[360][181];
  float   td_d[360][181];

  int32_t wm_c[360][181];
  int32_t wd_c[360][181];
  float   tm_c[360][181];
  float   td_c[360][181];

  float   stdta_rad[NUMBER_OF_BEAMS][RADIOMETER_POLARIZATIONS][360][181];
  float   stdta_cnd[NUMBER_OF_BEAMS][RADIOMETER_POLARIZATIONS];

  int32_t idata;

} rfi_data;


//      PLC structure
typedef struct plc_struct {
  //       ref. temp for FE losses and ND
  float Tref_3;
  float Tref_4;
  float Tref_5;

  //       ref. temp for FE losses and ND
  float Tref_CND;
  float Tref_ND;

  float Tref_ND_offset;
  float Tref_DL_offset;

  //       noise diodes
  float TCND_0[2][NUMBER_OF_BEAMS];
  float dTCND_dT4[2][NUMBER_OF_BEAMS];
  float dTCND_dT[2][NUMBER_OF_BEAMS];

  float TND_0[2][NUMBER_OF_BEAMS];
  float DTND_DT[2][NUMBER_OF_BEAMS];
  float TND_offset_0[2][NUMBER_OF_BEAMS];
  float DTND_offset_DT[2][NUMBER_OF_BEAMS];
  float TDL_offset_0[2][NUMBER_OF_BEAMS];
  float dTDL_offset_DT[2][NUMBER_OF_BEAMS];

  float CL1[2][NUMBER_OF_BEAMS];
  float CL2A[2][NUMBER_OF_BEAMS];
  float CL2B[2][NUMBER_OF_BEAMS];
  float CL3[2][NUMBER_OF_BEAMS];
  float CL4[2][NUMBER_OF_BEAMS];
  float dL4_dT4[2][NUMBER_OF_BEAMS];
  float CL5[2][NUMBER_OF_BEAMS];
  float dL5_dT5[2][NUMBER_OF_BEAMS];
  float CLMM[2][NUMBER_OF_BEAMS];

  //       non-linearity correction coefficients
  float Dnl20[4][NUMBER_OF_BEAMS];
  float Dnl21[4][NUMBER_OF_BEAMS];
  float Dnl22[4][NUMBER_OF_BEAMS];
  float Dnl30[4][NUMBER_OF_BEAMS];
  float Dnl31[4][NUMBER_OF_BEAMS];
  float Dnl32[4][NUMBER_OF_BEAMS];
  float Tref_nl[4][NUMBER_OF_BEAMS];
} plc_data;


// subroutine prototypes

int getGeoNavSun( Hdf::hdf5_Aquarius *l1afile, double *rpy_adj, 
		  float *cellon, float *cellat, float *celtht, float *celphi,
		  float *suntht, float *sunphi, float *sunglt, float *moonglt,
		  float *glxlon, float *glxlat,
		  double *zang, double *sun_zenith, 
		  double *sclon, double *sclat, double *scalt,  
		  float *cellonfoot, float *cellatfoot,
		  double *sund, double *sunr, double *moond, double *moonr, 
		  double *bore_sight,
		  double *Pos, double *Vel, double *rpy);

int getHKT( Hdf::hdf5_Aquarius *l1afile, uint8_t *acsMode);

int parseInput( int argc, char* argv[], instr*, uint32_t*, string*);
int initOptions( clo_optionList_t* list, instr *l2genInput);

int initialize_static_data( instr *l2genInput, static_data *staticData);

int getDeflectionRatio( char* incalfilelist, double granStart, 
			float* bestfitDR);

int comp_c_deltaTND( char* incalfilelist, float *deltaTaDR, float *c_delta);
int comp_c_deltaTND( char* incalfilelist, float *c_delta);


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
		 float *tagal_ref_GO, float *dtagal_ref);

int l_acc_corr( string filename, float *L_acc_raw);

extern "C" void count_to_ta_(char *, char *,
			     int32_t *, float *,  itemp *, gainoff *, 
			     rfi_data *, 
			     float *, float *, float *, float *,
			     float *, float *, double *, float *, float *,
			     float *, float *, float *, float *);

extern "C" void get_prelaunch_calibration_coeffs_(char *, char *, plc_data *);

extern "C" void  get_zone_temperatures_( int32_t *, float *, itemp *);

// extern "C" void fe_loss_fac_( int32_t *, plc_data *, itemp *, float *);

extern "C" void aq_rad_rtr_(unsigned short *, unsigned short *, 
			    unsigned short *, float *, float *, float *);

extern "C" void geolocation_(double *, double *, double *, double *, double *, 
			     float *, float *, float *, float *, float *, 
			     float *, float *, float *, float *, float *,
			     double *, double *, double *, double *, double *,
			     float *, float *,
			     double *, double *, double *, double *, double *);

extern "C" void fd_date_2000_(double *, double *, int *, int *, int *, int*,
			   double *); 

extern "C" void fd_ta_sun_(double *, double *, double *, float *, float *, 
			   float *, float *);

extern "C" void fd_ta_moon_(int32_t *, double *, double *, float *, float *, float *); 

extern "C" void fd_ta_galaxy_( int32_t *, double *, double *, float *, 
			       float *, float *, float *, double *, 
			       float *, float *, float *);

extern "C" void fd_dta_galaxy_( int32_t *, double *, double *, double *, 
				double *, float *);
extern "C" void fd_dta_symm_( int32_t *, double *, double *, double *, 
			      double *, float *);

extern "C" void fd_water_refl_( float *, float *, float *, float *, float *, 
				uint16_t *, uint16_t *, float *, float *);

extern "C" void fd_water_refl_exact_( float *, float *, float *, float *, 
				      float *, static_data *, float *);

extern "C" void	fd_dtb_roughness_( int32_t *, int32_t *, int32_t *, 
				   float *, float *, float *, float *,
				   float *, float*,
				   static_data *, int32_t *, float *,
                                   float *);

extern "C" void find_refl_tot_( int *, 
				float *, float *, float *, 
				float *, float *, float *, float *,
				float *, float *, float *, 
				static_data *, float *, float *);

extern "C" void fd_sun_backscatter_( int32_t *, float *, float *, double *, 
				     float *, float *);

extern "C" void adjust_tagal_ref_( int32_t *, float *, float *, float *, 
				   float *, float *, 
				   float *, float *, float *, float *);

extern "C" void fd_sss_sigma0_( float *, float *, float (*)[2],
				float *, int32_t *, float *);

extern "C" void vh_to_stokes_(float *);
extern "C" void stokes_to_vh_(float *);

extern "C" void get_ancillary_data_(int *, int *, int*,
				    float *, float *, float *,
				    float *, float *, float *, 
				    float *, float *, float *,
				    float *, float *, float *, 
				    float *, float*, float *,
				    int32_t *, double *, double *,
				    float *, float *, uint16_t *,
				    float *, float *, float *, float *,
				    float *, float *, float *,
				    float *, char *, char *, char *);

extern "C" void invert_3by3_( double (*)[3], double (*)[3], double *);

extern "C" void find_di_u_( int32_t *, float *, float *, float *);

extern "C" void interpolate_l2_error_(char *l2_unc_maps_file, 
                                      int32_t *idayjl, 
                                      float *clat, float *con,
                                      int32_t *iasc,
                                      float *err_ran, float *err_sys);

#ifdef IFORT_MODULE_NAME_MANGLING
   #define FD_SSS salinity_module_mp_fd_sss_
   #define FD_WSPD wind_speed_retrieval_module_mp_fd_wspd_
   #define MEISSNER_WENTZ_SALINITY salinity_module_mp_fdem0_meissner_wentz_salinity_
#else
   #define FD_SSS __salinity_module_MOD_fd_sss
   #define FD_WSPD __wind_speed_retrieval_module_MOD_fd_wspd
   #define MEISSNER_WENTZ_SALINITY __salinity_module_MOD_fdem0_meissner_wentz_salinity
#endif

extern "C" void FD_SSS( int32_t *, float *, float *, 
			float (*)[2], float *, float *, 
			int32_t *);

extern "C" void FD_WSPD( int *ibeam, 
			 float *phir, 
			 float *winspd, 
			 float *dew,
			 float *xsigma0_vv,  
			 float *xsigma0_hh,
			 float *ysst,
			 float *ysss,
			 static_data *staticData,
			 float *winspd_hh,  
			 float *chisq_scat, 
			 int *iflag_wspd_scat);

extern "C" void MEISSNER_WENTZ_SALINITY( float *freq_aq, float *thtadj, 
					 float *sst, 
					 float *sss_clim0, float *em0);


extern "C" void fd_sss_clm_( int32_t *isecyr,  float *clat, float *clon, 
			     static_data *staticData, float *sss_clim);

// Scatterometer
extern "C" void initialize_constants_();
extern "C" void l1b_read_limits_file_( char *);
extern "C" void l1b_read_params_( char *);
extern "C" void read_jpl_cut_antenna_pattern_relative_( const char *, 
							double *, double *);

int32_t prodOffset[MAXNAQPROD];
int8_t prodMultiplicity[MAXNAQPROD];
float *prodPtr[MAXNAQPROD];


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

