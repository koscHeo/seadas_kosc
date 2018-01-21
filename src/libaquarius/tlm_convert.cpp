#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#include "hdf5_Aquarius.h"
#include "tlm_convert.h"
#include "passthebuck.h"

#define FAIL -1

float crs_dspl_temp_10( uint16_t tel_16) {
  // Rule: CRS_DSPL_TEMP_10
  float b2f_10 = tel_16 * 250000 / (1 << 10);
  float tmpf = b2f_10 * CNST_CRS_DSPL_G + CNST_CRS_DSPL_0;
  tmpf = (tmpf * CNST_CRS_LRI) / (CNST_CRS_LRI - tmpf);
  return (1. / 
	  (CNST_NTCA + CNST_NTCB * logf(tmpf) + CNST_NTCC * pow(logf(tmpf),3)) 
	  - 273.15);
}

float crs_dspl_temp_12( uint16_t tel_16) {
  // Rule: CRS_DSPL_TEMP_12
  float b2f_12 = tel_16 * 250000 / (1 << 12);
  float tmpf = b2f_12 * CNST_CRS_DSPL_G + CNST_CRS_DSPL_0;
  tmpf = (tmpf * CNST_CRS_LRI) / (CNST_CRS_LRI - tmpf);
  return (1. / 
	  (CNST_NTCA + CNST_NTCB * logf(tmpf) + CNST_NTCC * pow(logf(tmpf),3)) 
	  - 273.15);
}

float pls_dspl_ohms_12( uint16_t tel_16) {
  // Rule: PLS_DSPL_OHMS_12
  float b2f_12 = tel_16 * 250000 / (1 << 12);
  return (b2f_12 * CNST_PLS_DSPL_G + CNST_PLS_DSPL_0);
}

float mns_dspl_ohms_12( uint16_t tel_16) {
  // Rule: MNS_DSPL_OHMS_12
  float b2f_12 = tel_16 * 250000 / (1 << 12);
  return (b2f_12 * CNST_MNS_DSPL_G + CNST_MNS_DSPL_0);
}

float crs_dspl_ohms_12( uint16_t tel_16) {
  // Rule: CRS_DSPL_OHMS_12
  float b2f_12 = tel_16 * 250000 / (1 << 12);
  return (b2f_12 * CNST_CRS_DSPL_G + CNST_CRS_DSPL_0);
}

float pls_dspl_temp_12( uint16_t tel_16) {
  // Rule: PLS_DSPL_TEMP_12
  float b2f_12 = tel_16 * 250000 / (1 << 12);
  float tmpf = b2f_12 * CNST_PLS_DSPL_G + CNST_PLS_DSPL_0;
  return (CNST_PRT_DSPL_A * pow(tmpf, 2) + CNST_PRT_DSPL_B * tmpf + 
	  CNST_PRT_DSPL_C);
}

float mns_dspl_temp_12( uint16_t tel_16) {
  // Rule: MNS_DSPL_TEMP_12
  float b2f_12 = tel_16 * 250000 / (1 << 12);
  float tmpf = b2f_12 * CNST_MNS_DSPL_G + CNST_MNS_DSPL_0;
  return (CNST_PRT_DSPL_A * pow(tmpf, 2) + CNST_PRT_DSPL_B * tmpf + 
	  CNST_PRT_DSPL_C);
}

float fnep_ohms( float tel_16, float *extram) {

  float gfinep = (CNST_RR1 - CNST_RR2) / (extram[0] - extram[2]);

  float ofinep = CNST_RR1 - gfinep * extram[0];
	
  float b2f_12 = tel_16 * 250000 / (1 << 12);

  return b2f_12 * gfinep + ofinep;
}

float fnen_ohms( float tel_16, float *extram) {

  float gfinen = (CNST_RR1 - CNST_RR2) / (extram[1] - extram[3]);

  float ofinen = CNST_RR1 - gfinen * extram[1];
	
  float b2f_12 = tel_16 * 250000 / (1 << 12);

  return b2f_12 * gfinen + ofinen;
}

float finep_temp_12( uint16_t tel_16, float *extram,
		     float const_a, float const_b, float const_r0) {

  float ohms = fnep_ohms( tel_16, extram);
  float tmpa = (sqrt(pow(const_a, 2) - 4 * const_b * (1 - ohms / const_r0)) 
		- const_a) / (2 * const_b);

  return tmpa + 0.045 * (tmpa/100) * (tmpa/100 - 1) * 
    (tmpa/419.58 - 1) * (tmpa/630.74 - 1);
}

float finen_temp_12( uint16_t tel_16, float *extram,
		     float const_a, float const_b, float const_r0) {

  float ohms = fnen_ohms( tel_16, extram);
  float tmpa = (sqrt(pow(const_a, 2) - 4 * const_b * (1 - ohms / const_r0)) 
		- const_a) / (2 * const_b);

  return tmpa + 0.045 * (tmpa/100) * (tmpa/100 - 1) * 
    (tmpa/419.58 - 1) * (tmpa/630.74 - 1);
}

float coarse_temp_12( uint16_t tel_16, float *extram) {

  float gcrs = (CNST_RR1 - CNST_RR3) / (extram[5] - extram[4]);

  float ocrs = CNST_RR1 - gcrs * extram[5];

  float b2f_12 = tel_16 * 250000 / (1 << 12);
  float crs_ohms = b2f_12 * gcrs + ocrs;

  float ntc_ohms = (crs_ohms * CNST_CRS_LRI) / (CNST_CRS_LRI - crs_ohms);

  return (1. / (CNST_NTCA + CNST_NTCB * logf(ntc_ohms) + CNST_NTCC * 
		pow(logf(ntc_ohms), 3)) - 273.15);
}

float coarse_temp_10( uint16_t tel_16, float *extram) {

  float gcrs = (CNST_RR1 - CNST_RR3) / (extram[5] - extram[4]);

  float ocrs = CNST_RR1 - gcrs * extram[5];

  float b2f_10 = tel_16 * 250000 / (1 << 10);
  float crs_ohms = b2f_10 * gcrs + ocrs;

  float ntc_ohms = (crs_ohms * CNST_CRS_LRI) / (CNST_CRS_LRI - crs_ohms);

  return (1. / (CNST_NTCA + CNST_NTCB * logf(ntc_ohms) + CNST_NTCC * 
		pow(logf(ntc_ohms), 3)) - 273.15);
}

float ie_volt( uint8_t tel_8) {
  // Rule: IE_VOLT
  return tel_8 * 256. * CNST_VLT_GAIN + CNST_VLT_OFFSET;
}

float resist_16bit_norm_disp( uint16_t tel_16) {
  return tel_16 * CNST_GAIN_NORM_NOM + CNST_OFFSET_NORM_NOM;
}

float resist_16bit_ext_disp( uint16_t tel_16) {
  return tel_16 * CNST_GAIN_EXT_NOM + CNST_OFFSET_EXT_NOM;
}

float temp_8bit_ext_fine( uint8_t tel_8,
			  float const_a, float const_b, float const_r0) {
  // Rule: TEMP_8BIT_EXT
  float resist_8bit_ext = 
    tel_8 * 256.* CNST_GAIN_EXT_NOM + CNST_OFFSET_EXT_NOM; 

  return (sqrt(pow(const_a, 2) - 4 * const_b * 
	       (1 - resist_8bit_ext / const_r0)) - const_a) / (2 * const_b);
}

float temp_8bit_norm_fine( uint8_t tel_8,
			   float const_a, float const_b, float const_r0) {
  // Rule: TEMP_8BIT_EXT
  float resist_8bit_norm = 
    tel_8 * 256.* CNST_GAIN_NORM_NOM + CNST_OFFSET_NORM_NOM; 

  return (sqrt(pow(const_a, 2) - 4 * const_b * 
	       (1 - resist_8bit_norm / const_r0)) - const_a) / (2 * const_b);
}

float temp_8bit_ext( uint8_t tel_8) {
  // Rule: TEMP_8BIT_EXT
  float resist_8bit_ext = 
    tel_8 * 256.* CNST_GAIN_EXT_NOM + CNST_OFFSET_EXT_NOM; 

  return CNST_POLY0 + (CNST_POLY1 + (CNST_POLY2 + CNST_POLY3*resist_8bit_ext)
		       * resist_8bit_ext) * resist_8bit_ext;
}

float temp_8bit_rad6k( uint8_t tel_8) {
  // Rule: TEMP_8BIT_RAD6K
  float resist_8bit_ext = 
    tel_8 * 256.* CNST_GAIN_EXT_NOM + CNST_OFFSET_EXT_NOM;

  return -64.8 + (-0.7885 + 1.84E-3*resist_8bit_ext) * resist_8bit_ext;
}

float temp_8bit_norm( uint8_t tel_8) {
  // Rule: TEMP_8BIT_NORM
  float resist_8bit_norm = 
    tel_8 * 256.* CNST_GAIN_NORM_NOM + CNST_OFFSET_NORM_NOM;

  return CNST_POLY0 + (CNST_POLY1 + (CNST_POLY2 + CNST_POLY3*resist_8bit_norm)
		       * resist_8bit_norm) * resist_8bit_norm;
}

float temp_apdu_8bit( uint8_t tel_8, float const_ap) {
  float resist_8bit_ext = 
    tel_8 * 256.* CNST_GAIN_EXT_NOM + CNST_OFFSET_EXT_NOM;
    
  return temp_resist( CNST_AP_A, CNST_AP_B, const_ap, resist_8bit_ext);
}


float temp_8bit_norm_abr0( uint8_t tel_8, 
			   float const_a, float const_b, float const_r0) {
  float resist_8bit_norm = 
    tel_8 * 256.* CNST_GAIN_NORM_NOM + CNST_OFFSET_NORM_NOM;

  return temp_resist( const_a, const_b, const_r0, resist_8bit_norm);
}

float temp_8bit_ext_abr0( uint8_t tel_8, 
			  float const_a, float const_b, float const_r0) {
  float resist_8bit_ext = 
    tel_8 * 256.* CNST_GAIN_EXT_NOM + CNST_OFFSET_EXT_NOM;

  return temp_resist( const_a, const_b, const_r0, resist_8bit_ext);
}


float temp_8bit_norm_scat( uint8_t tel_8, 
			   float const_a, float const_b, float const_r0) {
  float resist_8bit_norm = 
    tel_8 * 256.* CNST_GAIN_NORM_NOM + CNST_OFFSET_NORM_NOM;

  return temp_resist( const_a, const_b, const_r0, resist_8bit_norm);
}


float temp_16bit_norm_disp( uint16_t tel_16, float gain, float offset,
			    float const_a, float const_b, float const_r0) {
	
  float resist_16bit_norm = tel_16 * gain + offset;

  return temp_resist( const_a, const_b, const_r0, resist_16bit_norm);
}

float temp_16bit_ext_disp( uint16_t tel_16, float gain, float offset,
			   float const_a, float const_b, float const_r0,
			   float const_h_r) {

  float resist_16bit_ext = tel_16 * gain + offset;

  return temp_resist( const_a, const_b, const_r0, resist_16bit_ext-const_h_r);
}


float temp_16bit_norm_fine( uint16_t tel_16, uint16_t *extram_et,
			    float const_a, float const_b, float const_r0) {

  float gain = 1.0 * (CNST_CALRES2_NOM-CNST_CALRES3_NOM) / (extram_et[1] - extram_et[2]);
  float offset = CNST_CALRES2_NOM - extram_et[1] * gain;

  float resist_16bit_norm = tel_16 * gain + offset;

  return temp_resist( const_a, const_b, const_r0, resist_16bit_norm);
}


// extram = IE_ICDS_TLM_CAL_RES[1-4]_DN

float temp_16bit_ext_fine( uint16_t tel_16, uint16_t *extram_et,
			   float const_a, float const_b, float const_r0,
			   float const_h_r) {

  float gain = 1.0 * (CNST_CALRES1_NOM-CNST_CALRES4_NOM) / (extram_et[0] - extram_et[3]);
  float offset = CNST_CALRES1_NOM - extram_et[0] * gain;

  float resist_16bit_ext = tel_16 * gain + offset;

  return temp_resist( const_a, const_b, const_r0, resist_16bit_ext-const_h_r);
}


float temp_resist( float a, float b, float r0, float resist) {
  float A = 4 * b * (1 - (resist/r0)) / (a * a);

  //  if ( A < 0.2 && -A > -0.2)
  //  return -(a * (A/2 + A*A/8 + A*A*A/16)) / (2 * b);
  //else
  return (a * (sqrt(1 - A) - 1)) / (2 * b);
}


float dwnlk_fltg_pt( uint16_t tel_16) {
  int8_t presign = (tel_16 & 0x8000) >> 15;
  int8_t sign = 1 - 2 * presign;
  int16_t preexp = tel_16 - (presign << 15);
  int16_t ex_p = (preexp >> 10) - 15;
  float premant = ( (float) preexp) / (1 << 10);
  float mant = 1 + premant - ex_p - 15;

  return (sign * mant * pow(2.0, ex_p));
}

uint16_t in32_out16( uint32_t i32, uint16_t nbits, uint16_t mask) {
  if        ( nbits == 8  && mask == 0xf) {
    return (i32 & mask) * 16 + (i32 >> (16+12));
  } else if ( nbits == 10 && mask == 0x3) {
    return (i32 & mask) * 256 + (i32 >> (16+8));
  } else if ( nbits == 10 && mask == 0xf) {
    return (i32 & mask) * 64 + (i32 >> (16+10));
  } else if ( nbits == 10 && mask == 0x3f) {
    return (i32 & mask) * 16 + (i32 >> (16+12));
  } else if ( nbits == 10 && mask == 0xff) {
    return (i32 & mask) * 4 + (i32 >> (16+14));
  } else if ( nbits == 12 && mask == 0xf) {
    return (i32 & mask) * 256 + (i32 >> (16+8));
  } else if ( nbits == 12 && mask == 0xff) {
    return (i32 & mask) * 16 + (i32 >> (16+12));
  } else if ( nbits == 16 && mask == 0xff) {
    return (i32 & mask) * 256 + (i32 >> (16+8));
  }

  // uint16_t x = (uint16_t) (alog10(mask+1)/alog10(2))
  // uint16_t b = x + 16 - nbits
  // uint16_t a = (uint16_t) pow(2, nbits-x)
  // return (i32 & mask) * a + (i32 >> (16+b));
  return -1;
}


