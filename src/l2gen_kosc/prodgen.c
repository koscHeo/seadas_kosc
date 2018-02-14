/* =========================================================== */
/* Module prodgen.c                                            */
/*                                                             */
/* Function prodgen returns a pointer to the geophysical       */
/* product associated with the input product catalog entry and */
/* the input level-2 record.  The specific product will be     */
/* extracted from the level-2 record or computed using the     */
/* information in the level-2 record and knowledge of the      */
/* required algorithm to be called.  The output is stored in   */
/* a static array, and a pointer is returned to the caller.    */ 
/*                                                             */ 
/* Written By:                                                 */
/*     Bryan A. Franz, NASA/OBPG, August 2008.                 */
/* =========================================================== */

#include "l12_proto.h"

#include <stdint.h>
#include <inttypes.h>

static int32    numScans;
static int32    numPixels;
static int32    numBands;
static int32    numBandsIR;


/* ----------------------------------------------------------- */
/* extract_band() - extracts a product from a BIL array        */
/* ----------------------------------------------------------- */
VOIDP extract_band(float *fbuf, l2prodstr *p, int32 nbands)
{
    static  float32 *fbuf2 = NULL;

    int32 band = MIN(MAX(p->prod_ix,0),nbands-1);
    static int32 npix =0;
    int32 ip;

    if (fbuf2 == NULL) {
        npix = p->dim[1];
        if ( (fbuf2 = calloc(npix,sizeof(float32))) == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: Unable to allocate buffer space.\n",
                    __FILE__,__LINE__);
            exit(1);
        }
    } else if(npix != p->dim[1]) {
        npix = p->dim[1];
        free(fbuf2);
        // this is put here for l3gen which varies the number of pixels on each line

        if ((fbuf2 = (float32 *) calloc(npix, sizeof(float32)))
                == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: Unable to allocate buffer space.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
    }

    for (ip=0; ip<npix && ip<numPixels; ip++) {
        fbuf2[ip] = fbuf[ip*nbands + band];
    }

    return ((VOIDP) fbuf2);
}

/**
 * loop through a float array applying a multiplier
 *
 * @param in input array
 * @param out output array
 * @param count size of array
 * @param multiplier factor to apply
 */
void applyMultiplier(float* in, float* out, int count, float multiplier) {
    int i;
    for(i=0; i<count; i++) {
        if(in[i] == BAD_FLT)
            out[i] = BAD_FLT;
        else
            out[i] = in[i] * multiplier;
     }
}

/* ----------------------------------------------------------- */
/* prodgen() - returns pointer the the requested product       */
/* ----------------------------------------------------------- */
VOIDP prodgen( l2prodstr *p, l2str *l2rec, initstr *initrec) 
{
    static int firstCall = 1;

    static float32 *fbuf = NULL; 

    float *tmpFloat;
    VOIDP pbuf;

    uint32_t mask;


    if (firstCall) {
        firstCall = 0;
        numPixels  = l2rec->npix;
        numScans   = l2rec->nscans;
        numBands   = l2rec->nbands;
        numBandsIR = l2rec->nbandsir;
        if ( (fbuf = (float32 *) calloc(numPixels,sizeof(float32))) == NULL) {
            fprintf(stderr,
            "-E- %s line %d: Unable to allocate buffer space.\n",
            __FILE__,__LINE__);
            exit(1);
        }
    } else if(l2rec->npix != numPixels) {

        // this is put here for l3gen which varies the number of pixels on each line
        if (l2rec->npix > numPixels) {

            free(fbuf);
            if ((fbuf = (float32 *) calloc(l2rec->npix, sizeof(float32)))
                    == NULL) {
                fprintf(stderr,
                        "-E- %s line %d: Unable to allocate buffer space.\n",
                        __FILE__, __LINE__);
                exit(1);
            }

        }
        numPixels = l2rec->npix;
    }


    switch (p->cat_ix) {

      //
      // Band-dependent, precomputed products
      //
      case CAT_nLw:
          applyMultiplier(extract_band(l2rec->nLw,p,numBands), fbuf, numPixels, 10.0);
          pbuf = fbuf;
        break;
      case CAT_nLw_unc:
        applyMultiplier(extract_band(l2rec->nLw_unc,p,numBands), fbuf, numPixels, 10.0);
        pbuf = fbuf;
        break;
      case CAT_Lw:
        applyMultiplier(extract_band(l2rec->Lw,p,numBands), fbuf, numPixels, 10.0);
        pbuf = fbuf;
        break;      
      case CAT_Rrs:
        pbuf = extract_band(l2rec->Rrs,p,numBands);
        break;      
      case CAT_Rrs_unc:
        pbuf = extract_band(l2rec->Rrs_unc,p,numBands);
        break;      
      case CAT_Taua:
        pbuf = extract_band(l2rec->taua,p,numBands);
        break;      
      case CAT_Lr:
        applyMultiplier(extract_band(l2rec->Lr,p,numBands), fbuf, numPixels, 10.0);
        pbuf = fbuf;
        break;
      case CAT_L_q:
        applyMultiplier(extract_band(l2rec->L_q,p,numBands), fbuf, numPixels, 10.0);
        pbuf = fbuf;
        break;
      case CAT_L_u:
        applyMultiplier(extract_band(l2rec->L_u,p,numBands), fbuf, numPixels, 10.0);
        pbuf = fbuf;
        break;
      case CAT_polcor:
        pbuf = extract_band(l2rec->polcor,p,numBands);
        break;
      case CAT_La:
        applyMultiplier(extract_band(l2rec->La,p,numBands), fbuf, numPixels, 10.0);
        pbuf = fbuf;
        break;
      case CAT_TLg:
        applyMultiplier(extract_band(l2rec->TLg,p,numBands), fbuf, numPixels, 10.0);
        pbuf = fbuf;
        break;
      case CAT_tLf:
        applyMultiplier(extract_band(l2rec->tLf,p,numBands), fbuf, numPixels, 10.0);
        pbuf = fbuf;
        break;
      case CAT_brdf:
        pbuf = extract_band(l2rec->brdf,p,numBands);
        break;
      case CAT_Lt:
          if(p->prod_ix < numBands) {
              applyMultiplier(extract_band(l2rec->Lt,p,numBands), fbuf, numPixels, 10.0);
              pbuf = fbuf;
          } else {
              // first subtract the num of visible bands
              p->prod_ix -= numBands;
              applyMultiplier(extract_band(l2rec->Ltir,p,NBANDSIR), fbuf, numPixels, 10.0);
              pbuf = fbuf;
              p->prod_ix += numBands;
          }
          break;
      case CAT_Lt_unc:
        applyMultiplier(extract_band(l2rec->Lt_unc,p,numBands), fbuf, numPixels, 10.0);
        pbuf = fbuf;
        break;
      case CAT_BT:
        // first subtract the num of visible bands
        p->prod_ix -= numBands;
        pbuf = extract_band(l2rec->Bt,p,NBANDSIR);
        p->prod_ix += numBands;
        break;
      case CAT_rhos:
        pbuf = extract_band(l2rec->rhos,p,numBands);
        break;
      case CAT_nw:
        pbuf = extract_band(l2rec->sw_n,p,numBands);
        break;
      case CAT_aw:
        pbuf = extract_band(l2rec->sw_a,p,numBands);
        break;
      case CAT_bbw:
        pbuf = extract_band(l2rec->sw_bb,p,numBands);
        break;
      case CAT_bbws:
        get_bbws(l2rec,p,fbuf);
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_a:
        pbuf = extract_band(l2rec->a,p,numBands);
        break;
      case CAT_bb:
        pbuf = extract_band(l2rec->bb,p,numBands);
        break;
      case CAT_t_sol:
        pbuf = extract_band(l2rec->t_sol,p,numBands);
        break;
      case CAT_t_sen:
        pbuf = extract_band(l2rec->t_sen,p,numBands);
        break;
      case CAT_tg_sol:
        pbuf = extract_band(l2rec->tg_sol,p,numBands);
        break;
      case CAT_tg_sen:
        pbuf = extract_band(l2rec->tg_sen,p,numBands);
        break;
      case CAT_t_h2o:
        pbuf = extract_band(l2rec->t_h2o,p,numBands);
        break;
      case CAT_t_o2:
        pbuf = extract_band(l2rec->t_o2,p,numBands);
        break;
      case CAT_dpol:
        pbuf = extract_band(l2rec->dpol,p,numBands);
        break;
      case CAT_BT_39:
        pbuf = extract_band(l2rec->Bt,p,NBANDSIR);
        break;
      case CAT_BT_40:
        pbuf = extract_band(l2rec->Bt,p,NBANDSIR);
        break;
      case CAT_BT_11:
        pbuf = extract_band(l2rec->Bt,p,NBANDSIR);
        break;
      case CAT_BT_12:
        pbuf = extract_band(l2rec->Bt,p,NBANDSIR);
        break;

      //
      // Band independent, precomputed products
      //
      case CAT_epsilon:
        pbuf = (VOIDP) l2rec->eps;
        break;
      case CAT_solz:
        pbuf = (VOIDP) l2rec->solz;
        break;
      case CAT_sola:
        pbuf = (VOIDP) l2rec->sola;
        break;
      case CAT_senz:
        pbuf = (VOIDP) l2rec->senz;
        break;
      case CAT_sena:
        pbuf = (VOIDP) l2rec->sena;
        break;
      case CAT_relaz:
        pbuf = (VOIDP) l2rec->delphi;
        break;
      case CAT_alpha:
        pbuf = (VOIDP) l2rec->alpha;
        break;
      case CAT_scattang:
        pbuf = (VOIDP) l2rec->scattang;
        break;
      case CAT_ozone:
        pbuf = (VOIDP) l2rec->oz;
        break;
      case CAT_no2_tropo:
        pbuf = (VOIDP) l2rec->no2_tropo;
        break;
      case CAT_no2_strat:
        pbuf = (VOIDP) l2rec->no2_strat;
        break;
      case CAT_no2_frac:
        pbuf = (VOIDP) l2rec->no2_frac;
        break;
      case CAT_windspeed:
        pbuf = (VOIDP) l2rec->ws;
        break;
      case CAT_windangle:
        pbuf = (VOIDP) l2rec->wd;
        break;
      case CAT_zwind:
        pbuf = (VOIDP) l2rec->zw;
        break;
      case CAT_mwind:
        pbuf = (VOIDP) l2rec->mw;
        break;
      case CAT_pressure:
        pbuf = (VOIDP) l2rec->pr;
        break;
      case CAT_water_vapor:
        pbuf = (VOIDP) l2rec->wv;
        break;
      case CAT_humidity:
        pbuf = (VOIDP) l2rec->rh;
        break;
      case CAT_height:
        pbuf = (VOIDP) l2rec->height;
        break;
      case CAT_sstref:
        pbuf = (VOIDP) l2rec->sstref;
        break;
      case CAT_sssref:
        pbuf = (VOIDP) l2rec->sssref;
        break;
      case CAT_glint_coef:
        pbuf = (VOIDP) l2rec->glint_coef;
        break;
      case CAT_cloud_albedo:
        pbuf = (VOIDP) l2rec->cloud_albedo;
        break;
      case CAT_rho_cirrus:
        pbuf = (VOIDP) l2rec->rho_cirrus;
        break;
      case CAT_aerindex:
        pbuf = (VOIDP) l2rec->aerindex;
        break;
      case CAT_aer_ratio:
	if (p->prod_ix == 1) 
          pbuf = (VOIDP) l2rec->aerratio;
        else
          pbuf = (VOIDP) l2rec->aerratio2;
        break;

      //
      // Integer, precomputed  products
      //
      case CAT_l2_flags:
        pbuf = (VOIDP) l2rec->flags; 
        break;
      case CAT_num_iter:
        pbuf = (VOIDP) l2rec->num_iter; 
        break;
      case CAT_slot:
        pbuf = (VOIDP) l2rec->slot; 
        break;
      case CAT_aer_model:
	switch (p->prod_ix) {
	  case 1: 
            pbuf = l2rec->aermodmin;
	    break;
	  case 2: 
            pbuf = l2rec->aermodmax;
            break;
	  case 3: 
            pbuf = l2rec->aermodmin2;
	    break;
	  case 4: 
            pbuf = l2rec->aermodmax2;
            break;
	}
        break;

      //
      // 1-dimensional, precomputed products
      //
      case CAT_fsol:
        pbuf = (VOIDP) &l2rec->fsol;
        break;
      case CAT_pixnum:
        pbuf = (VOIDP) l2rec->pixnum;
        break;
      case CAT_detnum:
        pbuf = (VOIDP) &l2rec->detnum;
        break;
      case CAT_mside:
        pbuf = (VOIDP) &l2rec->mside;
        break;

      //
      // Derived products
      //
      case CAT_angstrom:
        get_angstrom(l2rec,p->prod_ix,fbuf);
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_ms_epsilon:
        get_ms_epsilon(l2rec,fbuf);
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_Es:
        get_es(l2rec,p->prod_ix,fbuf);
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_rhot:
        get_toa_refl(l2rec,p->prod_ix,fbuf); 
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_rhom:
        get_rho_mumm(l2rec,-1,p->prod_ix,fbuf); 
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_depth:
        get_depth(l2rec, fbuf);
        pbuf =  (VOIDP) fbuf;
        break;
      case CAT_fqy:
        get_fqy(l2rec, fbuf, initrec);
        pbuf =  (VOIDP) fbuf;
        break;
      case CAT_flh:
        get_flh(l2rec, fbuf);
        applyMultiplier(fbuf, fbuf, numPixels, 10.0);
        pbuf =  (VOIDP) fbuf;
        break;
      case CAT_fsat:
        get_fsat(l2rec, fbuf);
        applyMultiplier(fbuf, fbuf, numPixels, 10.0);
        pbuf =  (VOIDP) fbuf;
        break;
      case CAT_ipar:
        get_ipar(l2rec, fbuf, initrec);
        pbuf =  (VOIDP) fbuf;
        break;
      case CAT_BSi:
        get_bsi(l2rec, fbuf);
        pbuf =  (VOIDP) fbuf;
        break;
      case CAT_Kd_mueller:
      case CAT_Kd_532:
      case CAT_Kd_obpg:
      case CAT_Kd_lee:
      case CAT_Kd_morel:
      case CAT_Kd_KD2:
      case CAT_KPAR_morel:
      case CAT_KPAR_lee:
      case CAT_Zhl_morel:
      case CAT_Kd_jamet:
      case CAT_Kd_rhos:
        get_Kd(l2rec, p, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_Zphotic_lee:
      case CAT_Zeu_morel:
      case CAT_Zsd_morel:
      case CAT_Zsd_gbr:
        get_photic_depth(l2rec, p, fbuf); 
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_tindx_morel:
        tindx_morel(l2rec, -1, fbuf); 
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_tindx_shi:
        tindx_shi(l2rec, -1, fbuf); 
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_iCDOM_morel:
      case CAT_pCDOM_morel:
      case CAT_adg_morel:
      case CAT_chl_morel:
      case CAT_chl_cdomcorr_morel:
        get_cdom_morel(l2rec, p, fbuf); 
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_owt:
      case CAT_owtn:
      case CAT_chl_owterr:
        optical_water_type(l2rec, p, fbuf); 
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_owtd:
        optical_water_type(l2rec, p, fbuf); 
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_ndvi:
        get_ndvi(l2rec, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_evi:
        get_evi(l2rec, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_smoke:
        get_smoke(l2rec, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_par:
          get_par(l2rec,fbuf);
          pbuf = (VOIDP) fbuf;
          break;
      case CAT_chl_oc2:
      case CAT_chl_oc3:
      case CAT_chl_oc3c:
      case CAT_chl_oc4:
      case CAT_chl_hu:
      case CAT_chl_oci:
      case CAT_chl_oci2:
      case CAT_chl_cdr:
      case CAT_chl_abi:
        get_chl(l2rec, p->cat_ix, fbuf); 
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_chl_mgiop:
      case CAT_bbp_mgiop:
      case CAT_adg_mgiop:
      case CAT_aph_mgiop:
      case CAT_npix_mgiop:
      case CAT_crat_mgiop:
      case CAT_fitpar_mgiop:
        get_mgiop(l2rec, p, fbuf); 
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_a_pml:
      case CAT_aph_pml:
      case CAT_adg_pml:
      case CAT_bb_pml:
      case CAT_bbp_pml:
        get_pml(l2rec, p, fbuf); 
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_a_qaa:
      case CAT_b_qaa:
      case CAT_c_qaa:
      case CAT_aph_qaa:
      case CAT_adg_qaa:
      case CAT_bb_qaa:
      case CAT_bbp_qaa:
      case CAT_mod_rrs_qaa:
        get_qaa(l2rec, p, fbuf); 
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_flags_qaa:
        pbuf = (VOIDP) get_flags_qaa(l2rec); 
        break;
      case CAT_chl_carder:
      case CAT_a_carder:
      case CAT_bb_carder:
      case CAT_aph_carder:
      case CAT_adg_carder:
      case CAT_bbp_carder:
        get_carder(l2rec, p, fbuf); 
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_chl_carder_emp:
        chl_carder_empirical(l2rec, fbuf); 
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_flags_carder:
        pbuf = (VOIDP) get_flags_carder(l2rec); 
        break;
      case CAT_flags_giop:
        pbuf = (VOIDP) get_flags_giop(l2rec); 
        break;
      case CAT_iter_giop:
        pbuf = (VOIDP) get_iter_giop(l2rec); 
        break;
      case CAT_chl_giop:
      case CAT_a_giop:
      case CAT_bb_giop:
      case CAT_aph_giop:
      case CAT_adg_giop:
      case CAT_bbp_giop:
      case CAT_chl_unc_giop:
      case CAT_a_unc_giop:
      case CAT_bb_unc_giop:
      case CAT_aph_unc_giop:
      case CAT_adg_unc_giop:
      case CAT_bbp_unc_giop:
      case CAT_aphs_giop:
      case CAT_adgs_giop:
      case CAT_bbps_giop:
      case CAT_rrsdiff_giop:
      case CAT_mRrs_giop:
      case CAT_chisqr_giop:
      case CAT_fitpar_giop:
      case CAT_acdom_giop:
      case CAT_anap_giop:
      case CAT_bbph_giop:
      case CAT_bbnap_giop:
      case CAT_acdom_unc_giop:
      case CAT_anap_unc_giop:
      case CAT_bbph_unc_giop:
      case CAT_bbnap_unc_giop:
      case CAT_opt_siop_giop:
        get_giop(l2rec, p, fbuf); 
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_iter_gsm:
        pbuf = (VOIDP) get_iter_gsm(l2rec); 
        break;
      case CAT_chl_gsm:
      case CAT_a_gsm:
      case CAT_bb_gsm:
      case CAT_aph_gsm:
      case CAT_adg_gsm:
      case CAT_bbp_gsm:
        get_gsm(l2rec, p, fbuf); 
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_a_las:
      case CAT_b_las:
      case CAT_c_las:
      case CAT_bb_las:
      case CAT_bbp_las:
      case CAT_bbps_las:
        get_las(l2rec, p, fbuf); 
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_a_niwa:
      case CAT_bb_niwa:
        get_niwa(l2rec, p, fbuf); 
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_flags_niwa:
        pbuf = (VOIDP) get_flags_niwa(l2rec); 
        break;
      case CAT_chl_soa:
      case CAT_adg_soa:
      case CAT_bbp_soa:
      case CAT_pcentcdm_soa:
      case CAT_w0_soa:
      case CAT_v_soa:
        get_soa(l2rec,p->cat_ix, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_chl_sma:
      case CAT_adg_sma:
      case CAT_bbp_sma:
      case CAT_w0_sma:
      case CAT_dom_sma:
        get_sma(l2rec, p->cat_ix, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_poc_stramski_443:
      case CAT_poc_stramski_490:
        get_poc(l2rec, p->cat_ix, fbuf); 
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_ag_412_mlrc:
      case CAT_Sg_275_295_mlrc:
      case CAT_Sg_300_600_mlrc:
        cdom_mannino(l2rec, p->cat_ix, fbuf); 
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_calcite :
      case CAT_calcite_2b :
      case CAT_calcite_3b :
        calcite(l2rec, p, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_Rrs_vc :
	if(l2rec->input->band_shift_opt==1)
		bioOptBandShift(l2rec,p,fbuf);
	else
		virtual_constellation(l2rec, p, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_chl_vc :
        virtual_constellation(l2rec, p, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_sst:
        pbuf = (VOIDP) get_sst(l2rec);
        break;
      case CAT_sst4:
        pbuf = (VOIDP) get_sst4(l2rec);
        break;
      case CAT_sst_triple:
        pbuf = (VOIDP) get_sst_triple(l2rec);
        break;
      case CAT_bias_sst:
       pbuf = (VOIDP) get_bias_sst(l2rec);
        break;
      case CAT_bias_sst4:
        pbuf = (VOIDP) get_bias_sst4(l2rec);
        break;
      case CAT_bias_sst_triple:
        pbuf = (VOIDP) get_bias_sst_triple(l2rec);
        break;
      case CAT_bias_mean_sst:
        pbuf = (VOIDP) get_bias_mean_sst(l2rec);
        break;
      case CAT_bias_mean_sst4:
        pbuf = (VOIDP) get_bias_mean_sst4(l2rec);
        break;
      case CAT_bias_mean_sst_triple:
        pbuf = (VOIDP) get_bias_mean_sst_triple(l2rec);
        break;
      case CAT_counts_sst:
        pbuf = (VOIDP) get_counts_sst(l2rec);
        break;
      case CAT_counts_sst4:
        pbuf = (VOIDP) get_counts_sst4(l2rec);
        break;
      case CAT_counts_sst_triple:
        pbuf = (VOIDP) get_counts_sst_triple(l2rec);
        break;
      case CAT_stdv_sst:
        pbuf = (VOIDP) get_stdv_sst(l2rec);
        break;
      case CAT_stdv_sst4:
        pbuf = (VOIDP) get_stdv_sst4(l2rec);
        break;
      case CAT_stdv_sst_triple:
        pbuf = (VOIDP) get_stdv_sst_triple(l2rec);
        break;
      case CAT_flags_sst:
        pbuf = (VOIDP) get_flags_sst(l2rec); 
        break;
      case CAT_flags_sst4:
        pbuf = (VOIDP) get_flags_sst4(l2rec); 
        break;
      case CAT_flags_sst_triple:
        pbuf = (VOIDP) get_flags_sst_triple(l2rec); 
        break;
      case CAT_qual_sst:
        pbuf = (VOIDP) get_qual_sst(l2rec); 
        break;
      case CAT_qual_sst4:
        pbuf = (VOIDP) get_qual_sst4(l2rec); 
        break;
      case CAT_qual_sst_triple:
        pbuf = (VOIDP) get_qual_sst_triple(l2rec); 
        break;
      case CAT_sst_treesum:
        pbuf = (VOIDP) get_sst_treesum(l2rec); 
        break;
      case CAT_vLt:
      case CAT_vtLw:
      case CAT_vLw:
      case CAT_vnLw:
          vcal(l2rec, p, fbuf);
          applyMultiplier(fbuf, fbuf, numPixels, 10.0);
          pbuf = (VOIDP) fbuf;
          break;
      case CAT_vgain:
      case CAT_vbsat:
      case CAT_vbtgt:
        vcal(l2rec, p, fbuf); 
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_ice_frac:
        get_ice_frac(l2rec, fbuf);
        pbuf =  (VOIDP) fbuf;
        break;
      case CAT_class_ward_owmc:
        pbuf = get_class_ward_owmc(l2rec);
        break;
      case CAT_class_k_owmc:
        pbuf = get_class_k_owmc(l2rec);
        break;
      case CAT_class_34k_w_owmc:
        pbuf = get_class_34k_w_owmc(l2rec);
        break;
      case CAT_a_swim:
      case CAT_bb_swim:
      case CAT_adg_swim:
      case CAT_aph_swim:
      case CAT_bbp_swim:
      case CAT_Kd_swim:
      case CAT_edz_swim:
      case CAT_chl_swim:
      case CAT_tsm_swim:
        get_swim(l2rec, p, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_elev:
        pbuf = (VOIDP) l2rec->elev;
        break;
      case CAT_microplankton_hirata:
      case CAT_diatoms_hirata:
      case CAT_greenalgae_hirata:
      case CAT_picoplankton_hirata:
      case CAT_prokaryotes_hirata:
      case CAT_prochlorococcus_hirata:
      case CAT_dinoflagellates_hirata:
      case CAT_nanoplankton_hirata:
      case CAT_picoeukaryotes_hirata:
      case CAT_prymnesiophytes_hirata:
        get_pft_hirata(l2rec, p, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_microplankton_uitz:
      case CAT_nanoplankton_uitz:
      case CAT_picoplankton_uitz:
              get_pft_uitz(l2rec, p, fbuf);
              pbuf = (VOIDP) fbuf;
              break;
      case CAT_opp_befa:
      case CAT_opp_eppley:
      case CAT_opp_cbpm2:
      case CAT_opp_mld:
      case CAT_opp_zno3:
      case CAT_opp_bbp:
      case CAT_opp_par:
        get_opp(l2rec, p->cat_ix, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_CI_stumpf:
      case CAT_CI_cyano:
      case CAT_CI_noncyano:
      case CAT_MCI_stumpf:
        get_habs_ci(l2rec, p, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_MPH_chl:
        get_habs_mph(l2rec, p, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_MPH_flags:
        get_habs_mph_flags(l2rec, p, fbuf);
        pbuf = (VOIDP) fbuf;
        break;
      case CAT_flags_habs:
          get_habs_cldmask(l2rec, fbuf);
        pbuf = (VOIDP) fbuf;
        break;

      case CAT_microplankton_abundanceksm:
      case CAT_nanoplankton_abundanceksm:
      case CAT_picoplankton_abundanceksm:
      case CAT_microplankton_volumeksm:
      case CAT_nanoplankton_volumeksm:
      case CAT_picoplankton_volumeksm:
      case CAT_microplankton_ratioksm:
      case CAT_nanoplankton_ratioksm:
      case CAT_picoplankton_ratioksm:
      	get_psd_ksm(l2rec, p, fbuf);
      	pbuf = (VOIDP) fbuf;
      	break;
     default:
	fprintf(stderr, "-E- %s Line %d: Unknown product catalogue ID %d.\n",
            __FILE__,__LINE__,p->cat_ix);
        exit(1);
    }

    return(pbuf);
}

