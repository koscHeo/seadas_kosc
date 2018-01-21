#ifndef _L2PROD_H
#define _L2PROD_H

#define MAXPROD             1000  /* Max number of individual products    */
#define PRODSTRLEN          2048  /* String length limit for product spec */

#define CAT_Lt                 0
#define CAT_Lr                 1
#define CAT_La                 2
#define CAT_Lw                 3
#define CAT_nLw                4
#define CAT_tLf                5
#define CAT_Taua               6
#define CAT_Es                 7
#define CAT_t_sol              8
#define CAT_t_sen              9

#define CAT_tg_sol            10
#define CAT_tg_sen            11
#define CAT_solz              12
#define CAT_sola              13
#define CAT_sena              14
#define CAT_senz              15
#define CAT_ozone             16
#define CAT_windspeed         17
#define CAT_pressure          18
#define CAT_humidity          19
#define CAT_water_vapor       20

#define CAT_no2_tropo         21
#define CAT_epsilon           22
#define CAT_aer_model         23
#define CAT_aer_ratio         24
#define CAT_l2_flags          25
#define CAT_chl_oc2           26
#define CAT_depth             27
#define CAT_TLg               28
#define CAT_par               29

#define CAT_angstrom          30
#define CAT_Kd_mueller        31
#define CAT_ndvi              32
#define CAT_glint_coef        33
#define CAT_num_iter          34
#define CAT_brdf              35
#define CAT_aerindex          36
#define CAT_rhos              37
#define CAT_evi               38
#define CAT_smoke             39

#define CAT_no2_strat         40
#define CAT_chl_oc4           41
#define CAT_windangle         42
#define CAT_mwind             43
#define CAT_zwind             44
#define CAT_cloud_albedo      45
#define CAT_t_o2              46
#define CAT_fsol              47
#define CAT_rhot              48
#define CAT_height            49

#define CAT_calcite_3b        50
#define CAT_flh               51
#define CAT_Rrs               52
#define CAT_sst               53
#define CAT_chl_gsm           54
#define CAT_adg_gsm           55
#define CAT_bbp_gsm           56
#define CAT_t_h2o             57
#define CAT_sstref            58
#define CAT_L_q               59

#define CAT_L_u               60
#define CAT_polcor            61
#define CAT_chl_oc3           62
/* #define CAT_tsm_clark         63  open for use */
/* #define CAT_poc_clark         64  open for use */
#define CAT_slot              65
#define CAT_pixnum            66
#define CAT_detnum            67
#define CAT_mside             68
#define CAT_alpha             69

#define CAT_no2_frac          70
#define CAT_fsat              71
#define CAT_a_carder          72
#define CAT_bb_carder         73
#define CAT_aph_carder        74
#define CAT_adg_carder        75
#define CAT_chl_carder        76
#define CAT_flags_carder      77
#define CAT_dpol              78
#define CAT_calcite_2b        79

#define CAT_a_gsm             80
#define CAT_bb_gsm            81
#define CAT_aph_gsm           82
/* #define CAT_fqy2              83  open for use */ 
#define CAT_ipar              84
#define CAT_fqy               85
#define CAT_bbp_carder        86
#define CAT_a                 87
#define CAT_bb                88
#define CAT_a_qaa             89

#define CAT_bb_qaa            90
#define CAT_aph_qaa           91
#define CAT_adg_qaa           92
#define CAT_bbp_qaa           93
#define CAT_Kd_lee            94
#define CAT_Kd_obpg           95
#define CAT_iter_gsm          96
#define CAT_sst4              97
#define CAT_chl_soa           98
#define CAT_bbp_soa           99

#define CAT_adg_soa          100
#define CAT_pcentcdm_soa     101
#define CAT_w0_soa           102
#define CAT_v_soa            103
#define CAT_sssref           104
#define CAT_flags_sst        105
#define CAT_flags_sst4       106
#define CAT_qual_sst         107
#define CAT_qual_sst4        108
#define CAT_calcite          109

#define CAT_myprod1          110
#define CAT_myprod2          111
#define CAT_myprod3          112
#define CAT_myprod4          113
#define CAT_myprod5          114
#define CAT_myprod6          115
#define CAT_myprod7          116
#define CAT_myprod8          117
#define CAT_myprod9          118
#define CAT_myprod10         119

#define CAT_bias_sst         120
#define CAT_bias_sst4        121
#define CAT_stdv_sst         122
#define CAT_stdv_sst4        123
#define CAT_rhom             124
#define CAT_Kd_morel         125
#define CAT_tindx_shi        126
#define CAT_KPAR_morel       127
#define CAT_Zhl_morel        128
#define CAT_Zeu_morel        129

#define CAT_Zsd_morel        130
#define CAT_tindx_morel      131
#define CAT_Kd_KD2           132
#define CAT_vgain            133
#define CAT_vLt              134
#define CAT_vtLw             135
#define CAT_vLw              136
#define CAT_vnLw             137
#define CAT_vbsat            138
#define CAT_vbtgt            139

#define CAT_chl_carder_emp   140
#define CAT_Zphotic_lee      141
#define CAT_b_qaa            142
#define CAT_c_qaa            143
#define CAT_Kd_532           144
#define CAT_KPAR_lee         145
#define CAT_BT               146
#define CAT_BT_39            147  /* phase-out */
#define CAT_BT_40            148  /* phase-out */
#define CAT_BT_11            149  /* phase-out */

#define CAT_BT_12            150  /* phase-out */
#define CAT_Ltir             151
#define CAT_poc_stramski_443 152
#define CAT_poc_stramski_490 153
#define CAT_chl_sma          154
#define CAT_bbp_sma          155
#define CAT_adg_sma          156
#define CAT_w0_sma           157
#define CAT_dom_sma          158
#define CAT_a_pml            159

#define CAT_bb_pml           160
#define CAT_bbp_pml          161
#define CAT_aph_pml          162
#define CAT_adg_pml          163
#define CAT_mod_rrs_qaa      164
#define CAT_a_las            165
#define CAT_b_las            166
#define CAT_c_las            167
#define CAT_bb_las           168
#define CAT_bbp_las          169

#define CAT_a_giop           170
#define CAT_bb_giop          171
#define CAT_bbp_giop         172
#define CAT_aph_giop         173
#define CAT_adg_giop         174
#define CAT_chl_giop         175
#define CAT_a_unc_giop       176
#define CAT_bb_unc_giop      177
#define CAT_bbp_unc_giop     178
#define CAT_aph_unc_giop     179
#define CAT_adg_unc_giop     180
#define CAT_chl_unc_giop     181
#define CAT_aphs_giop        182
#define CAT_adgs_giop        183
#define CAT_bbps_giop        184
#define CAT_iter_giop        185
#define CAT_rrsdiff_giop     186
#define CAT_chisqr_giop      187
#define CAT_fitpar_giop      188

#define CAT_mRrs_giop        190
#define CAT_flags_giop       191
#define CAT_relaz            192
#define CAT_flags_qaa        193
#define CAT_bbps_las         194
#define CAT_a_niwa           195
#define CAT_bb_niwa          196
#define CAT_flags_niwa       197
#define CAT_rho_cirrus       198
/* #define CAT_chl1_meris       199    open for use */
/* #define CAT_chl2_meris       200    open for use */
/* #define CAT_adg_meris        201    open for use */
/* #define CAT_tsm_meris        202    open for use */
/* #define CAT_par_meris        203    open for use */
/* #define CAT_angstrom_meris   204    open for use */
/* #define CAT_Taua_meris       205    open for use */
#define CAT_iCDOM_morel      206
#define CAT_pCDOM_morel      207
#define CAT_chl_morel        208
#define CAT_adg_morel        209
#define CAT_scattang         210

#define CAT_ms_epsilon       211
#define CAT_ice_frac         212
/* #define CAT_ipar2            213    open for use */
/* #define CAT_par2             214    open for use */
#define CAT_owt              215
#define CAT_owtn             216
#define CAT_owtd             217
#define CAT_chl_owterr       218
#define CAT_class_ward_owmc  219

#define CAT_class_k_owmc     220
#define CAT_class_34k_w_owmc 221
/* #define CAT_fdiff            222     open for use */
#define CAT_Zsd_gbr          223
#define CAT_chl_cdomcorr_morel   224
#define CAT_chl_hu           225
#define CAT_Lt_unc           226
#define CAT_nLw_unc          227
#define CAT_Rrs_unc          228
#define CAT_chl_oci          229

#define CAT_chl_oc3c         230
#define CAT_chl_oci2         231

#define CAT_Rrs_vc           233
#define CAT_chl_vc           234
#define CAT_aw               235
#define CAT_bbw              236
#define CAT_nw               237

#define CAT_chl_mgiop        238
#define CAT_bbp_mgiop        239
#define CAT_adg_mgiop        240
#define CAT_aph_mgiop        241
#define CAT_npix_mgiop       242
#define CAT_crat_mgiop       243
#define CAT_fitpar_mgiop     244

#define CAT_BSi              245
#define CAT_bbws             246

#define CAT_a_swim           247
#define CAT_bb_swim          248
#define CAT_adg_swim         249
#define CAT_aph_swim         250
#define CAT_bbp_swim         251

#define CAT_elev             252
#define CAT_Kd_jamet         253
#define CAT_chl_cdr          254

#define CAT_Kd_swim          255
#define CAT_edz_swim         256
#define CAT_chl_swim         257
#define CAT_tsm_swim         258

#define CAT_microplankton_hirata    259
#define CAT_diatoms_hirata          260
#define CAT_greenalgae_hirata       261
#define CAT_picoplankton_hirata     262
#define CAT_prokaryotes_hirata      263
#define CAT_prochlorococcus_hirata  264
#define CAT_dinoflagellates_hirata  265
#define CAT_nanoplankton_hirata     266
#define CAT_picoeukaryotes_hirata   267
#define CAT_prymnesiophytes_hirata  268

#define CAT_microplankton_uitz      270
#define CAT_nanoplankton_uitz       271
#define CAT_picoplankton_uitz       272

/* additional SST SSES products */
#define CAT_bias_mean_sst    273
#define CAT_bias_mean_sst4   274
#define CAT_counts_sst       275
#define CAT_counts_sst4      276

#define CAT_ag_412_mlrc             277
#define CAT_Sg_275_295_mlrc         278
#define CAT_Sg_300_600_mlrc         279

#define CAT_opp_befa                280
#define CAT_opp_eppley              281
#define CAT_opp_cbpm2               282

#define CAT_chl_abi                 283

#define CAT_Kd_rhos                 284
#define CAT_CI_stumpf               285
#define CAT_MCI_stumpf              286
#define CAT_MPH_chl                 287
#define CAT_MPH_flags               288

#define CAT_sst_triple             289
#define CAT_flags_sst_triple       290
#define CAT_qual_sst_triple        291
#define CAT_bias_sst_triple        292
#define CAT_stdv_sst_triple        293
#define CAT_bias_mean_sst_triple   294
#define CAT_counts_sst_triple      295

#define CAT_microplankton_abundanceksm     296
#define CAT_nanoplankton_abundanceksm      297
#define CAT_picoplankton_abundanceksm      298

#define CAT_microplankton_volumeksm        299
#define CAT_nanoplankton_volumeksm         300
#define CAT_picoplankton_volumeksm         301

#define CAT_microplankton_ratioksm         302
#define CAT_nanoplankton_ratioksm          303
#define CAT_picoplankton_ratioksm          304

#define CAT_flags_habs             305

#define CAT_opp_mld                 306
#define CAT_opp_zno3                307
#define CAT_opp_par                 308
#define CAT_opp_bbp                 309

/*additional giop products for aLMI*/
#define CAT_acdom_giop              310
#define CAT_anap_giop               311
#define CAT_bbph_giop               312
#define CAT_bbnap_giop              313
#define CAT_acdom_unc_giop          314
#define CAT_anap_unc_giop           315
#define CAT_bbph_unc_giop           316
#define CAT_bbnap_unc_giop          317
#define CAT_opt_siop_giop           318

/* Additional Cyano Index products */
#define CAT_CI_cyano                319
#define CAT_CI_noncyano             320

#define CAT_sst_treesum             321
#endif
