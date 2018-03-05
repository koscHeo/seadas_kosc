cmake_minimum_required(VERSION 2.6)

project(OCSSW)

#set up versioning info
set (L2GEN_WC_REVISION 0)
# check if this is a subversion repo and set working copy revision number
# appropriately
if(USE_SVNVERSION)
    find_package(Subversion)
    if(SUBVERSION_FOUND)
        Subversion_WC_INFO(${PROJECT_SOURCE_DIR} L2GEN)
        message("Current subversion revision is ${L2GEN_WC_REVISION}")
    endif(SUBVERSION_FOUND)
endif(USE_SVNVERSION)

set (L2GEN_VERSION_MAJOR 8)
set (L2GEN_VERSION_MINOR 10)
set (L2GEN_VERSION_PATCH_LEVEL 3)

configure_file (
  "${CMAKE_CURRENT_LIST_DIR}/version.h.in"
  "${CMAKE_CURRENT_LIST_DIR}/version.h"
  )

if (BUILD_VIIRS_L1)
    set (VIIRS_L1_INCLUDES 
        ${CMAKE_SOURCE_DIR}/include/ViirsL1BCal 
        ${CMAKE_SOURCE_DIR}/include/VL1control
    )
    set (VIIRS_L1_SRCS calibrate_viirs.c l1b_viirs_nc.c l1_viirs_nc.c)
    set (VIIRS_L1_LIBS ViirsL1BCal)
else (BUILD_VIIRS_L1)
    set (VIIRS_L1_INCLUDES )
    set (VIIRS_L1_SRCS calibrate_viirs.c l1b_viirs_nc.c l1_viirs_nc_stub.c)
    set (VIIRS_L1_LIBS )
endif (BUILD_VIIRS_L1)

add_executable(l2gen_kosc
  main_l2gen.c
  l1subpix.c
  get_l2prod_index.c
  msl12_input.c
  filter.c
  getl1rec.c
  loadl1.c
  xcal.c
  cpl1rec.c
  l1_io.c
  scene_meta.c
  target_io.c
  l1_hdf_generic_read.c
  l1_nc_generic_read.c
  l1_generic_write.c
  l1_mos_hdf.c
  l1_octs_hdf.c
  l1_czcs_hdf.c
  ll2vec.c
  l1_hmodis_hdf.c
  ${VIIRS_L1_SRCS}
  l1_xcal_hdf.c
  etbpsub.f
  ftrim.c
  hmf8.f
  lenstr.f
  raygetpol.f
  simpsn.f
  airmass_avhrr.f
  avhrrsub5h.f
  l1_aci_hdf.c
  l1_meris_N1.c
  l1_meris_CC.c
  l1_orca.c
  l1_aviris.c
  l1_olci.c
  l1_prism.c
  l1_viirs_h5.c
  l1_hico_h5.c
  l1_goci.c
  jplaeriallib.c
  viirs_utls.c
  bioOptBandShift.c
  prodgen.c
  prodlist.c
  l2_generic.c
  l1a_seawifs.c
  l1a_osmi.c
  l1_ocm_hdf.c
  l1_ocm2_hdf.c
  l1_ocmdb_hdf.c
  l1_oli.c
  intpos.f
  l12_seawifs.c
  getformat.c
  alloc_l1.c
  alloc_l2.c
  alloc_target.c
  alloc_vcal.c
  convl12.c
  convl21.c
  atmocor1.c
  atmocor1_land.c
  whitecaps.c
  rayleigh.c
  polcor.c
  setanc.c
  setflags.c
  cloud_flag.c
  get_chl.c
  get_habs.c
  get_es.c
  get_Kd.c
  get_ndvi.c
  get_smoke.c
  get_depth.c
  get_poc.c
  get_par.c
  get_toa_refl.c
  get_tricho.c
  get_opp.c
  photic_depth.c
  cdom_morel.c
  calc_par.c
  cdom_mannino.c
  par_utils.c
  ipar.c
  fluorescence.c
  water_vapor.c
  airmass.c
  aerosol.c
  atmocor2.c
  glint.c
  b128_msk_get.c
  b128_wd_bit.c
  b128_box_num.c
  read_mask.c
  getglint.f
  fresnel.c
  get_rhown_nir.c
  mumm.c
  nlw_outband.c
  brdf.c
  filehandle_init.c
  get_dem_height.c
  get_nc_height.c
  get_height.c
  get_rhos.c
  calcite.c
  ice_mask.c
  get_ice_frac.c
  smi_climatology.c
  bin_climatology.c
  aer_io.c
  alloc_aer.c
  amoeba.c
  gsm.c
  giop.c
  mgiop.c
  carder.c
  las_iop.c
  pml.c
  get_pml.c
  pml_iop_calculate.c
  pml_iop_config.c
  pml_iop_tables.c
  qaa.c
  get_qaa.c
  flags_iop.c
  niwa_iop.c
  get_niwa_iop.c
  aph.c
  water.c
  get_f0.c
  windex.c
  brightness.c
  sst.c
  sstref.c
  sssref.c
  soa_sma.c
  atmcor_soa.f
  soa_sma_utils.f
  gas_trans.c
  vcal.c
  myprod.c
  dtran_brdf.f
  turbid.c
  init_l2.c
  init_l1.c
  h5io.c
  convert_band.c
  fuzzy_func_v3.c
  gammln.c
  gcf.c
  gser.c
  nrutil.c
  sprsax.c
  sprsin.c
  covariance_inversion.c
  lubksb.c
  ludcmp.c
  owt.c
  get_owmc.c
  lonlat2pixline.c
  seawater.c
  seawater_get.c
  viirs_pxcvt.c
  virtual_constellation.c
  get_bsi.c
  swim.c
  elev.c
  read_l3bin.cpp
  anc_acq.c
  met_cvt.c
  ncio.c
  get_pft_hirata.c  
  get_pft_uitz.c
  get_mld.c
  get_zno3.c
  get_atrem_corl1.c
  rdatreminfo.c
  atrem_app_refl_plus_gas_removal_l2.f90
  atrem_cubeio.f90
  atrem_bndprms.f
  atrem_tpvmr_init.f 
  numerical.c
  raman.c
  get_psd_ksm.c
  alloc_l1_n.c
  alloc_l2_n.c
  alloc_aer_n.c
  alloc_target_n.c
)

include_directories (
  ../lib24to8
  ${CMAKE_SOURCE_DIR}/include/osmi
  ${CMAKE_SOURCE_DIR}/include/meris
  ${VIIRS_L1_INCLUDES}
  )

target_link_libraries(l2gen_kosc
  bin++
  l2
  anc
  czcs
  osmi
  meris
  goci
  seawifs
  oli
  nav
  levmar
  lapack
  blas
  gsl
  gslcblas
  geotiff
  tiff
  netcdff
  dfutils
  fann
  ${VIIRS_L1_LIBS}
  ${IFORT_CXX_LIBS}
  ${EXTRA_OS_LIBS}
  )

#set (CMAKE_C_FLAGS_RELEASE "-O1") 
# Use -Wpadded -Wpacked to warn about align-double. Use  '#pragma GCC diagnostic ignored "-Wpadded"' in file if you don't want to see messages
# -malign-double needed for accessing fortran common blocks from C structs

set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -malign-double -mpreferred-stack-boundary=8")
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -malign-double -mpreferred-stack-boundary=8")

if (${Fortran_COMPILER_NAME} MATCHES "ifort")
  set_target_properties(l3gen PROPERTIES LINKER_LANGUAGE Fortran)
elseif (${Fortran_COMPILER_NAME} MATCHES "gfortran")
    set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -malign-double")
    set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} -malign-double")
endif (${Fortran_COMPILER_NAME} MATCHES "ifort")

# Include the test definitions files
include(test_l2gen_config.cmake OPTIONAL)
include(test_l1bgen_generic_config.cmake OPTIONAL)
include(test_l1brsgen_config.cmake OPTIONAL)
include(test_l1mapgen_config.cmake OPTIONAL)
include(test_l1info_config.cmake OPTIONAL)
include(test_lonlat2pixline_config.cmake OPTIONAL)
include(test_l3gen_config.cmake OPTIONAL)

install (TARGETS l2gen_kosc DESTINATION bin/${OCSSW_ARCH})

