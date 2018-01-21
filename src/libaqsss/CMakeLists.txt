cmake_minimum_required(VERSION 2.6)

PROJECT(OCSSW)

ADD_LIBRARY(aqsss STATIC
  get_ancillary_data.f
  fd_ta_sun.f
  stokes_converters.f
  fd_ta_galaxy.f
  fd_dta_symm.f
  fd_dta_galaxy.f
  fd_dta_symm.f
  fd_water_refl.f
  find_refl_tot.f
  fd_seaice_refl.f
  fd_climate_land_refl.f
  fd_climate_sss.f90
  fd_sun_backscatter.f
  adjust_tagal_ref.f
  fd_land.f
  fd_tb_toa_land_corr.f
  fd_land_refl.f
  fd_tausq.f
  fd_ta_moon.f
  phir_routines.f
  fd_dtb_roughness.f90
  static_data_module.f 
  wind_speed_retrieval_module.f90
  salinity_module.f90
  fd_sss_sigma0_1par.f
  fd_DI_U.f90
  interpolate_L2_error.f90
)

target_link_libraries(aqsss
  aquarius
  hdf5utils
)


