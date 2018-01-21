cmake_minimum_required(VERSION 2.6)

PROJECT(OCSSW)

ADD_LIBRARY(seawifs STATIC
  geonav.f
  geonav_lonlat.f
  get_cal_swf.c
  get_cal.c
  get_cal_misc.c
  calibrate_l1a.c
  l1a.c
  lac_st.c
  gac_st.c
  out_band_corr.c
  st_lt.c
  GetStationInfo.c
)

target_link_libraries(seawifs
  nav
  hdf4utils
)
