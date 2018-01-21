cmake_minimum_required(VERSION 2.6)

PROJECT(OCSSW)

INCLUDE_DIRECTORIES (
  ${CMAKE_SOURCE_DIR}/include/osmi
)

ADD_LIBRARY(osmi STATIC
  get_cal_osmi.c
  get_cal_misc_osmi.c
  calibrate_l1a_osmi.c
  calday.c
  obliq.c
  sunpos.c
  ctogd.c
  gmha.c
  tconv.c
  ctotc.c
  julian.c
  time-utils.c
  eanom.c
  locate.c
)

target_link_libraries(osmi
  seawifs
  hdf4utils
)
