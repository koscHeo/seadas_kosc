cmake_minimum_required(VERSION 2.6)

project(GRIB)

add_executable(wgrib
  wgrib.c
)

install (TARGETS wgrib DESTINATION bin3/$ENV{OCSSW_ARCH})

