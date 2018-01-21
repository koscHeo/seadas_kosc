cmake_minimum_required(VERSION 2.6)

PROJECT(OCSSW)

ADD_LIBRARY(aquarius STATIC
  b128_msk_get.c
  FileStager.cpp
  hdf5_Aquarius.cpp
  invert_3by3.f
  l1_io.cpp
  l2_io.cpp
  sedgesort.c
  trig_deg.f
  tlm_convert.cpp
  tai_utc_unix.cpp
)

target_link_libraries(aquarius
  hdf5utils
  ${IFORT_CXX_LIBS}
)
