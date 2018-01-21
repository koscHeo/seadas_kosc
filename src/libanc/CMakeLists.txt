cmake_minimum_required(VERSION 2.6)

PROJECT(OCSSW)

ADD_LIBRARY(anc STATIC
  getanc.c
  HDFroutines.c
  dataintp.f
  julian.f
)

target_link_libraries(anc
  mfhdf
  df
  jpeg
  z
)
