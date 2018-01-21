cmake_minimum_required(VERSION 2.6)

PROJECT(OCSSW)

ADD_LIBRARY(genutils STATIC
  alloc_2d.c
  endianess.c
  parse_file_name.c
  swapc_bytes.c
  fread_swap.c
  genutils_globals.c
  clo.c
  table_io.cpp
  nr_spline.c
  lspline.c
  lenstr.f
  filesize.c
  getlut_file.c
  lowcase.c
  upcase.c
  allocateMemory.c
  trimBlanks.c
  xmlUtils.c
  isValidInt.c
  fileFormatUtils.c
)

target_link_libraries(genutils
  roxml
  pthread
  timeutils
  m
)
