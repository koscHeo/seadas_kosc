cmake_minimum_required (VERSION 2.6)

project (HEALPIX)

ADD_LIBRARY(chealpix STATIC
chealpix.c
  )

target_link_libraries(chealpix libcfitsio)

install (TARGETS chealpix ARCHIVE DESTINATION $ENV{LIB3_DIR}/lib)
