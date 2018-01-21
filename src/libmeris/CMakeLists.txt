cmake_minimum_required(VERSION 2.6)

PROJECT(OCSSW)

INCLUDE_DIRECTORIES (
  ${CMAKE_SOURCE_DIR}/include/meris
)

ADD_LIBRARY(meris STATIC
  epr_api.c
  epr_band.c
  epr_bitmask.c
  epr_core.c
  epr_dataset.c
  epr_dddb.c
  epr_dsd.c
  epr_dump.c
  epr_field.c
  epr_msph.c
  epr_param.c
  epr_product.c
  epr_ptrarray.c
  epr_record.c
  epr_string.c
  epr_swap.c
  epr_typconv.c
)

target_link_libraries(meris
  m
)
