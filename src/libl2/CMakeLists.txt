cmake_minimum_required(VERSION 2.6)

PROJECT(OCSSW)

ADD_LIBRARY(l2 STATIC
  readL2scan.c
  setupflags.c
  get_product_table.c
)

target_link_libraries(l2
  seawifs
  dfutils
)
