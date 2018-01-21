cmake_minimum_required(VERSION 2.6)

PROJECT(OCSSW)

ADD_LIBRARY(bin STATIC
  read_l3b_meta.c
  getl3b_misc.c
  l3b_misc.c
  put_l3b.c
  get_beg_ext.c
  wr_vdata.c
  bin_csub.c
  write_l3b_meta.c
)

target_link_libraries(bin
  dfutils
)
