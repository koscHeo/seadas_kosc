cmake_minimum_required(VERSION 2.6)

PROJECT(GRIB)

ADD_LIBRARY(g2 STATIC
  g2grids.f
  params_ecmwf.f
  gribmod.f
  realloc.f
  gridtemplates.f
  pdstemplates.f
  drstemplates.f
  params.f
  addfield.F
  addgrid.f
  addlocal.f
  cmplxpack.f
  compack.f
  comunpack.f
  gb_info.f
  gbytesc.f
  gdt2gds.f
  getdim.f
  getfield.f
  getg2i.f
  getg2ir.f
  getgb2l.f
  getgb2.f
  getgb2p.f
  getgb2r.f
  getgb2rp.f
  getgb2s.f
  getidx.f
  getlocal.f
  getpoly.f
  gettemplates.f
  gf_free.f
  gf_getfld.f
  gf_unpack1.f
  gf_unpack2.f
  gf_unpack3.f
  gf_unpack4.f
  gf_unpack5.f
  gf_unpack6.f
  gf_unpack7.F
  gribcreate.f
  gribend.f
  gribinfo.f
  ixgb2.f
  jpcpack.f
  jpcunpack.f
  misspack.f
  mkieee.f
  pack_gp.f
  pngpack.f
  pngunpack.f
  putgb2.f
  rdieee.f
  reduce.f
  simpack.f
  simunpack.f
  skgb.f
  specpack.f
  specunpack.f
 )

target_link_libraries(g2
  w3
  jpeg
  m
)

install (TARGETS g2 ARCHIVE DESTINATION $ENV{LIB3_DIR}/lib)

