cmake_minimum_required(VERSION 2.6)

project(GRIB)

add_executable(wgrib2
  Alias.c
  CodeTable.c
  Code_Values.c
  Data.c
  Dump.c
  Ensemble.c
  File.c
  FlagTable.c
  Flush.c
  Grads.c
  Grib.c
  Help.c
  Inv.c
  Latlon.c
  Level.c
  Lola.c
  Macro.c
  Match.c
  Misc.c
  Netcdf.c
  Range.c
  Scan.c
  Sec0.c
  Sec1.c
  Sec3.c
  Sec4.c
  Sec5.c
  Sec6.c
  Sec_len.c
  Tigge.c
  Undefine.c
  VerfTime.c
  addtime.c
  cname.c
  earth.c
  fatal_error.c
  flt2ieee.c
  fnlist.c
  geo.c
  gribtab.c
  ieee2flt.c
  int8.c
  intpower.c
  missing.c
  rd_grib2_msg.c
  seekgrib2.c
  wgrib2.c
  wrt_grib.c
  wrtieee.c
  g2_free.c
  g2_getfld.c
  gbits.c
  g2_unpack1.c
  g2_unpack2.c
  g2_unpack3.c
  g2_unpack4.c
  g2_unpack5.c
  g2_unpack6.c
  g2_unpack7.c
  gridtemplates.c
  pdstemplates.c
  rdieee.c
  drstemplates.c
  simunpack.c
  comunpack.c
  specunpack.c
  int_power.c
  jpcunpack.c
  dec_jpeg2000.c
  jas_stream.c
  jas_image.c
  jas_seq.c
  jpc_dec.c
  jas_malloc.c
  jas_cm.c
  jas_string.c
  jpc_t1dec.c
  jas_debug.c
  jpc_tsfb.c
  jpc_mct.c
  jpc_cs.c
  jpc_t2cod.c
  jas_tvp.c
  jpc_t1cod.c
  jpc_bs.c
  jpc_mqdec.c
  jpc_tagtree.c
  jpc_t2dec.c
  jas_icc.c
  jpc_qmfb.c
  jpc_mqcod.c
  jas_iccdata.c
  jpc_math.c
)

add_definitions(-DUSE_TIGGE -DUSE_JPEG2000 -DUSE_REGEX)

include_directories (
  .
)

target_link_libraries(wgrib2
  m
)

install (TARGETS wgrib2 DESTINATION bin3/$ENV{OCSSW_ARCH})

