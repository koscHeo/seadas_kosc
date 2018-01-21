cmake_minimum_required(VERSION 2.6)

project(OCSSW)

add_executable(l1agen_seawifs
  main_l1agen.c
  swl0_utils.c
  swl1_hdf.c
  getl0indx.c
  printindx.c
  temporal_anomaly.c
  getl0scene.c
  getl0scene_nav.c
  printscene.c
  getorbdata.c
  getnavdata.c
  printnav.c
  getl1rec.c
  getorbnum.c
  instlm.c
  mkmeta.c
  getmetrics.c
)

add_executable(l0info_seawifs
  main_l0info.c
  getl0indx.c
  printindx.c
  getorbnum.c
  getl0scene.c
  printscene.c
  swl0_utils.c
)

add_executable(l0merge_seawifs
  main_l0merge.c
  getl0indx.c
  printindx.c
  swl0_utils.c
)

target_link_libraries(l1agen_seawifs
  seawifs
  swfnav
  nav
  hdfutils
  genutils
)

target_link_libraries(l0info_seawifs
  genutils
)

target_link_libraries(l0merge_seawifs
  genutils
)

if (${Fortran_COMPILER_NAME} MATCHES "ifort")
    set_target_properties(l1agen_seawifs PROPERTIES LINKER_LANGUAGE Fortran)
endif (${Fortran_COMPILER_NAME} MATCHES "ifort")

include(test_l1agen_seawifs_config.cmake OPTIONAL)

install (TARGETS l1agen_seawifs l0info_seawifs l0merge_seawifs 
  DESTINATION bin/${OCSSW_ARCH})
