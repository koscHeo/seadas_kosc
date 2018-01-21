cmake_minimum_required(VERSION 2.6)

project(OCSSW)

add_executable(l1bextract_meris
  l1bextract_meris.cpp
  EnvsatFile.cpp
  EnvsatMPH.cpp
  EnvsatSPH.cpp
  EnvsatDSD.cpp
  EnvsatDSR.cpp
  EnvsatUtil.cpp
  MerisSPH.cpp
  TiepointDSR.cpp
  MeasurementDSR.cpp
  RadianceDSR.cpp
  FlagsDSR.cpp
)

target_link_libraries(l1bextract_meris
  genutils
)

install (TARGETS l1bextract_meris DESTINATION bin/${OCSSW_ARCH})
