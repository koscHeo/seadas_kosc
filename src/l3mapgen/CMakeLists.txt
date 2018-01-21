cmake_minimum_required(VERSION 2.6)

project(OCSSW)

set (L3MAPGEN_WC_REVISION 0)

if(USE_SVNVERSION)
  find_package(Subversion)
  if(SUBVERSION_FOUND)
      Subversion_WC_INFO(${PROJECT_SOURCE_DIR} L3MAPGEN)
      message("Current subversion revision is ${L3MAPGEN_WC_REVISION}")
  endif(SUBVERSION_FOUND)
endif(USE_SVNVERSION)

set (L3MAPGEN_VERSION_MAJOR 1)
set (L3MAPGEN_VERSION_MINOR 0)
set (L3MAPGEN_VERSION_PATCH_LEVEL 1)

configure_file (
  "${CMAKE_CURRENT_LIST_DIR}/version.h.in"
  "${CMAKE_CURRENT_LIST_DIR}/version.h"
  )

add_executable(l3mapgen
  OutFile.cpp
  l3mapgen.cpp
  l3mapgen_input.cpp
  )

target_link_libraries(l3mapgen
  bin++
  l2
  png
  geotiff
  tiff
  proj
  )

# Include the test definitions files
include(test_l3mapgen_config.cmake OPTIONAL)

install (TARGETS l3mapgen DESTINATION bin/${OCSSW_ARCH})
