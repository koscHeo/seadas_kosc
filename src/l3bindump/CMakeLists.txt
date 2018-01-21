cmake_minimum_required(VERSION 2.6)

project(OCSSW)

set (L3BINDUMP_WC_REVISION 0)

if(USE_SVNVERSION)
    find_package(Subversion)
    if(SUBVERSION_FOUND)
        Subversion_WC_INFO(${PROJECT_SOURCE_DIR} L3BINDUMP)
        message("Current subversion revision is ${L3BINDUMP_WC_REVISION}")
    endif(SUBVERSION_FOUND)
endif(USE_SVNVERSION)

set (L3BINDUMP_VERSION_MAJOR 3)
set (L3BINDUMP_VERSION_MINOR 0)
set (L3BINDUMP_VERSION_PATCH_LEVEL 0)

configure_file (
  "${CMAKE_CURRENT_LIST_DIR}/version.h.in"
  "${CMAKE_CURRENT_LIST_DIR}/version.h"
  )


add_executable(l3bindump
  l3bindump_input.c
  compareObj.cpp
  main_l3bindump.cpp
)

target_link_libraries(l3bindump
  bin++
  ${EXTRA_OS_LIBS}
)


# Include the test definitions files
include(test_l3bindump_config.cmake OPTIONAL)

install (TARGETS l3bindump DESTINATION bin/${OCSSW_ARCH})

