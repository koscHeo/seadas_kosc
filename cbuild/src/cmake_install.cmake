# Install script for directory: /home/jaemoo/seadas-7.4/ocssw/build/src

# Set the install prefix
IF(NOT DEFINED CMAKE_INSTALL_PREFIX)
  SET(CMAKE_INSTALL_PREFIX "/home/jaemoo/seadas-7.4/ocssw/run")
ENDIF(NOT DEFINED CMAKE_INSTALL_PREFIX)
STRING(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
IF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  IF(BUILD_TYPE)
    STRING(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  ELSE(BUILD_TYPE)
    SET(CMAKE_INSTALL_CONFIG_NAME "Release")
  ENDIF(BUILD_TYPE)
  MESSAGE(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
ENDIF(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)

# Set the component getting installed.
IF(NOT CMAKE_INSTALL_COMPONENT)
  IF(COMPONENT)
    MESSAGE(STATUS "Install component: \"${COMPONENT}\"")
    SET(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  ELSE(COMPONENT)
    SET(CMAKE_INSTALL_COMPONENT)
  ENDIF(COMPONENT)
ENDIF(NOT CMAKE_INSTALL_COMPONENT)

# Install shared libraries without execute permission?
IF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  SET(CMAKE_INSTALL_SO_NO_EXE "0")
ENDIF(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)

IF(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/lib24to8/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/libanc/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/libnav/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/libczcs/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/libgenutils/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/libhdf4utils/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/libhdf5utils/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/libnetcdfutils/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/libtimeutils/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/libbin/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/libbin++/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/libl2/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/libmap/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/libmeris/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/libgoci/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/libseawifs/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/libswfnav/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/libosmi/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/liboli/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/libpiutils/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/libdfutils/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l2gen/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l2gen_kosc/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l2brsgen/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l2bin/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l3bin/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l1aextract_modis/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l1aextract_seawifs/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l1bextract_meris/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/smigen/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l2extract/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l2mapgen/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/smitoppm/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l3bindump/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/get_product_info/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l3mapgen/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l1agen_modis/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l1bgen_modist/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l1bgen_modisa/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/geogen_modis/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l0fix_modis/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l0cnst_write_modis/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l0chunk_modis/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/l1asubset_modis/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/atteph_info_modis/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/geocheck_modis/cmake_install.cmake")
  INCLUDE("/home/jaemoo/seadas-7.4/ocssw/build/cbuild/src/pdsmerge/cmake_install.cmake")

ENDIF(NOT CMAKE_INSTALL_LOCAL_ONLY)

