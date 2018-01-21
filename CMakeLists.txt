cmake_minimum_required (VERSION 2.6)

project (OCSSW)

# search for test dierctory - if found, enable testing
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(OCSSWTest)

if (OCSSWTestDir_FOUND)
  message("Test Directory found! ...enabling testing...")
  enable_testing()
endif(OCSSWTestDir_FOUND)

# if Doxygen exists, create the API documentation 
find_package(Doxygen)
if(DOXYGEN_FOUND AND IS_DIRECTORY "${$ENV{OCSSWROOT}/doc}")
  configure_file($ENV{OCSSWROOT}/doc/Doxyfile
  ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(doc
      ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
      WORKING_DIRECTORY $ENV{OCSSWROOT}
      COMMENT "Generating API documentation with Doxygen" VERBATIM
  )
endif(DOXYGEN_FOUND AND IS_DIRECTORY "${$ENV{OCSSWROOT}/doc}")

# set this to Release or Debug
# note: rerun cmake if you change the environment var OCSSW_DEBUG
if ($ENV{OCSSW_DEBUG})
  set (CMAKE_BUILD_TYPE Debug)
else($ENV{OCSSW_DEBUG})
  set (CMAKE_BUILD_TYPE Release)
endif($ENV{OCSSW_DEBUG})
OPTION(USE_SVNVERSION "Build code with SVN revision tag in version" OFF)

# Finer-grained control...
OPTION(BUILD_AQUARIUS "Build the Aquarius code" OFF)
OPTION(BUILD_JPL "Build the JPL 'proprietary' code" OFF)
OPTION(BUILD_HISTORICAL "Build the historical code (e.g. old seawifs binaries
no longer needed)" OFF)
OPTION(BUILD_MISC "Build non-essential code" OFF)
OPTION(BUILD_VIIRS_L1 "Build VIIRS L1 code" OFF)
OPTION(BUILD_AHMAD_FRASER "Build the AF radiative transfer code" OFF)
OPTION(BUILD_ALL "Build all including  Aquarius, historical, misc, VIIRS L1 and Ahmad-Fraser code" OFF)

if(BUILD_JPL)
  set(BUILD_AQUARIUS ON)
endif(BUILD_JPL)

if(BUILD_ALL)
   set(BUILD_AQUARIUS ON)
   set(BUILD_JPL ON)
   set(BUILD_HISTORICAL ON)
   set(BUILD_MISC ON)
   set(BUILD_VIIRS_L1 ON)
   set(BUILD_AHMAD_FRASER ON)
endif(BUILD_ALL)

# set this to print out the compile lines
# can also use "make VERBOSE=1" when building
#set (CMAKE_VERBOSE_MAKEFILE TRUE)

enable_language (C)
enable_language (CXX)
enable_language (Fortran)

# make the build static
set (CMAKE_EXE_LINKER_FLAGS "-static -mcmodel=large")
set (CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")
set (CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "")
set (CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS "")

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

# extra libs needed by some OSs
set (EXTRA_OS_LIBS "")

# extra libs needed for C++ when using ifort as the linker
set (IFORT_CXX_LIBS "")

# init for the most normal cases
set (EXTRA_PGS_LIB "")
set (HDF4_LIBRARIES mfhdf df jpeg z)
set (HDF5_LIBRARIES hdf5_hl hdf5 z m)

# Set up C and CXX flags
set (CMAKE_C_FLAGS_RELEASE "-O3 -Wno-unused-result")
if (DEFINED ENV{GCC_TUNE})
  set (CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE}  $ENV{GCC_TUNE}")
endif(DEFINED ENV{GCC_TUNE}) 
set (CMAKE_C_FLAGS_DEBUG "-g3 -gdwarf-2 -Wall -Wno-unused-variable -Wno-unused-label")
set (CMAKE_CXX_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -Wno-unused-local-typedefs")
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wno-unused-local-typedefs")
set (CMAKE_C_FLAGS "-DLINUX -Dlinux")
set (CMAKE_CXX_FLAGS ${CMAKE_C_FLAGS})

# Set up FORTRAN flags
if (${Fortran_COMPILER_NAME} MATCHES "gfortran")
    set (CMAKE_Fortran_FLAGS "-DLINUX -Dlinux -ffixed-line-length-none -fconvert=big-endian -fno-automatic -fno-range-check -cpp")
    set (CMAKE_Fortran_FLAGS_RELEASE "-O3 -funroll-all-loops")
    set (CMAKE_Fortran_FLAGS_DEBUG "-g3 -gdwarf-2 -Wall -Wno-unused-variable -Wno-unused-label -Wno-unused-dummy-argument -Wno-conversion")
elseif (${Fortran_COMPILER_NAME} MATCHES "ifort")
    set (CMAKE_Fortran_FLAGS "-DLINUX -Dlinux -extend_source -convert big_endian -assume byterecl -save -cpp")
    set (CMAKE_Fortran_FLAGS_RELEASE "-O3")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
    set (CMAKE_Fortran_LINK_FLAGS "-static-intel -nofor_main")
    set (IFORT_CXX_LIBS gcc_eh stdc++)
else (${Fortran_COMPILER_NAME} MATCHES "gfortran")
    # other FORTRAN...*might* work...
    message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
    message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
    message ("No optimized Fortran compiler flags are known, we just try -O2...")
    set (CMAKE_Fortran_FLAGS "-DLINUX -Dlinux")
    set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
    set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (${Fortran_COMPILER_NAME} MATCHES "gfortran")

# set up OS-specific bits
if(${CMAKE_SYSTEM_NAME} MATCHES "Linux" )
  # standard Linux bits...
  set (HDF4_LIBRARIES ${HDF4_LIBRARIES} xdr)
  set (OCSSW_ARCH linux)
  set (EOS_ARCH linux)
  if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")
  # 64bit bits...
    #set (HDF5_LIBRARIES ${HDF5_LIBRARIES} dl)
    set (OCSSW_ARCH linux_64)
    set (EOS_ARCH linux64)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m64")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64")
    if (${Fortran_COMPILER_NAME} MATCHES "ifort")
      set (CMAKE_Fortran_FLAGS "-m64 ${CMAKE_Fortran_FLAGS}")
    endif(${Fortran_COMPILER_NAME} MATCHES "ifort")
  endif(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
# Mac bits
  set (OCSSW_ARCH macosx_intel)
  set (EOS_ARCH macintel64)
  set (CMAKE_C_FLAGS "-DMACINTOSH ${CMAKE_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "-DMACINTOSH ${CMAKE_CXX_FLAGS}")
  set (CMAKE_Fortran_FLAGS "-DMACINTOSH ${CMAKE_Fortran_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "-Wl,-no_compact_unwind -static-libstdc++ -static-libgcc -mmacosx-version-min=10.7")
#  set (EXTRA_OS_LIBS "gcc_eh") 
#  set (EXTRA_PGS_LIB /usr/lib/libiconv.dylib)
else(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
# Miscellaneous - not sure what we are...
  message ("CMAKE_SYSTEM_NAME: " ${CMAKE_SYSTEM_NAME})
  message ("System not recognized, using other")
  set (OCSSW_ARCH linux)
  set (EOS_ARCH linux)
endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

# stuff 'fer MODIS code
set (PGS_INC ${CMAKE_SOURCE_DIR}/lib3/EOS/include)
set (PGS_LIB ${CMAKE_SOURCE_DIR}/lib3/EOS/lib/${EOS_ARCH}/libPGSTK.a ${EXTRA_PGS_LIB})
set (HDFEOS_LIB ${CMAKE_SOURCE_DIR}/lib3/EOS/lib/${EOS_ARCH}/libhdfeos.a)
set (HDFEOS5_LIB ${CMAKE_SOURCE_DIR}/lib3/EOS/lib/${EOS_ARCH}/libhe5_hdfeos.a)
set (GCTP_LIB ${CMAKE_SOURCE_DIR}/lib3/EOS/lib/${EOS_ARCH}/libGctp.a)

# where, oh where are your FORTRAN mod files...
set (CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/fortran_modules)

# Need path to third party libraries
link_directories(
  ${CMAKE_SOURCE_DIR}/lib3/lib
)

include_directories (
  ${CMAKE_SOURCE_DIR}/lib3/include
  ${CMAKE_SOURCE_DIR}/include
)

# where we put the binaries
set (CMAKE_INSTALL_PREFIX $ENV{OCSSWROOT}/run)

add_subdirectory(src)

