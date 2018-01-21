cmake_minimum_required (VERSION 2.6)

project (HEALPIX)

# set this to Release or Debug
# note: rerun cmake if you change the environment var OCSSW_DEBUG
if ($ENV{OCSSW_DEBUG})
  set (CMAKE_BUILD_TYPE Debug)
else($ENV{OCSSW_DEBUG})
  set (CMAKE_BUILD_TYPE Release)
endif($ENV{OCSSW_DEBUG})

# set this to print out the compile lines
# can also use "make VERBOSE=1" when building
#set (CMAKE_VERBOSE_MAKEFILE TRUE)

enable_language (C)
enable_language (CXX)
enable_language (Fortran)

# make the build static
set (CMAKE_EXE_LINKER_FLAGS "-static")
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

# Set up C and CXX flags
set (CMAKE_C_FLAGS_RELEASE "-O3")
if (DEFINED ENV{GCC_TUNE})
  set (CMAKE_C_FLAGS_RELEASE "-O3 $ENV{GCC_TUNE}")
endif(DEFINED ENV{GCC_TUNE}) 
set (CMAKE_C_FLAGS_DEBUG "-g3 -gdwarf-2 -Wall -Wno-unused-variable -Wno-unused-label")
set (CMAKE_CXX_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE})
set (CMAKE_CXX_FLAGS_DEBUG ${CMAKE_C_FLAGS_DEBUG})
set (CMAKE_C_FLAGS "-DLINUX -Dlinux -Df2cFortran -DSYSV -ffloat-store")
set (CMAKE_CXX_FLAGS ${CMAKE_C_FLAGS})

# Set up FORTRAN flags
if (${Fortran_COMPILER_NAME} MATCHES "gfortran")
  set (CMAKE_Fortran_FLAGS "${CMAKE_C_FLAGS} -DGFORTRAN -ffixed-line-length-none -fconvert=big-endian -fno-automatic -fno-range-check -cpp")
  set (CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -funroll-all-loops")
  set (CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG}")
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
  if(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")
    # 64bit bits...
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -m64")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64")
    if (${Fortran_COMPILER_NAME} MATCHES "ifort")
      set (CMAKE_Fortran_FLAGS "-m64 ${CMAKE_Fortran_FLAGS}")
    endif(${Fortran_COMPILER_NAME} MATCHES "ifort")
  endif(${CMAKE_SYSTEM_PROCESSOR} MATCHES "x86_64")
elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  # Mac bits
  set (CMAKE_C_FLAGS "-DMACINTOSH ${CMAKE_C_FLAGS}")
  set (CMAKE_CXX_FLAGS "-DMACINTOSH ${CMAKE_CXX_FLAGS}")
  set (CMAKE_Fortran_FLAGS "-DMACINTOSH ${CMAKE_Fortran_FLAGS}")
  set (CMAKE_EXE_LINKER_FLAGS "-Wl,-no_compact_unwind")
  set (EXTRA_OS_LIBS gcc_eh)
# set (EXTRA_PGS_LIB /usr/lib/libiconv.dylib)
else(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  # Miscellaneous - not sure what we are...
  message ("CMAKE_SYSTEM_NAME: " ${CMAKE_SYSTEM_NAME})
  message ("System not recognized, using other")
endif(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

# where, oh where are your FORTRAN mod files...
set (CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/fortran_modules)

# Need path to third party libraries
link_directories    ($ENV{LIB3_LIB} $ENV{PGSLIB})
include_directories ($ENV{LIB3_INC} $ENV{PGSINC})

# Next build directory
message(STATUS ${CMAKE_SOURCE_DIR})
add_subdirectory($ENV{srcdir}/src/C/subs)
add_subdirectory($ENV{srcdir}/src/f90/mod)
add_subdirectory($ENV{srcdir}/src/f90/sharp)
#add_subdirectory($ENV{srcdir}/src/f90/lib)
