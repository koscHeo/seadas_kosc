cmake_minimum_required(VERSION 2.6)

project(OCSSW)

#add_executable(val_extract val_extract.c)
#target_link_libraries(val_extract hdfutils genutils)


#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include)

#set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_library(phash phash.c)
add_library(shash shash.c)
add_library(pqueue pqueue-ll.c)

add_library(argpar argpar.c argpar-help.c)

add_library(valextract val_extract.c)
target_link_libraries(valextract argpar pqueue m netcdf ${HDF5_LIBRARIES})

add_library(sstvalextract sstval_extract.c)
target_link_libraries(sstvalextract valextract argpar)

add_executable(val_extract val_extract-main.c)
target_link_libraries(val_extract valextract argpar ${HDF5_LIBRARIES})

add_executable(sstval_extract sstval_extract-main.c)
target_link_libraries(sstval_extract sstvalextract argpar)

install(TARGETS val_extract DESTINATION bin/${OCSSW_ARCH})
install(TARGETS sstval_extract DESTINATION bin/${OCSSW_ARCH})

