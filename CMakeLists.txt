cmake_minimum_required(VERSION 2.8.12.2)

project(xFitter CXX Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

#Install paths
set(DEST_BIN bin)
set(DEST_LIB lib)
set(DEST_MODULE lib/xfitter)

message("$ENV{CXXFLAGS}")
message("${CMAKE_CXX_FLAGS}")
#Use c++11
if(CMAKE_VERSION VERSION_LESS "3.1")
  #workaround for older versions of cmake
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
else()
  set(CMAKE_CXX_STANDARD 11)
endif()

#Use c preprocessor with fortran
set(CMAKE_Fortran_FLAGS "-cpp")

if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffixed-line-length-132")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -132")
else()
  message(FATAL_ERROR "Unknown Fortran compiler ${CMAKE_Fortran_COMPILER_ID}")
endif()

find_package(PkgConfig)

find_package(LHAPDF)

add_subdirectory(src)
#add_subdirectory(pdfparams/HERAPDF)
