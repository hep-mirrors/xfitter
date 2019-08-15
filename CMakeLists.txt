cmake_minimum_required(VERSION 2.8.12.2)

project(xFitter CXX C Fortran)

set(xFitter_VERSION_MAJOR 2)
set(xFitter_VERSION_MINOR 2)
set(xFitter_VERSION_PATCH 1)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

#Install paths
set(DEST_BIN bin)
set(DEST_LIB lib)
set(DEST_MODULE lib/xfitter)

#Find packages
#see cmake/FindXXX.cmake for find scripts
find_package(LAPACK REQUIRED)
find_package(QCDNUM REQUIRED)
find_package(yaml-cpp REQUIRED)
find_package(LHAPDF)

#Use c++11
if(CMAKE_VERSION VERSION_LESS "3.1")
  #workaround for older versions of cmake
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
else()
  set(CMAKE_CXX_STANDARD 11)
endif()

#Use c preprocessor with fortran
set(CMAKE_Fortran_FLAGS "-cpp")

#For Fortran, set max line length to 132
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffixed-line-length-132")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -132")
else()
  message(FATAL_ERROR "Unknown Fortran compiler ${CMAKE_Fortran_COMPILER_ID}")
endif()

add_subdirectory(src)
#add_subdirectory(pdfparams/HERAPDF)
