cmake_minimum_required(VERSION 2.8.12.2)

project(xFitter CXX C Fortran)

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(SEND_ERROR "In-source builds are not allowed. The recommended way to use cmake is:\nmkdir build; cd build; cmake .. -DCMAKE_INSTALL_PREFIX=../install")
endif()

set(xFitter_VERSION_MAJOR 2)
set(xFitter_VERSION_MINOR 2)
set(xFitter_VERSION_PATCH 1)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

#Install paths, relative to install prefix
set(DEST_BIN bin) #executables
set(DEST_LIB lib) #main library
set(DEST_MODULE lib/xfitter) #dynamically-loaded libraries
set(DEST_YAML share/xfitter) #YAML configuration files

#Find packages
#see cmake/FindXXX.cmake for find scripts
find_package(LAPACK REQUIRED)
find_package(QCDNUM REQUIRED)
find_package(yaml-cpp REQUIRED)
find_package(LHAPDF)
find_package(APFEL)
find_package(Ceres)

#Use c++11
if(CMAKE_VERSION VERSION_LESS "3.1")
  #workaround for older versions of cmake
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++11")
else()
  set(CMAKE_CXX_STANDARD 11)
endif()

#Use c preprocessor with fortran
set(CMAKE_Fortran_FLAGS "-cpp")

#For Fortran
#set max line length to 132
#and treat all local variables as static by default (as if SAVE was specified) (without this hf_errlog will crash)
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffixed-line-length-132 -fno-automatic")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -132 -noauto")
else()
  message(FATAL_ERROR "Unknown Fortran compiler ${CMAKE_Fortran_COMPILER_ID}")
endif()

add_subdirectory(src)

#==Dynamically-loaded modules==

#Minimizers
add_subdirectory(minimizers/MINUIT)
if(CERES_FOUND)
  add_subdirectory(minimizers/CERES)
else()
  message(STATUS "Skipping optional CERES minimizer because Ceres library was not found
(If Ceres is actually installed, set Ceres_DIR environment variable the its install prefix)")
endif()

#Evolutions
add_subdirectory(evolutions/QCDNUM)
if(APFEL_FOUND)
  add_subdirectory(evolutions/APFEL)
else()
  message(STATUS "Skipping optional APFEL evolution because APFEL library was not found")
endif()

#Reactions
add_subdirectory(reactions/BaseDISNC)
add_subdirectory(reactions/FFABM_DISNC)
