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
find_package(APPLgrid)
find_package(APFEL)
find_package(APFELxx)
find_package(Ceres)

#Use c++11
if(CMAKE_VERSION VERSION_LESS "3.1")
  #workaround for older versions of cmake
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
  set(CMAKE_CXX_STANDARD 11)
  set(CMAKE_CXX_STANDARD_REQUIRED YES)
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
#ABM is a shared library used by the FFABM_DISNC and FFABM_DISCC reactions
add_subdirectory(src/ABM)

#==Dynamically-loaded modules==

#Minimizers
add_subdirectory(minimizers/MINUIT)
add_subdirectory(minimizers/CERES)

#Parameterisations
add_subdirectory(pdfparams/HERAPDF)
add_subdirectory(pdfparams/Factor)
add_subdirectory(pdfparams/NegativeGluon)
add_subdirectory(pdfparams/ABMPgluon)
add_subdirectory(pdfparams/ABMPsea)
add_subdirectory(pdfparams/ABMPvalence)
add_subdirectory(pdfparams/Expression)
add_subdirectory(pdfparams/PolySqrt)

#Decompositions
add_subdirectory(pdfdecomps/UvDvUbarDbarS)
add_subdirectory(pdfdecomps/LHAPDF)
add_subdirectory(pdfdecomps/SU3_Pion)
add_subdirectory(pdfdecomps/UvDvUbarDbarSSbar)

#Evolutions
add_subdirectory(evolutions/QCDNUM)
add_subdirectory(evolutions/APFEL)
add_subdirectory(evolutions/APFELxx)
add_subdirectory(evolutions/FlipCharge)
add_subdirectory(evolutions/FlipUD)
add_subdirectory(evolutions/LHAPDF)

#Reactions
add_subdirectory(reactions/BaseDISNC)
add_subdirectory(reactions/BaseDISCC)
add_subdirectory(reactions/FFABM_DISNC)
add_subdirectory(reactions/RT_DISNC)
add_subdirectory(reactions/FONLL_DISNC)
add_subdirectory(reactions/FONLL_DISCC)
add_subdirectory(reactions/FFABM_DISCC)
add_subdirectory(reactions/AFB)
add_subdirectory(reactions/APPLgrid)
add_subdirectory(reactions/Fractal_DISNC)
