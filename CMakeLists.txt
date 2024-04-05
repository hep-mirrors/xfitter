cmake_minimum_required(VERSION 2.8.12.2)
set(CMAKE_MACOSX_RPATH 1)

project(xFitter CXX C Fortran)

if(PROJECT_SOURCE_DIR STREQUAL PROJECT_BINARY_DIR)
  message(SEND_ERROR "In-source builds are not allowed. The recommended way to use cmake is:\nmkdir build; cd build; cmake .. -DCMAKE_INSTALL_PREFIX=../install")
endif()

set(xFitter_VERSION_MAJOR 2)
set(xFitter_VERSION_MINOR 2)
set(xFitter_VERSION_PATCH 1)
set(xFitter_VERSION "${xFitter_VERSION_MAJOR}.${xFitter_VERSION_MINOR}.${xFitter_VERSION_PATCH}")
message(STATUS "Configuring xFitter ${xFitter_VERSION}")

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

#Install paths, relative to install prefix
set(DEST_BIN bin) #executables
set(DEST_LIB lib) #main library
set(DEST_MODULE lib/xfitter) #dynamically-loaded libraries
set(DEST_YAML share/xfitter) #YAML configuration files

#Find packages
#see cmake/FindXXX.cmake for find scripts
#TODO: rewrite find modules to use cmake caching to speed up subsequent configuration
find_package(GSL REQUIRED)
find_package(LAPACK REQUIRED)
find_package(yaml-cpp REQUIRED)
find_package(ROOT)
find_package(QCDNUM)
find_package(LHAPDF)
find_package(APPLgrid)
find_package(APFELgrid)
find_package(APFEL)
find_package(APFELxx)
find_package(Ceres)
find_package(Hathor)
find_package(DYTurbo)
find_package(yaml)
find_package(ZLIB)
find_package(PineAPPL)
find_package(Ploughshare)

set(CMAKE_CXX_FLAGS "-O ${CMAKE_CXX_FLAGS}")
set(CMAKE_Fortran_FLAGS "-O ${CMAKE_Fortran_FLAGS}")
set(CMAKE_C_FLAGS "-O ${CMAKE_C_FLAGS}")

if(OPENMP)
  set(CMAKE_CXX_FLAGS "-fopenmp ${CMAKE_CXX_FLAGS}")
  set(CMAKE_Fortran_FLAGS "-fopenmp ${CMAKE_Fortran_FLAGS}")
  set(CMAKE_C_FLAGS "-fopenmp ${CMAKE_C_FLAGS}")
endif()

if(ROOT_FOUND)
  message("ROOT FOUND")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS}")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${ROOT_CC_FLAGS}")
else()

  #Use c++11
  if(CMAKE_VERSION VERSION_LESS "3.1")
    #workaround for older versions of cmake
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  else()
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED YES)
  endif()
endif()

#Use c preprocessor with fortran

if(UNIX AND NOT APPLE)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -Wno-argument-mismatch")
endif()

if(APPLE)
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -fallow-argument-mismatch")
endif()


if (CMAKE_MAJOR_VERSION VERSION_GREATER_EQUAL "3")
if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "10")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -fallow-argument-mismatch")
endif()
endif()


#For Fortran
#set max line length to 132
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffixed-line-length-132")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -132")
else()
  message(FATAL_ERROR "Unknown Fortran compiler ${CMAKE_Fortran_COMPILER_ID}")
endif()

#For Debug builds, enable checks of array bounds and trap floating-point errors
if(CMAKE_BUILD_TYPE STREQUAL "Debug")
  if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fbounds-check -ffpe-trap='invalid,overflow,zero'")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ggdb")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ggdb")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ggdb")
  elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -check bounds -fpe0")
  endif()
endif()

# For semaphore in multithreading
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -lrt")

#When linking, set rpath to install directories
#This eliminates the need to add xfitter libraries to LD_LIBRARY_PATH
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${DEST_LIB};${CMAKE_INSTALL_PREFIX}/${DEST_MODULE}")
set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)

# add the to the install RPATH directories for external libraries which are outside the build tree
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

#Suppress install messages
#Only works since cmake 3.1
#LAZY   to show "-- Installing", but not "-- Up-to-date"
#ALWAYS to show both
#NEVER  to hide both
set(CMAKE_INSTALL_MESSAGE NEVER)

include_directories(${PROJECT_SOURCE_DIR}/include)

add_subdirectory(src)
#ABM is a shared library used by the FFABM_DISNC and FFABM_DISCC reactions
add_subdirectory(src/ABM)

#==Dynamically-loaded modules==

#Minimizers
add_subdirectory(minimizers/MINUIT)
add_subdirectory(minimizers/CERES)

#Parameterisations
add_subdirectory(pdfparams/Pion_FF)
add_subdirectory(pdfparams/HERAPDF)
add_subdirectory(pdfparams/Normalized)
add_subdirectory(pdfparams/Factor)
add_subdirectory(pdfparams/NegativeGluon)
add_subdirectory(pdfparams/ABMPgluon)
add_subdirectory(pdfparams/ABMPsea)
add_subdirectory(pdfparams/ABMPvalence)
add_subdirectory(pdfparams/Expression)
add_subdirectory(pdfparams/PolySqrt)
add_subdirectory(pdfparams/Chebyschev)

#Decompositions
add_subdirectory(pdfdecomps/Pion_FF_B_C)
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
add_subdirectory(reactions/NC_SIA)
add_subdirectory(reactions/BaseDISNC)
add_subdirectory(reactions/BaseDISCC)
add_subdirectory(reactions/FFABM_DISNC)
add_subdirectory(reactions/RT_DISNC)
add_subdirectory(reactions/FONLL_DISNC)
add_subdirectory(reactions/FONLL_DISCC)
add_subdirectory(reactions/FFABM_DISCC)
add_subdirectory(reactions/AFB)
add_subdirectory(reactions/APPLgrid)
add_subdirectory(reactions/APFELgrid)
add_subdirectory(reactions/Fractal_DISNC)
add_subdirectory(reactions/KFactor)
add_subdirectory(reactions/KMatrix)
add_subdirectory(reactions/TensorPomeron)
add_subdirectory(reactions/BaseHVQMNR)
add_subdirectory(reactions/cbdiff)
add_subdirectory(reactions/HVQMNR_LHCb_7TeV_beauty)
add_subdirectory(reactions/HVQMNR_LHCb_7TeV_charm)
add_subdirectory(reactions/Hathor)
add_subdirectory(reactions/HathorMSR)
add_subdirectory(reactions/HathorSingleTop)
add_subdirectory(reactions/fastNLO)
add_subdirectory(reactions/CIJET)
add_subdirectory(reactions/DYTurbo)
add_subdirectory(reactions/PineAPPL)

#Tools
add_subdirectory(tools/draw/src)
add_subdirectory(pdf2yaml)
add_subdirectory(tools/process)
