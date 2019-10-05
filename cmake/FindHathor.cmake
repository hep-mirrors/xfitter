#We recognize two versions of Hathor: 1.5 and 2.0
#In 1.5:
#  both headers and the library are in $HATHOR_DIR
#  the library is libHathor.a
#  does not require LHAPDF
#In 2.0:
#  headers are in $HATHOR_DIR/include, libraries are in $HATHOR_DIR/lib
#  libraries are libHathor.a and libff.a, need to link both
#  requires LHAPDF

#Begin by checking HATHOR_DIR environment variable, it should point to Hathor installation prefix
set(HATHOR_FOUND 0)
if(EXISTS "$ENV{HATHOR_DIR}")
  set(HATHOR_DIR "$ENV{HATHOR_DIR}")
  if(EXISTS "${HATHOR_DIR}/libHathor.a")
    #no lib/ directory, therefore this is 1.5
    set(HATHOR_FOUND 1)
    set(HATHOR_VERSION "1.5")
    add_library(HATHOR STATIC IMPORTED)
    set_target_properties(HATHOR PROPERTIES IMPORTED_LOCATION "${HATHOR_DIR}/libHathor.a")
    set_target_properties(HATHOR PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${HATHOR_DIR}")
  elseif(EXISTS "${HATHOR_DIR}/lib/libHathor.a")
    #found libHathor.a in lib/ directory, therefore this is 2.0
    set(HATHOR_FOUND 1)
    set(HATHOR_VERSION "2.0")
    add_library(HATHOR STATIC IMPORTED)
    set_target_properties(HATHOR PROPERTIES IMPORTED_LOCATION "${HATHOR_DIR}/lib/libHathor.a")
    set_target_properties(HATHOR PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${HATHOR_DIR}/include")
  endif()
endif()
#if HATHOR_DIR is not provided, or if we failed to find Hathor there, search for libHathor.a in standard paths
if(NOT HATHOR_FOUND)
  find_library(HATHOR_IMPORTED_LOCATION "libHathor.a")
  if(EXISTS ${HATHOR_IMPORTED_LOCATION})
    get_filename_component(HATHOR_LIBDIR "${HATHOR_IMPORTED_LOCATION}" DIRECTORY)
    if(EXISTS "${HATHOR_LIBDIR}/libff.a")
      #this is probably Hathor-2.0
      get_filename_component(HATHOR_DIR "${HATHOR_LIBDIR}" DIRECTORY) #HATHOR_DIR=$HATHOR_LIBDIR/..
      #Search for Hathor.h
      find_path(HATHOR_HEADER "Hathor.h" HINTS "${HATHOR_DIR}/include")
      if(EXISTS ${HATHOR_HEADER})
        set(HATHOR_FOUND 1)
        set(HATHOR_VERSION "2.0")
        add_library(HATHOR STATIC IMPORTED)
        set_target_properties(HATHOR PROPERTIES IMPORTED_LOCATION "${HATHOR_IMPORTED_LOCATION}")
        set_target_properties(HATHOR PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${HATHOR_DIR}/include")
      endif()
    else()
      #this is probably Hathor-1.5
      #Search for Hathor.h
      get_filename_component(HATHOR_DIR "${HATHOR_IMPORTED_LOCATION}" DIRECTORY)
      find_path(HATHOR_HEADER "Hathor.h" HINTS "${HATHOR_DIR}")
      if(EXISTS ${HATHOR_HEADER})
        set(HATHOR_FOUND 1)
        set(HATHOR_VERSION "1.5")
        add_library(HATHOR STATIC IMPORTED)
        set_target_properties(HATHOR PROPERTIES IMPORTED_LOCATION "${HATHOR_IMPORTED_LOCATION}")
        set_target_properties(HATHOR PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${HATHOR_DIR}")
      endif()
    endif()
  endif()
endif()

if(NOT HATHOR_FOUND AND HATHOR_FIND_REQUIRED)
  message(FATAL_ERROR "Hathor not found, set HATHOR_DIR to its install prefix")
endif()

if(NOT HATHOR_FIND_QUIETLY)
  if(HATHOR_FOUND)
    message(STATUS "Found Hathor ${HATHOR_VERSION}: ${HATHOR_DIR}")
  else()
    message(STATUS "Hathor not found, set HATHOR_DIR to its install prefix")
  endif()
endif()

#if found Hathor 2.0, make sure LHAPDF is linked
#if found Hathor 2.0 and LHAPDF is not found, disable Hathor
if(HATHOR_FOUND AND (HATHOR_VERSION STREQUAL "2.0"))
  if(NOT DEFINED LHAPDF_FOUND)
    find_package(LHAPDF)
  endif()
  if(LHAPDF_FOUND)
    set_target_properties(HATHOR PROPERTIES INTERFACE_LINK_LIBRARIES "${HATHOR_DIR}/lib/libff.a;LHAPDF")
  else()
    set(HATHOR_FOUND 0)
    message(STATUS "Not using Hathor because LHAPDF, which Hathor-2.0 requires, is not found")
  endif()
endif()
