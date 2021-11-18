if(CMAKE_DISABLE_FIND_PACKAGE_ROOT)
  set(ROOT_FOUND 0)
  message(STATUS "Skipping disabled package ROOT")
  return()
else()

set(ROOT_HINTS)
if(EXISTS $ENV{ROOT_DIR})
  set(ROOT_DIR $ENV{ROOT_DIR})
  list(APPEND ROOT_HINTS ${ROOT_DIR})
endif()
set(root-config_HINTS)
if(EXISTS ${ROOT_DIR})
  #if ROOT_DIR was given, search for root-config there too
  list(APPEND ${root-config_HINTS} ${ROOT_DIR} ${ROOT_DIR}/bin)
endif()
#if root-config is available, run "root-config --prefix" and use it to find ROOTConfig.cmake
find_program(root-config root-config HINTS ${root-config_HINTS})
if(EXISTS ${root-config})
  execute_process(COMMAND ${root-config} --prefix OUTPUT_VARIABLE ROOT_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE)
  list(APPEND ROOT_HINTS ${ROOT_PREFIX})
  # also store C++
  execute_process(COMMAND ${root-config} --cflags OUTPUT_VARIABLE ROOT_CFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()

if(NOT CMAKE_VERSION VERSION_LESS 3.0.0)
  find_package(ROOT QUIET HINTS ${ROOT_HINTS})
endif()
if(ROOT_FOUND)
  #Use variables provided by ROOTConfig.cmake
  #Remove some extra flags from ROOT_DEFINITIONS
  list(REMOVE_ITEM ROOT_DEFINITIONS "-fPIC" "-std=c++11" "-W" "-Wall" "-Woverloaded-virtual")
elseif(EXISTS ${root-config})
  #Fall back to root-config
  set(ROOT_FOUND 1)
  execute_process(COMMAND ${root-config} --version OUTPUT_VARIABLE ROOT_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${root-config} --prefix OUTPUT_VARIABLE ROOT_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${root-config} --libs OUTPUT_VARIABLE ROOT_LIBFLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${root-config} --libdir OUTPUT_VARIABLE ROOT_LIBRARY_DIRS OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${root-config} --incdir OUTPUT_VARIABLE ROOT_INCLUDE_DIRS OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${root-config} --cflags OUTPUT_VARIABLE ROOT_DEFINITIONS OUTPUT_STRIP_TRAILING_WHITESPACE)
  #parse libraries from flags returned by root-config --libs
  #convert lib flags to a cmake list
  separate_arguments(ROOT_LIBFLAGS)
  set(ROOT_LIBRARIES)
  set(CHAR_l "l") #Some hack to prevent variable expansion if a variable with name "l" is defined somewhere
  foreach(FLAG IN LISTS ROOT_LIBFLAGS)
    string(SUBSTRING "${FLAG}" 1 1 FLAGTYPE) #get second character, for libs like "-lCore" second character is "l"
    if(FLAGTYPE STREQUAL CHAR_l) #skip flags that are not libraries
      string(SUBSTRING "${FLAG}" 2 -1 LIBNAME)
      find_library(ROOT_${LIBNAME}_PATH "${LIBNAME}" HINTS ${ROOT_LIBRARY_DIRS})
      if(EXISTS "${ROOT_${LIBNAME}_PATH}")
        list(APPEND ROOT_LIBRARIES "${ROOT_${LIBNAME}_PATH}")
      else()
        message(WARNING "Failed to find ROOT component library ${LIBNAME}, this might cause problems when linking ROOT")
      endif()
    endif()
  endforeach()
  #Remove some extra flags from ROOT_DEFINITIONS
  separate_arguments(ROOT_DEFINITIONS)
  list(REMOVE_ITEM ROOT_DEFINITIONS "-std=c++11" "-I${ROOT_INCLUDE_DIRS}")
endif()
endif()

if(ROOT_FOUND)
  #Create a ROOT imported target
  add_library(ROOT SHARED IMPORTED)
  list(GET ROOT_LIBRARIES 0 ROOT_IMPORTED_LOCATION)
  set_target_properties(ROOT PROPERTIES IMPORTED_LOCATION "${ROOT_IMPORTED_LOCATION}")
  set_target_properties(ROOT PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${ROOT_INCLUDE_DIRS}")
  set_target_properties(ROOT PROPERTIES INTERFACE_COMPILE_OPTIONS "${ROOT_DEFINITIONS}")
  set_target_properties(ROOT PROPERTIES INTERFACE_LINK_LIBRARIES "${ROOT_LIBRARIES}")
  set_target_properties(ROOT PROPERTIES INTERFACE_LINK_DIRECTORIES "${ROOT_LIBRARY_DIRS}")
  message(STATUS "Found ROOT ${ROOT_VERSION}: ${ROOT_DIR}")
else()
  message(STATUS "ROOT not found
(If ROOT is actualy installed, make sure root-config is in PATH, or, alternatively, set the ROOT_DIR environment variable to ROOT's install prefix)")
endif()
