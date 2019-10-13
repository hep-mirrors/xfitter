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
endif()

find_package(ROOT HINTS ${ROOT_HINTS})
if(ROOT_FOUND)
  #Create a ROOT imported target using variables provided by ROOTConfig.cmake

  #Remove some extra flags from ROOT_DEFINITIONS
  list(REMOVE_ITEM ROOT_DEFINITIONS "-fPIC" "-std=c++11" "-W" "-Wall" "-Woverloaded-virtual")
  add_library(ROOT SHARED IMPORTED)
  list(GET ROOT_LIBRARIES 0 ROOT_IMPORTED_LOCATION)
  set_target_properties(ROOT PROPERTIES IMPORTED_LOCATION "${ROOT_IMPORTED_LOCATION}")
  set_target_properties(ROOT PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${ROOT_INCLUDE_DIRS}")
  set_target_properties(ROOT PROPERTIES INTERFACE_COMPILE_OPTIONS "${ROOT_DEFINITIONS}")
  set_target_properties(ROOT PROPERTIES INTERFACE_LINK_LIBRARIES "${ROOT_LIBRARIES}")
endif()

if(ROOT_FOUND)
  message(STATUS "Found ROOT ${ROOT_VERSION}: ${ROOT_DIR}")
else()
  message(STATUS "ROOT not found
(If ROOT is actualy installed, make sure root-config is in PATH, or, alternatively, set the ROOT_DIR environment variable to ROOT's install prefix)")
endif()
