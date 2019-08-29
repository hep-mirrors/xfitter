#Create a ROOT imported target using variables provided by find_package(ROOT)
if(NOT ROOT_FOUND)
  message(STATUS "ROOT not found")
  return()
endif()
message(STATUS "Found ROOT ${ROOT_VERSION}")

#Remove some extra flags from ROOT_DEFINITIONS
list(REMOVE_ITEM ROOT_DEFINITIONS "-fPIC" "-std=c++11" "-W" "-Wall" "-Woverloaded-virtual")

add_library(ROOT SHARED IMPORTED)
list(GET ROOT_LIBRARIES 0 ROOT_IMPORTED_LOCATION)
set_target_properties(ROOT PROPERTIES IMPORTED_LOCATION "${ROOT_IMPORTED_LOCATION}")
set_target_properties(ROOT PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${ROOT_INCLUDE_DIRS}")
set_target_properties(ROOT PROPERTIES INTERFACE_COMPILE_OPTIONS "${ROOT_DEFINITIONS}")
set_target_properties(ROOT PROPERTIES INTERFACE_LINK_LIBRARIES "${ROOT_LIBRARIES}")
