#if(APFELgrid_FOUND)
  #return()
#endif()
find_program(apfelgrid-config apfelgrid-config PATHS ${CMAKE_SOURCE_DIR}/deps/applgrid/bin)
if(EXISTS "${apfelgrid-config}")
  set(APFELgrid_FOUND 1)
  execute_process(COMMAND ${apfelgrid-config} --version OUTPUT_VARIABLE APFELgrid_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${apfelgrid-config} --prefix OUTPUT_VARIABLE APFELgrid_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${apfelgrid-config} --libdir OUTPUT_VARIABLE APFELgrid_LIBRARY_DIRS OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${apfelgrid-config} --cxxflags OUTPUT_VARIABLE APFELgrid_COMPILE_OPTIONS OUTPUT_STRIP_TRAILING_WHITESPACE)
  separate_arguments(APFELgrid_COMPILE_OPTIONS UNIX_COMMAND ${APFELgrid_COMPILE_OPTIONS}) #convert space-separated list of compiler flags into a ';'-separated cmake list
  add_library(APFELgrid SHARED IMPORTED)
  if(APPLE)
    set_target_properties(APFELgrid PROPERTIES IMPORTED_LOCATION "${APFELgrid_LIBRARY_DIRS}/libAPFELgrid.dylib")
  endif()
  if(UNIX AND NOT APPLE)
    set_target_properties(APFELgrid PROPERTIES IMPORTED_LOCATION "${APFELgrid_LIBRARY_DIRS}/libAPFELgrid.so")
  endif()
  set_target_properties(APFELgrid PROPERTIES INTERFACE_COMPILE_OPTIONS "${APFELgrid_COMPILE_OPTIONS}")
else()
  set(APFELgrid_FOUND 0)
endif()
if(NOT APFELgrid_FIND_QUIETLY)
  if(APFELgrid_FOUND)
    message(STATUS "Found APFELgrid ${APFELgrid_VERSION}: ${APFELgrid_PREFIX}")
  else()
    message(STATUS "APFELgrid not found")
  endif()
endif()