find_program(hoppet-config hoppet-config PATHS ${CMAKE_SOURCE_DIR}/deps/hoppet/bin)
if(EXISTS "${hoppet-config}")
  set(HOPPET_FOUND 1)
  execute_process(COMMAND ${hoppet-config} --prefix OUTPUT_VARIABLE HOPPET_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${hoppet-config} --libs   OUTPUT_VARIABLE HOPPET_LIBRARY OUTPUT_STRIP_TRAILING_WHITESPACE)
  add_library(HOPPET SHARED IMPORTED)

  set_target_properties(HOPPET PROPERTIES IMPORTED_LOCATION "${HOPPET_PREFIX}/lib/libhoppet_v1.a")

  set_target_properties(HOPPET PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${HOPPET_PREFIX}/include/")

  #libhoppet_v1.a is a static Fortran archive, so the Fortran runtime it depends
  #on (libgfortran and friends) is not pulled in automatically when the modules
  #that use it are linked by the C++ compiler. hoppet-config --libs does not
  #report those either, so take them from CMake's Fortran compiler introspection.
  set_target_properties(HOPPET PROPERTIES
    INTERFACE_LINK_LIBRARIES  "${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES}"
    INTERFACE_LINK_DIRECTORIES "${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}")
else()
  set(HOPPET_FOUND 0)
endif()
if(NOT HOPPET_FOUND AND HOPPET_FIND_REQUIRED)
  message(FATAL_ERROR "HOPPET not found")
endif()
if(NOT HOPPET_FIND_QUIETLY)
  if(HOPPET_FOUND)
    message(STATUS "Found HOPPET ${HOPPET_VERSION}: ${HOPPET_PREFIX}")
  else()
    message(STATUS "HOPPET not found")
  endif()
endif()
