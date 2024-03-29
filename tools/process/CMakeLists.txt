if(NOT yaml_FOUND)
  message(STATUS "Skipping optional profiling tools since yaml (C-version) is not found")
  return()
endif()

set(xfitter-process_SOURCES
  xfitter-process.c
  utils.c
  profile.c
  rotate.c
  custom.c
  symmetrize.c
  reweight.c
  scale90to68.c
  )


if (LHAPDF_FOUND)
  list(APPEND xfitter-process_SOURCES interpolation.cc combine.c)
endif()

if (ROOT_FOUND)
  list(APPEND xfitter-process_SOURCES rootplot.cc)
else()
  list(APPEND xfitter-process_SOURCES rootplot_dummy.cc)
endif()

add_executable(xfitter-process
  ${xfitter-process_SOURCES}
  )

target_link_libraries(xfitter-process PRIVATE xfitter pdf2yaml yaml)
  
if(LHAPDF_FOUND)
  target_link_libraries(xfitter-process PRIVATE LHAPDF)
endif()

install(TARGETS xfitter-process DESTINATION ${DEST_BIN})

