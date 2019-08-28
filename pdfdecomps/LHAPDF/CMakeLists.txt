if(NOT LHAPDF_FOUND)
  message(STATUS "Skipping optional LHAPDF decomposition because LHAPDF library was not found")
  return()
endif()

set(TARGET pdfdecompLHAPDF)
add_library(${TARGET} MODULE LHAPDFDecomposition.cc)
target_link_libraries(${TARGET} PRIVATE xfitter LHAPDF)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
