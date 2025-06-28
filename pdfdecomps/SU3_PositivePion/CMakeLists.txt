set(TARGET pdfdecompSU3_PositivePion)
add_library(${TARGET} MODULE SU3_PositivePionPdfDecomposition.cc)
target_link_libraries(${TARGET} PRIVATE xfitter)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
