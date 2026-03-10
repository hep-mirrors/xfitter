set(TARGET pdfdecompSU3_PositivePion_gluon)
add_library(${TARGET} MODULE SU3_PositivePion_gluonPdfDecomposition.cc)
target_link_libraries(${TARGET} PRIVATE xfitter)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
