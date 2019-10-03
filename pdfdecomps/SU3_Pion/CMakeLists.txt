set(TARGET pdfdecompSU3_Pion)
add_library(${TARGET} MODULE SU3_PionPdfDecomposition.cc)
target_link_libraries(${TARGET} PRIVATE xfitter)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
