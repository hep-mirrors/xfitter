set(TARGET pdfdecompMSHT20)
add_library(${TARGET} MODULE MSHT20PdfDecomposition.cc)
target_link_libraries(${TARGET} PRIVATE xfitter)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
