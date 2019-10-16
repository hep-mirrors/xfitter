set(TARGET pdfdecompUvDvUbarDbarSSbar)
add_library(${TARGET} MODULE UvDvUbarDbarSSbarPdfDecomposition.cc)
target_link_libraries(${TARGET} PRIVATE xfitter)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
