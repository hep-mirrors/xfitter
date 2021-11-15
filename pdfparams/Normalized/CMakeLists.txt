set(TARGET pdfparamNormalized)
add_library(${TARGET} MODULE Normalized_PdfParam.cc)
target_link_libraries(${TARGET} PRIVATE xfitter)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
