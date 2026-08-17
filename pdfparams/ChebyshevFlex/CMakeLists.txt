set(TARGET pdfparamChebyshevFlex)
add_library(${TARGET} MODULE ChebyshevFlex_PdfParam.cc)
target_link_libraries(${TARGET} PRIVATE xfitter)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})