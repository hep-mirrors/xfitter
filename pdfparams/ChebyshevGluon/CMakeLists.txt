set(TARGET pdfparamChebyshevGluon)
add_library(${TARGET} MODULE ChebyshevGluon_PdfParam.cc)
target_link_libraries(${TARGET} PRIVATE xfitter)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})