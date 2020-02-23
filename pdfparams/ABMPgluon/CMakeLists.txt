set(TARGET pdfparamABMPgluon)
add_library(${TARGET} MODULE ABMPgluonPdfParam.cc)
target_link_libraries(${TARGET} PRIVATE xfitter)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
