set(TARGET pdfparamNegativeGluon)
add_library(${TARGET} MODULE NegativeGluonPdfParam.cc)
target_link_libraries(${TARGET} PRIVATE xfitter)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
