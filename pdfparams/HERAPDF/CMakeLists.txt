set(TARGET pdfparamHERAPDF)
add_library(${TARGET} MODULE HERAPDF_PdfParam.cc)
target_link_libraries(${TARGET} PRIVATE xfitter)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
