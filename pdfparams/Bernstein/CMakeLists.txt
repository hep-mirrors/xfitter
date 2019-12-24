set(TARGET pdfparamBernstein)
add_library(${TARGET} MODULE Bernstein_PdfParam.cc)
target_link_libraries(${TARGET} PRIVATE xfitter)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
