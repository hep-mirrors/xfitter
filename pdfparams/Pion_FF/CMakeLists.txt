set(TARGET pdfparamPion_FF)
add_library(${TARGET} MODULE Pion_FF_PdfParam.cc)
target_link_libraries(${TARGET} PRIVATE xfitter)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
