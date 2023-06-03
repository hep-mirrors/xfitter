set(TARGET pdfparamChebyschev)
add_library(${TARGET} MODULE Chebyschev_PdfParam.cc)
target_link_libraries(${TARGET} PRIVATE xfitter)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
