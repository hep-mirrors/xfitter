set(TARGET pdfparamABMPvalence)
add_library(${TARGET} MODULE ABMPvalencePdfParam.cc)
target_link_libraries(${TARGET} PRIVATE xfitter)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
