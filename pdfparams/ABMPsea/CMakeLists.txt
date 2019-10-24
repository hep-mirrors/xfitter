set(TARGET pdfparamABMPsea)
add_library(${TARGET} MODULE ABMPseaPdfParam.cc)
target_link_libraries(${TARGET} PRIVATE xfitter)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
