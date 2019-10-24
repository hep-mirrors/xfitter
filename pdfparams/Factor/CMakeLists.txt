set(TARGET pdfparamFactor)
add_library(${TARGET} MODULE FactorPdfParam.cc)
target_link_libraries(${TARGET} PRIVATE xfitter)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
