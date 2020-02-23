set(TARGET pdfparamPolySqrt)
add_library(${TARGET} MODULE PolySqrtPdfParam.cc)
target_link_libraries(${TARGET} PRIVATE xfitter)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
