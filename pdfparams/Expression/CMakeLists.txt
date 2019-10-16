set(TARGET pdfparamExpression)
add_library(${TARGET} MODULE ExpressionPdfParam.cc)
target_link_libraries(${TARGET} PRIVATE xfitter)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
