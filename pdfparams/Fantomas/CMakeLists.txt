set(TARGET pdfparamFantomas)
add_library(${TARGET} MODULE Fantomas_PdfParam.cc metamorphCollection.cc CToWrapper.cc)
target_link_libraries(${TARGET} PRIVATE xfitter)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
