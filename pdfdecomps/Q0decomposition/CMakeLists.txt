set(TARGET pdfdecompQ0decomposition)
add_library(${TARGET} MODULE Q0decomposition.cc)
target_link_libraries(${TARGET} PRIVATE xfitter)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
