set(TARGET pdfdecompUvDvUbarDbarS)
add_library(${TARGET} MODULE UvDvUbarDbarS.cc)
target_link_libraries(${TARGET} PRIVATE xfitter)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
