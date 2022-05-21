set(TARGET pdfdecompPion_FF_BC)
add_library(${TARGET} MODULE Pion_FF_BC.cc)
target_link_libraries(${TARGET} PRIVATE xfitter)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
