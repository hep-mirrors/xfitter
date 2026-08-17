set(TARGET pdfdecompMSHT)
add_library(${TARGET} MODULE MSHT_Decomposition.cc)
target_link_libraries(${TARGET} PRIVATE xfitter)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})