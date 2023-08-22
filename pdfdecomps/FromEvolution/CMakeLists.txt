set(TARGET pdfdecompFromEvolution)
add_library(${TARGET} MODULE FromEvolutionDecomposition.cc)
target_link_libraries(${TARGET} PRIVATE xfitter)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
