message(STATUS "Adding OPENQCDRAD")

set(TARGET evolutionOPENQCDRAD)
add_library(${TARGET} MODULE EvolutionOPENQCDRAD.cc)
target_link_libraries(${TARGET} PRIVATE xfitter ABM)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
install(FILES OPENQCDRAD.yaml DESTINATION ${DEST_YAML}/evolutions)
