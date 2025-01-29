if(NOT HOPPET_FOUND)
  message(STATUS "Skipping optional evolution HOPPET since HOPPET library was not found")
  return()
endif()

message(STATUS "Adding HOPPET")

set(TARGET evolutionHOPPET)
add_library(${TARGET} MODULE EvolutionHOPPET.cc)
target_link_libraries(${TARGET} PRIVATE xfitter HOPPET)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
install(FILES HOPPET.yaml DESTINATION ${DEST_YAML}/evolutions)
