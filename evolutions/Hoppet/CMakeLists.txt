if(NOT Hoppet_FOUND)
  message(STATUS "Skipping optional evolution Hoppet since Hoppet library was not found")
  return()
endif()

message(STATUS "Adding Hoppet")

set(TARGET evolutionHoppet)
add_library(${TARGET} MODULE EvolutionHoppet.cc)
target_link_libraries(${TARGET} PRIVATE xfitter Hoppet)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
install(FILES Hoppet.yaml DESTINATION ${DEST_YAML}/evolutions)
