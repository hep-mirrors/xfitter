if(NOT HOPPET_FOUND)
  message(STATUS "Skipping optional evolution HOPPET since HOPPET library was not found")
  return()
endif()

set(TARGET reactionHOPPET_DISCC)
add_library(${TARGET} MODULE ReactionHOPPET_DISCC.cc)
target_link_libraries(${TARGET} PRIVATE xfitter HOPPET reactionBaseDISCC)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
install(FILES HOPPET_DISCC.yaml DESTINATION ${DEST_YAML}/reactions)
