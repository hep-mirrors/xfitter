if(NOT Hoppet_FOUND)
  message(STATUS "Skipping optional evolution Hoppet since Hoppet library was not found")
  return()
endif()

set(TARGET Reaction_DISNC_Hoppet)
add_library(${TARGET} MODULE Reaction_DISNC_Hoppet.cc)
target_link_libraries(${TARGET} PRIVATE xfitter Hoppet reactionBaseDISNC)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
