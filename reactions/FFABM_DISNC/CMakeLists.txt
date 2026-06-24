#if(NOT QCDNUM_FOUND)
#  message(STATUS "Skipping optional reaction FFABM_DISNC because QCDNUM library was not found")
#  return()
#endif()
# not sure if QCDNUM is needed

set(TARGET reactionFFABM_DISNC)
add_library(${TARGET} MODULE ReactionFFABM_DISNC.cc)
target_link_libraries(${TARGET} PRIVATE xfitter reactionBaseFFABM reactionBaseDISNC)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
install(FILES FFABM_DISNC.yaml DESTINATION ${DEST_YAML}/reactions)
