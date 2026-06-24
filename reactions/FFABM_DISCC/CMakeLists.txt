#if(NOT QCDNUM_FOUND)
#  message(STATUS "Skipping optional reaction FFABM_DISNC because QCDNUM library was not found")
#  return()
#endif()
# not sure if QCDNUM is needed

set(TARGET reactionFFABM_DISCC)
add_library(${TARGET} MODULE ReactionFFABM_DISCC.cc)
target_link_libraries(${TARGET} PRIVATE xfitter reactionBaseFFABM reactionBaseDISCC)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
install(FILES FFABM_DISCC.yaml DESTINATION ${DEST_YAML}/reactions)
