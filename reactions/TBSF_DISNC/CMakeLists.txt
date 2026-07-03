#if(NOT QCDNUM_FOUND)
#  message(STATUS "Skipping optional reaction TBSF_DISNC because QCDNUM library was not found")
#  return()
#endif()
# not sure if QCDNUM is needed

set(TARGET reactionTBSF_DISNC)
add_library(${TARGET} MODULE ReactionTBSF_DISNC.cc)
target_link_libraries(${TARGET} PRIVATE xfitter reactionBaseTBSF reactionBaseFFABM reactionBaseDISNC)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
install(FILES TBSF_DISNC.yaml DESTINATION ${DEST_YAML}/reactions)
