#if(NOT QCDNUM_FOUND)
#  message(STATUS "Skipping optional reaction FFABM_DISNC because QCDNUM library was not found")
#  return()
#endif()

#set(TARGET reactionBaseFFABM)
#add_library(${TARGET} SHARED ReactionBaseFFABM.cc)
#target_link_libraries(${TARGET} PRIVATE xfitter ABM DIS_HT DIS_TMC DIS_NUKE reactionBaseDISNC reactionBaseDISCC)
#target_include_directories(${TARGET} PRIVATE ${CMAKE_SOURCE_DIR}/reactions/BaseDISNC ${CMAKE_SOURCE_DIR}/reactions/BaseDISCC)
#target_include_directories(${TARGET} PUBLIC .)
#install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
##install(FILES BaseFFABM.yaml DESTINATION ${DEST_YAML}/reactions)

#add_library(reactionBaseFFABM STATIC ReactionBaseFFABM.cc)