if(NOT ROOT_FOUND)
  message(STATUS "Skipping optional reaction ReactionHVQMNR_LHCb_7TeV_charm because ROOT library was not found")
  return()
endif()

set(TARGET reactionHVQMNR_LHCb_7TeV_charm)
add_library(${TARGET} MODULE ReactionHVQMNR_LHCb_7TeV_charm.cc)
target_link_libraries(${TARGET} PRIVATE reactionBaseHVQMNR)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
