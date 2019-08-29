if(NOT ROOT_FOUND)
  message(STATUS "Skipping optional reaction cbdiff because ROOT library was not found")
  return()
endif()

set(TARGET reactioncbdiff)
add_library(${TARGET} MODULE Reactioncbdiff.cc)
target_link_libraries(${TARGET} PRIVATE reactionBaseHVQMNR)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
