set(TARGET reactionKRunning)
add_library(${TARGET} MODULE ReactionKRunning.cc)
target_link_libraries(${TARGET} PRIVATE xfitter)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
