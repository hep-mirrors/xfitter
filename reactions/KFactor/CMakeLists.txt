set(TARGET reactionKFactor)
add_library(${TARGET} MODULE ReactionKFactor.cc)
target_link_libraries(${TARGET} PRIVATE xfitter)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
