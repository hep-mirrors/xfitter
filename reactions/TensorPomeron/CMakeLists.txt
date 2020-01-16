set(TARGET reactionTensorPomeron)
add_library(${TARGET} MODULE ReactionTensorPomeron.cc)
target_link_libraries(${TARGET} PRIVATE xfitter)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
