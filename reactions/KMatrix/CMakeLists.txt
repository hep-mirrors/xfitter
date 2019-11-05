set(TARGET reactionKMatrix)
add_library(${TARGET} MODULE ReactionKMatrix.cc)
target_link_libraries(${TARGET} PRIVATE xfitter)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
