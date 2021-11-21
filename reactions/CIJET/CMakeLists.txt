set(TARGET reactionCIJET)
add_library(${TARGET} MODULE
  ReactionCIJET.cc
  CIJETReader.cc
  cijet.f
)
target_link_libraries(${TARGET} PRIVATE xfitter)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
install(FILES CIJET.yaml DESTINATION ${DEST_YAML}/reactions)
