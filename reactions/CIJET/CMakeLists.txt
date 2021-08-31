set(TARGET reactionCIJET)
add_library(${TARGET} MODULE
  src/ReactionCIJET.cc
  src/CIJETReader.cc
  src/cijet.f
)
target_link_libraries(${TARGET} PRIVATE xfitter)
target_include_directories(${TARGET} PUBLIC include)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
install(FILES yaml/parameters.yaml DESTINATION ${DEST_YAML}/reactions RENAME ReactionCIJET.yaml)

