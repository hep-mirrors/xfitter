set(TARGET reactionEFT)
add_library(${TARGET} MODULE 
ReactionEFT.cc
EFTReader.cc
Vec.cc
)
target_link_libraries(${TARGET} PRIVATE xfitter)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
install(FILES EFT.yaml DESTINATION ${DEST_YAML}/reactions)
