set(TARGET reactionBaseTBSF)
add_library(${TARGET} SHARED
  ReactionBaseTBSF.cc
)
target_link_libraries(${TARGET} PUBLIC 
  xfitter
  reactionBaseFFABM
  reactionBaseDISNC
  reactionBaseDISCC
  #cuba
)

target_include_directories(${TARGET} PUBLIC .)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
install(FILES BaseTBSF.yaml DESTINATION ${DEST_YAML}/reactions)
