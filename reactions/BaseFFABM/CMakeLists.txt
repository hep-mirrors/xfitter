# not sure if QCDNUM is needed, but...
if(NOT QCDNUM_FOUND)
  message(STATUS "Skipping optional reaction BaseFFABM because QCDNUM library was not found")
  return()
endif()

set(TARGET reactionBaseFFABM)
add_library(${TARGET} SHARED
  ReactionBaseFFABM.cc
  numufcalflux1.F
)
target_link_libraries(${TARGET} PUBLIC 
  xfitter
  ABM
  DIS_HT
  DIS_TMC
  DIS_NUKE
  reactionBaseDISNC
  reactionBaseDISCC
  #cuba
)
target_include_directories(${TARGET} PUBLIC .)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
install(FILES BaseFFABM.yaml DESTINATION ${DEST_YAML}/reactions)
