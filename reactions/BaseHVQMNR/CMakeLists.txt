if(NOT ROOT_FOUND)
  message(STATUS "Skipping optional reaction BaseHVQMNR because ROOT library was not found")
  return()
endif()

set(TARGET reactionBaseHVQMNR)
add_library(${TARGET} SHARED
  hvqcrsx.f
  MNR.cc
  MNRFrag.cc
  MNRGrid.cc
  ReactionBaseHVQMNR.cc
)
target_link_libraries(${TARGET} PUBLIC xfitter ROOT)
target_include_directories(${TARGET} PUBLIC .)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
