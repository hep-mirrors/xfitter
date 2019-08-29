if(NOT ROOT_FOUND)
  message(STATUS "Skipping optional reaction BaseHVQMNR because ROOT library was not found")
  return()
endif()

set(TARGET reactionBaseHVQMNR)
add_library(${TARGET} SHARED
  src/hvqcrsx.f
  src/Makefile.am
  src/MNR.cc
  src/MNRFrag.cc
  src/MNRGrid.cc
  src/ReactionBaseHVQMNR.cc
)
target_link_libraries(${TARGET} PUBLIC xfitter ROOT)
target_include_directories(${TARGET} PUBLIC include)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
