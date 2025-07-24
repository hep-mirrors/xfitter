if(NOT APPLgrid_FOUND)
  message(STATUS "Skipping optional reaction APPLgrid because APPLgrid library was not found")
  return()
endif()

if(NOT ROOT_FOUND)
  message(STATUS "Skipping optional reaction APPLgrid because ROOT library was not found")
  return()
endif()

if (PLOUGHSHARE_FOUND)
  message(STATUS "Found ploughshare, including in APPLgrid install ")
  add_definitions(-DPLOUGHSHARE_FOUND)
  add_definitions(-I${PLOUGHSHARE_INCLUDE_DIRS})
endif()

# Add fapplwrap subdirectory
add_subdirectory(fapplwrap)

set(TARGET reactionAPPLgrid)
add_library(${TARGET} MODULE ReactionAPPLgrid.cc)
if (PLOUGHSHARE_FOUND)
  target_link_libraries(${TARGET} PRIVATE xfitter APPLgrid PLOUGHSHARE fapplwrap)
else()
  target_link_libraries(${TARGET} PRIVATE xfitter APPLgrid fapplwrap)
endif()

# Add fapplwrap include directory
target_include_directories(${TARGET} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/fapplwrap)

install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
install(FILES APPLgrid.yaml DESTINATION ${DEST_YAML}/reactions)
