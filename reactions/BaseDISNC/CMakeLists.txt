if(NOT QCDNUM_FOUND)
  message(STATUS "Skipping optional reaction BaseDISNC because QCDNUM library was not found")
  return()
endif()

set(TARGET reactionBaseDISNC)
add_library(${TARGET} SHARED ReactionBaseDISNC.cc)
target_link_libraries(${TARGET} PRIVATE xfitter evolutionQCDNUM)
target_include_directories(${TARGET} PUBLIC .)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
