if(NOT PineAPPL_FOUND)
  message(STATUS "Skipping optional reaction PineAPPL because PineAPPL library was not found")
  return()
endif()

set(TARGET reactionPineAPPL)
add_library(${TARGET} MODULE ReactionPineAPPL.cc)
target_link_libraries(${TARGET} PRIVATE xfitter PineAPPL)
#target_include_directories(${TARGET} PRIVATE /home/makelat/work/xfitter-dev/deps/pineappl/include/)
#-L/home/makelat/work/xfitter-dev/deps/pineappl/lib -lpineappl_capi
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
install(FILES PineAPPL.yaml DESTINATION ${DEST_YAML}/reactions)
