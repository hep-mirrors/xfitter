if(PineAPPL_FOUND)
  message(STATUS "PineAPPL library found; EFT reaction will support PineAPPL grids.")
  option(WITH_PINEAPPL "Include pineappl library support" ON)
else()
  message(STATUS "PineAPPL library not found, you can still compile and use the EFT reaction w/o PineAPPL support.")
  option(WITH_PINEAPPL "Include pineappl library support" OFF)
endif()

set(TARGET reactionEFT)
add_library(${TARGET} MODULE 
  ReactionEFT.cc
  EFTTerm.cc
  Vec.cc
)

if(PineAPPL_FOUND)
  target_compile_definitions(${TARGET} PRIVATE WITH_PINEAPPL)
endif()

target_link_libraries(${TARGET} PRIVATE xfitter)

# Conditionally link PineAPPL
if(PineAPPL_FOUND)
    target_link_libraries(${TARGET} PRIVATE xfitter PineAPPL)
endif()

target_link_libraries(${TARGET} PRIVATE xfitter APPLgrid)

install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
install(FILES EFT.yaml DESTINATION ${DEST_YAML}/reactions)
