if(PineAPPL_FOUND)
  message(STATUS "PineAPPL library found; EFT reaction will support PineAPPL grids.")
  option(WITH_PINEAPPL "add support for PineAPPL" ON)
else()
  message(STATUS "PineAPPL library not found, you can still compile and use the EFT reaction w/o PineAPPL support.")
  option(WITH_PINEAPPL "add support for PineAPPL" OFF)
endif()

if(APPLgrid_FOUND)
  option(WITH_APPLGRID "add support for APPLgrid" ON)
else()
  option(WITH_APPLGRID "add support for APPLgrid" OFF)
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

if(APPLgrid_FOUND)
  target_compile_definitions(${TARGET} PRIVATE WITH_APPLGRID)
endif()

target_link_libraries(${TARGET} PRIVATE xfitter)

# Conditionally link PineAPPL
if(PineAPPL_FOUND)
  target_link_libraries(${TARGET} PRIVATE xfitter PineAPPL)
endif()

if(APPLgrid_FOUND)
  target_link_libraries(${TARGET} PRIVATE xfitter APPLgrid)
endif()

install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
install(FILES EFT.yaml DESTINATION ${DEST_YAML}/reactions)
