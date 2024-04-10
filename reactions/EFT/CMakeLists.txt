option(WITH_PINEAPPL "Include pineappl library support" ON)

set(TARGET reactionEFT)
add_library(${TARGET} MODULE 
  ReactionEFT.cc
  EFTTerm.cc
  Vec.cc
)

# Add a definition based on WITH_PINEAPPL option
if(WITH_PINEAPPL)
    target_compile_definitions(${TARGET} PRIVATE WITH_PINEAPPL)
endif()

target_link_libraries(${TARGET} PRIVATE xfitter)

# Conditionally link PineAPPL based on WITH_PINEAPPL option
if(WITH_PINEAPPL)
    target_link_libraries(${TARGET} PRIVATE xfitter PineAPPL)
endif()

target_link_libraries(${TARGET} PRIVATE xfitter APPLgrid)

install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
install(FILES EFT.yaml DESTINATION ${DEST_YAML}/reactions)
