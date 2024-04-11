option(WITH_PINEAPPL "Include pineappl library support" ON)
# if the PineAPPL library is not installed together with xFitter, you can still compile and use the EFT reaction partially with:
# option(WITH_PINEAPPL "Include pineappl library support" OFF)

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
# todo: fastNLO

install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
install(FILES EFT.yaml DESTINATION ${DEST_YAML}/reactions)
