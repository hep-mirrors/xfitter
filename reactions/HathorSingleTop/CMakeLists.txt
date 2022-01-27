if(NOT HATHOR_FOUND)
  message(STATUS "Skipping optional reaction HathorSingleTop because Hathor library was not found")
  return()
endif()

set(TARGET reactionHathorSingleTop)
add_library(${TARGET} MODULE 
  ReactionHathorSingleTop.cc 
  HathorPdfxFitter.cc
)
target_link_libraries(${TARGET} PRIVATE xfitter HATHOR)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
install(FILES HathorSingleTop.yaml DESTINATION ${DEST_YAML}/reactions)
