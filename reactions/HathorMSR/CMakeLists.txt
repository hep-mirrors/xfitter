if(NOT HATHOR_FOUND)
  message(STATUS "Skipping optional reaction HathorMSR because Hathor library was not found")
  return()
endif()

set(TARGET reactionHathorMSR)
add_library(${TARGET} MODULE 
  ReactionHathorMSR.cc 
  HathorPdfxFitter.cc
)
target_link_libraries(${TARGET} PRIVATE xfitter HATHOR)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
install(FILES HathorMSR.yaml DESTINATION ${DEST_YAML}/reactions)
