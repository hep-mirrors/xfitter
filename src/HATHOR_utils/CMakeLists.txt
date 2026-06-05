if(NOT HATHOR_FOUND)
  message(STATUS "Skipping optional library HATHOR_utils because Hathor library was not found")
  return()
endif()
set(TARGET HATHOR_utils)
add_library(HATHOR_utils SHARED
  HathorPdfxFitter.cc
  HathorGenericIntegrator.cc
)
target_link_libraries(HATHOR_utils xfitter HATHOR)
install(TARGETS HATHOR_utils DESTINATION ${DEST_MODULE})
