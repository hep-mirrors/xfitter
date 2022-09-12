if(NOT Boost_FOUND)
  message(STATUS "Skipping optional parametrisation Chebyschev because BOOST library was not found")
  return()
endif()



set(TARGET pdfparamChebyschev)
include_directories(${Boost_INCLUDE_DIRS})
add_library(${TARGET} MODULE Chebyschev_PdfParam.cc)
target_link_libraries(${TARGET} PRIVATE xfitter)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
