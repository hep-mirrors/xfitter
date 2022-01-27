if(NOT yaml_FOUND)
  message(STATUS "Skipping optional profiling tools since yaml (C-version) is not found")
  return()
endif()

set(TARGET pdf2yaml)
add_library(${TARGET} SHARED c2yaml.c  list.c  pdf2yaml.c)
target_link_libraries(${TARGET} PRIVATE xfitter yaml)
install(TARGETS ${TARGET} DESTINATION ${DEST_LIB})

# also test executable
add_executable(test_pdf2yaml tests.c)
target_link_libraries(test_pdf2yaml pdf2yaml xfitter yaml)
# Note that test_pdf2yaml executable can be found in ./build/pdf2yaml/ directory"
