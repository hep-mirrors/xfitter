if(NOT ROOT_FOUND)
  message(STATUS "Skipping optional tool xfitter-draw because ROOT library was not found")
  return()
endif()

# Find ROOT version
execute_process(COMMAND root-config --version OUTPUT_VARIABLE ROOT_VERSION)

# Remove newline from the end of the output
string(STRIP ${ROOT_VERSION} ROOT_VERSION)

# Parse the ROOT version
string(REGEX REPLACE "([0-9]+)\\.([0-9]+)/([0-9]+)" "\\1;\\2;\\3" ROOT_VERSION_LIST ${ROOT_VERSION})
list(GET ROOT_VERSION_LIST 0 ROOT_MAJOR)
list(GET ROOT_VERSION_LIST 1 ROOT_MINOR)
list(GET ROOT_VERSION_LIST 2 ROOT_PATCH)

add_executable(xfitter-draw
  Chi2.cc
  Chi2scanData.cc
  Chi2scanGauss.cc
  Chi2scanPainter.cc
  Chi2scanUnc.cc
  CommandParser.cc
  DataPainter.cc
  Dataset.cc
  DrawLogo.cc
  FitPainter.cc
  Outdir.cc
  Par.cc
  ParPainter.cc
  PdfData.cc
  pdferrors.cc
  PdfsPainter.cc
  ShiftPainter.cc
  xfitter-draw.cc
)
target_include_directories(xfitter-draw PRIVATE
  ../include
)
target_compile_definitions(xfitter-draw PRIVATE
  -DVERSION="${xFitter_VERSION}"
)

# Compare version (assuming ROOT version format is always 3 numbers separated by periods or slashes)
if((ROOT_MAJOR GREATER 6) OR (ROOT_MAJOR EQUAL 6 AND ROOT_MINOR GREATER 26) OR (ROOT_MAJOR EQUAL 6 AND ROOT_MINOR EQUAL 26 AND ROOT_PATCH GREATER_EQUAL 10))
  add_definitions(-DBUGFIXEDROOT)
  message (STATUS "Use system fonts in xFitter logo")
else()
  message (STATUS "Use internal fonts in xFitter logo ${ROOT_MAJOR}.${ROOT_MINOR}.${ROOT_PATCH} < 6.26.10")  
endif()

target_link_libraries(xfitter-draw xfitter ROOT)
#Disable warning "deprecated conversion from string constant to ‘char*’" for one source file
#set_property(SOURCE DrawLogo.cc PROPERTY COMPILE_FLAGS "-Wnowrite-strings")
install(TARGETS xfitter-draw DESTINATION ${DEST_BIN})
