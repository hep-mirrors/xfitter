if(NOT ROOT_FOUND)
  message(STATUS "Skipping optional tool xfitter-draw because ROOT library was not found")
  return()
endif()

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
target_link_libraries(xfitter-draw xfitter ROOT)
#Disable warning "deprecated conversion from string constant to ‘char*’" for one source file
#set_property(SOURCE DrawLogo.cc PROPERTY COMPILE_FLAGS "-Wnowrite-strings")
install(TARGETS xfitter-draw DESTINATION ${DEST_BIN})
