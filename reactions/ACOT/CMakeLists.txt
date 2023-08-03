if(NOT QCDNUM_FOUND)
  message(STATUS "Skipping optional reaction ACOT because QCDNUM library was not found")
  return()
endif()

set(TARGET reactionACOT)
# USE MSTW JUST FOR CROSS CHECKING
add_library(${TARGET} MODULE ReactionACOT.cc ext_pdfs.cc mstw2008.f acotNC_wrap.f  sf_acot_wrap.f  alphas70.f c2terms.f cc116.f dinteg.f fgen122qcdnum.f ingo.f ingoUtil.f qcdnumPDF.f  tot130.f unified121.f utilities113.f UseZmvnsScheme.f ACOTindexRESET.f )
target_link_libraries(${TARGET} PRIVATE xfitter reactionBaseDISNC)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
install(FILES ACOT.yaml DESTINATION ${DEST_YAML}/reactions)
