set(TARGET reactionfastNLO)
add_library(${TARGET} MODULE fastNLOCoeffAddBase.cc fastNLOCoeffAddFlex.cc fastNLOCoeffData.cc fastNLOCoeffMult.cc fastNLOPDFLinearCombinations.cc fastNLOTable.cc speaker.cc fastNLOCoeffAddFix.cc fastNLOCoeffBase.cc fastNLOCoefficients.cc fastNLOReader.cc fastNLOTools.cc ReactionfastNLO.cc)
target_link_libraries(${TARGET} PRIVATE xfitter)
target_include_directories(${TARGET} PRIVATE .)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
