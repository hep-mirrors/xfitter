set(TARGET evolutionFlipCharge)
add_library(${TARGET} MODULE FlipCharge.cc)
target_link_libraries(${TARGET} PRIVATE xfitter)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
