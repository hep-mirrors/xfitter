set(TARGET evolutionFlipUD)
add_library(${TARGET} MODULE FlipUD.cc)
target_link_libraries(${TARGET} PRIVATE xfitter)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
