add_library(minimizerMINUIT MODULE MINUITMinimizer.cc)
target_link_libraries(minimizerMINUIT PRIVATE xfitter)
install(TARGETS minimizerMINUIT DESTINATION ${DEST_MODULE})
