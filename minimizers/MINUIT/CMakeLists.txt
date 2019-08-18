add_library(minimizerMINUIT MODULE MINUITMinimizer.cc)
target_link_libraries(minimizerMINUIT PRIVATE xfitter yaml-cpp)
install(TARGETS minimizerMINUIT DESTINATION ${DEST_MODULE})
