add_library(DIS_TMC SHARED
DIS_TMC.cc
)
target_link_libraries(DIS_TMC xfitter gsl)
install(TARGETS DIS_TMC DESTINATION ${DEST_MODULE})
