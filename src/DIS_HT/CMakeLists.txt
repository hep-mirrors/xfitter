add_library(DIS_HT SHARED
  DIS_HT.cc
)
target_link_libraries(DIS_HT xfitter)
install(TARGETS DIS_HT DESTINATION ${DEST_MODULE})
