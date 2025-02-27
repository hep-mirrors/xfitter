add_library(DIS_NUKE SHARED
DIS_NUKE.cc
nuke_fast.F
)
target_link_libraries(DIS_NUKE xfitter)
install(TARGETS DIS_NUKE DESTINATION ${DEST_MODULE})
