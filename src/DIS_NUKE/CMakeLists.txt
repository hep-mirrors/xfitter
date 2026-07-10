add_library(DIS_NUKE SHARED
DIS_NUKE.cc
nuke_fast.F
)
target_compile_definitions(DIS_NUKE PRIVATE NUKE_DATA_DIR="${CMAKE_CURRENT_SOURCE_DIR}/petti14")
target_compile_options(DIS_NUKE PRIVATE $<$<COMPILE_LANGUAGE:Fortran>:-ffixed-line-length-none>)
target_link_libraries(DIS_NUKE xfitter)
install(TARGETS DIS_NUKE DESTINATION ${DEST_MODULE})
