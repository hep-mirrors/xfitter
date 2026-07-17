set(TARGET reactionKFactor)
#add_library(${TARGET} MODULE ReactionKFactor.cc)
# https://www.ttp.kit.edu/_media/preprints/2017/ttp17-011.pdf?utm_source=chatgpt.com
add_library(${TARGET} MODULE ReactionKFactor.cc CRunDec.cpp)
target_link_libraries(${TARGET} PRIVATE xfitter)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
