set(TARGET reactionBaseTBSF)
add_library(${TARGET} SHARED
  ReactionBaseTBSF.cc
)
set(STRUCTFUNCS_DIR /home/zenaiev/soft/structfuncs/build/src)
target_link_libraries(${TARGET} PUBLIC 
  xfitter
  #structfuncs
  ${STRUCTFUNCS_DIR}/structfuncs/libstructfuncs.a
  ${STRUCTFUNCS_DIR}/coefffuncs/libcoefffuncs.a
  ${STRUCTFUNCS_DIR}/constants/libconstants.a
  ${STRUCTFUNCS_DIR}/splitfuncs/libsplitfuncs.a
  ${STRUCTFUNCS_DIR}/logarithms/libHPLs.a
  ${STRUCTFUNCS_DIR}/logarithms/libli2.a
  ${STRUCTFUNCS_DIR}/logarithms/libLi2.a
  ${STRUCTFUNCS_DIR}/logarithms/libhplog5.a
  ${STRUCTFUNCS_DIR}/utils/libutils.a
  ${STRUCTFUNCS_DIR}/n3loblumlein/libn3loblumlein.a
  ${STRUCTFUNCS_DIR}/n3loblumlein/libn3loblumleinfortran.a
  ${STRUCTFUNCS_DIR}/WILS3HQF2/libheavyasymp.a
  ${STRUCTFUNCS_DIR}/WILS3HQF2/libWILS3HQF2.a
  reactionBaseFFABM
  reactionBaseDISNC
  reactionBaseDISCC
  #cuba
)

target_include_directories(${TARGET} PUBLIC .)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
install(FILES BaseTBSF.yaml DESTINATION ${DEST_YAML}/reactions)
