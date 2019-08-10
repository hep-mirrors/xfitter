add_library(xfitter SHARED
  chi2scan.cc
  #c_interface.f
  dataset_tools.f
  dependent_pars.cc
  dis_sigma_ht.f
  error_bands_pumplin.f
  error_logging.f
  expression_utils.cc
  ExtraParConstr.cc
  fcn.f
  fortran_interface.cc
  fortran_io.f
  ftheor_eval.cc
  GetChisquare.f
  GetDummyXsection.f
  GetPDFUncType.f
  g_offset.cc
  #hfbanner.cc
  init_theory.f
  IntegrateDIS.cc
  lhapdf6_output.cc
  lhapdferrors.cc
  #main.f
  mc_errors.f
  minuit_ini.f
  offset_fns.f
  prep_corr.f
  Profiler.cc
  QCDNUM_Manager.cc
  ReactionTheory.cc
  read_data.f
  read_steer.f
  rediagonalize.f
  store_output.f
  systematics.f
  TermData.cc
  TheorEval.cc
  theory_dispatcher.f
  tinyexpr.c
  xfitter_cpp_base.cc
  xfitter_pars.cc
  xfitter_steer.cc
)
target_include_directories(xfitter PUBLIC ${CMAKE_SOURCE_DIR}/include)
target_compile_definitions(xfitter PRIVATE -DDLOPEN_PATH="${CMAKE_INSTALL_PREFIX}/${DEST_MODULE}")
install(TARGETS xfitter DESTINATION ${DEST_LIB})

#Executable
add_executable(xfitter-bin main.cc)
set_target_properties(xfitter-bin PROPERTIES OUTPUT_NAME xfitter)
target_link_libraries(xfitter-bin xfitter)
install(TARGETS xfitter-bin DESTINATION ${DEST_BIN})
