#Determine sources
set(xfitter_SOURCES
  BaseMinimizer.cc
  BasePdfParam.cc
  chi2scan.cc
  #c_interface.f
  dataset_tools.f
  dependent_pars.cc
  #dis_sigma_ht.f
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
  g_offset.cc
  #hfbanner.cc
  init_theory.f
  IntegrateDIS.cc
  lhapdf6_output.cc
  #lhapdferrors.cc
  main.f
  mc_errors.f
  minuit_ini.f
  #offset_fns.f
  prep_corr.f
  Profiler.cc
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
  num_utils/ranlux.f
  num_utils/rnpssn.f
  num_utils/dinv.f
  num_utils/dfact.f
  num_utils/deqinv.f
  num_utils/ranmar.f
  num_utils/rnormx.f
  num_utils/deqn.f
  num_utils/dfeqn.f
  num_utils/dfinv.f
  num_utils/rzero.f
  num_utils/rnorml.f
  num_utils/datime.f #placeholder for datetime
  num_utils/ddilog.f
  num_utils/snleq64.f
  Xstring.cc
  FitPars_base.cc
  FTNFitPars.cc
  ../tools/draw/src/pdferrors.cc
)
#Conditionally compiled sources
if(QCDNUM_FOUND)
  list(APPEND xfitter_SOURCES QCDNUM_Manager.cc)
endif()
if(LHAPDF_FOUND)
  list(APPEND xfitter_SOURCES CheckForPDF.cxx)
else()
  list(APPEND xfitter_SOURCES lhapdf_dum.f)
endif()

add_library(xfitter SHARED ${xfitter_SOURCES})
target_include_directories(xfitter PUBLIC
  ${CMAKE_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}/include/linalg
  ${CMAKE_SOURCE_DIR}/tools/draw/include/
)
target_compile_definitions(xfitter PRIVATE
  XFITTER_DEFAULT_MODULE_PATH="${CMAKE_INSTALL_PREFIX}/${DEST_MODULE}"
  XFITTER_YAML_PATH="${CMAKE_INSTALL_PREFIX}/${DEST_YAML}/"
  PREFIX="${CMAKE_INSTALL_PREFIX}"
)
target_link_libraries(xfitter PRIVATE
  MINUIT
  ${CMAKE_DL_LIBS} #for dlopen
  ${LAPACK_LIBRARIES} ZLIB::ZLIB z
)
target_link_libraries(xfitter PUBLIC yaml-cpp)

#Optional libraries
if(QCDNUM_FOUND)
  target_link_libraries(xfitter PRIVATE QCDNUM)
  target_compile_definitions(xfitter PUBLIC QCDNUM_FOUND="1")
endif()
if(ROOT_FOUND)
  if(ROOT_VERSION VERSION_GREATER 6.30.05)
    target_link_libraries(xfitter PUBLIC ROOT::Graf)
  else()
    target_link_libraries(xfitter PUBLIC ROOT)
  endif()
  target_compile_definitions(xfitter PUBLIC ROOT_FOUND="1")
endif()
if(APPLgrid_FOUND)
  target_compile_definitions(xfitter PUBLIC APPLgrid_FOUND="1")
endif()
if(LHAPDF_FOUND)
  target_link_libraries(xfitter PRIVATE LHAPDF)
  target_compile_definitions(xfitter PUBLIC LHAPDF_FOUND="1")
endif()

#Install location
install(TARGETS xfitter DESTINATION ${DEST_LIB})

#Executable
add_executable(xfitter-bin main.cc)
set_target_properties(xfitter-bin PROPERTIES OUTPUT_NAME xfitter)
target_link_libraries(xfitter-bin xfitter ZLIB::ZLIB z)
install(TARGETS xfitter-bin DESTINATION ${DEST_BIN})

add_subdirectory(minuit)

#Use -fno-automatic flag for fortran
#It means to treat all local variables as static by default (as if SAVE was specified)
#Without this flag hf_errlog will crash, and some other code will also work incorrectly
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fno-automatic")
elseif(CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -noauto")
else()
  message(FATAL_ERROR "Unknown Fortran compiler ${CMAKE_Fortran_COMPILER_ID}")
endif()
