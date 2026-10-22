# GaussNewton minimizer: Levenberg-Marquardt on the profiled chi2 with a
# theory Jacobian. See README.md in this directory.
if(CMAKE_VERSION VERSION_LESS 3.19)
  find_package(Eigen3 3.3 QUIET NO_MODULE)
else()
  find_package(Eigen3 3.3...5 QUIET NO_MODULE)
endif()
if(NOT TARGET Eigen3::Eigen)
  message(STATUS "Skipping optional GaussNewton minimizer because Eigen3 was not found")
  return()
endif()

# core/: pure numerics (Eigen + std), shared by the plugin and the unit tests
set(GN_CORE_SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/core/ExtraSystRescale.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/core/LevenbergMarquardt.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/core/ObservedHessian.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/core/PointModel.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/core/ProfiledNuisances.cc
)

set(TARGET minimizerGaussNewton)
add_library(${TARGET} MODULE
  GaussNewtonMinimizer.cc
  ProfiledChi2.cc
  DerivativeCheck.cc
  bridge/xfitter_gn_bridge.f
  bridge/FortranBridge.cc
  bridge/TheoryEvaluator.cc
  jacobian/FiniteDifferenceJacobian.cc
  util/ForkPool.cc
  util/Report.cc
  ${GN_CORE_SOURCES}
)
target_compile_features(${TARGET} PRIVATE cxx_std_17)
target_include_directories(${TARGET} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(${TARGET} PRIVATE xfitter Eigen3::Eigen)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})

option(XFITTER_GN_TESTS "Build the GaussNewton core unit tests (gn_core_tests)" OFF)
if(XFITTER_GN_TESTS)
  add_subdirectory(test)
endif()
