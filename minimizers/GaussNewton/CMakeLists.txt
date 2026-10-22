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

set(TARGET minimizerGaussNewton)
add_library(${TARGET} MODULE
  GaussNewtonMinimizer.cc
)
target_compile_features(${TARGET} PRIVATE cxx_std_17)
target_include_directories(${TARGET} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(${TARGET} PRIVATE xfitter Eigen3::Eigen)
install(TARGETS ${TARGET} DESTINATION ${DEST_MODULE})
