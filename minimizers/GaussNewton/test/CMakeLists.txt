# Unit tests of core/ (no xFitter needed). Run: build/minimizers/GaussNewton/test/gn_core_tests
add_executable(gn_core_tests test_core.cc ${GN_CORE_SOURCES})
target_compile_features(gn_core_tests PRIVATE cxx_std_17)
target_include_directories(gn_core_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_link_libraries(gn_core_tests PRIVATE Eigen3::Eigen)
