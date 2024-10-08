find_package(PkgConfig REQUIRED)
#save PKG_CONFIG_PATH from environment
set(PKG_CONFIG_PATH_SAVED "$ENV{PKG_CONFIG_PATH}")
if(EXISTS "$ENV{yaml_cpp_DIR}")
  #if yaml_cpp_DIR environment variable is set, make pkg-config search there too
  set(ENV{PKG_CONFIG_PATH} "$ENV{yaml_cpp_DIR}/share/pkgconfig:$ENV{yaml_cpp_DIR}/lib/pkgconfig:$ENV{yaml_cpp_DIR}/lib64/pkgconfig:$ENV{yaml_cpp_DIR}/pkgconfig:$ENV{yaml_cpp_DIR}:$ENV{PKG_CONFIG_PATH}")
endif()
pkg_check_modules(yaml-cpp QUIET yaml-cpp)
if(yaml-cpp_FOUND)
  if(APPLE)
    set(yaml-cpp_SO "lib${yaml-cpp_LIBRARIES}.dylib")
  endif()
  if(UNIX AND NOT APPLE)
    set(yaml-cpp_SO "lib${yaml-cpp_LIBRARIES}.so")
  endif()
  find_library(yaml-cpp_IMPORTED_LOCATION ${yaml-cpp_SO} HINTS "${yaml-cpp_LIBRARY_DIRS}")
  if(EXISTS ${yaml-cpp_IMPORTED_LOCATION})
    add_library(yaml-cpp SHARED IMPORTED)
    set_target_properties(yaml-cpp PROPERTIES IMPORTED_LOCATION "${yaml-cpp_IMPORTED_LOCATION}")
    set_target_properties(yaml-cpp PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${yaml-cpp_INCLUDE_DIRS}")
  else()
    message(WARNING "yaml-cpp found by pkg-config, but the actual library file \"${yaml-cpp_SO}\" could not be found")
    set(yaml-cpp_FOUND 0)
  endif()
endif()
if(NOT yaml-cpp_FOUND AND yaml-cpp_FIND_REQUIRED)
  message(FATAL_ERROR "yaml-cpp not found
(If it is actually installed, set the environment variable yaml_cpp_DIR to its install prefix)")
endif()
if(NOT yaml-cpp_FIND_QUIETLY)
  message(STATUS "Found yaml-cpp ${yaml-cpp_VERSION}: ${yaml-cpp_IMPORTED_LOCATION}")
endif()
#restore PKG_CONFIG_PATH
set(ENV{PKG_CONFIG_PATH} "${PKG_CONFIG_PATH_SAVED}")
