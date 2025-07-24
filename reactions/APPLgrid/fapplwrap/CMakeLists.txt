find_path(APPLgrid_INCLUDE_DIR
    appl_grid/appl_file.h
    PATHS /data/atlas/atlasdata/giuli/xFitter_fapplwrap/deps/applgrid/include
    NO_DEFAULT_PATH
)
find_library(APPLgrid_LIBRARY
    NAMES libAPPLgrid.so
    PATHS /data/atlas/atlasdata/giuli/xFitter_fapplwrap/deps/applgrid/lib
    NO_DEFAULT_PATH
)
if(APPLgrid_INCLUDE_DIR AND APPLgrid_LIBRARY)
    set(APPLgrid_FOUND TRUE)
    set(APPLgrid_INCLUDE_DIRS ${APPLgrid_INCLUDE_DIR})
    set(APPLgrid_LIBRARIES ${APPLgrid_LIBRARY})
else()
  message(FATAL_ERROR "APPLgrid not found. Please check your paths.")
endif()
  
add_library(fapplwrap STATIC
    applwrap.cxx
    fapplwrap.cxx
    smooth.cxx
)

set_target_properties(fapplwrap PROPERTIES POSITION_INDEPENDENT_CODE ON)

# Link APPLgrid and ROOT if available
if(APPLgrid_FOUND)
    target_include_directories(fapplwrap PUBLIC ${APPLgrid_INCLUDE_DIRS})
    target_link_libraries(fapplwrap PUBLIC ${APPLgrid_LIBRARIES})
endif()

if(ROOT_FOUND)
    target_include_directories(fapplwrap PUBLIC ${ROOT_INCLUDE_DIRS})
    target_link_libraries(fapplwrap PUBLIC ${ROOT_LIBRARIES})
endif()

target_include_directories(fapplwrap PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)
