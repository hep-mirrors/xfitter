# Only find_package if target not already created
if(NOT TARGET APPLgrid)
    find_package(APPLgrid REQUIRED)
endif()

add_library(fapplwrap STATIC
    applwrap.cxx
    fapplwrap.cxx
    smooth.cxx
)

set_target_properties(fapplwrap PROPERTIES POSITION_INDEPENDENT_CODE ON)

target_include_directories(fapplwrap PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Link to APPLgrid if found
if(APPLgrid_FOUND)
    target_link_libraries(fapplwrap PUBLIC APPLgrid)
endif()

# Link to ROOT if found
if(ROOT_FOUND)
    target_include_directories(fapplwrap PUBLIC ${ROOT_INCLUDE_DIRS})
    target_link_libraries(fapplwrap PUBLIC ${ROOT_LIBRARIES})
endif()
