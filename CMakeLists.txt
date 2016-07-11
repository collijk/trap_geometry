cmake_minimum_required(VERSION 3.5)
project(trap_geometry)

find_package(CGAL REQUIRED)

include(${CGAL_USE_FILE})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES mesh_electrodes.cpp)
add_executable(mesh_electrodes ${SOURCE_FILES})

add_to_cached_list( CGAL_EXECUTABLE_TARGETS mesh_electrodes )

target_link_libraries( mesh_electrodes ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})