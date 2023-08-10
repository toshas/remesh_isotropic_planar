cmake_minimum_required(VERSION 3.25)
project(remesh_isotropic_planar)

set(CMAKE_CXX_STANDARD 17)

list(INSERT CMAKE_MODULE_PATH 0 /usr/local/share/cmake/Modules)

find_package(Eigen3 REQUIRED)
message(STATUS "Eigen version: ${EIGEN3_VERSION_STRING}")
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(CGAL REQUIRED)
if(CGAL_FOUND)
    message(STATUS "CGAL library found.")
else()
    message(FATAL_ERROR "Cannot find CGAL.")
endif()
include_directories(${CGAL_INCLUDE_DIRS})

find_package(Boost REQUIRED COMPONENTS program_options)

include_directories(include)

link_directories(${CGAL_LIBRARY_DIRS})

add_executable(remesh_isotropic_planar
        src/remesh_isotropic_planar.cpp
        src/utils.cpp
        src/main.cpp
)
target_link_libraries(remesh_isotropic_planar
        ${CGAL_LIBRARIES}
        Boost::program_options
)
