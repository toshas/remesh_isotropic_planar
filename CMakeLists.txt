cmake_minimum_required(VERSION 3.25)

project(remesh_isotropic_planar)

set(CMAKE_CXX_STANDARD 17)

find_package(Eigen3 REQUIRED)
find_package(CGAL REQUIRED)
find_package(Boost REQUIRED COMPONENTS program_options)

add_executable(remesh_isotropic_planar
    src/remesh_isotropic_planar.cpp
    src/utils.cpp
    src/main.cpp
)

target_include_directories(remesh_isotropic_planar PRIVATE
    ${EIGEN3_INCLUDE_DIR}
    ${CGAL_INCLUDE_DIRS}
    include
)

target_link_libraries(remesh_isotropic_planar PRIVATE
    ${CGAL_LIBRARIES}
    Boost::program_options
)
