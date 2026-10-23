cmake_minimum_required(VERSION 3.16)
project(VolumeConductionFEM LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 3.3 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
    find_path(EIGEN3_INCLUDE_DIR Eigen/Core PATH_SUFFIXES eigen3)
    if(NOT EIGEN3_INCLUDE_DIR)
        message(FATAL_ERROR "Eigen3 manquant : installer libeigen3-dev ou fournir -DEIGEN3_INCLUDE_DIR=...")
    endif()
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    set_target_properties(Eigen3::Eigen PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${EIGEN3_INCLUDE_DIR}")
endif()

add_executable(conduction
    src/main.cpp
    src/maillage.cpp
    src/solveur.cpp
    src/output.cpp
)
target_include_directories(conduction PRIVATE include)
target_link_libraries(conduction PRIVATE Eigen3::Eigen)
