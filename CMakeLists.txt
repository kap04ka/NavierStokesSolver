cmake_minimum_required(VERSION 3.20)
project(cfd_pipe LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_library(core
    src/core/Geometry.cpp)
add_library(math  src/math/PoissonSolver.cpp)
add_library(solvers
    src/solvers/base/Solver.cpp
    src/solvers/velocity_pressure/VelocityPressureSolver.cpp)

add_executable(test_obstacle tests/test_obstacle.cpp)

foreach(t core math solvers)
    target_include_directories(${t} PUBLIC src)
endforeach()

target_link_libraries(math    PUBLIC core)
target_link_libraries(solvers PUBLIC core math)
target_link_libraries(test_obstacle PRIVATE solvers)

find_package(OpenMP QUIET)
if(OpenMP_CXX_FOUND)
  foreach(t core math solvers)
    target_link_libraries(${t} PUBLIC OpenMP::OpenMP_CXX)
    target_compile_definitions(${t} PUBLIC USE_OPENMP)
  endforeach()
endif()