cmake_minimum_required(VERSION 3.20)
project(cfd_pipe LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(ENABLE_OMP "Enable OpenMP parallelization (requires OpenMP support from compiler)" ON) 

# ─── Core library ────────────────────────────────────────────────────────────
add_library(core
    src/core/Geometry.cpp
)
target_include_directories(core PUBLIC src)

# ─── Math (Poisson) library ──────────────────────────────────────────────────
add_library(math
    src/math/PoissonSolver.cpp
    src/math/PoissonSolverDyn.cpp
    src/math/sequential/JacobiSolverSeq.cpp
    src/math/sequential/SORSolverSeq.cpp
    src/math/omp/JacobiSolverOMP.cpp
    src/math/omp/SORSolverRedBlackOMP.cpp
)
target_link_libraries(math PUBLIC core)
target_include_directories(math PUBLIC src)

# ─── Solvers library ─────────────────────────────────────────────────────────
add_library(solvers
    src/solvers/base/Solver.cpp
    src/solvers/velocity_pressure/VelocityPressureSolver.cpp
    src/solvers/vorticity_streamfunction/VorticityStreamfunctionSolver.cpp
    src/turbulence/base/TurbulenceModel.cpp
    src/turbulence/k_epsilon/KEpsilonModel.cpp
)
target_link_libraries(solvers PUBLIC core math)
target_include_directories(solvers PUBLIC src)

# ─── Test executable ─────────────────────────────────────────────────────────
add_executable(test_obstacle
    tests/test_obstacle.cpp
)
target_include_directories(test_obstacle PRIVATE src)
target_link_libraries(test_obstacle PRIVATE solvers)

# ─── Visualizer executable ───────────────────────────────────────────────────
add_executable(cfd_visualizer
    src/visualizer/main.cpp
    src/visualizer/FlowVisualizer.cpp
    external/imgui/imgui.cpp
    external/imgui/imgui_draw.cpp
    external/imgui/imgui_widgets.cpp
    external/imgui/imgui_tables.cpp
    external/imgui/backends/imgui_impl_glfw.cpp
    external/imgui/backends/imgui_impl_opengl3.cpp
)
target_include_directories(cfd_visualizer PRIVATE
    src
    external/imgui
    external/imgui/backends
)
find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)
find_package(glfw3 REQUIRED)
target_link_libraries(cfd_visualizer PRIVATE
    core
    math
    solvers
    glfw
    GLEW::GLEW
    OpenGL::GL
)

# ─── OpenMP support ───────────────────────────────────────────────────────────
if(ENABLE_OMP) # Проверяем опцию пользователя
    find_package(OpenMP QUIET) # Ищем OpenMP
    if(OpenMP_CXX_FOUND)
        message(STATUS "OpenMP support enabled and OpenMP found.")
        foreach(t core math solvers cfd_visualizer) # Применяем ко всем нужным целям
            target_link_libraries(${t} PUBLIC OpenMP::OpenMP_CXX)
            target_compile_definitions(${t} PUBLIC USE_OPENMP)
        endforeach()
    else()
        message(WARNING "OpenMP support was requested (ENABLE_OMP=ON), but OpenMP was not found. Proceeding without OpenMP.")
    endif()
else()
    message(STATUS "OpenMP support disabled by user option (ENABLE_OMP=OFF).")
endif()