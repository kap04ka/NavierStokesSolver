cmake_minimum_required(VERSION 3.20)
project(cfd_pipe LANGUAGES CXX CUDA)

if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_CUDA_ARCHITECTURES 61 CACHE STRING "Compute Capability 6.1 for GTX 1060")
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(ENABLE_OMP "Enable OpenMP parallelization (requires OpenMP support from compiler)" ON) 
option(ENABLE_CUDA "Enable CUDA parallelization" ON) 

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
    src/math/cuda/JacobiSolverCUDA.cu
    src/math/cuda/SORSolverRedBlackCUDA.cu
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

# ─── CUDA support ───────────────────────────────────────────────────────────

if(ENABLE_CUDA)
    find_package(CUDAToolkit QUIET) 
    if(CUDAToolkit_FOUND)
        message(STATUS "CUDA Toolkit found, enabling CUDA support specifically for 'math' library.")

        if(MATH_CUDA_SOURCES) 
            target_sources(math PRIVATE ${MATH_CUDA_SOURCES}) 
        endif()
        
        target_link_libraries(math PUBLIC CUDA::cudart) 
        target_compile_definitions(math PUBLIC USE_CUDA) 
    else()
        message(WARNING "ENABLE_CUDA is ON, but CUDA Toolkit not found. 'math' library built without CUDA specifics.")
    endif()
endif()

add_custom_command(
    TARGET cfd_visualizer POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
        "${CMAKE_CURRENT_SOURCE_DIR}/src/resources/ofont.ru_ISOCPEUR.ttf"  # <--- Откуда копировать
        "$<TARGET_FILE_DIR:cfd_visualizer>"                       # <--- Куда копировать
    COMMENT "Copying font to executable directory"
)