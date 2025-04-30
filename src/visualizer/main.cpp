#include "visualizer/FlowVisualizer.hpp"
#include "core/Geometry.hpp"
#include "solvers/velocity_pressure/VelocityPressureSolver.hpp"

#include <GLFW/glfw3.h>
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <thread>
#include <mutex>
#include <atomic>
#include <iostream>

struct SimulationConfig {
    int NX = 120;
    int NY = 60;
    double Lx = 1.0;
    double Ly = 0.5;
    bool useObstacle = true;
    int obs_i0 = 40, obs_j0 = 20;
    int obs_i1 = 80, obs_j1 = 40;
    double rho = 1000.0;
    double nu = 1e-3;
    double umax = 1.0;
    double cfl = 0.4;
    double omega = 1.7;
    int ptype = static_cast<int>(cfd::PoissonType::SOR);
};

int main() {
    SimulationConfig cfg;
    bool started = false;

    // ------------------- CONFIGURATION WINDOW -------------------
    if (!glfwInit()) return -1;
    GLFWwindow* cfgWindow = glfwCreateWindow(600, 400, "Simulation Setup", nullptr, nullptr);
    if (!cfgWindow) {
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(cfgWindow);
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui_ImplGlfw_InitForOpenGL(cfgWindow, true);
    ImGui_ImplOpenGL3_Init("#version 130");
    ImGui::StyleColorsDark();

    while (!started && !glfwWindowShouldClose(cfgWindow)) {
        glfwPollEvents();
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        ImGui::Begin("Simulation Configuration");
        ImGui::InputInt("NX", &cfg.NX);
        ImGui::InputInt("NY", &cfg.NY);
        ImGui::InputDouble("Lx", &cfg.Lx);
        ImGui::InputDouble("Ly", &cfg.Ly);
        ImGui::Checkbox("Use Obstacle", &cfg.useObstacle);
        if (cfg.useObstacle) {
            ImGui::InputInt2("Obs i0,j0", &cfg.obs_i0);
            ImGui::InputInt2("Obs i1,j1", &cfg.obs_i1);
        }
        ImGui::InputDouble("rho", &cfg.rho);
        ImGui::InputDouble("nu", &cfg.nu);
        ImGui::InputDouble("umax", &cfg.umax);
        ImGui::InputDouble("CFL", &cfg.cfl);
        ImGui::InputDouble("omega", &cfg.omega);
        static const char* poissonItems[] = {"Jacobi","SOR"};
        ImGui::Combo("Poisson Type", &cfg.ptype, poissonItems, IM_ARRAYSIZE(poissonItems));
        if (ImGui::Button("Start Simulation")) started = true;
        ImGui::End();

        ImGui::Render();
        int w, h;
        glfwGetFramebufferSize(cfgWindow, &w, &h);
        glViewport(0, 0, w, h);
        glClearColor(0.1f,0.1f,0.1f,1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(cfgWindow);
    }

    // Cleanup configuration UI
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(cfgWindow);
    glfwTerminate();

    // ------------------- SIMULATION SETUP -------------------
    cfd::Geometry geom(cfg.NX, cfg.NY, cfg.Lx, cfg.Ly);
    if (cfg.useObstacle) {
        geom.add_rectangle(cfg.obs_i0, cfg.obs_j0,
                           cfg.obs_i1, cfg.obs_j1);
    }
    cfd::VelocityPressureSolver solver(
        geom,
        cfg.rho, cfg.nu,
        static_cast<cfd::PoissonType>(cfg.ptype),
        cfg.cfl, cfg.omega
    );
    solver.set_inlet_parabola(cfg.umax);

    // Prepare shared buffers and threading
    auto u_buffer = solver.u();
    auto v_buffer = solver.v();
    std::mutex data_mutex;
    std::atomic<bool> stopFlag{false};
    double dt = 0.001;

    std::thread worker([&]() {
        while (!stopFlag) {
            solver.step(dt);
            {
                std::lock_guard<std::mutex> lock(data_mutex);
                u_buffer = solver.u();
                v_buffer = solver.v();
            }
        }
    });

    // ------------------- VISUALIZATION -------------------
    FlowVisualizer vis(geom, u_buffer, v_buffer, geom.tags());
    vis.init();
    while (vis.shouldRun()) {
        {
            std::lock_guard<std::mutex> lock(data_mutex);
            // share latest u_buffer, v_buffer
        }
        vis.render();
    }

    stopFlag = true;
    worker.join();
    vis.shutdown();
    return 0;
}