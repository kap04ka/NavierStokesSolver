#include "visualizer/FlowVisualizer.hpp"
#include "core/Geometry.hpp"
#include "solvers/base/Solver.hpp" // Базовый класс Solver
#include "solvers/velocity_pressure/VelocityPressureSolver.hpp"
#include "solvers/vorticity_streamfunction/VorticityStreamfunctionSolver.hpp" // Подключаем новый решатель

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <thread>
#include <mutex>
#include <atomic>
#include <iostream>
#include <string>
#include <iomanip> // Для std::fixed, std::setprecision
#include <memory>  // Для std::unique_ptr
#include <variant>
#include <vector>

// Тип решателя
enum class SolverChoice { VelocityPressure, VorticityStreamfunction };

enum class ObstacleType { Rectangle, Circle };

const char* ObstacleTypeToString(ObstacleType type) {
    switch (type) {
        case ObstacleType::Rectangle: return "Rectangle";
        case ObstacleType::Circle:    return "Circle";
        default:                      return "Unknown";
    }
}

struct RectangleObstacleParams {
    int i0 = 10, j0 = 10;
    int i1 = 20, j1 = 20;
};

struct CircleObstacleParams {
    int center_i = 0;      // Индекс i центра круга
    int center_j = 0;      // Индекс j центра круга
    double radius_phys = 0.0;  // Радиус в физ. единицах
};;

struct ObstacleConfig {
    ObstacleType type;
    std::variant<RectangleObstacleParams, CircleObstacleParams> params;
    bool enabled = true;
    std::string name;

    int imgui_id; 
    static int next_imgui_id;

    ObstacleConfig(ObstacleType t, decltype(params) p, std::string n = "") 
        : type(t), params(std::move(p)), enabled(true), name(std::move(n)), imgui_id(next_imgui_id++) {}
};
int ObstacleConfig::next_imgui_id = 0; // Статический член для уникальных ID


// Обновленная конфигурация симуляции
struct SimulationConfig {
    // Общие параметры сетки и домена
    int NX = 60;
    int NY = 30;
    double Lx = 1.0;
    double Ly = 0.5;

    // Препятствие
    std::vector<ObstacleConfig> obstacles;

    // Физические свойства
    double rho = 1.0;
    double nu = 1e-3; // Кинематическая вязкость
    double umax = 1.0; // Макс. скорость на входе (для профиля и Q)

    // Общие параметры симуляции
    double sim_duration = 10.0;
    double dt_user = 0.001;

    // Выбор решателя
    SolverChoice solverType = SolverChoice::VelocityPressure;

    // Параметры, которые могут быть специфичны или иметь разные оптимальные значения
    double cfl = 0.4; // Общий CFL, будет интерпретирован решателем (0.4 для VP, 0.5 для VS по умолчанию)

    // Параметры для решателя типа Пуассона (давление в VP, функция тока в VS)
    int poisson_solver_type = static_cast<int>(cfd::PoissonType::SOR);
    double poisson_sor_omega = 1.7;
    unsigned poisson_max_iter = 500; 
    double poisson_tol = 1e-5;

    // Параметры турбулентности (общие для обоих, если они ее поддерживают)
    int turbulenceChoice = static_cast<int>(cfd::TurbulenceModelType::None); // 0=None, 1=KEpsilon
    double inletTurbIntensity = 0.05;
    double inletLengthScaleFactor = 0.07;

    SimulationConfig() {
        ObstacleConfig::next_imgui_id = 0; // Сбрасываем счетчик ID при создании новой конфигурации
    }
};


int main() {
    SimulationConfig cfg;
    bool started = false;
    double simulationTime = 0.0;

    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return -1;
    }
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    ImGui::StyleColorsDark();

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
    GLFWwindow* cfgWindow = glfwCreateWindow(700, 650, "Simulation Setup", nullptr, nullptr);
    if (!cfgWindow) { 
        std::cerr << "Failed to create GLFW config window" << std::endl;
        ImGui::DestroyContext();
        glfwTerminate();
        return -1; 
    }
    glfwMakeContextCurrent(cfgWindow);
    glfwSwapInterval(1); // Enable vsync
    ImGui_ImplGlfw_InitForOpenGL(cfgWindow, true);
    ImGui_ImplOpenGL3_Init("#version 130");

    while (!started && !glfwWindowShouldClose(cfgWindow)) {
        glfwPollEvents();
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        ImGui::Begin("Simulation Configuration");
        SolverChoice previousSolverType = cfg.solverType; // Сохраняем предыдущий тип
        ImGui::Text("Solver Type:");
        ImGui::RadioButton("Velocity-Pressure", reinterpret_cast<int*>(&cfg.solverType), static_cast<int>(SolverChoice::VelocityPressure)); ImGui::SameLine();
        ImGui::RadioButton("Vorticity-Streamfunction", reinterpret_cast<int*>(&cfg.solverType), static_cast<int>(SolverChoice::VorticityStreamfunction));
        
        ImGui::Separator();

        ImGui::InputInt("NX (Grid Cells X)", &cfg.NX);
        ImGui::InputInt("NY (Grid Cells Y)", &cfg.NY);
        ImGui::InputDouble("Lx (Domain Length)", &cfg.Lx);
        ImGui::InputDouble("Ly (Domain Height)", &cfg.Ly);
        ImGui::Separator();
        ImGui::InputDouble("rho (Density)", &cfg.rho);
        ImGui::InputDouble("nu (Kinematic Viscosity)", &cfg.nu, 0.0, 0.0, "%.1e");
        ImGui::InputDouble("Umax (Inlet Profile)", &cfg.umax);
        ImGui::Separator();
        ImGui::InputDouble("User dt (Max)", &cfg.dt_user, 0.0, 0.0, "%.1e");
        ImGui::InputDouble("Simulation Duration (s)", &cfg.sim_duration);
        ImGui::Separator();

        ImGui::InputDouble("CFL Factor", &cfg.cfl, 0.0, 0.0, "%.2f");
        ImGui::Separator();
        ImGui::Text("Poisson-like Solver Settings (for P or Psi):");
        static const char* poissonItems[] = {"Jacobi","SOR"};
        ImGui::Combo("Solver Type##Poisson", &cfg.poisson_solver_type, poissonItems, IM_ARRAYSIZE(poissonItems));
        if (cfg.poisson_solver_type == static_cast<int>(cfd::PoissonType::SOR)) {
            ImGui::InputDouble("SOR Omega", &cfg.poisson_sor_omega, 0.0, 0.0, "%.2f");
        }
        ImGui::InputInt("Max Iterations##Poisson", (int*)&cfg.poisson_max_iter);
        ImGui::InputDouble("Tolerance##Poisson", &cfg.poisson_tol, 0.0, 0.0, "%.1e");
        ImGui::Separator();

        ImGui::Text("Flow Model (Common):");
        ImGui::RadioButton("Laminar##flow", &cfg.turbulenceChoice, static_cast<int>(cfd::TurbulenceModelType::None)); ImGui::SameLine();
        ImGui::RadioButton("k-epsilon##flow", &cfg.turbulenceChoice, static_cast<int>(cfd::TurbulenceModelType::KEpsilon));
        if (cfg.turbulenceChoice != static_cast<int>(cfd::TurbulenceModelType::None)) {
            ImGui::Indent();
            ImGui::InputDouble("Inlet Turb Intensity", &cfg.inletTurbIntensity, 0.0,0.0, "%.3f");
            cfg.inletTurbIntensity = std::max(0.0, std::min(1.0, cfg.inletTurbIntensity));
            ImGui::InputDouble("Inlet L Scale Factor", &cfg.inletLengthScaleFactor, 0.0,0.0, "%.3f");
            cfg.inletLengthScaleFactor = std::max(0.001, cfg.inletLengthScaleFactor);
            ImGui::Unindent();
        }

        ImGui::Separator();
        ImGui::Text("Obstacles Configuration");

        // Кнопка для открытия модального окна добавления препятствия
        if (ImGui::Button("Add New Obstacle")) {
            ImGui::OpenPopup("AddObstacleTypePopup");
        }

        // Модальное окно для выбора типа нового препятствия
        if (ImGui::BeginPopupModal("AddObstacleTypePopup", NULL, ImGuiWindowFlags_AlwaysAutoResize)) {
            ImGui::Text("Select obstacle type to add:");
            // Статическая переменная для хранения выбора в модальном окне
            static ObstacleType selected_new_obstacle_type = ObstacleType::Rectangle; 
            
            ImGui::RadioButton("Rectangle##AddType", reinterpret_cast<int*>(&selected_new_obstacle_type), static_cast<int>(ObstacleType::Rectangle));
            ImGui::RadioButton("Circle##AddType",    reinterpret_cast<int*>(&selected_new_obstacle_type), static_cast<int>(ObstacleType::Circle));
            // Сюда можно будет добавить другие типы препятствий в будущем

            ImGui::Separator();
            if (ImGui::Button("Add This Type", ImVec2(150, 0))) {
                std::string new_obs_name;
                if (selected_new_obstacle_type == ObstacleType::Rectangle) {
                    new_obs_name = "Rectangle " + std::to_string(cfg.obstacles.size() + 1);
                    cfg.obstacles.emplace_back(
                        ObstacleType::Rectangle, 
                        RectangleObstacleParams{cfg.NX/4, cfg.NY/4, cfg.NX/4+10, cfg.NY/4+10}, 
                        new_obs_name
                    );
                } else if (selected_new_obstacle_type == ObstacleType::Circle) {
                    new_obs_name = "Circle " + std::to_string(cfg.obstacles.size() + 1);
                    // Параметры по умолчанию для круга: центр в (NX/2, NY/2), радиус = Ly/10
                    cfg.obstacles.emplace_back(
                        ObstacleType::Circle, 
                        CircleObstacleParams{cfg.NX/2, cfg.NY/2, cfg.Ly/10.0}, 
                        new_obs_name
                    );
                }
                ImGui::CloseCurrentPopup();
            }
            ImGui::SameLine();
            if (ImGui::Button("Cancel", ImVec2(120, 0))) {
                ImGui::CloseCurrentPopup();
            }
            ImGui::EndPopup();
        }

        ImGui::Spacing();

                // Отображаем список препятствий и даем возможность удалить любое
                int obstacle_to_remove_idx = -1; // Индекс препятствия для удаления
                for (int k = 0; k < static_cast<int>(cfg.obstacles.size()); ++k) { // Используем int для k для сравнения с obstacle_to_remove_idx
                    ObstacleConfig& obs_cfg_ref = cfg.obstacles[k]; // Берем ссылку для модификации имени
        
                    ImGui::PushID(obs_cfg_ref.imgui_id); // Используем уникальный ID для каждого препятствия
        
                    // Имя препятствия (можно редактировать)
                    char name_buf[128]; // Буфер для имени
                    strncpy(name_buf, obs_cfg_ref.name.c_str(), sizeof(name_buf) - 1);
                    name_buf[sizeof(name_buf) - 1] = '\0'; // Гарантируем нуль-терминацию
                    if (ImGui::InputText("##ObsName", name_buf, sizeof(name_buf))) {
                        obs_cfg_ref.name = name_buf;
                    }
        
                    ImGui::SameLine();
                    ImGui::Text("(%s)", ObstacleTypeToString(obs_cfg_ref.type)); // Показываем тип
        
                    ImGui::SameLine(ImGui::GetWindowWidth() - 100); // Кнопка удаления справа
                    if (ImGui::Button("Remove")) {
                        obstacle_to_remove_idx = k;
                    }
                    
                    ImGui::Checkbox("Enabled", &obs_cfg_ref.enabled);
        
                    if (obs_cfg_ref.enabled) {
                        // Используем std::visit для отображения параметров в зависимости от типа
                        std::visit([&](auto& params) { // params здесь auto&, т.к. мы можем их менять через UI
                            using T_params = std::decay_t<decltype(params)>; // Получаем чистый тип
        
                            if constexpr (std::is_same_v<T_params, RectangleObstacleParams>) {
                                ImGui::InputInt("i0", &params.i0); ImGui::InputInt("j0", &params.j0);
                                ImGui::InputInt("i1", &params.i1); ImGui::InputInt("j1", &params.j1);
                            } else if constexpr (std::is_same_v<T_params, CircleObstacleParams>) {
                                // Параметры для круга: центр в индексах, радиус физический
                                ImGui::InputInt("Center i", &params.center_i); 
                                ImGui::InputInt("Center j", &params.center_j);
                                ImGui::InputDouble("Radius (physical units)", &params.radius_phys, 0.0, 0.0, "%.4f");
                            }
                            // Сюда можно будет добавить else if для других типов препятствий
                        }, obs_cfg_ref.params);
                    }
                    ImGui::Separator();
                    ImGui::PopID();
                }
        
                if (obstacle_to_remove_idx != -1) {
                    cfg.obstacles.erase(cfg.obstacles.begin() + obstacle_to_remove_idx);
                }

        ImGui::Separator();
        if (ImGui::Button("Start Simulation")) { started = true; }
        ImGui::End();

        ImGui::Render();
        int w_cfg, h_cfg; glfwGetFramebufferSize(cfgWindow, &w_cfg, &h_cfg); glViewport(0, 0, w_cfg, h_cfg);
        glClearColor(0.1f, 0.1f, 0.1f, 1.0f); glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(cfgWindow);
    }

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    glfwDestroyWindow(cfgWindow);

    if (!started) { 
        ImGui::DestroyContext(); 
        glfwTerminate(); 
        return 0; 
    }

    std::cout << "Simulation Setup:" << std::endl;
    std::cout << " Grid: " << cfg.NX << " x " << cfg.NY << ", Domain: " << cfg.Lx << " x " << cfg.Ly << std::endl;
    std::cout << " Solver: " << (cfg.solverType == SolverChoice::VelocityPressure ? "Velocity-Pressure" : "Vorticity-Streamfunction") << std::endl;
    std::cout << " Flow: " << (cfg.turbulenceChoice == static_cast<int>(cfd::TurbulenceModelType::None) ? "Laminar" : "k-epsilon Turbulent") << std::endl;
    std::cout << " Re = " << (cfg.umax * cfg.Ly / cfg.nu) << " (based on L=Ly, U=umax)" << std::endl; // cfg.Ly более характерный размер
    
    cfd::Geometry geom(cfg.NX, cfg.NY, cfg.Lx, cfg.Ly);
    std::cout << "Processing " << cfg.obstacles.size() << " defined obstacle(s)..." << std::endl;
    for (const auto& obs_config_item : cfg.obstacles) { // Изменил имя переменной цикла для ясности
        if (obs_config_item.enabled) {
            // Используем obs_config_item.user_name, как определено в ObstacleConfig
            std::cout << "  Obstacle Name: " << (obs_config_item.name.empty() ? "(unnamed)" : obs_config_item.name.c_str()) 
                      << ", Type: " << ObstacleTypeToString(obs_config_item.type) << std::endl;

            std::visit([&](const auto& params) { 
                using T_params = std::decay_t<decltype(params)>;
                if constexpr (std::is_same_v<T_params, RectangleObstacleParams>) {
                    // Валидация координат прямоугольника
                    int i0 = std::max(0, std::min(cfg.NX - 1, params.i0));
                    int i1 = std::max(i0, std::min(cfg.NX - 1, params.i1)); 
                    int j0 = std::max(0, std::min(cfg.NY - 1, params.j0));
                    int j1 = std::max(j0, std::min(cfg.NY - 1, params.j1)); 

                    // geom.add_rectangle ожидает i0,j0 (включительно) и i1,j1 (включительно)
                    // Убедимся, что i0 <= i1 и j0 <= j1 (для add_rectangle это нормально, если i0=i1 - линия)
                    // Но для физического препятствия лучше i0 < i1 и j0 < j1
                    if (i0 < i1 && j0 < j1) { // Препятствие должно иметь хотя бы 1x1 внутреннюю ячейку (т.е. 2x2 по индексам)
                                              // Или, если i0=i1, то это вертикальная линия толщиной в одну ячейку.
                                              // geom.add_rectangle(i0,j0,i1,j1) обработает это.
                         geom.add_rectangle(static_cast<std::size_t>(i0), 
                                            static_cast<std::size_t>(j0), 
                                            static_cast<std::size_t>(i1), 
                                            static_cast<std::size_t>(j1));
                         std::cout << "    Added Rectangle: i=[" << i0 << "," << i1 
                                   << "], j=[" << j0 << "," << j1 << "]" << std::endl;
                    } else {
                        std::cout << "    Warning: Rectangle with invalid/degenerate dimensions (i0="<<i0<<", i1="<<i1<<", j0="<<j0<<", j1="<<j1<<") not added or will be a line/point." << std::endl;
                    }
                } else if constexpr (std::is_same_v<T_params, CircleObstacleParams>) {
                    if (params.radius_phys > 1e-9) { // Проверка на положительный радиус
                        // Проверка, что центр круга (в индексах) внутри сетки
                        if (params.center_i >= 0 && params.center_i < cfg.NX &&
                            params.center_j >= 0 && params.center_j < cfg.NY) {
                            geom.add_circle(params.center_i, 
                                            params.center_j, 
                                            params.radius_phys);
                            std::cout << "    Added Circle: center_cell_idx=(" << params.center_i 
                                      << "," << params.center_j << "), R_phys=" << params.radius_phys << std::endl;
                        } else {
                            std::cout << "    Warning: Circle center indices (" << params.center_i << "," << params.center_j 
                                      << ") are outside mesh bounds [" << cfg.NX << "," << cfg.NY << "]. Circle not added." << std::endl;
                        }
                    } else {
                         std::cout << "    Warning: Circle with zero/negative radius (" << params.radius_phys << ") not added." << std::endl;
                    }
                }
            }, obs_config_item.params);
        } else {
            std::cout << "  Obstacle Name: " << (obs_config_item.name.empty() ? "(unnamed)" : obs_config_item.name.c_str()) 
                      << " (Type: " << ObstacleTypeToString(obs_config_item.type) << ") - Disabled." << std::endl;
        }
    }

    std::unique_ptr<cfd::Solver> solver_ptr;
    cfd::VelocityPressureSolver* vp_solver_raw_ptr = nullptr;
    cfd::VorticityStreamfunctionSolver* vs_solver_raw_ptr = nullptr;

    if (cfg.solverType == SolverChoice::VelocityPressure) {
        solver_ptr = std::make_unique<cfd::VelocityPressureSolver>(
            geom, cfg.rho, cfg.nu,
            static_cast<cfd::TurbulenceModelType>(cfg.turbulenceChoice),
            cfg.umax, cfg.inletTurbIntensity, cfg.inletLengthScaleFactor,
            static_cast<cfd::PoissonType>(cfg.poisson_solver_type),
            cfg.cfl, cfg.poisson_sor_omega, cfg.poisson_max_iter, cfg.poisson_tol
        );
        vp_solver_raw_ptr = static_cast<cfd::VelocityPressureSolver*>(solver_ptr.get());
        vp_solver_raw_ptr->set_inlet_parabola(cfg.umax);
    } else { 
        solver_ptr = std::make_unique<cfd::VorticityStreamfunctionSolver>(
            geom, cfg.rho, cfg.nu,
            static_cast<cfd::TurbulenceModelType>(cfg.turbulenceChoice),
            cfg.umax, cfg.inletTurbIntensity, cfg.inletLengthScaleFactor,
            static_cast<cfd::PoissonType>(cfg.poisson_solver_type),
            cfg.cfl, // VS_Solver будет использовать этот CFL для вихря
            cfg.poisson_sor_omega, cfg.poisson_max_iter, cfg.poisson_tol
        );
        vs_solver_raw_ptr = static_cast<cfd::VorticityStreamfunctionSolver*>(solver_ptr.get());
    }

    cfd::Field2D<double> u_buffer(geom.mesh().nx(), geom.mesh().ny());
    cfd::Field2D<double> v_buffer(geom.mesh().nx(), geom.mesh().ny());
    cfd::Field2D<double> scalar_buffer_for_vis(geom.mesh().nx(), geom.mesh().ny(), 0.0); // Для VP, или dummy для VS

    if (vp_solver_raw_ptr) {
        u_buffer = vp_solver_raw_ptr->u();
        v_buffer = vp_solver_raw_ptr->v();
        scalar_buffer_for_vis = vp_solver_raw_ptr->p();
    } else if (vs_solver_raw_ptr) {
        u_buffer = vs_solver_raw_ptr->u_velocity_from_psi();
        v_buffer = vs_solver_raw_ptr->v_velocity_from_psi();
        scalar_buffer_for_vis = vs_solver_raw_ptr->streamfunction();
    }
    
    cfd::SolverMonitorInfo      shared_vp_monitor_info; // Для VP
    cfd::VS_SolverMonitorInfo   shared_vs_monitor_info; // Для VS
    double last_actual_dt_shared = 0.0;                 // Общий actualDt для потоков

    std::mutex data_mutex;
    std::atomic<bool> stopFlag{false};

    std::thread worker([&]() {
        simulationTime = 0.0;
        int stepCount = 0;
        while (!stopFlag.load() && simulationTime < cfg.sim_duration) {
            if (!solver_ptr) break; // Защита, если решатель не создан
            solver_ptr->step(cfg.dt_user); 

            { 
                std::lock_guard<std::mutex> lock(data_mutex);
                if (vp_solver_raw_ptr) {
                    u_buffer = vp_solver_raw_ptr->u();
                    v_buffer = vp_solver_raw_ptr->v();
                    scalar_buffer_for_vis = vp_solver_raw_ptr->p();
                    shared_vp_monitor_info = vp_solver_raw_ptr->getMonitorInfo();
                    last_actual_dt_shared = shared_vp_monitor_info.actualDt;
                } else if (vs_solver_raw_ptr) {
                    u_buffer = vs_solver_raw_ptr->u_velocity_from_psi();
                    v_buffer = vs_solver_raw_ptr->v_velocity_from_psi();
                    scalar_buffer_for_vis = vs_solver_raw_ptr->streamfunction();
                    shared_vs_monitor_info = vs_solver_raw_ptr->getMonitorInfo();
                    last_actual_dt_shared = shared_vs_monitor_info.actualDt;
                }
                if (last_actual_dt_shared > 0) { 
                    simulationTime += last_actual_dt_shared;
                }
                stepCount++;
            }
             if (stepCount % 100 == 0) {
                 std::cout << "Sim Time: " << std::fixed << std::setprecision(3) << simulationTime 
                           << " s, Step: " << stepCount 
                           << ", dt: " << std::scientific << std::setprecision(3) << last_actual_dt_shared;
                 if (vp_solver_raw_ptr) {
                     std::cout << ", P.Iter: " << shared_vp_monitor_info.pressureIterations
                               << ", P.Res: " << shared_vp_monitor_info.pressureResidual;
                 } else if (vs_solver_raw_ptr) {
                     std::cout << ", Psi.Iter: " << shared_vs_monitor_info.streamfunctionIterations
                               << ", Psi.Res: " << shared_vs_monitor_info.streamfunctionResidual;
                 }
                 std::cout << std::endl;
             }
        }
        stopFlag = true;
        std::cout << "Simulation thread finished. Total time: " << simulationTime << "s, Steps: " << stepCount << std::endl;
    });

    FlowVisualizer vis(geom, u_buffer, v_buffer, scalar_buffer_for_vis, geom.tags());

    if (cfg.solverType == SolverChoice::VorticityStreamfunction) {
        vis.showPressure_ = false;    // Давление не рисуем
        vis.showStreamlines_ = false;  // Линии тока включаем по умолчанию для VS
        vis.showVelocity_ = true;    // Скорости выключаем, чтобы не мешали линиям тока
    } else { // VelocityPressure
        vis.showPressure_ = true;     // Давление включаем по умолчанию для VP
        vis.showStreamlines_ = false; // Линии тока не рисуем (пока)
        vis.showVelocity_ = true;     // Скорости включаем
    }

    int winWidth = std::max(600, std::min(1600, cfg.NX * 8)); 
    int winHeight = std::max(400, std::min(1000, cfg.NY * 8));
    if (!vis.init(winWidth, winHeight, "CFD Visualization")) {
        std::cerr << "Failed to initialize FlowVisualizer" << std::endl;
        stopFlag = true; if(worker.joinable()) worker.join();
        ImGui::DestroyContext(); glfwTerminate(); return -1; 
    }

    glfwMakeContextCurrent(vis.getWindowHandle()); 
    ImGui_ImplGlfw_InitForOpenGL(vis.getWindowHandle(), true);
    ImGui_ImplOpenGL3_Init("#version 130");

    cfd::SolverMonitorInfo      display_vp_monitor_info; 
    cfd::VS_SolverMonitorInfo   display_vs_monitor_info;
    double display_actual_dt = 0.0;

    while (vis.shouldRun()) { // Убрал stopFlag из условия цикла окна, чтобы можно было смотреть результат после остановки потока
        glfwPollEvents();
        
        bool current_stop_flag_val = stopFlag.load(std::memory_order_relaxed);

        { 
            std::lock_guard<std::mutex> lock(data_mutex);
            // Копируем данные для отображения даже если поток остановлен, чтобы видеть последнее состояние
            if (cfg.solverType == SolverChoice::VelocityPressure) {
                display_vp_monitor_info = shared_vp_monitor_info;
            } else {
                display_vs_monitor_info = shared_vs_monitor_info;
            }
            display_actual_dt = last_actual_dt_shared;
            // u_buffer, v_buffer, p_buffer обновляются в рабочем потоке,
            // FlowVisualizer использует их по ссылкам.
        }

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        ImGui::Begin("Solver Monitor");
        ImGui::Text("Simulation Time: %.3f s / %.1f s", simulationTime, cfg.sim_duration);
        ImGui::Text("Actual dt: %.3e s", display_actual_dt);
        ImGui::Separator();
        const char* poissonSolverName = cfg.poisson_solver_type == static_cast<int>(cfd::PoissonType::Jacobi) ? "Jacobi" : "SOR";
        if (cfg.solverType == SolverChoice::VelocityPressure) {
            ImGui::Text("Pressure Solver (%s):", poissonSolverName);
            ImGui::Text(" Last Iterations: %u / %u", display_vp_monitor_info.pressureIterations, cfg.poisson_max_iter);
            ImGui::Text(" Last Residual: %.3e (Tol: %.1e)", display_vp_monitor_info.pressureResidual, cfg.poisson_tol);
            if (display_vp_monitor_info.pressureIterations >= cfg.poisson_max_iter && 
                display_vp_monitor_info.pressureResidual > cfg.poisson_tol && // Проверяем оба условия
                simulationTime > 0) { // Не показываем в самом начале
                 ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.3f, 0.3f, 1.0f));
                 ImGui::Text("PRESSURE SOLVER DID NOT CONVERGE!");
                 ImGui::PopStyleColor();
            }
        } else { 
            ImGui::Text("Streamfunction Solver (%s):", poissonSolverName);
            ImGui::Text(" Last Iterations: %u / %u", display_vs_monitor_info.streamfunctionIterations, cfg.poisson_max_iter);
            ImGui::Text(" Last Residual: %.3e (Tol: %.1e)", display_vs_monitor_info.streamfunctionResidual, cfg.poisson_tol);
             if (display_vs_monitor_info.streamfunctionIterations >= cfg.poisson_max_iter && 
                 display_vs_monitor_info.streamfunctionResidual > cfg.poisson_tol && // Проверяем оба условия
                 simulationTime > 0) { // Не показываем в самом начале
                 ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.3f, 0.3f, 1.0f));
                 ImGui::Text("STREAMFUNCTION SOLVER DID NOT CONVERGE!");
                 ImGui::PopStyleColor();
            }
        }
        ImGui::Separator();
        if (current_stop_flag_val && simulationTime >= cfg.sim_duration) {
            ImGui::TextColored(ImVec4(0.2f, 1.0f, 0.2f, 1.0f), "Simulation Target Time Reached.");
        } else if (current_stop_flag_val) {
             ImGui::TextColored(ImVec4(1.0f, 1.0f, 0.2f, 1.0f), "Simulation Thread Finished/Stopped.");
        }
        ImGui::End();

        ImGui::Begin("Visualization Settings");
        ImGui::Checkbox("Show Velocity Vectors", &vis.showVelocity_);

        if (cfg.solverType == SolverChoice::VelocityPressure) {
            ImGui::Checkbox("Show Pressure Field", &vis.showPressure_);
        } else {
            ImGui::Checkbox("Show Streamfunction lines", &vis.showStreamlines_);
        }
        
        if (vis.showVelocity_) {
            ImGui::SliderInt("Velocity Samples (Y)", &vis.maxSamplesY_, 5, 100);
            ImGui::SliderFloat("Arrow Scale", &vis.arrowScale_, 0.001f, 0.2f, "%.4f");
            ImGui::SliderFloat("Head Length Factor", &vis.headLengthFactor_, 0.001f, 0.2f, "%.4f");
            ImGui::SliderFloat("Head Width Factor",  &vis.headWidthFactor_,  0.001f, 0.2f, "%.4f");
            ImGui::Separator();
        }
        if (cfg.solverType == SolverChoice::VelocityPressure && vis.showPressure_) {
             ImGui::Text("Pressure Min: %.3f", vis.minScalarValue_); 
             ImGui::Text("Pressure Max: %.3f", vis.maxScalarValue_);
        } else if (cfg.solverType == SolverChoice::VorticityStreamfunction && vis.showStreamlines_) {
            ImGui::SliderInt("Num Streamline Levels", &vis.numStreamlineLevels_, 3, 50);
            ImGui::SliderFloat("Streamline Thickness", &vis.streamlineThickness_, 0.5f, 5.0f);
        }

        ImGui::End();

        int display_w, display_h;
        glfwGetFramebufferSize(vis.getWindowHandle(), &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        vis.renderOpenGLScene(); 
        
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(vis.getWindowHandle());

        if (current_stop_flag_val && glfwWindowShouldClose(vis.getWindowHandle())) {
             break; // Если поток завершен и окно закрывается
        }
    }

    if (!stopFlag.load()) { // Если вышли из цикла визуализации, а поток еще работает
      stopFlag = true; // Сигнализируем потоку остановиться
    }
    if (worker.joinable()) { 
        worker.join(); 
    }

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    vis.shutdown(); 

    std::cout << "Application terminated cleanly." << std::endl;
    return 0;
}