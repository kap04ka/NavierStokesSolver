#include "visualizer/FlowVisualizer.hpp"
#include "core/Geometry.hpp"
#include "solvers/base/Solver.hpp" // Базовый класс Solver
#include "solvers/velocity_pressure/VelocityPressureSolver.hpp"
#include "solvers/vorticity_streamfunction/VorticityStreamfunctionSolver.hpp" // Подключаем новый решатель
#include "utils/export_utils.hpp" 

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
#include <chrono>

// Тип решателя
enum class SolverChoice { VelocityPressure, VorticityStreamfunction };

enum class ObstacleType { Rectangle, Circle };

enum class ParallelizationChoice { Sequential, OpenMP, CUDA, MPI };

const char* ObstacleTypeToString(ObstacleType type) {
    switch (type) {
        case ObstacleType::Rectangle: return "Rectangle";
        case ObstacleType::Circle:    return "Circle";
        default:                      return "Unknown";
    }
}

const char* ParallelizationChoiceToString(ParallelizationChoice mode) {
    switch (mode) {
        case ParallelizationChoice::Sequential: return "Sequential";
        case ParallelizationChoice::OpenMP:     return "OpenMP";
        case ParallelizationChoice::CUDA:       return "CUDA";
        default:                                return "Unknown";
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
    int NX = 120;
    int NY = 61;
    double Lx = 2.0;
    double Ly = 1.0;

    // Препятствие
    std::vector<ObstacleConfig> obstacles;

    // Физические свойства
    double rho = 1.0;
    double nu = 1e-3; // Кинематическая вязкость
    double umax = 1.0; // Макс. скорость на входе (для профиля и Q)

    // Общие параметры симуляции
    double sim_duration = 3.0;
    double dt_user = 0.001;

    // Выбор решателя
    SolverChoice solverType = SolverChoice::VelocityPressure;

    // Параметры, которые могут быть специфичны или иметь разные оптимальные значения
    double cfl = 0.4; // Общий CFL, будет интерпретирован решателем (0.4 для VP, 0.5 для VS по умолчанию)

    // Параметры для решателя типа Пуассона (давление в VP, функция тока в VS)
    int poisson_solver_type = static_cast<int>(cfd::PoissonType::SOR);
    double poisson_sor_omega = 1.7;
    unsigned poisson_max_iter = 500; 
    double poisson_tol = 1e-3;

    // Параметры турбулентности (общие для обоих, если они ее поддерживают)
    int turbulenceChoice = static_cast<int>(cfd::TurbulenceModelType::None); // 0=None, 1=KEpsilon
    double inletTurbIntensity = 0.05;
    double inletLengthScaleFactor = 0.07;

    // Настройки параллельности
    ParallelizationChoice parallelModePoisson = ParallelizationChoice::Sequential; 

    SimulationConfig() {
        ObstacleConfig::next_imgui_id = 0;
        parallelModePoisson = ParallelizationChoice::Sequential;
    }
};


int main() {
    SimulationConfig cfg;
    bool started = false;
    double simulationTime = 0.0;

    // Переменные для измерения реального времени работы worker-а
    std::chrono::high_resolution_clock::time_point wall_time_start;
    std::chrono::high_resolution_clock::time_point wall_time_end;
    bool wall_time_measured = false;
    double total_real_calc_time_sec = 0.0;
    bool results_saved_this_session = false;

    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return -1;
    }

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    ImGui::StyleColorsLight();

    const char* font_path = "ofont.ru_ISOCPEUR.ttf";
    float font_size = 16.0f;
    const ImWchar* ranges = io.Fonts->GetGlyphRangesCyrillic();
    io.Fonts->AddFontFromFileTTF(font_path, font_size, nullptr, ranges);

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
    GLFWwindow* cfgWindow = glfwCreateWindow(1920, 1080, "Настройка Симуляции", nullptr, nullptr);
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

        ImGui::Begin("Конфигурация симуляции", nullptr, ImGuiWindowFlags_NoCollapse);
        SolverChoice previousSolverType = cfg.solverType; // Сохраняем предыдущий тип

        // Секция 1: Тип солвера
        ImGui::TextColored(ImVec4(0.2f, 0.5f, 0.8f, 1.0f), "Тип солвера");
        ImGui::Separator();
        ImGui::RadioButton("Скорость-Давление", reinterpret_cast<int*>(&cfg.solverType), static_cast<int>(SolverChoice::VelocityPressure)); 
        ImGui::SameLine();
        ImGui::RadioButton("Вихрь-Функция тока", reinterpret_cast<int*>(&cfg.solverType), static_cast<int>(SolverChoice::VorticityStreamfunction));
        ImGui::Dummy(ImVec2(0, 10));
        
        // Секция 2: Расчетная область
        ImGui::TextColored(ImVec4(0.2f, 0.5f, 0.8f, 1.0f), "Расчетная область");
        ImGui::Separator();
        ImGui::InputInt("NX (ячеек по X)", &cfg.NX);
        ImGui::InputInt("NY (ячеек по Y)", &cfg.NY);
        ImGui::InputDouble("Lx (длина области)", &cfg.Lx);
        ImGui::InputDouble("Ly (высота области)", &cfg.Ly);
        ImGui::Dummy(ImVec2(0, 10));

        // Секция 3: Параметры жидкости
        ImGui::TextColored(ImVec4(0.2f, 0.5f, 0.8f, 1.0f), "Параметры жидкости");
        ImGui::Separator();
        ImGui::InputDouble("Плотность (rho)", &cfg.rho);
        ImGui::InputDouble("Кинематическая вязкость (nu)", &cfg.nu, 0.0, 0.0, "%.1e");
        ImGui::InputDouble("Umax (макс. скорость)", &cfg.umax);
        ImGui::Dummy(ImVec2(0, 10));

        // Секция 4: Настройки симуляции
        ImGui::TextColored(ImVec4(0.2f, 0.5f, 0.8f, 1.0f), "Настройки симуляции");
        ImGui::Separator();
        ImGui::InputDouble("Шаг по времени (dt)", &cfg.dt_user, 0.0, 0.0, "%.1e");
        ImGui::InputDouble("Длительность симуляции", &cfg.sim_duration);
        ImGui::InputDouble("CFL число", &cfg.cfl, 0.0, 0.0, "%.2f");
        ImGui::Dummy(ImVec2(0, 10));

        // Секция 5: Решатель Пуассона
        ImGui::TextColored(ImVec4(0.2f, 0.5f, 0.8f, 1.0f), "Решатель Пуассона");
        ImGui::Separator();
        static const char* poissonItems[] = {"Якоби", "SOR"};
        ImGui::Combo("Тип решателя##Poisson", &cfg.poisson_solver_type, poissonItems, IM_ARRAYSIZE(poissonItems));
        if (cfg.poisson_solver_type == static_cast<int>(cfd::PoissonType::SOR)) {
            ImGui::InputDouble("Параметр релаксации SOR", &cfg.poisson_sor_omega, 0.0, 0.0, "%.2f");
        }
        ImGui::InputInt("Макс. итераций##Poisson", (int*)&cfg.poisson_max_iter);
        ImGui::InputDouble("Допуск сходимости##Poisson", &cfg.poisson_tol, 0.0, 0.0, "%.1e");
        ImGui::Dummy(ImVec2(0, 10));

        // Секция 6: Турбулентность
        ImGui::TextColored(ImVec4(0.2f, 0.5f, 0.8f, 1.0f), "Модель турбулентности");
        ImGui::Separator();
        ImGui::RadioButton("Ламинарный поток##flow", &cfg.turbulenceChoice, static_cast<int>(cfd::TurbulenceModelType::None)); 
        ImGui::SameLine();
        ImGui::RadioButton("k-epsilon##flow", &cfg.turbulenceChoice, static_cast<int>(cfd::TurbulenceModelType::KEpsilon));
        if (cfg.turbulenceChoice != static_cast<int>(cfd::TurbulenceModelType::None)) {
            ImGui::Indent();
            ImGui::InputDouble("Турбулентность на входе", &cfg.inletTurbIntensity, 0.0,0.0, "%.3f");
            cfg.inletTurbIntensity = std::max(0.0, std::min(1.0, cfg.inletTurbIntensity));
            ImGui::InputDouble("Масштаб длины на входе", &cfg.inletLengthScaleFactor, 0.0,0.0, "%.3f");
            cfg.inletLengthScaleFactor = std::max(0.001, cfg.inletLengthScaleFactor);
            ImGui::Unindent();
        }
        ImGui::Dummy(ImVec2(0, 10));
        
        // Секция 7: Препятствия
        ImGui::TextColored(ImVec4(0.2f, 0.5f, 0.8f, 1.0f), "Препятствия");
        ImGui::Separator();
        if (ImGui::Button("Добавить препятствие", ImVec2(-1, 30))) {
            ImGui::OpenPopup("AddObstacleTypePopup");
        }       

        // Модальное окно для выбора типа нового препятствия
        if (ImGui::BeginPopupModal("AddObstacleTypePopup", NULL, ImGuiWindowFlags_AlwaysAutoResize)) {
            ImGui::Text("Выберите какое препятствие добавить:");
            // Статическая переменная для хранения выбора в модальном окне
            static ObstacleType selected_new_obstacle_type = ObstacleType::Rectangle; 
            
            ImGui::RadioButton("Прямоугольник", reinterpret_cast<int*>(&selected_new_obstacle_type), static_cast<int>(ObstacleType::Rectangle));
            ImGui::RadioButton("Круг",    reinterpret_cast<int*>(&selected_new_obstacle_type), static_cast<int>(ObstacleType::Circle));
            // Сюда можно будет добавить другие типы препятствий в будущем

            ImGui::Separator();
            if (ImGui::Button("Добавить", ImVec2(150, 0))) {
                std::string new_obs_name;
                if (selected_new_obstacle_type == ObstacleType::Rectangle) {
                    new_obs_name = "Прямоугольник " + std::to_string(cfg.obstacles.size() + 1);
                    cfg.obstacles.emplace_back(
                        ObstacleType::Rectangle, 
                        RectangleObstacleParams{cfg.NX/4, cfg.NY/4, cfg.NX/4+10, cfg.NY/4+10}, 
                        new_obs_name
                    );
                } else if (selected_new_obstacle_type == ObstacleType::Circle) {
                    new_obs_name = "Круг " + std::to_string(cfg.obstacles.size() + 1);
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
                    if (ImGui::Button("Удалить")) {
                        obstacle_to_remove_idx = k;
                    }
                    
                    ImGui::Checkbox("Включено", &obs_cfg_ref.enabled);
        
                    if (obs_cfg_ref.enabled) {
                        // Используем std::visit для отображения параметров в зависимости от типа
                        std::visit([&](auto& params) { // params здесь auto&, т.к. мы можем их менять через UI
                            using T_params = std::decay_t<decltype(params)>; // Получаем чистый тип
        
                            if constexpr (std::is_same_v<T_params, RectangleObstacleParams>) {
                                ImGui::InputInt("i0", &params.i0); ImGui::InputInt("j0", &params.j0);
                                ImGui::InputInt("i1", &params.i1); ImGui::InputInt("j1", &params.j1);
                            } else if constexpr (std::is_same_v<T_params, CircleObstacleParams>) {
                                // Параметры для круга: центр в индексах, радиус физический
                                ImGui::InputInt("Центр i", &params.center_i); 
                                ImGui::InputInt("Центр j", &params.center_j);
                                ImGui::InputDouble("Радиус (физический)", &params.radius_phys, 0.0, 0.0, "%.4f");
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

        // Секция 8: Параллелизация
        ImGui::TextColored(ImVec4(0.2f, 0.5f, 0.8f, 1.0f), "Параллелизация");
        ImGui::Separator();
        ImGui::RadioButton("Последовательный", reinterpret_cast<int*>(&cfg.parallelModePoisson), static_cast<int>(ParallelizationChoice::Sequential)); 
        ImGui::SameLine();
        ImGui::RadioButton("OpenMP", reinterpret_cast<int*>(&cfg.parallelModePoisson), static_cast<int>(ParallelizationChoice::OpenMP));
        ImGui::SameLine();
        ImGui::RadioButton("CUDA", reinterpret_cast<int*>(&cfg.parallelModePoisson), static_cast<int>(ParallelizationChoice::CUDA));
        ImGui::Dummy(ImVec2(0, 15));

        // Кнопка запуска
        ImGui::Dummy(ImVec2(0, 20));
        if (ImGui::Button("НАЧАТЬ РАСЧЕТ", ImVec2(-1, 50))) {
            started = true;
        }
        ImGui::SetItemDefaultFocus();
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

    cfd::ParallelizationMode pmode_for_solvers;

    if (cfg.parallelModePoisson == ParallelizationChoice::Sequential) pmode_for_solvers = cfd::ParallelizationMode::Sequential;
    else if (cfg.parallelModePoisson == ParallelizationChoice::OpenMP) pmode_for_solvers = cfd::ParallelizationMode::OpenMP;
    else if (cfg.parallelModePoisson == ParallelizationChoice::CUDA) pmode_for_solvers = cfd::ParallelizationMode::CUDA;
    else pmode_for_solvers = cfd::ParallelizationMode::Sequential;

    std::unique_ptr<cfd::Solver> solver_ptr;
    cfd::VelocityPressureSolver* vp_solver_raw_ptr = nullptr;
    cfd::VorticityStreamfunctionSolver* vs_solver_raw_ptr = nullptr;

    if (cfg.solverType == SolverChoice::VelocityPressure) {
        solver_ptr = std::make_unique<cfd::VelocityPressureSolver>(
            geom, cfg.rho, cfg.nu,
            static_cast<cfd::TurbulenceModelType>(cfg.turbulenceChoice),
            cfg.umax, cfg.inletTurbIntensity, cfg.inletLengthScaleFactor,
            static_cast<cfd::PoissonType>(cfg.poisson_solver_type), pmode_for_solvers,
            cfg.cfl, cfg.poisson_sor_omega, cfg.poisson_max_iter, cfg.poisson_tol
        );
        vp_solver_raw_ptr = static_cast<cfd::VelocityPressureSolver*>(solver_ptr.get());
        vp_solver_raw_ptr->set_inlet_parabola(cfg.umax);
    } else { 
        solver_ptr = std::make_unique<cfd::VorticityStreamfunctionSolver>(
            geom, cfg.rho, cfg.nu,
            static_cast<cfd::TurbulenceModelType>(cfg.turbulenceChoice),
            cfg.umax, cfg.inletTurbIntensity, cfg.inletLengthScaleFactor,
            static_cast<cfd::PoissonType>(cfg.poisson_solver_type), pmode_for_solvers,
            cfg.cfl, 
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

    wall_time_start = std::chrono::high_resolution_clock::now(); 

    std::thread worker([&]() {
        simulationTime = 0.0;
        int stepCount = 0;

        {
            std::lock_guard<std::mutex> lock(data_mutex);
            shared_vp_monitor_info = {};
            shared_vs_monitor_info = {};
            last_actual_dt_shared = 0.0;
        }

        while (!stopFlag.load(std::memory_order_acquire) && simulationTime < cfg.sim_duration) {
            if (!solver_ptr) { stopFlag.store(true, std::memory_order_release); break; }

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
        std::cout << "Simulation thread finished. Target/Actual SimTime: " 
                  << std::fixed << std::setprecision(4) << cfg.sim_duration << " / " << simulationTime
                  << " s, Steps: " << stepCount << std::endl; });

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

    int winWidth = std::max(1920, std::min(1600, cfg.NX * 8)); 
    int winHeight = std::max(1080, std::min(1000, cfg.NY * 8));
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
        
        bool current_stop_flag_val = stopFlag.load(std::memory_order_acquire);

        { 
            std::lock_guard<std::mutex> lock(data_mutex);
            // Копируем данные для отображения даже если поток остановлен, чтобы видеть последнее состояние
            if (cfg.solverType == SolverChoice::VelocityPressure) {
                display_vp_monitor_info = shared_vp_monitor_info;
            } else {
                display_vs_monitor_info = shared_vs_monitor_info;
            }
            display_actual_dt = last_actual_dt_shared;
        }

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        ImGui::Begin("Показатели расчета");
        ImGui::Text("Время симуляции: %.3f s / %.1f s", simulationTime, cfg.sim_duration);
        ImGui::Text("Актуальный шаг по времени (dt): %.3e s", display_actual_dt);
        ImGui::Separator();
        const char* poissonSolverName = cfg.poisson_solver_type == static_cast<int>(cfd::PoissonType::Jacobi) ? "Jacobi" : "SOR";
        if (cfg.solverType == SolverChoice::VelocityPressure) {
            ImGui::Text("Решатель давление-скорость (%s):", poissonSolverName);
            ImGui::Text(" Последняя итерация: %u / %u", display_vp_monitor_info.pressureIterations, cfg.poisson_max_iter);
            ImGui::Text(" Последняя точность: %.3e (Tol: %.1e)", display_vp_monitor_info.pressureResidual, cfg.poisson_tol);
            if (display_vp_monitor_info.pressureIterations >= cfg.poisson_max_iter && 
                display_vp_monitor_info.pressureResidual > cfg.poisson_tol && // Проверяем оба условия
                simulationTime > 0) { // Не показываем в самом начале
                 ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.3f, 0.3f, 1.0f));
                 ImGui::Text("Решение не сошлось!");
                 ImGui::PopStyleColor();
            }
        } else { 
            ImGui::Text("Решатель вихрь-функция тока (%s):", poissonSolverName);
            ImGui::Text(" Последняя итерация: %u / %u", display_vs_monitor_info.streamfunctionIterations, cfg.poisson_max_iter);
            ImGui::Text(" Последняя точность: %.3e (Tol: %.1e)", display_vs_monitor_info.streamfunctionResidual, cfg.poisson_tol);
             if (display_vs_monitor_info.streamfunctionIterations >= cfg.poisson_max_iter && 
                 display_vs_monitor_info.streamfunctionResidual > cfg.poisson_tol && // Проверяем оба условия
                 simulationTime > 0) { // Не показываем в самом начале
                 ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.3f, 0.3f, 1.0f));
                 ImGui::Text("Решение не сошлось!");
                 ImGui::PopStyleColor();
            }
        }
        ImGui::Separator();
        
        if (current_stop_flag_val) { // Если рабочий поток завершился
            if (!wall_time_measured) { // Вычисляем реальное время один раз
                wall_time_end = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> wall_time_duration = wall_time_end - wall_time_start;
                total_real_calc_time_sec = wall_time_duration.count();
                wall_time_measured = true;
            }
            // Отображаем сообщение о завершении
            if (simulationTime >= cfg.sim_duration - 1e-3*cfg.sim_duration) { // Добавил небольшой допуск
                ImGui::TextColored(ImVec4(0.2f, 1.0f, 0.2f, 1.0f), "Досигнуто время симуляции.");
            } else {
                 ImGui::TextColored(ImVec4(1.0f, 1.0f, 0.2f, 1.0f), "Расчет завершен.");
            }
            // Отображаем реальное время расчета
            if (wall_time_measured) {
                 ImGui::Text("Общее время расчета %.3f секунд", total_real_calc_time_sec);
            }

            ImGui::Dummy(ImVec2(0.0f, 10.0f)); // Отступ
            ImGui::Separator();
            ImGui::Text("Сохранение результатов:");
            if (!results_saved_this_session) { 
                if (ImGui::Button("Сохранить финальные результаты в CSV", ImVec2(-1, 30))) { 
                    // Формируем базовое имя файла
                    std::ostringstream time_ss_final;
                    time_ss_final << std::fixed << std::setprecision(3) << simulationTime;
                    std::string time_str_final = time_ss_final.str();
                    std::replace(time_str_final.begin(), time_str_final.end(), '.', '_');
                    
                    std::string base_filename = "cfd_output_t" + time_str_final; 
                    std::cout << "Saving final results to files with prefix: " << base_filename << "..." << std::endl;

                    // Экспорт полей скоростей (они всегда есть)
                    cfd::exportFieldToCSV(u_buffer, base_filename + "_U.csv");
                    cfd::exportFieldToCSV(v_buffer, base_filename + "_V.csv");

                    // Экспорт полей в зависимости от типа решателя
                    if (cfg.solverType == SolverChoice::VelocityPressure) {
                        cfd::exportFieldToCSV(scalar_buffer_for_vis, base_filename + "_Pressure.csv");
                    } else { // VorticityStreamfunction
                        cfd::exportFieldToCSV(scalar_buffer_for_vis, base_filename + "_Psi.csv");
                        if (vs_solver_raw_ptr) { // Убедимся, что указатель валидный
                            cfd::exportFieldToCSV(vs_solver_raw_ptr->vorticity(), base_filename + "_Omega.csv");
                        }
                    }

                    // Экспорт полей турбулентности (если режим турбулентный)
                    if (static_cast<cfd::TurbulenceModelType>(cfg.turbulenceChoice) != cfd::TurbulenceModelType::None && 
                        solver_ptr && solver_ptr->getTurbulenceModel()) 
                    {
                        auto turb_model = solver_ptr->getTurbulenceModel();
                        // Проверяем, что указатели на поля не нулевые
                        if (turb_model->k()) cfd::exportFieldToCSV(*(turb_model->k()), base_filename + "_k.csv");
                        if (turb_model->epsilon()) cfd::exportFieldToCSV(*(turb_model->epsilon()), base_filename + "_epsilon.csv");
                        if (turb_model->nu_t()) cfd::exportFieldToCSV(*(turb_model->nu_t()), base_filename + "_nut.csv");
                    }

                    // Экспорт информации о сетке в отдельный файл
                    std::ofstream mesh_file(base_filename + "_mesh_info.txt");
                    if (mesh_file.is_open()) {
                        mesh_file << "NX: " << geom.mesh().nx() << std::endl;
                        mesh_file << "NY: " << geom.mesh().ny() << std::endl;
                        mesh_file << "Lx: " << geom.mesh().Lx() << std::endl;
                        mesh_file << "Ly: " << geom.mesh().Ly() << std::endl;
                        mesh_file << "dx: " << geom.mesh().dx() << std::endl;
                        mesh_file << "dy: " << geom.mesh().dy() << std::endl;
                        mesh_file.close();
                    }
                    results_saved_this_session = true;
                    std::cout << "Final results saved." << std::endl;
                }
            } else { 
                ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(0.2f, 1.0f, 0.2f, 1.0f)); // Зеленый цвет для подтверждения
                ImGui::Text("Результаты успешно сохранены!");
                ImGui::PopStyleColor();
            }
        } else { // Если поток еще работает
            ImGui::Text("Идет симуляция");
        }
        ImGui::End(); 

        ImGui::Begin("Настройки визуализации");
        ImGui::Checkbox("Показывать поле скоростей", &vis.showVelocity_);

        if (cfg.solverType == SolverChoice::VelocityPressure) {
            ImGui::Checkbox("Показывать поле давления", &vis.showPressure_);
        } else {
            ImGui::Checkbox("Показывать изолинии функции тока", &vis.showStreamlines_);
        }
        
        if (vis.showVelocity_) {
            ImGui::SliderInt("Максимальное кол-во стрелок по Y", &vis.maxSamplesY_, 5, 100);
            ImGui::SliderFloat("Масштаб стрелок", &vis.arrowScale_, 0.001f, 0.2f, "%.4f");
            ImGui::SliderFloat("Множитель длины наконечника", &vis.headLengthFactor_, 0.001f, 0.2f, "%.4f");
            ImGui::SliderFloat("Множитель ширины наконечника",  &vis.headWidthFactor_,  0.001f, 0.2f, "%.4f");
            ImGui::Separator();
        }
        if (cfg.solverType == SolverChoice::VelocityPressure && vis.showPressure_) {
             ImGui::Text("Минимальное давление: %.3f", vis.minScalarValue_); 
             ImGui::Text("Максимальное давление: %.3f", vis.maxScalarValue_);
        } else if (cfg.solverType == SolverChoice::VorticityStreamfunction && vis.showStreamlines_) {
            ImGui::SliderInt("Количество изолиний уровня тока", &vis.numStreamlineLevels_, 3, 50);
            ImGui::SliderFloat("Толщина изолиний уровня тока", &vis.streamlineThickness_, 0.5f, 5.0f);
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

    if (!stopFlag.load(std::memory_order_acquire)) {
        stopFlag.store(true, std::memory_order_release); 
    }

    if (worker.joinable()) { 
        worker.join(); 

        if (!wall_time_measured && started) { 
            wall_time_end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> wall_time_duration = wall_time_end - wall_time_start;
            total_real_calc_time_sec = wall_time_duration.count();
            std::cout << "Final Real Calculation Time (measured after join): " << total_real_calc_time_sec << " seconds" << std::endl;
        }
    }

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    vis.shutdown(); 

    std::cout << "Application terminated cleanly." << std::endl;
    return 0;
}