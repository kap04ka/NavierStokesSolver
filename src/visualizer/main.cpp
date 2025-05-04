#include "visualizer/FlowVisualizer.hpp"
#include "core/Geometry.hpp"
#include "solvers/velocity_pressure/VelocityPressureSolver.hpp" // Для SolverMonitorInfo

// GLEW нужно подключать ДО GLFW
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

// Конфигурация симуляции
struct SimulationConfig {
    int NX = 120;
    int NY = 60;
    double Lx = 1.0;
    double Ly = 0.5;
    bool useObstacle = true;
    int obs_i0 = 40, obs_j0 = 20; // Координаты препятствия (индексы ячеек)
    int obs_i1 = 80, obs_j1 = 39; // Верхняя граница j на 1 меньше NY/2
    double rho = 1.0;      // Плотность (часто 1 для простоты)
    double nu = 1e-3;    // Вязкость
    double umax = 1.0;     // Макс. скорость на входе
    double cfl = 0.4;      // Число CFL
    double omega = 1.7;    // Параметр SOR
    int ptype = static_cast<int>(cfd::PoissonType::SOR); // Тип решателя Пуассона
    unsigned p_max_iter = 500; // Макс. итераций Пуассона
    double p_tol = 1e-5;       // Точность Пуассона
    double sim_duration = 10.0; // Желаемая длительность симуляции (в секундах)
    double dt_user = 0.001;   // Желаемый шаг по времени
    int turbulenceChoice = static_cast<int>(cfd::TurbulenceModelType::None); // 0=None, 1=KEpsilon
    double inletTurbIntensity = 0.05; // 5%
    double inletLengthScaleFactor = 0.07; // 7% от Ly
};

// Структура для передачи данных решателя (определена в хедере решателя)
// namespace cfd { struct SolverMonitorInfo { ... }; }

int main() {
    SimulationConfig cfg;
    bool started = false;
    double simulationTime = 0.0; // Текущее время симуляции

    // --- Первичная инициализация GLFW и ImGui ---
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW" << std::endl;
        return -1;
    }
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    ImGui::StyleColorsDark();

    // ------------------- ОКНО КОНФИГУРАЦИИ -------------------
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
    GLFWwindow* cfgWindow = glfwCreateWindow(600, 600, "Simulation Setup", nullptr, nullptr); // Увеличим высоту
    if (!cfgWindow) {
        std::cerr << "Failed to create GLFW config window" << std::endl;
        ImGui::DestroyContext();
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(cfgWindow);
    glfwSwapInterval(1);
    ImGui_ImplGlfw_InitForOpenGL(cfgWindow, true);
    ImGui_ImplOpenGL3_Init("#version 130");

    while (!started && !glfwWindowShouldClose(cfgWindow)) {
        glfwPollEvents();
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        ImGui::Begin("Simulation Configuration");
        ImGui::InputInt("NX (Grid Cells X)", &cfg.NX);
        ImGui::InputInt("NY (Grid Cells Y)", &cfg.NY);
        ImGui::InputDouble("Lx (Domain Length)", &cfg.Lx);
        ImGui::InputDouble("Ly (Domain Height)", &cfg.Ly);
        ImGui::Separator();
        ImGui::Checkbox("Use Obstacle", &cfg.useObstacle);
        if (cfg.useObstacle) {
            ImGui::InputInt("Obstacle i0 (Start X Index)", &cfg.obs_i0);
            ImGui::InputInt("Obstacle j0 (Start Y Index)", &cfg.obs_j0);
            ImGui::InputInt("Obstacle i1 (End X Index)", &cfg.obs_i1);
            ImGui::InputInt("Obstacle j1 (End Y Index)", &cfg.obs_j1);
        }
        ImGui::Separator();
        ImGui::InputDouble("rho (Density)", &cfg.rho);
        ImGui::InputDouble("nu (Viscosity)", &cfg.nu, 0.0, 0.0, "%.1e"); // Формат для вязкости
        ImGui::InputDouble("Umax (Inlet Profile)", &cfg.umax);
        ImGui::Separator();
        ImGui::Text("Solver Settings:");
        ImGui::InputDouble("CFL Safety Factor", &cfg.cfl, 0.0, 0.0, "%.2f");
        ImGui::InputDouble("User dt (Max)", &cfg.dt_user, 0.0, 0.0, "%.1e");
        ImGui::InputDouble("Simulation Duration (s)", &cfg.sim_duration);
        ImGui::Separator();
        ImGui::Text("Pressure Solver:");
        static const char* poissonItems[] = {"Jacobi","SOR"};
        ImGui::Combo("Type", &cfg.ptype, poissonItems, IM_ARRAYSIZE(poissonItems));
        if (cfg.ptype == static_cast<int>(cfd::PoissonType::SOR)) {
            ImGui::InputDouble("SOR Omega", &cfg.omega, 0.0, 0.0, "%.2f");
        }
        ImGui::InputInt("Max Iterations", (int*)&cfg.p_max_iter);
        ImGui::InputDouble("Tolerance", &cfg.p_tol, 0.0, 0.0, "%.1e");
        ImGui::Separator();
        ImGui::Text("Flow Model:");
        // Используем RadioButton для выбора режима
        ImGui::RadioButton("Laminar", &cfg.turbulenceChoice, static_cast<int>(cfd::TurbulenceModelType::None)); ImGui::SameLine();
        ImGui::RadioButton("k-epsilon", &cfg.turbulenceChoice, static_cast<int>(cfd::TurbulenceModelType::KEpsilon));
        if (cfg.turbulenceChoice != static_cast<int>(cfd::TurbulenceModelType::None)) {
            ImGui::Indent();
            ImGui::InputDouble("Inlet Intensity (0-1)", &cfg.inletTurbIntensity, 0.001, 0.01, "%.3f");
            cfg.inletTurbIntensity = std::max(0.0, std::min(1.0, cfg.inletTurbIntensity)); // Ограничиваем 0..1
            ImGui::InputDouble("Inlet L Scale Factor (vs Ly)", &cfg.inletLengthScaleFactor, 0.001, 0.01, "%.3f");
            cfg.inletLengthScaleFactor = std::max(0.001, cfg.inletLengthScaleFactor); // Ограничиваем снизу
            ImGui::Unindent();
        }
        ImGui::Separator();
        if (ImGui::Button("Start Simulation")) { started = true; }
        ImGui::End();

        ImGui::Render();
        int w, h;
        glfwGetFramebufferSize(cfgWindow, &w, &h);
        glViewport(0, 0, w, h);
        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(cfgWindow);
    }

    // --- Очистка окна конфигурации ---
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    glfwDestroyWindow(cfgWindow);
    // НЕ УДАЛЯЕМ КОНТЕКСТ ImGui и НЕ ЗАВЕРШАЕМ GLFW

    if (!started) { // Выход, если не начали симуляцию
        ImGui::DestroyContext();
        glfwTerminate();
        return 0;
    }

    // ------------------- НАСТРОЙКА СИМУЛЯЦИИ -------------------
    std::cout << "Simulation Setup:" << std::endl;
    std::cout << " Grid: " << cfg.NX << " x " << cfg.NY << std::endl;
    std::cout << " Domain: " << cfg.Lx << " x " << cfg.Ly << std::endl;
    std::cout << " Re = " << (cfg.umax * cfg.Ly / cfg.nu) << " (based on L=Ly, U=umax)" << std::endl;

    cfd::Geometry geom(cfg.NX, cfg.NY, cfg.Lx, cfg.Ly);
    if (cfg.useObstacle) {
        // Проверка индексов препятствия
        cfg.obs_i0 = std::max(0, std::min(cfg.NX - 1, cfg.obs_i0));
        cfg.obs_i1 = std::max(0, std::min(cfg.NX - 1, cfg.obs_i1));
        cfg.obs_j0 = std::max(0, std::min(cfg.NY - 1, cfg.obs_j0));
        cfg.obs_j1 = std::max(0, std::min(cfg.NY - 1, cfg.obs_j1));
        if (cfg.obs_i0 <= cfg.obs_i1 && cfg.obs_j0 <= cfg.obs_j1) {
             geom.add_rectangle(cfg.obs_i0, cfg.obs_j0, cfg.obs_i1, cfg.obs_j1);
             std::cout << " Obstacle added: i=[" << cfg.obs_i0 << "," << cfg.obs_i1 << "], j=[" << cfg.obs_j0 << "," << cfg.obs_j1 << "]" << std::endl;
        } else {
            std::cerr << "Warning: Invalid obstacle indices, obstacle not added." << std::endl;
        }
    }
    cfd::VelocityPressureSolver solver(
        geom,
        cfg.rho, 
        cfg.nu, 
        static_cast<cfd::TurbulenceModelType>(cfg.turbulenceChoice),
        cfg.umax,
        cfg.inletTurbIntensity,
        cfg.inletLengthScaleFactor,
        static_cast<cfd::PoissonType>(cfg.ptype),
        cfg.cfl, 
        cfg.omega, 
        cfg.p_max_iter, 
        cfg.p_tol
    );
    solver.set_inlet_parabola(cfg.umax);

    // --- Подготовка буферов и потока ---
    auto u_buffer = solver.u(); // Копируем начальное состояние
    auto v_buffer = solver.v();
    auto p_buffer = solver.p(); // Копируем начальное давление
    cfd::SolverMonitorInfo shared_monitor_info; // Структура для обмена
    std::mutex data_mutex;                      // Мьютекс для защиты общих данных
    std::atomic<bool> stopFlag{false};          // Флаг для остановки потока

    std::thread worker([&]() {
        simulationTime = 0.0; // Сбрасываем время симуляции
        int stepCount = 0;
        while (!stopFlag.load() && simulationTime < cfg.sim_duration) {
            // Выполняем шаг решателя
            solver.step(cfg.dt_user); // Передаем макс. желаемый dt

            // Копируем данные под мьютексом
            {
                std::lock_guard<std::mutex> lock(data_mutex);
                u_buffer = solver.u();
                v_buffer = solver.v();
                p_buffer = solver.p();
                shared_monitor_info = solver.getMonitorInfo(); // Копируем инфо мониторинга
                simulationTime += shared_monitor_info.actualDt; // Увеличиваем время на фактический dt
                stepCount++;
            }
             // Опционально выводим прогресс в консоль
             // if (stepCount % 100 == 0) {
             //     std::cout << "Sim Time: " << simulationTime << " s, Step: " << stepCount << ", dt: " << shared_monitor_info.actualDt << std::endl;
             // }
        }
        stopFlag = true; // Устанавливаем флаг, если вышли по времени
        std::cout << "Simulation finished. Total time: " << simulationTime << "s, Steps: " << stepCount << std::endl;
    });

    // ------------------- ОКНО ВИЗУАЛИЗАЦИИ и Цикл -------------------
    FlowVisualizer vis(geom, u_buffer, v_buffer, p_buffer, geom.tags()); // Передаем p_buffer
    // Адаптируем размер окна под сетку, но не слишком маленькое/большое
    int winWidth = std::max(600, std::min(1600, cfg.NX * 6));
    int winHeight = std::max(400, std::min(1000, cfg.NY * 6));
    if (!vis.init(winWidth, winHeight, "CFD Visualization")) {
        std::cerr << "Failed to initialize FlowVisualizer" << std::endl;
        stopFlag = true;
        if (worker.joinable()) worker.join();
        ImGui::DestroyContext();
        glfwTerminate();
        return -1;
    }

    // Инициализация бэкендов ImGui для окна визуализации
    glfwMakeContextCurrent(vis.getWindowHandle());
    ImGui_ImplGlfw_InitForOpenGL(vis.getWindowHandle(), true);
    ImGui_ImplOpenGL3_Init("#version 130");

    bool simulation_finished_message_shown = false;

    // Основной цикл визуализации
    while (vis.shouldRun()) { // Выходим также, если симуляция закончилась
        glfwPollEvents();

        // Получаем последние данные из потока симуляции
        cfd::SolverMonitorInfo current_info_for_display;
        // Копируем поля (если FlowVisualizer их не копирует сам)
        // Для простоты предполагаем, что FlowVisualizer использует ссылки,
        // а мьютекс гарантирует, что данные не меняются во время рендеринга
        // (это не совсем верно, рендеринг может занять время, лучше копировать)
        // Безопаснее: сделать копии полей ПОСЛЕ блокировки мьютекса
        {
            std::lock_guard<std::mutex> lock(data_mutex);
            current_info_for_display = shared_monitor_info;
            // Если бы FlowVisualizer работал с копиями:
            // u_vis_copy = u_buffer; v_vis_copy = v_buffer; p_vis_copy = p_buffer;
        }

        // Начинаем кадр ImGui
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // --- UI Окно Мониторинга ---
        ImGui::Begin("Solver Monitor");
        ImGui::Text("Simulation Time: %.3f s / %.1f s", simulationTime, cfg.sim_duration);
        ImGui::Text("Actual dt: %.3e s", current_info_for_display.actualDt);
        ImGui::Separator();
        ImGui::Text("Pressure Solver (%s):", cfg.ptype == 0 ? "Jacobi" : "SOR");
        ImGui::Text(" Last Iterations: %u / %u", current_info_for_display.pressureIterations, cfg.p_max_iter);
        ImGui::Text(" Last Residual: %.3e (Tol: %.1e)", current_info_for_display.pressureResidual, cfg.p_tol);
        if (current_info_for_display.pressureIterations >= cfg.p_max_iter) {
             ImGui::PushStyleColor(ImGuiCol_Text, ImVec4(1.0f, 0.3f, 0.3f, 1.0f));
             ImGui::Text(" POISSON DID NOT CONVERGE!");
             ImGui::PopStyleColor();
        }
        ImGui::Separator();
        if (stopFlag.load() && !simulation_finished_message_shown) { // Используем load() для atomic
            ImGui::Separator();
            ImGui::TextColored(ImVec4(0.8f, 0.8f, 0.8f, 1.0f), "Simulation Finished.");
            // simulation_finished_message_shown = true; // Раскомментировать, если нужно показать только один раз
        }
        ImGui::End();

        // --- UI Окно Настроек Визуализации ---
        ImGui::Begin("Visualization Settings");
        ImGui::Checkbox("Show Velocity Vectors", &vis.showVelocity_);
        ImGui::Checkbox("Show Pressure Field", &vis.showPressure_);
        if (vis.showVelocity_) {
            ImGui::SliderInt("Velocity Samples (Y)", &vis.maxSamplesY_, 5, 100);
            ImGui::SliderFloat("Arrow Scale", &vis.arrowScale_, 0.001f, 0.2f, "%.4f");
            ImGui::SliderFloat("Head Length Factor", &vis.headLengthFactor_, 0.001f, 0.2f, "%.4f");
            ImGui::SliderFloat("Head Width Factor",  &vis.headWidthFactor_,  0.001f, 0.2f, "%.4f");
            ImGui::Separator();
        }
        if (vis.showPressure_) {
            //  if (ImGui::Button("Recalculate Pressure Range")) {
            //       vis.resetPressureRange(); // Используем метод для сброса флага
            //  }
             // Отображаем диапазон, делаем их public или через геттеры в FlowVisualizer
             ImGui::Text("Pressure Min: %.3f", vis.minPressure_);
             ImGui::Text("Pressure Max: %.3f", vis.maxPressure_);
             ImGui::Separator();
        }
        ImGui::End();

        // --- Рендеринг OpenGL ---
        int display_w, display_h;
        glfwGetFramebufferSize(vis.getWindowHandle(), &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        // 1. Рисуем сцену CFD
        vis.renderOpenGLScene();

        // 2. Рисуем ImGui поверх
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        // 3. Обмен буферов
        glfwSwapBuffers(vis.getWindowHandle());
    }

    // --- Очистка после цикла визуализации ---
    stopFlag = true; // Сигнал рабочему потоку (на случай, если вышли по закрытию окна)
    if (worker.joinable()) {
        worker.join();
    }  // Ждем завершения

    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext(); // Уничтожаем контекст ImGui

    vis.shutdown(); // Уничтожает окно и вызывает glfwTerminate()

    std::cout << "Application terminated cleanly." << std::endl;
    return 0;
}