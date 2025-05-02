#pragma once

#include "core/Geometry.hpp"
#include "core/Field2D.hpp"
#include "core/CellTag.hpp"
#include <vector>
#include <string>

// Подключаем заголовки GLFW и GLEW здесь или в cpp
#include <GL/glew.h> // GLEW должен быть до GLFW
#include <GLFW/glfw3.h>

namespace cfd {
// Пусто, т.к. предварительные объявления не нужны
}

class FlowVisualizer {
public:
    // Конструктор теперь принимает и поле давления p
    FlowVisualizer(const cfd::Geometry& geom,
                   const cfd::Field2D<double>& u,
                   const cfd::Field2D<double>& v,
                   const cfd::Field2D<double>& p, // Поле давления
                   const cfd::Field2D<cfd::CellTag>& tag);
    ~FlowVisualizer();

    // Инициализация окна и OpenGL
    bool init(int width, int height, const char* title);
    // Проверка, должно ли окно продолжать работать
    bool shouldRun() const;
    // Метод для отрисовки только CFD сцены (без ImGui, без очистки/смены буферов)
    void renderOpenGLScene();
    // Завершение работы (уничтожение окна, завершение GLFW)
    void shutdown();

    // Геттер для хендла окна GLFW
    [[nodiscard]] GLFWwindow* getWindowHandle() const { return window_; }

    // Метод для сброса рассчитанного диапазона давления (вызывать из ImGui, если нужно)
    //void resetPressureRange() { pressureRangeCalculated_ = false; }

    // --- Публичные параметры для контроля из ImGui в main ---
    bool showVelocity_{true};
    bool showPressure_{true}; // Включаем по умолчанию для проверки
    int maxSamplesY_{20};
    float arrowScale_{0.05f};
    float headLengthFactor_{0.05f};
    float headWidthFactor_{0.05f};
    // Для отображения диапазона давления в ImGui
    double minPressure_{0.0};
    double maxPressure_{1.0};
    // --- Конец публичных параметров ---


private:
    // Настройка проекции OpenGL
    void setupProjection(int width, int height);
    // Обновление точек для отрисовки векторов скорости
    void updateSampling();
    // Методы отрисовки отдельных компонентов сцены
    void drawBoundary();
    void drawObstacles();
    void drawVelocityField();
    void drawPressureField(); // Отдельный метод для давления
    // Вспомогательная функция для цветовой карты
    void getColorMap(float value, float& r, float& g, float& b);
    // Вспомогательная функция для расчета диапазона давления
    //void calculatePressureRangeIfNeeded();
    void calculatePressureRange();

    // Ссылки на данные CFD
    const cfd::Geometry& geom_;
    const cfd::Field2D<double>& u_;
    const cfd::Field2D<double>& v_;
    const cfd::Field2D<double>& p_; // Ссылка на давление
    const cfd::Field2D<cfd::CellTag>& tag_;

    // Внутренние параметры визуализации
    std::vector<int> xs_;
    std::vector<int> ys_;
    //bool pressureRangeCalculated_{false}; // Флаг, что диапазон посчитан

    GLFWwindow* window_{nullptr}; // Хендл окна GLFW
};