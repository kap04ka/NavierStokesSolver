#pragma once

#include "core/Geometry.hpp"
#include "core/Field2D.hpp"
#include "core/CellTag.hpp"
#include <vector>
#include <string>

// Подключаем заголовки GLFW и GLEW здесь или в cpp
#include <GL/glew.h> // GLEW должен быть до GLFW
#include <GLFW/glfw3.h>

class FlowVisualizer {
public:
    FlowVisualizer(const cfd::Geometry& geom,
                   const cfd::Field2D<double>& u,
                   const cfd::Field2D<double>& v,
                   const cfd::Field2D<double>& p_or_psi,
                   const cfd::Field2D<cfd::CellTag>& tag);
    ~FlowVisualizer();

    bool init(int width, int height, const char* title);
    bool shouldRun() const;
    void renderOpenGLScene();
    void shutdown();

    // Геттер для хендла окна GLFW
    [[nodiscard]] GLFWwindow* getWindowHandle() const { return window_; }

    bool showVelocity_{true};
    bool showPressure_{true}; 
    bool showStreamlines_{true};

    int maxSamplesY_{20};
    float arrowScale_{0.05f};
    float headLengthFactor_{0.05f};
    float headWidthFactor_{0.05f};

    int numStreamlineLevels_ = 20;  // Количество уровней для линий тока
    float streamlineColorR_ = 0.7f;
    float streamlineColorG_ = 0.7f;
    float streamlineColorB_ = 0.7f;
    float streamlineThickness_ = 1.0f;
    
    // Для отображения диапазона давления в ImGui
    double minScalarValue_{0.0};
    double maxScalarValue_{1.0};


private:
    // Настройка проекции OpenGL
    void setupProjection(int width, int height);
    void updateSampling(); // Обновление точек для отрисовки векторов скорости

    // Методы отрисовки отдельных компонентов сцены
    void drawBoundary();
    void drawObstacles();
    void drawVelocityField();
    void drawPressureField(); // Отдельный метод для давления
    void drawStreamlines();

    // Вспомогательная функция для расчета диапазона давления
    void getColorMap(float value, float& r, float& g, float& b); // Вспомогательная функция для цветовой карты
    void calculateScalarFieldRange();

    std::pair<double, double> interpolate_isoline_point(
        double x1, double y1, double val1,
        double x2, double y2, double val2,
        double level_val) const;

    double get_value_at_node(std::size_t i_node, std::size_t j_node) const;

    // Ссылки на данные CFD
    const cfd::Geometry& geom_;
    const cfd::Field2D<double>& u_;
    const cfd::Field2D<double>& v_;
    const cfd::Field2D<double>& scalar_field_; // Ссылка на давление
    const cfd::Field2D<cfd::CellTag>& tag_;

    // Внутренние параметры визуализации
    std::vector<int> xs_;
    std::vector<int> ys_;

    GLFWwindow* window_{nullptr}; // Хендл окна GLFW
};