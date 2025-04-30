#pragma once

#include "core/Geometry.hpp"
#include "core/Field2D.hpp"
#include "core/CellTag.hpp"
#include <vector>

struct GLFWwindow;

namespace cfd {
class VelocityPressureSolver;
}

class FlowVisualizer {
public:
    FlowVisualizer(const cfd::Geometry& geom,
                   const cfd::Field2D<double>& u,
                   const cfd::Field2D<double>& v,
                   const cfd::Field2D<cfd::CellTag>& tag);
    ~FlowVisualizer();

    bool init(int width = 800, int height = 600,
              const char* title = "CFD Visualizer");
    bool shouldRun() const;
    void render();
    void shutdown();

private:
    void setupProjection();
    void setupImGui();
    void cleanupImGui();
    void updateSampling();
    void drawBoundary();
    void drawObstacles();
    void drawVelocityField();

    const cfd::Geometry& geom_;
    const cfd::Field2D<double>& u_;
    const cfd::Field2D<double>& v_;
    const cfd::Field2D<cfd::CellTag>& tag_;

    std::vector<int> xs_;
    std::vector<int> ys_;
    int maxSamplesY_{20};
    float arrowScale_{0.05f};

    // Фиксированные факторы размера наконечника (относительно arrowScale_)
    float headLengthFactor_{0.05f};  // длина наконечника
    float headWidthFactor_{0.05f};   // ширина наконечника



    GLFWwindow* window_{nullptr};
};