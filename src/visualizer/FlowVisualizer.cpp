#include "visualizer/FlowVisualizer.hpp"
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <imgui.h>
#include <backends/imgui_impl_glfw.h>
#include <backends/imgui_impl_opengl3.h>
#include <cmath>
#include <algorithm>

FlowVisualizer::FlowVisualizer(
    const cfd::Geometry& geom,
    const cfd::Field2D<double>& u,
    const cfd::Field2D<double>& v,
    const cfd::Field2D<cfd::CellTag>& tag)
    : geom_(geom), u_(u), v_(v), tag_(tag) {}

FlowVisualizer::~FlowVisualizer() = default;

bool FlowVisualizer::init(int width, int height, const char* title) {
    if (!glfwInit()) return false;
    window_ = glfwCreateWindow(width, height, title, nullptr, nullptr);
    if (!window_) { glfwTerminate(); return false; }
    glfwMakeContextCurrent(window_);
    if (glewInit() != GLEW_OK) return false;

    setupProjection();
    setupImGui();
    updateSampling();
    return true;
}

bool FlowVisualizer::shouldRun() const {
    return window_ && !glfwWindowShouldClose(window_);
}

void FlowVisualizer::render() {
    glfwPollEvents();
    int w, h;
    glfwGetFramebufferSize(window_, &w, &h);
    glViewport(0, 0, w, h);
    glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);

    drawBoundary();
    drawObstacles();
    drawVelocityField();

    ImGui_ImplOpenGL3_NewFrame();
    ImGui_ImplGlfw_NewFrame();
    ImGui::NewFrame();
    ImGui::Begin("Settings");
    ImGui::SliderInt("Max Samples Y", &maxSamplesY_, 5, 50);
    ImGui::SliderFloat("Arrow Scale", &arrowScale_, 0.001f, 0.1f);
    ImGui::SliderFloat("Head Length Factor", &headLengthFactor_, 0.001f, 0.15f);
    ImGui::SliderFloat("Head Width Factor",  &headWidthFactor_,  0.001f, 0.15f);
    ImGui::End();
    ImGui::Render();
    ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
    glfwSwapBuffers(window_);
    updateSampling();
}

void FlowVisualizer::shutdown() {
    cleanupImGui();
    if (window_) glfwDestroyWindow(window_);
    glfwTerminate();
}

void FlowVisualizer::setupProjection() {
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    double Lx = geom_.mesh().dx() * geom_.mesh().nx();
    double Ly = geom_.mesh().dy() * geom_.mesh().ny();
    glOrtho(0.0, Lx, 0.0, Ly, -1.0, 1.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

void FlowVisualizer::setupImGui() {
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGui_ImplGlfw_InitForOpenGL(window_, true);
    ImGui_ImplOpenGL3_Init("#version 130");
}

void FlowVisualizer::cleanupImGui() {
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
}

void FlowVisualizer::updateSampling() {
    int ny = static_cast<int>(geom_.mesh().ny());
    int nx = static_cast<int>(geom_.mesh().nx());
    int stepY = std::max(1, ny / maxSamplesY_);
    int cntY = (ny + stepY - 1) / stepY;
    int cntX = static_cast<int>(std::round((double)cntY * nx / ny));
    int stepX = std::max(1, nx / cntX);
    ys_.clear(); xs_.clear();

    ys_.push_back(0);
    for (int j = stepY/2; j < ny; j += stepY) 
        ys_.push_back(j);
    ys_.push_back(ny - 1);

    xs_.push_back(0);
    for (int i = stepX/2; i < nx; i += stepX) 
        xs_.push_back(i);
    xs_.push_back(nx-1);
}

void FlowVisualizer::drawBoundary() {
    double Lx = geom_.mesh().dx() * geom_.mesh().nx();
    double Ly = geom_.mesh().dy() * geom_.mesh().ny();
    glLineWidth(4.0f);

    glColor3f(1.0f, 1.0f, 1.0f);
    glBegin(GL_LINE_LOOP);
        glVertex2d(0.0, 0.0);
        glVertex2d(Lx, 0.0);
    glEnd();

    glLineWidth(4.0f);
    glColor3f(1.0f, 1.0f, 1.0f);
    glBegin(GL_LINE_LOOP);
        glVertex2d(0.0, Ly);
        glVertex2d(Lx, Ly);
    glEnd();
    glLineWidth(1.0f);
}

void FlowVisualizer::drawObstacles() {
    glColor3f(0.5f, 0.5f, 0.5f);
    glBegin(GL_QUADS);
    for (int j = 0; j < geom_.mesh().ny(); ++j) {
        for (int i = 0; i < geom_.mesh().nx(); ++i) {
            if (tag_(i, j) == cfd::CellTag::SOLID) {
                double x = i * geom_.mesh().dx();
                double y = j * geom_.mesh().dy();
                glVertex2d(x, y);
                glVertex2d(x + geom_.mesh().dx(), y);
                glVertex2d(x + geom_.mesh().dx(), y + geom_.mesh().dy());
                glVertex2d(x, y + geom_.mesh().dy());
            }
        }
    }
    glEnd();
}

void FlowVisualizer::drawVelocityField() {
    // Определяем наибольшую скорость для нормировки стрелок
    double uMax = 0.0;
    for (std::size_t j = 0; j < geom_.mesh().ny(); ++j) {
        for (std::size_t i = 0; i < geom_.mesh().nx(); ++i) {
            if (tag_(i,j) == cfd::CellTag::SOLID) continue;
            double uz = u_(i,j);
            double vz = v_(i,j);
            uMax = std::max(uMax, std::hypot(uz, vz));
        }
    }
    if (uMax < 1e-8) uMax = 1.0;

    // Отрисовка стрелок скорости с учетом относительной длины
    glColor3f(0.0f, 0.8f, 0.2f);
    for (int j : ys_) {
        for (int i : xs_) {
            if (tag_(i, j) == cfd::CellTag::SOLID) continue;
            auto [cx, cy] = geom_.mesh().centre(i, j);
            double uz = u_(i, j);
            double vz = v_(i, j);
            double mag = std::hypot(uz, vz);
            if (mag < 1e-6) continue;
            // Направление вектора
            double dx = uz / mag;
            double dy = vz / mag;
            // длина стрелки пропорциональна mag/uMax
            double len = arrowScale_ * (mag / uMax);
            double ex = cx + dx * len;
            double ey = cy + dy * len;

            // Рисуем стержень стрелки
            glBegin(GL_LINES);
                glVertex2d(cx, cy);
                glVertex2d(ex, ey);
            glEnd();

            // Наконечник
            double headLen  = arrowScale_ * headLengthFactor_;
            double headWidth = arrowScale_ * headWidthFactor_;
            double px = -dy;
            double py = dx;
            double bx = ex - dx * headLen;
            double by = ey - dy * headLen;
            double lx = bx + px * headWidth;
            double ly = by + py * headWidth;
            double rx = bx - px * headWidth;
            double ry = by - py * headWidth;

            glBegin(GL_TRIANGLES);
                glVertex2d(ex,  ey);
                glVertex2d(lx,  ly);
                glVertex2d(rx,  ry);
            glEnd();
        }
    }
}