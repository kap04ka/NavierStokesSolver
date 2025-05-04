#include "visualizer/FlowVisualizer.hpp"
// Заголовки OpenGL/GLFW уже включены в .hpp
#include <cmath>
#include <algorithm>
#include <iostream>
#include <limits>

FlowVisualizer::FlowVisualizer(
    const cfd::Geometry& geom,
    const cfd::Field2D<double>& u,
    const cfd::Field2D<double>& v,
    const cfd::Field2D<double>& p, // Принимаем давление
    const cfd::Field2D<cfd::CellTag>& tag)
    : geom_(geom), u_(u), v_(v), p_(p), tag_(tag) {} // Сохраняем ссылку

FlowVisualizer::~FlowVisualizer() {
    // Если бы были ресурсы, созданные в классе (кроме окна),
    // их нужно было бы очистить здесь.
}

bool FlowVisualizer::init(int width, int height, const char* title) {
    // --- Создание окна GLFW ---
    // Подсказки для нужной версии OpenGL (если требуются)
    // glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    // glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    // glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE); // Core profile

    window_ = glfwCreateWindow(width, height, title, nullptr, nullptr);
    if (!window_) {
        std::cerr << "FlowVisualizer Error: Failed to create GLFW window." << std::endl;
        return false;
    }
    glfwMakeContextCurrent(window_);
    glfwSwapInterval(1); // V-Sync

    // --- Инициализация GLEW ---
    // GLEW должен быть инициализирован ПОСЛЕ создания контекста OpenGL
    if (glewInit() != GLEW_OK) {
        std::cerr << "FlowVisualizer Error: Failed to initialize GLEW." << std::endl;
        glfwDestroyWindow(window_);
        window_ = nullptr;
        return false;
     }
    std::cout << "GLEW Initialized: Using OpenGL version " << glGetString(GL_VERSION) << std::endl;

    // --- Настройка OpenGL ---
    setupProjection(width, height);

    // --- Первоначальное обновление сэмплирования ---
    updateSampling();
    return true;
}

// Проверка окна остается без изменений
bool FlowVisualizer::shouldRun() const {
    return window_ && !glfwWindowShouldClose(window_);
}

// Основной метод отрисовки ТОЛЬКО сцены OpenGL
void FlowVisualizer::renderOpenGLScene() {
    if (!window_) return; // Защита

    // Обновляем точки сэмплирования, если параметры изменились извне
    updateSampling();

    // Рассчитываем диапазон давления, если он еще не посчитан или сброшен
    //calculatePressureRangeIfNeeded();
    if(showPressure_) {
        calculatePressureRange();
    }

    // Рисуем компоненты сцены
    drawBoundary();
    drawObstacles();

    if (showPressure_) {
        drawPressureField(); // << Вызываем отрисовку давления
    }
    if (showVelocity_) {
        drawVelocityField(); // << Вызываем отрисовку скорости
    }
    // ВАЖНО: Здесь НЕТ ImGui::Render, glClear, glfwSwapBuffers
}

// Завершение работы
void FlowVisualizer::shutdown() {
    if (window_) {
        glfwDestroyWindow(window_); // Уничтожаем окно
        window_ = nullptr;
    }
    // Завершаем GLFW, т.к. это конец работы визуализатора
    glfwTerminate();
}

// Настройка проекции OpenGL
void FlowVisualizer::setupProjection(int /*width*/, int /*height*/) {
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    double Lx = geom_.mesh().dx() * geom_.mesh().nx();
    double Ly = geom_.mesh().dy() * geom_.mesh().ny();
    glOrtho(0.0, Lx, 0.0, Ly, -1.0, 1.0); // Простая ортографическая проекция
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

// Обновление точек сэмплирования
void FlowVisualizer::updateSampling() {
    // Код остается тем же, использует публичный член maxSamplesY_
    int ny = static_cast<int>(geom_.mesh().ny());
    int nx = static_cast<int>(geom_.mesh().nx());
    if (ny <= 0 || nx <= 0) { // Защита
        xs_.clear();
        ys_.clear();
        return;
    }

    int stepY = std::max(1, ny / maxSamplesY_);
    int cntY = (stepY > 0) ? (ny + stepY - 1) / stepY : ny; // Избегаем деления на 0
    int cntX = (ny > 0 && cntY > 0) ? static_cast<int>(std::round((double)cntY * nx / ny)) : nx;
    int stepX = (cntX > 0) ? std::max(1, nx / cntX) : 1; // Избегаем деления на 0

    ys_.clear(); xs_.clear();

    // Добавляем граничные точки и точки внутри с шагом
    ys_.push_back(0);
    for (int j = stepY / 2; j < ny; j += stepY)
        if(j > 0) ys_.push_back(j); // Не добавляем 0 дважды
    if (ny - 1 > 0) ys_.push_back(ny - 1);

    xs_.push_back(0);
    for (int i = stepX / 2; i < nx; i += stepX)
        if (i > 0) xs_.push_back(i); // Не добавляем 0 дважды
    if (nx - 1 > 0) xs_.push_back(nx - 1);
}

// Отрисовка границ расчетной области
void FlowVisualizer::drawBoundary() {
    const std::size_t nx = geom_.mesh().nx();
    const std::size_t ny = geom_.mesh().ny();
    if (ny == 0) return; // Нечего рисовать, если нет высоты

    const double dx = geom_.mesh().dx();
    const double dy = geom_.mesh().dy();

    // Задаем цвет для стенок (можно такой же, как у препятствий, или другой)
    glColor3f(0.5f, 0.5f, 0.5f); // Серый цвет, как у препятствий

    glBegin(GL_QUADS); // Начинаем рисовать квадраты

    // --- Рисуем НИЖНЮЮ стенку (все ячейки в строке j=0) ---
    std::size_t j_bottom = 0;
    double y_bottom = j_bottom * dy; // y координата нижнего края ячеек j=0
    for (std::size_t i = 0; i < nx; ++i) {
        double x = i * dx;
        // Вершины квадрата для ячейки (i, 0)
        glVertex2d(x,      y_bottom);      // Нижний левый угол
        glVertex2d(x + dx, y_bottom);      // Нижний правый угол
        glVertex2d(x + dx, y_bottom + dy); // Верхний правый угол
        glVertex2d(x,      y_bottom + dy); // Верхний левый угол
    }

    // --- Рисуем ВЕРХНЮЮ стенку (все ячейки в строке j=ny-1) ---
    std::size_t j_top = ny - 1;
    double y_top = j_top * dy; // y координата нижнего края ячеек j=ny-1
    for (std::size_t i = 0; i < nx; ++i) {
        double x = i * dx;
        // Вершины квадрата для ячейки (i, ny-1)
        glVertex2d(x,      y_top);      // Нижний левый угол
        glVertex2d(x + dx, y_top);      // Нижний правый угол
        glVertex2d(x + dx, y_top + dy); // Верхний правый угол (y = Ly)
        glVertex2d(x,      y_top + dy); // Верхний левый угол (y = Ly)
    }

    glEnd(); // Завершаем рисование квадратов
}

// Отрисовка препятствий
void FlowVisualizer::drawObstacles() {
    glColor3f(0.4f, 0.4f, 0.4f); // Цвет препятствий
    glBegin(GL_QUADS);
    for (std::size_t j = 0; j < geom_.mesh().ny(); ++j) {
        for (std::size_t i = 0; i < geom_.mesh().nx(); ++i) {
            if (tag_(i, j) == cfd::CellTag::SOLID) {
                double x = i * geom_.mesh().dx();
                double y = j * geom_.mesh().dy();
                double dx = geom_.mesh().dx();
                double dy = geom_.mesh().dy();
                glVertex2d(x, y);
                glVertex2d(x + dx, y);
                glVertex2d(x + dx, y + dy);
                glVertex2d(x, y + dy);
            }
        }
    }
    glEnd();
}

// Отрисовка поля скорости (вектора)
void FlowVisualizer::drawVelocityField() {
    // Определяем максимальную скорость для нормализации (опционально)
    double maxVelMag = 0.0;
    for (std::size_t j = 0; j < geom_.mesh().ny(); ++j) {
        for (std::size_t i = 0; i < geom_.mesh().nx(); ++i) {
             if (tag_(i,j) == cfd::CellTag::FLUID) {
                 maxVelMag = std::max(maxVelMag, std::hypot(u_(i,j), v_(i,j)));
             }
        }
    }
    if (maxVelMag < 1e-9) maxVelMag = 1.0; // Избегаем деления на ноль

    if (showPressure_) {
        glColor3f(0.0f, 0.0f, 0.0f); // Черный цвет, если давление показывается
    } else {
        glColor3f(0.2f, 0.9f, 0.2f); // Зеленый цвет, если давление скрыто
    }

    for (int j : ys_) {
        for (int i : xs_) {
            if (tag_(i, j) == cfd::CellTag::SOLID) continue;

            auto [cx, cy] = geom_.mesh().centre(i, j);
            double u = u_(i, j);
            double v = v_(i, j);
            double mag = std::hypot(u, v);

            if (mag < 1e-5) continue; // Не рисуем слишком маленькие вектора

            double scale = arrowScale_ * (mag / maxVelMag); // Масштабируем по относительной величине
            // Альтернатива: использовать фиксированный масштаб arrowScale_
            // double scale = arrowScale_;

            double dx_n = u / mag; // Нормализованное направление
            double dy_n = v / mag;
            double end_x = cx + dx_n * scale;
            double end_y = cy + dy_n * scale;

            // Стержень стрелки
            glBegin(GL_LINES);
                glVertex2d(cx, cy);
                glVertex2d(end_x, end_y);
            glEnd();

            // Наконечник
            double headLen  = arrowScale_ * headLengthFactor_;
            double headWidth = arrowScale_ * headWidthFactor_;
            double px = -dy_n;
            double py = dx_n;
            double bx = end_x - dx_n * headLen;
            double by = end_y - dy_n * headLen;
            double lx = bx + px * (headWidth / 2.0); // Сдвиг влево на ПОЛОВИНУ headWidth
            double ly = by + py * (headWidth / 2.0);
            double rx = bx - px * (headWidth / 2.0); // Сдвиг вправо на ПОЛОВИНУ headWidth
            double ry = by - py * (headWidth / 2.0);

            glBegin(GL_TRIANGLES);
                glVertex2d(end_x,  end_y);
                glVertex2d(lx,  ly);
                glVertex2d(rx,  ry);
            glEnd();
        }
    }
}


// ОТДЕЛЬНЫЙ Метод отрисовки поля давления (цветовая карта)
void FlowVisualizer::drawPressureField() {
    // Диапазон давления уже должен быть рассчитан в calculatePressureRangeIfNeeded()

    if (std::abs(maxPressure_ - minPressure_) < 1e-9) return; // Не рисуем, если диапазон нулевой

    double pressureDiffInv = 1.0 / (maxPressure_ - minPressure_); // Обратный диапазон для нормализации

    glBegin(GL_QUADS);
    for (std::size_t j = 1; j < geom_.mesh().ny() - 1; ++j) {
        for (std::size_t i = 0; i < geom_.mesh().nx(); ++i) {
            if (tag_(i, j) == cfd::CellTag::SOLID) continue; // Пропускаем твердые тела

            double pressure = p_(i, j);
            // Нормализуем давление к диапазону [0, 1]
            float t = static_cast<float>((pressure - minPressure_) * pressureDiffInv);
            t = std::max(0.0f, std::min(1.0f, t)); // Ограничиваем [0, 1]

            float r, g, b;
            getColorMap(t, r, g, b); // Получаем цвет из карты
            glColor3f(r, g, b);

            // Рисуем квадрат ячейки
            double x = i * geom_.mesh().dx();
            double y = j * geom_.mesh().dy();
            double dx = geom_.mesh().dx();
            double dy = geom_.mesh().dy();

            glVertex2d(x, y);
            glVertex2d(x + dx, y);
            glVertex2d(x + dx, y + dy);
            glVertex2d(x, y + dy);
        }
    }
    glEnd();
}

// Расчет диапазона давления (вызывается перед отрисовкой, если нужно)
//void FlowVisualizer::calculatePressureRangeIfNeeded() {
void FlowVisualizer::calculatePressureRange() {
    //if (pressureRangeCalculated_) return;

    minPressure_ = std::numeric_limits<double>::max();
    maxPressure_ = std::numeric_limits<double>::lowest();
    bool found_fluid = false;

    for (std::size_t j = 0; j < geom_.mesh().ny(); ++j) {
        for (std::size_t i = 0; i < geom_.mesh().nx(); ++i) {
             if (tag_(i,j) == cfd::CellTag::FLUID) {
                minPressure_ = std::min(minPressure_, p_(i,j));
                maxPressure_ = std::max(maxPressure_, p_(i,j));
                found_fluid = true;
             }
        }
    }

    if (!found_fluid) { // Если нет жидкости, ставим диапазон по умолчанию
         minPressure_ = 0.0;
         maxPressure_ = 1.0;
    } else if (std::abs(maxPressure_ - minPressure_) < 1e-9) {
        // Если диапазон очень мал, немного расширяем его
        maxPressure_ += 0.5;
        minPressure_ -= 0.5;
        // Убедимся, что min не стал > max
        if (minPressure_ > maxPressure_) std::swap(minPressure_, maxPressure_);
    }
    //pressureRangeCalculated_ = true; // Устанавливаем флаг
    // std::cout << "Pressure range calculated: [" << minPressure_ << ", " << maxPressure_ << "]" << std::endl; // Отладка
}


// Простая функция для получения цвета из карты (Синий -> Зеленый -> Красный)
// (Реализация как была в предыдущем ответе)
void FlowVisualizer::getColorMap(float value, float& r, float& g, float& b) {
    // value находится в диапазоне [0, 1]
    if (value < 0.0f) value = 0.0f; // Доп. защита
    if (value > 1.0f) value = 1.0f;

    if (value < 0.5f) {
        // От синего (0,0,1) к зеленому (0,1,0)
        r = 0.0f;
        g = 2.0f * value;
        b = 1.0f - 2.0f * value;
    } else {
        // От зеленого (0,1,0) к красному (1,0,0)
        r = 2.0f * (value - 0.5f);
        g = 1.0f - 2.0f * (value - 0.5f);
        b = 0.0f;
    }
}