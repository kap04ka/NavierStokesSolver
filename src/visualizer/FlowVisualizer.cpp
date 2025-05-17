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
    const cfd::Field2D<double>& p_or_psi, // Принимаем давление
    const cfd::Field2D<cfd::CellTag>& tag)
    : geom_(geom), u_(u), v_(v), scalar_field_(p_or_psi), tag_(tag) {} // Сохраняем ссылку

FlowVisualizer::~FlowVisualizer() {}

bool FlowVisualizer::init(int width, int height, const char* title) {
    // --- Создание окна GLFW ---

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

    updateSampling();

    // Рассчитываем диапазон давления, если он еще не посчитан или сброшен
    if(showPressure_ || showStreamlines_) {
        calculateScalarFieldRange();
    }

    // Рисуем компоненты сцены
    drawBoundary();
    drawObstacles();

    if (showPressure_) {
        drawPressureField(); // << Вызываем отрисовку давления
    }
    if (showStreamlines_) {
        drawStreamlines();  // << Вызываем отрисовку изолиний функции тока
        drawBoundary();
    }
    if (showVelocity_) {
        drawVelocityField(); // << Вызываем отрисовку скорости
    }

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
    } else if (showStreamlines_) {
        glColor3f(0.9f, 0.1f, 0.1f);
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

    if (std::abs(maxScalarValue_ - minScalarValue_) < 1e-9) return; // Не рисуем, если диапазон нулевой

    double scalarDiffInv = 1.0 / (maxScalarValue_ - minScalarValue_); // Обратный диапазон для нормализации

    glBegin(GL_QUADS);
    for (std::size_t j = 1; j < geom_.mesh().ny() - 1; ++j) {
        for (std::size_t i = 0; i < geom_.mesh().nx(); ++i) {
            if (tag_(i, j) == cfd::CellTag::SOLID) continue; // Пропускаем твердые тела

            double pressure = scalar_field_(i, j);
            // Нормализуем давление к диапазону [0, 1]
            float t = static_cast<float>((pressure - minScalarValue_) * scalarDiffInv);
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
void FlowVisualizer::calculateScalarFieldRange() {
    minScalarValue_ = std::numeric_limits<double>::max();
    maxScalarValue_ = std::numeric_limits<double>::lowest();
    bool found_fluid = false;

    for (std::size_t j = 0; j < geom_.mesh().ny(); ++j) {
        for (std::size_t i = 0; i < geom_.mesh().nx(); ++i) {
             if (tag_(i,j) == cfd::CellTag::FLUID) {
                minScalarValue_ = std::min(minScalarValue_, scalar_field_(i,j));
                maxScalarValue_ = std::max(maxScalarValue_, scalar_field_(i,j));
                found_fluid = true;
             }
        }
    }

    if (!found_fluid) { // Если нет жидкости, ставим диапазон по умолчанию
        minScalarValue_ = 0.0;
        maxScalarValue_ = 1.0;
    } else if (std::abs(maxScalarValue_ - minScalarValue_) < 1e-9) {
        // Если диапазон очень мал, немного расширяем его
        maxScalarValue_ += 0.5;
        maxScalarValue_ -= 0.5;
        // Убедимся, что min не стал > max
        if (minScalarValue_ > maxScalarValue_) std::swap(minScalarValue_, maxScalarValue_);
    }
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

double FlowVisualizer::get_value_at_node(std::size_t i_node, std::size_t j_node) const {
    const std::size_t nx_c = geom_.mesh().nx(); // Количество ячеек по X
    const std::size_t ny_c = geom_.mesh().ny(); // Количество ячеек по Y

    // Если nx_c или ny_c равны 0, это ошибка конфигурации, но нужно обработать
    if (nx_c == 0 || ny_c == 0) {
        // std::cerr << "Warning: get_value_at_node called on an empty mesh." << std::endl;
        return 0.0;
    }

    // --- Обработка узлов, лежащих на внешних границах расчетной области ---
    // Эти узлы имеют одну или две координаты, совпадающие с границами сетки узлов (0, nx_c, 0, ny_c)

    bool on_left_boundary = (i_node == 0);
    bool on_right_boundary = (i_node == nx_c);
    bool on_bottom_boundary = (j_node == 0);
    bool on_top_boundary = (j_node == ny_c);

    // 1. Угловые узлы: берем значение из ближайшего центра ячейки
    if (on_left_boundary && on_bottom_boundary) return scalar_field_(0, 0);
    if (on_right_boundary && on_bottom_boundary) return scalar_field_(nx_c - 1, 0);
    if (on_left_boundary && on_top_boundary) return scalar_field_(0, ny_c - 1);
    if (on_right_boundary && on_top_boundary) return scalar_field_(nx_c - 1, ny_c - 1);

    // 2. Узлы на ребрах (не угловые)
    if (on_bottom_boundary) { // Нижняя граница (j_node = 0, 0 < i_node < nx_c)
        // Узел (i_node, 0) лежит между центрами ячеек (i_node-1, 0) и (i_node, 0)
        return 0.5 * (scalar_field_(i_node - 1, 0) + scalar_field_(i_node, 0));
    }
    if (on_top_boundary) { // Верхняя граница (j_node = ny_c, 0 < i_node < nx_c)
        // Узел (i_node, ny_c) лежит между центрами ячеек (i_node-1, ny_c-1) и (i_node, ny_c-1)
        return 0.5 * (scalar_field_(i_node - 1, ny_c - 1) + scalar_field_(i_node, ny_c - 1));
    }
    if (on_left_boundary) { // Левая граница (i_node = 0, 0 < j_node < ny_c)
        // Узел (0, j_node) лежит между центрами ячеек (0, j_node-1) и (0, j_node)
        return 0.5 * (scalar_field_(0, j_node - 1) + scalar_field_(0, j_node));
    }
    if (on_right_boundary) { // Правая граница (i_node = nx_c, 0 < j_node < ny_c)
        // Узел (nx_c, j_node) лежит между центрами ячеек (nx_c-1, j_node-1) и (nx_c-1, j_node)
        return 0.5 * (scalar_field_(nx_c - 1, j_node - 1) + scalar_field_(nx_c - 1, j_node));
    }

    // 3. Внутренний узел: Узел (i_node, j_node) окружен центрами 4-х ячеек:
    // Ячейка BL: (i_node - 1, j_node - 1)
    // Ячейка BR: (i_node,     j_node - 1)
    // Ячейка TL: (i_node - 1, j_node)
    // Ячейка TR: (i_node,     j_node)
    // (Это при условии, что i_node и j_node > 0 и < nx_cells/ny_cells соответственно)
    // Такая нумерация узлов соответствует тому, что узел (i,j) - это верхний правый угол ячейки поля (i-1,j-1)
    // Или, если узел (i_node, j_node) - это левый нижний угол ячейки поля (i_node, j_node)
    // и он окружен ячейками (i_node-1,j_node-1), (i_node,j_node-1), (i_node-1,j_node), (i_node,j_node)
    // то i_node должно быть от 1 до nx_c-1, j_node от 1 до ny_c-1 для этого усреднения.

    // Перепроверим: если i_node от 0 до nx_c, j_node от 0 до ny_c
    // Узел (i_node, j_node) является общим углом для ячеек:
    // (i_node-1, j_node-1) - если i_node>0, j_node>0
    // (i_node,   j_node-1) - если i_node<nx_c, j_node>0
    // (i_node-1, j_node)   - если i_node>0, j_node<ny_c
    // (i_node,   j_node)   - если i_node<nx_c, j_node<ny_c

    // Эта проверка гарантирует, что все 4 индекса ячеек поля существуют:
    if (i_node > 0 && i_node < nx_c && j_node > 0 && j_node < ny_c) {
         return 0.25 * (scalar_field_(i_node - 1, j_node - 1) + 
                         scalar_field_(i_node,     j_node - 1) +
                         scalar_field_(i_node - 1, j_node)     + 
                         scalar_field_(i_node,     j_node));
    }
    
    // Если мы дошли сюда, значит, узел лежит на границе, но не был покрыт предыдущими случаями,
    // или сетка очень маленькая (1xN, Nx1), и часть предыдущих условий для ребер не сработала.
    // Этого быть не должно при nx_c, ny_c >= 1.
    // Для безопасности, вернем значение из ближайшей допустимой ячейки.
    std::size_t i_clamped = std::max(static_cast<std::size_t>(0), std::min(i_node, nx_c - 1));
     if (i_node == nx_c && nx_c > 0) i_clamped = nx_c - 1; // Если узел на правой границе, берем из последней ячейки
     else if (i_node > 0) i_clamped = i_node -1; // Если узел не на левой границе, но не покрыт выше

    std::size_t j_clamped = std::max(static_cast<std::size_t>(0), std::min(j_node, ny_c - 1));
    if (j_node == ny_c && ny_c > 0) j_clamped = ny_c - 1;
    else if (j_node > 0) j_clamped = j_node -1;

    std::cerr << "Warning: get_value_at_node reached fallback for node (" << i_node << "," << j_node 
              << "). Using value from cell (" << i_clamped << "," << j_clamped << ")." << std::endl;
    return scalar_field_(i_clamped, j_clamped); 
}

std::pair<double, double> FlowVisualizer::interpolate_isoline_point(
    double x1, double y1, double val1,
    double x2, double y2, double val2,
    double level_val) const
{
    // Проверка на случай, если val1 и val2 слишком близки (избегаем деления на ноль)
    // или если level_val вне диапазона (хотя это должно отсекаться выбором ребра)
    if (std::abs(val1 - val2) < 1e-9) { 
        return { x1, y1 }; // Или середина, или одна из точек
    }
    // Проверка, что level_val между val1 и val2
    if (!((val1 <= level_val && level_val <= val2) || (val2 <= level_val && level_val <= val1))) {
        // Это не должно происходить, если ребро действительно пересекается изолинией
        // Возвращаем одну из точек как заглушку
        return {x1,y1}; 
    }

    double t = (level_val - val1) / (val2 - val1);
    return { x1 + t * (x2 - x1), y1 + t * (y2 - y1) };
}

void FlowVisualizer::drawStreamlines() {
    if (numStreamlineLevels_ <= 0) return;
    if (std::abs(maxScalarValue_ - minScalarValue_) < 1e-9) return; // Диапазон нулевой

    const std::size_t nx_cells = geom_.mesh().nx();
    const std::size_t ny_cells = geom_.mesh().ny();
    const double dx = geom_.mesh().dx();
    const double dy = geom_.mesh().dy();

    // Устанавливаем цвет и толщину линий
    glColor3f(streamlineColorR_, streamlineColorG_, streamlineColorB_);
    glLineWidth(streamlineThickness_);

    // Итерация по "ячейкам Marching Squares", которые образованы узлами сетки.
    // Узлов по x: nx_cells+1, по y: ny_cells+1.
    // Таких ячеек будет nx_cells * ny_cells.
    // Ячейка Marching Squares с левым нижним узлом (i_node, j_node)
    for (std::size_t j_node = 0; j_node < ny_cells; ++j_node) { // от 0 до ny_cells-1
        for (std::size_t i_node = 0; i_node < nx_cells; ++i_node) { // от 0 до nx_cells-1
            
            // Пропускаем рисование, если ячейка (или ее центр) находится внутри SOLID препятствия
            // Это грубая проверка, можно улучшить, проверяя все 4 ячейки поля, образующие этот квадрат.
            // Центр квадрата Marching Squares: ( (i_node+0.5)*dx, (j_node+0.5)*dy )
            // Соответствует центру ячейки поля (i_node, j_node)
            if (i_node < nx_cells && j_node < ny_cells) { // Проверяем, что ячейка поля существует
                if (tag_(i_node, j_node) == cfd::CellTag::SOLID) {
                    continue; // Пропускаем рисование изолиний для SOLID ячеек поля
                }
            } else {
                continue; // Узел на самой правой или верхней границе сетки узлов, нет соответствующей ячейки поля справа/сверху
            }

            // Значения psi в 4-х узлах текущего квадрата Marching Squares
            double node_val[4];
            node_val[0] = get_value_at_node(i_node,     j_node);     // Bottom-left (BL)
            node_val[1] = get_value_at_node(i_node + 1, j_node);     // Bottom-right (BR)
            node_val[2] = get_value_at_node(i_node + 1, j_node + 1); // Top-right (TR)
            node_val[3] = get_value_at_node(i_node,     j_node + 1); // Top-left (TL)

            // Физические координаты узлов
            double node_x[] = { i_node * dx, (i_node + 1) * dx, (i_node + 1) * dx, i_node * dx };
            double node_y[] = { j_node * dy, j_node * dy, (j_node + 1) * dy, (j_node + 1) * dy };

            for (int k_level = 0; k_level < numStreamlineLevels_; ++k_level) {
                double psi_level = minScalarValue_ + 
                                   static_cast<double>(k_level + 0.5) * // Берем середину интервала для уровня
                                   (maxScalarValue_ - minScalarValue_) / static_cast<double>(numStreamlineLevels_);
                if (numStreamlineLevels_ == 1 && k_level == 0) { // Для одной линии берем середину диапазона
                     psi_level = (minScalarValue_ + maxScalarValue_) / 2.0;
                }


                int square_index = 0;
                if (node_val[0] > psi_level) square_index |= 1; // BL
                if (node_val[1] > psi_level) square_index |= 2; // BR
                if (node_val[2] > psi_level) square_index |= 4; // TR
                if (node_val[3] > psi_level) square_index |= 8; // TL

                // Таблица ребер для Marching Squares (16 случаев)
                // Ребра: 0: BL-BR (низ), 1: BR-TR (право), 2: TR-TL (верх), 3: TL-BL (лево)
                // Каждая запись: {v1_idx, v2_idx, v3_idx, v4_idx, -1} (пары индексов ребер)
                // Или {p1_edge_idx, p2_edge_idx, -1} для одного сегмента
                //    {p1_edge_idx, p2_edge_idx, p3_edge_idx, p4_edge_idx, -1} для двух сегментов (случаи 5 и 10)

                // Более простая таблица, указывающая, какие ребра пересекаются
                // (узел0-узел1), (узел1-узел2), (узел2-узел3), (узел3-узел0)
                //    Ребро 0: между узлом 0 и 1
                //    Ребро 1: между узлом 1 и 2
                //    Ребро 2: между узлом 2 и 3
                //    Ребро 3: между узлом 3 и 0
                
                std::vector<std::pair<double,double>> intersections;

                // Ребро 0 (низ: узел 0 -> узел 1)
                if (((square_index & 1) && !(square_index & 2)) || (!(square_index & 1) && (square_index & 2)))
                    intersections.push_back(interpolate_isoline_point(node_x[0], node_y[0], node_val[0], node_x[1], node_y[1], node_val[1], psi_level));
                // Ребро 1 (право: узел 1 -> узел 2)
                if (((square_index & 2) && !(square_index & 4)) || (!(square_index & 2) && (square_index & 4)))
                    intersections.push_back(interpolate_isoline_point(node_x[1], node_y[1], node_val[1], node_x[2], node_y[2], node_val[2], psi_level));
                // Ребро 2 (верх: узел 2 -> узел 3)
                if (((square_index & 4) && !(square_index & 8)) || (!(square_index & 4) && (square_index & 8)))
                    intersections.push_back(interpolate_isoline_point(node_x[2], node_y[2], node_val[2], node_x[3], node_y[3], node_val[3], psi_level));
                // Ребро 3 (лево: узел 3 -> узел 0)
                if (((square_index & 8) && !(square_index & 1)) || (!(square_index & 8) && (square_index & 1)))
                    intersections.push_back(interpolate_isoline_point(node_x[3], node_y[3], node_val[3], node_x[0], node_y[0], node_val[0], psi_level));

                // Рисуем отрезки
                if (intersections.size() == 2) {
                    glBegin(GL_LINES);
                    glVertex2d(intersections[0].first, intersections[0].second);
                    glVertex2d(intersections[1].first, intersections[1].second);
                    glEnd();
                } else if (intersections.size() == 4) {
                    glBegin(GL_LINES);
                    glVertex2d(intersections[0].first, intersections[0].second);
                    glVertex2d(intersections[1].first, intersections[1].second); // Может быть не тот сосед
                    glEnd();
                    glBegin(GL_LINES);
                    glVertex2d(intersections[2].first, intersections[2].second);
                    glVertex2d(intersections[3].first, intersections[3].second); // Может быть не тот сосед
                    glEnd();
                }
            } // конец цикла по уровням psi
        } // конец цикла по i_node
    } // конец цикла по j_node
}