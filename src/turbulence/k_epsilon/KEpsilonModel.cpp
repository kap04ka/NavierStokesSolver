#include "turbulence/k_epsilon/KEpsilonModel.hpp"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <iostream> 

namespace cfd {

    KEpsilonModel::KEpsilonModel(const Geometry& geom, double rho, double nu_molecular,
        double u_max_for_init,
        double inlet_turb_intensity, // Принимаем параметры
        double inlet_length_scale_factor, // Принимаем фактор для L
        const KEpsilonConstants& constants)
        : TurbulenceModel(geom, rho, nu_molecular),
        constants_(constants),
        inlet_turbulence_intensity_(inlet_turb_intensity),
        inlet_length_scale_(inlet_length_scale_factor * geom.mesh().Ly())
{
    std::cout << "Initializing Standard k-epsilon model with Wall Functions..." << std::endl;
    const auto nx = geom_.mesh().nx();
    const auto ny = geom_.mesh().ny();
    
    // Используем константы из структуры для инициализации
    k_           = std::make_unique<Field2D<double>>(nx, ny, constants_.k_min);
    epsilon_     = std::make_unique<Field2D<double>>(nx, ny, constants_.epsilon_min);
    nu_t_        = std::make_unique<Field2D<double>>(nx, ny, 0.0);
    Pk_          = std::make_unique<Field2D<double>>(nx, ny, 0.0);
    k_old_       = std::make_unique<Field2D<double>>(nx, ny, 0.0);
    epsilon_old_ = std::make_unique<Field2D<double>>(nx, ny, 0.0);
    tau_wx_      = std::make_unique<Field2D<double>>(nx, ny, 0.0);
    tau_wy_      = std::make_unique<Field2D<double>>(nx, ny, 0.0);

    double U_avg_guess = 0.5 * u_max_for_init; // Примерная средняя скорость для оценки (можно взять umax)
    double k0 = 1.5 * std::pow(U_avg_guess * inlet_turbulence_intensity_, 2);
    double eps0 = std::pow(constants_.Cmu, 0.75) * std::pow(k0, 1.5) / inlet_length_scale_;
    k_->fill(std::max(constants_.k_min, k0));
    epsilon_->fill(std::max(constants_.epsilon_min, eps0));
    update_nu_t(); // Обновляем nu_t

    std::cout << "k-epsilon model initialized." << std::endl;
}

// --- Метод применения всех ГУ ---
void KEpsilonModel::apply_bc(const Field2D<double>& u, const Field2D<double>& v) {
    // Сначала обнуляем поля напряжений (они будут вычислены в WF)
    tau_wx_->fill(0.0);
    tau_wy_->fill(0.0);

    // Применяем ГУ для входа/выхода
    apply_inlet_outlet_bc(u, v); // Передаем u, v на случай, если они нужны для ГУ

    // Применяем пристенные функции (они перезапишут k/epsilon у стенок и вычислят tau_w)
    apply_wall_functions(u, v);
}

// --- Реализация ГУ на входе/выходе ---
void KEpsilonModel::apply_inlet_outlet_bc(const Field2D<double>& u, const Field2D<double>& v) {
    const auto nx = geom_.mesh().nx();
    const auto ny = geom_.mesh().ny();
    if (nx < 2) return;

    // --- ВХОД (i=0) ---
    // Задаем k и epsilon на входе.

    for (std::size_t j = 0; j < ny; ++j) {
        // Используем скорость из поля u (которое было передано) как референсную U_ref
        // Если вход не вертикальный, нужно брать модуль скорости |U|
        double U_ref = std::max(1e-6, std::fabs(u(0, j))); // Берем скорость на входе, избегаем нуля

        double k_in = 1.5 * std::pow(U_ref * inlet_turbulence_intensity_, 2);

        double effective_length_scale = std::max(1e-6, inlet_length_scale_);
        double epsilon_in = std::pow(constants_.Cmu, 0.75) * std::pow(k_in, 1.5) / effective_length_scale;

        // Применяем минимальные значения
        k_->operator()(0, j) = std::max(constants_.k_min, k_in);
        epsilon_->operator()(0, j) = std::max(constants_.epsilon_min, epsilon_in);

        // nu_t на входе будет вычислено позже в update_nu_t()
    }

    // --- ВЫХОД (i=nx-1) ---
    // Условие Неймана (нулевой градиент по нормали - здесь по x)
    for (std::size_t j = 0; j < ny; ++j) {
        k_->operator()(nx - 1, j) = k_->operator()(nx - 2, j);
        epsilon_->operator()(nx - 1, j) = epsilon_->operator()(nx - 2, j);
        // nu_t на выходе также установится Нейманом после update_nu_t()
    }
}

void KEpsilonModel::apply_wall_functions(const Field2D<double>& u, const Field2D<double>& v) {
    const auto nx = geom_.mesh().nx();
    const auto ny = geom_.mesh().ny();
    const double dx = geom_.mesh().dx();
    const double dy = geom_.mesh().dy();
    const auto& tags = geom_.tags();

    // Цикл по ВСЕМ ячейкам, чтобы найти те, что РЯДОМ со стеной
    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {

            // Пропускаем не-жидкостные ячейки
            if (tags(i, j) != CellTag::FLUID) continue;

            // Проверяем соседей (Север, Юг, Запад, Восток)
            // --- Южная стенка (Граница j=0 или SOLID снизу) ---
            if (j == 0 || (j > 0 && tags(i, j - 1) == CellTag::SOLID)) {
                double yp = geom_.mesh().centre(i, j).second; // Расстояние до стенки y=0
                if (j > 0 && tags(i, j - 1) == CellTag::SOLID) { // Если SOLID снизу
                    yp = 0.5 * dy; // Расстояние до центра грани
                }
                double Up = u(i, j); // Скорость параллельная стенке (предполагаем u)
                double u_tau = calculate_u_tau(Up, yp);
                double tau_w = rho_ * u_tau * u_tau * std::copysign(1.0, Up); // Учитываем знак Up для tau_w
                double kP = std::max(constants_.k_min, u_tau * u_tau / std::sqrt(constants_.Cmu));
                double epsP = std::max(constants_.epsilon_min, std::pow(std::abs(u_tau), 3.0) / (constants_.kappa * yp));

                k_->operator()(i, j) = kP;
                epsilon_->operator()(i, j) = epsP;
                tau_wx_->operator()(i, j) = tau_w; // Сохраняем x-компоненту tau в ячейке P
                tau_wy_->operator()(i, j) = 0.0;   // Y-компонента для горизонтальной стенки = 0
            }
            // --- Северная стенка (Граница j=ny-1 или SOLID сверху) ---
            else if (j == ny - 1 || (j < ny - 1 && tags(i, j + 1) == CellTag::SOLID)) {
                double yp = geom_.mesh().Ly() - geom_.mesh().centre(i, j).second; // Расстояние до стенки y=Ly
                if (j < ny - 1 && tags(i, j + 1) == CellTag::SOLID) { // Если SOLID сверху
                    yp = 0.5 * dy; // Расстояние до центра грани
                }
                double Up = u(i, j); // Скорость параллельная стенке
                double u_tau = calculate_u_tau(Up, yp);
                double tau_w = rho_ * u_tau * u_tau * std::copysign(1.0, Up);
                double kP = std::max(constants_.k_min, u_tau * u_tau / std::sqrt(constants_.Cmu));
                double epsP = std::max(constants_.epsilon_min, std::pow(std::abs(u_tau), 3.0) / (constants_.kappa * yp));

                k_->operator()(i, j) = kP;
                epsilon_->operator()(i, j) = epsP;
                tau_wx_->operator()(i, j) = tau_w; // X-компонента tau
                tau_wy_->operator()(i, j) = 0.0;   // Y-компонента = 0
            }

            // --- Западная стенка (SOLID слева, т.к. i=0 это вход) ---
            if (i > 0 && tags(i - 1, j) == CellTag::SOLID) {
                double xp = 0.5 * dx; // Расстояние до центра грани
                double Vp = v(i, j); // Скорость параллельная стенке (теперь v)
                double u_tau = calculate_u_tau(Vp, xp); // Используем ту же функцию!
                double tau_w = rho_ * u_tau * u_tau * std::copysign(1.0, Vp);
                double kP = std::max(constants_.k_min, u_tau * u_tau / std::sqrt(constants_.Cmu));
                double epsP = std::max(constants_.epsilon_min, std::pow(std::abs(u_tau), 3.0) / (constants_.kappa * xp));

                k_->operator()(i, j) = kP;
                epsilon_->operator()(i, j) = epsP;
                tau_wx_->operator()(i, j) = 0.0;   // X-компонента = 0
                tau_wy_->operator()(i, j) = tau_w; // Y-компонента tau
            }
             // --- Восточная стенка (SOLID справа, т.к. i=nx-1 это выход) ---
            else if (i < nx - 1 && tags(i + 1, j) == CellTag::SOLID) {
                double xp = 0.5 * dx; // Расстояние до центра грани
                double Vp = v(i, j); // Скорость параллельная стенке
                double u_tau = calculate_u_tau(Vp, xp);
                double tau_w = rho_ * u_tau * u_tau * std::copysign(1.0, Vp);
                double kP = std::max(constants_.k_min, u_tau * u_tau / std::sqrt(constants_.Cmu));
                double epsP = std::max(constants_.epsilon_min, std::pow(std::abs(u_tau), 3.0) / (constants_.kappa * xp));

                k_->operator()(i, j) = kP;
                epsilon_->operator()(i, j) = epsP;
                tau_wx_->operator()(i, j) = 0.0;   // X-компонента = 0
                tau_wy_->operator()(i, j) = tau_w; // Y-компонента tau
            }
        }
    }
    // После этого цикла ячейки у стен обновлены, и поля tau_w* содержат напряжения
}

double KEpsilonModel::calculate_u_tau(double Up, double yp) const {
    if (std::abs(Up) < 1e-9 || yp < 1e-12) {
        return 0.0; // Если скорость или расстояние нулевые
    }

    const double tol = 1e-6;
    const int max_iter = 20;
    double u_tau_old = std::max(1e-4, std::sqrt(constants_.Cmu) * std::abs(Up)); // Начальное приближение
    double u_tau_new = u_tau_old;

    for(int iter = 0; iter < max_iter; ++iter) {
        double y_plus = u_tau_old * yp / nu_molecular_;
        if (y_plus < 11.225) { // Порог логарифмической зоны (прибл.)
            // Используем линейный закон u+ = y+ => u_tau = Up / y+ = Up * nu / (u_tau_old * yp) - нелинейно
            // Проще использовать mu * Up / yp как приближение tau_w в вязкой зоне
             u_tau_new = std::sqrt(nu_molecular_ * std::abs(Up) / yp); // Приближение для вязкого подслоя
             // Или использовать сшивку? Пока просто выходим со значением для вязкого подслоя
             // std::cout << "Warning: y+ < 11.225 in Wall Function. y+ = " << y_plus << std::endl;
             u_tau_old = u_tau_new; // Используем его для расчета k, eps
             break; // Выходим из итераций
        }

        // Логарифмический закон
        double u_plus_calc = (1.0 / constants_.kappa) * std::log(constants_.E_log * y_plus);
        // Обновляем u_tau (простая подстановка с релаксацией)
        double u_tau_calc = Up / u_plus_calc;
        u_tau_new = 0.7 * u_tau_old + 0.3 * u_tau_calc; // Релаксация для стабильности

        if (std::abs(u_tau_new - u_tau_old) < tol * u_tau_old) {
            u_tau_old = u_tau_new;
            break;
        }
        u_tau_old = u_tau_new;
        if (iter == max_iter - 1) {
            std::cout << "Warning: calculate_u_tau did not converge." << std::endl;
        }
    }
    return std::max(0.0, u_tau_old); // Возвращаем неотрицательное значение
}

std::pair<double, double> KEpsilonModel::get_wall_shear_stress(std::size_t i, std::size_t j, BoundarySide side) const {
    if (!tau_wx_ || !tau_wy_) return {0.0, 0.0}; // Если поля не созданы

    // Проверяем, что ячейка (i,j) - жидкость (иначе tau не имеет смысла)
    if (geom_.tags()(i,j) != CellTag::FLUID) return {0.0, 0.0};

    // Возвращаем сохраненное значение для этой ячейки
    // Примечание: Эта реализация не идеальна, т.к. не различает, к какой из 4х
    // возможных стенок относится напряжение, сохраненное в ячейке (i,j).
    // Но для простого случая канала с препятствием может сработать.
    // Правильнее было бы хранить tau на гранях.
    // Пока возвращаем оба сохраненных значения:
    double tau_x = tau_wx_->operator()(i,j);
    double tau_y = tau_wy_->operator()(i,j);

    // Можно добавить логику, чтобы обнулить ненужную компоненту, если знаем side
     switch(side) {
         case BoundarySide::BOTTOM:
         case BoundarySide::TOP:
             return {tau_x, 0.0}; // Для горизонтальных стенок важна tau_x
         case BoundarySide::LEFT:
         case BoundarySide::RIGHT:
              return {0.0, tau_y}; // Для вертикальных стенок важна tau_y
         default:
              return {0.0, 0.0};
     }
}

void KEpsilonModel::solve_k_equation(const Field2D<double>& u, const Field2D<double>& v, double dt) {
    if (!k_ || !epsilon_old_ || !nu_t_ || !Pk_ || !k_old_) {
        throw std::runtime_error("Required fields not initialized for solve_k_equation!");
    }

    const auto nx = geom_.mesh().nx();
    const auto ny = geom_.mesh().ny();
    const double dx = geom_.mesh().dx();
    const double dy = geom_.mesh().dy();
    const auto& tags = geom_.tags();

    // k_old_ содержит значения k на старом временном слое (после swap в solve_step)
    // Результат будем записывать в k_

    // Цикл по внутренним ячейкам
#ifdef USE_OPENMP
    #pragma omp parallel for collapse(2)
#endif
    for (std::size_t j = 1; j < ny - 1; ++j) {
        for (std::size_t i = 1; i < nx - 1; ++i) {

            // Пропускаем SOLID ячейки
            if (tags(i, j) == CellTag::SOLID) {
                 k_->operator()(i,j) = 0.0; // Убедимся, что в солидах 0
                 continue;
             }

            // Проверяем, не является ли ячейка пристенной (где k задается через WF)
            // Простой способ: проверяем соседей. Если хоть один - стена (SOLID или граница j=0/ny-1),
            // то пропускаем расчет здесь, т.к. apply_wall_functions уже задал значение в k_.
            // (Более строгий способ - использовать отдельный флаг для WF ячеек).
             bool is_wall_adjacent = (j == 1 || tags(i, j - 1) == CellTag::SOLID ||
                                      j == ny - 2 || tags(i, j + 1) == CellTag::SOLID ||
                                      i == 1 || tags(i - 1, j) == CellTag::SOLID || // Не учитываем вход/выход i=0, i=nx-1 как стенки
                                      i == nx - 2 || tags(i + 1, j) == CellTag::SOLID);

             if (is_wall_adjacent) {
                 // Значение k в этой ячейке уже установлено через apply_wall_functions.
                 // Но так как мы сделали swap в начале solve_step, текущее значение
                 // в k_ - это старое значение. Нужно скопировать из k_old_ обратно.
                 // (Это немного костыльно из-за swap'а, но проще, чем усложнять solve_step)
                 k_->operator()(i,j) = k_old_->operator()(i,j); // Просто оставляем значение от WF
                 continue;
             }


            // --- 1. Адвекция (Upwind 1-го порядка) ---
            double k_old_ij = k_old_->operator()(i, j); // k на старом шаге
            double u_ij = u(i, j);
            double v_ij = v(i, j);

            // Градиенты для адвекции (используем k_old_)
            double dkdx = (u_ij > 0) ? (k_old_ij - k_old_->operator()(i - 1, j)) / dx
                                     : (k_old_->operator()(i + 1, j) - k_old_ij) / dx;
            double dkdy = (v_ij > 0) ? (k_old_ij - k_old_->operator()(i, j - 1)) / dy
                                     : (k_old_->operator()(i, j + 1) - k_old_ij) / dy;
            double advection_term = -(u_ij * dkdx + v_ij * dkdy);

            // --- 2. Диффузия (Центр. разности, явная схема) ---
            // Эффективная диффузивность для k
            double nu_t_ij = std::max(0.0, nu_t_->operator()(i, j)); // Берем nu_t из текущего поля
            double nu_eff_k = nu_molecular_ + nu_t_ij / constants_.sigk;

            // Лапласиан от k_old_
            double laplacian_k = (k_old_->operator()(i + 1, j) + k_old_->operator()(i - 1, j) - 2.0 * k_old_ij) / (dx * dx)
                               + (k_old_->operator()(i, j + 1) + k_old_->operator()(i, j - 1) - 2.0 * k_old_ij) / (dy * dy);
            double diffusion_term = nu_eff_k * laplacian_k;

            // --- 3. Источники/Стоки (Pk - epsilon) ---
            // Берем Pk и epsilon со старого временного слоя для явной схемы
            double Pk_ij = std::max(0.0, Pk_->operator()(i, j)); // Pk уже вычислен на основе u, v, nu_t
            double epsilon_old_ij = std::max(constants_.epsilon_min, epsilon_old_->operator()(i, j)); // Используем старое epsilon
            double source_term = Pk_ij - epsilon_old_ij;

            // --- 4. Обновление (Явный Эйлер) ---
            // k_ новое = k_ старое + dt * (Адвекция + Диффузия + Источник)
            // Результат пишем в k_ (т.к. в начале solve_step был swap k_ и k_old_)
            k_->operator()(i, j) = k_old_ij + dt * (advection_term + diffusion_term + source_term);

            // Ограничение снизу (хотя лучше делать это после всех шагов)
            // k_->operator()(i, j) = std::max(constants_.k_min, k_->operator()(i, j));

        } // end for i
    } // end for j

    // Граничные условия для k уже применены в apply_bc
}

void KEpsilonModel::solve_epsilon_equation(const Field2D<double>& u, const Field2D<double>& v, double dt) {
    if (!k_ || !epsilon_ || !nu_t_ || !Pk_ || !k_old_ || !epsilon_old_) {
        throw std::runtime_error("Required fields not initialized for solve_epsilon_equation!");
    }

    const auto nx = geom_.mesh().nx();
    const auto ny = geom_.mesh().ny();
    const double dx = geom_.mesh().dx();
    const double dy = geom_.mesh().dy();
    const auto& tags = geom_.tags();

    // epsilon_old_ содержит значения epsilon на старом временном слое
    // k_old_ содержит значения k на старом временном слое
    // Результат будем записывать в epsilon_

#ifdef USE_OPENMP
    #pragma omp parallel for collapse(2)
#endif
    for (std::size_t j = 1; j < ny - 1; ++j) {
        for (std::size_t i = 1; i < nx - 1; ++i) {
            if (tags(i, j) == CellTag::SOLID) {
                epsilon_->operator()(i,j) = 0.0; // Или epsilon_min?
                continue;
            }

            // Пропускаем ячейки у стенок (где epsilon задается WF)
            bool is_wall_adjacent = (j == 1 || tags(i, j - 1) == CellTag::SOLID ||
                                     j == ny - 2 || tags(i, j + 1) == CellTag::SOLID ||
                                     i == 1 || tags(i - 1, j) == CellTag::SOLID ||
                                     i == nx - 2 || tags(i + 1, j) == CellTag::SOLID);

            if (is_wall_adjacent) {
                // Копируем значение от WF из epsilon_old_ обратно в epsilon_
                epsilon_->operator()(i,j) = epsilon_old_->operator()(i,j);
                continue;
            }

            // --- Значения на старом шаге ---
            double eps_old_ij = epsilon_old_->operator()(i, j);
            // k со старого шага, но не меньше k_min (для деления в источнике)
            double k_old_ij_clipped = std::max(constants_.k_min, k_old_->operator()(i, j));
            double u_ij = u(i, j);
            double v_ij = v(i, j);

            // --- 1. Адвекция epsilon (Upwind 1-го порядка) ---
            double depsdx = (u_ij > 0) ? (eps_old_ij - epsilon_old_->operator()(i - 1, j)) / dx
                                       : (epsilon_old_->operator()(i + 1, j) - eps_old_ij) / dx;
            double depsdy = (v_ij > 0) ? (eps_old_ij - epsilon_old_->operator()(i, j - 1)) / dy
                                       : (epsilon_old_->operator()(i, j + 1) - eps_old_ij) / dy;
            double advection_term = -(u_ij * depsdx + v_ij * depsdy);

            // --- 2. Диффузия epsilon (Центр. разности, явная схема) ---
            double nu_t_ij = std::max(0.0, nu_t_->operator()(i, j));
            // Эффективная диффузивность для epsilon
            double nu_eff_eps = nu_molecular_ + nu_t_ij / constants_.sige; // Используем sigma_epsilon

            // Лапласиан от epsilon_old_
            double laplacian_eps = (epsilon_old_->operator()(i + 1, j) + epsilon_old_->operator()(i - 1, j) - 2.0 * eps_old_ij) / (dx * dx)
                                 + (epsilon_old_->operator()(i, j + 1) + epsilon_old_->operator()(i, j - 1) - 2.0 * eps_old_ij) / (dy * dy);
            double diffusion_term = nu_eff_eps * laplacian_eps;

            // --- 3. Источники/Стоки S_eps = (C1e * Pk - C2e * epsilon) * epsilon / k ---
            double Pk_ij = std::max(0.0, Pk_->operator()(i, j)); // Используем уже вычисленный Pk
            // Гарантируем, что epsilon_old_ij >= epsilon_min для использования в члене C2e
            double eps_old_ij_clipped = std::max(constants_.epsilon_min, eps_old_ij);

            // Вычисляем источниковый член S_epsilon, используя значения со старого шага
            double source_term = (constants_.C1e * Pk_ij - constants_.C2e * eps_old_ij_clipped) * eps_old_ij_clipped / k_old_ij_clipped;

            // --- ВАЖНОЕ ЗАМЕЧАНИЕ О СТАБИЛЬНОСТИ ---
            // Явная обработка источникового члена S_epsilon, особенно стока (-C2e * epsilon^2 / k),
            // часто приводит к численным проблемам и нестабильности, особенно если dt большое.
            // Более робастные схемы используют полу-неявную обработку стокового члена,
            // например: S_eps ~ (C1e * Pk * eps_old / k_old) - (C2e * eps_old / k_old) * epsilon_NEW
            // Это требует модификации явной схемы Эйлера (например, решение простого лин. уравнения для epsilon_NEW).
            // Пока оставляем полностью явную схему для простоты.

            // --- 4. Обновление (Явный Эйлер) ---
            // epsilon_ новое = epsilon_ старое + dt * (Адвекция + Диффузия + Источник)
            // Результат пишем в epsilon_
            epsilon_->operator()(i, j) = eps_old_ij + dt * (advection_term + diffusion_term + source_term);

             // Ограничение снизу (опять же, лучше делать после всех шагов)
             // epsilon_->operator()(i, j) = std::max(constants_.epsilon_min, epsilon_->operator()(i, j));

        } // end for i
    } // end for j
}

void KEpsilonModel::solve_step(const Field2D<double>& u, const Field2D<double>& v, double dt) 
{
    if (!k_ || !epsilon_ || !nu_t_ || !Pk_ || !k_old_ || !epsilon_old_) {
        throw std::runtime_error("k-epsilon fields not properly initialized!");
    }
    if (!k_ || !epsilon_ || !nu_t_) {
         throw std::runtime_error("k-epsilon fields not initialized!");
    }

    k_old_->swap(*k_);
    epsilon_old_->swap(*epsilon_);

    apply_bc(u, v); // Обновит k, eps у стенок и tau_w
    calculate_Pk(u, v); // Вычислит Pk_
    solve_k_equation(u, v, dt); // Решит для k_
    solve_epsilon_equation(u, v, dt); // Решит для epsilon_

    // Обеспечение положительности
    const auto nx = geom_.mesh().nx();
    const auto ny = geom_.mesh().ny();
    const auto& cell_tags = geom_.tags();
    for (std::size_t j = 0; j < ny; ++j) {
       for (std::size_t i = 0; i < nx; ++i) {
           if (cell_tags(i, j) == CellTag::FLUID) {
                k_->operator()(i,j) = std::max(constants_.k_min, k_->operator()(i,j)); // Используем константу
                epsilon_->operator()(i,j) = std::max(constants_.epsilon_min, epsilon_->operator()(i,j)); // Используем константу
           } else {
                k_->operator()(i,j) = 0.0;
                epsilon_->operator()(i,j) = 0.0; // Или другое ГУ
           }
       }
   }

   update_nu_t(); // Обновит nu_t_ на основе новых k_ и epsilon_
}

void KEpsilonModel::update_nu_t() {
    const auto nx = geom_.mesh().nx();
    const auto ny = geom_.mesh().ny();
    const auto& cell_tags = geom_.tags();

    for(std::size_t j=0; j<ny; ++j) {
       for(std::size_t i=0; i<nx; ++i) {
           if(cell_tags(i,j) == CellTag::FLUID) {
               // Используем константу из структуры
               nu_t_->operator()(i,j) = constants_.Cmu * (k_->operator()(i,j) * k_->operator()(i,j)) / epsilon_->operator()(i,j);
               nu_t_->operator()(i,j) = std::max(0.0, nu_t_->operator()(i,j)); // Неотрицательность nu_t
           } else {
               nu_t_->operator()(i,j) = 0.0;
           }
       }
   }
}

const Field2D<double>* KEpsilonModel::nu_t() const {
    return nu_t_.get();
}

void KEpsilonModel::calculate_Pk(const Field2D<double>& u, const Field2D<double>& v) {
    if (!Pk_ || !nu_t_) { // Проверяем, что нужные поля инициализированы
        throw std::runtime_error("Pk or nu_t field not initialized in KEpsilonModel!");
    }

    const auto nx = geom_.mesh().nx();
    const auto ny = geom_.mesh().ny();
    const double dx = geom_.mesh().dx();
    const double dy = geom_.mesh().dy();
    const auto& tags = geom_.tags(); // Получаем доступ к тегам ячеек

    // Предварительно вычисляем инвертированные удвоенные шаги для эффективности
    // Добавляем проверку на нулевой шаг
    const double inv_2dx = (std::abs(dx) > 1e-12) ? 1.0 / (2.0 * dx) : 0.0;
    const double inv_2dy = (std::abs(dy) > 1e-12) ? 1.0 / (2.0 * dy) : 0.0;

    // Обнуляем поле Pk перед расчетом (важно для граничных/твердых ячеек)
    Pk_->fill(0.0);

    // Используем OpenMP, если включен (итерации независимы)
#ifdef USE_OPENMP
    #pragma omp parallel for collapse(2)
#endif
    // Цикл только по ВНУТРЕННИМ ячейкам, где можно безопасно взять соседей для центр. разности
    for (std::size_t j = 1; j < ny - 1; ++j) {
        for (std::size_t i = 1; i < nx - 1; ++i) {

            // Пропускаем не-жидкостные ячейки
            if (tags(i, j) != CellTag::FLUID) {
                continue;
            }

            // Вычисляем градиенты скорости центральными разностями
            // ВАЖНО: Этот простой вариант не учитывает SOLID соседей.
            // Для большей робастности можно было бы использовать односторонние разности
            // у границ SOLID или модифицировать центральные разности.
            // Пока оставляем так для простоты.
            double du_dx = (u(i + 1, j) - u(i - 1, j)) * inv_2dx;
            double dv_dy = (v(i, j + 1) - v(i, j - 1)) * inv_2dy;
            double du_dy = (u(i, j + 1) - u(i, j - 1)) * inv_2dy;
            double dv_dx = (v(i + 1, j) - v(i - 1, j)) * inv_2dx;

            // Вычисляем S^2 = 2*(S11^2 + S22^2) + 4*S12^2 = 2*( (du/dx)^2 + (dv/dy)^2 ) + (du/dy + dv/dx)^2
            double S_sq = 2.0 * (du_dx * du_dx + dv_dy * dv_dy) + (du_dy + dv_dx) * (du_dy + dv_dx);

            // Гарантируем неотрицательность S^2 (хотя по формуле она должна быть >=0)
            S_sq = std::max(0.0, S_sq);

            // Получаем nu_t для текущей ячейки (гарантируем неотрицательность)
            double nu_t_ij = std::max(0.0, nu_t_->operator()(i, j));

            // Вычисляем генерацию Pk = nu_t * S^2
            double Pk_ij = nu_t_ij * S_sq;

            // Сохраняем Pk (тоже гарантируем неотрицательность)
            Pk_->operator()(i, j) = std::max(0.0, Pk_ij);
        }
    }
    // Примечание: Для ячеек на границах (i=0, nx-1, j=0, ny-1) и в SOLID
    // значение Pk_ останется равным 0.0 после Pk_->fill(0.0).
    // Это обычно приемлемо, т.к. уравнения k/epsilon там не решаются или
    // перезаписываются граничными условиями / WF.
}
} // namespace cfd