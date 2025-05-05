#include "solvers/velocity_pressure/VelocityPressureSolver.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#ifdef USE_OPENMP
 #include <omp.h>
#endif

namespace cfd {

//------------------------------------------------------------------------
VelocityPressureSolver::VelocityPressureSolver(
    const Geometry& geom,
    double rho,
    double nu,
    PoissonType ptype,
    double cfl, 
    double omega,
    unsigned max_p_iter,
    double p_tol)
    : 
    Solver(geom, nu),
    u_(geom.mesh().nx(), geom.mesh().ny(), 0.0),
    v_(u_), p_(u_),
    u_star_(u_), v_star_(u_), rhs_(u_),
    tag_(geom.tags()),
    rho_(rho), cfl_(cfl),
    max_pressure_iter_(max_p_iter), pressure_tol_(p_tol),
    poisson_(ptype, omega) {}

//---------------------------------------------------------------- inlet --
void VelocityPressureSolver::set_inlet_parabola(double umax)
{
    const std::size_t ny = geom_.mesh().ny();
    const double Ly = geom_.mesh().Ly(); // <-- Используем полную высоту канала
    const double dy = geom_.mesh().dy(); // <-- Шаг по y

    // Проверка на случай нулевой высоты канала
    if (std::abs(Ly) < 1e-12) {
        for (std::size_t j = 0; j < ny; ++j) {
            u_(0, j) = 0.0;
            v_(0, j) = 0.0;
        }
        return; // Выходим, если высота нулевая
    }

    const double Ly_sq = Ly * Ly; // <-- Предварительно вычисляем H^2 = Ly^2

    for (std::size_t j = 0; j < ny; ++j) {
        // Вычисляем y-координату ЦЕНТРА ячейки j
        double y_center = (j + 0.5) * dy;

        u_(0, j) = 4.0 * umax * y_center * (Ly - y_center) / Ly_sq;

        // v-компонента на входе равна нулю
        v_(0, j) = 0.0;
    }
}

//---------------------------------------------------------------- bc -----
void VelocityPressureSolver::apply_bc()
{
    const std::size_t nx = geom_.mesh().nx();
    const std::size_t ny = geom_.mesh().ny();

    for (std::size_t i = 0; i < nx; ++i) {
        u_(i, 0) = v_(i, 0) = 0.0;
        u_(i, ny - 1) = v_(i, ny - 1) = 0.0;
    }
    for (std::size_t j = 0; j < ny; ++j) {
        u_(nx - 1, j) = u_(nx - 2, j);
        v_(nx - 1, j) = v_(nx - 2, j);
    }
    for (std::size_t j = 1; j < ny - 1; ++j)
        for (std::size_t i = 1; i < nx - 1; ++i) {
            if (tag_(i,j)==CellTag::SOLID) { u_(i,j)=v_(i,j)=0.0; continue; }
            if (tag_(i-1,j)==CellTag::SOLID||tag_(i+1,j)==CellTag::SOLID||
                tag_(i,j-1)==CellTag::SOLID||tag_(i,j+1)==CellTag::SOLID)
                u_(i,j)=v_(i,j)=0.0;
        }
}

//------------------------------ вспом. -----------------------------------
double VelocityPressureSolver::compute_cfl_dt(double safety) const
{
    double umax=0.0, vmax=0.0;
    #ifdef USE_OPENMP
       #pragma omp parallel for reduction(max:umax, vmax) collapse(2)
    #endif
    for(std::size_t j=0;j<geom_.mesh().ny();++j)
        for(std::size_t i=0;i<geom_.mesh().nx();++i){
            umax = std::max(umax, std::fabs(u_(i,j)));
            vmax = std::max(vmax, std::fabs(v_(i,j)));
        }
    double dt_x = umax>0? geom_.mesh().dx()/umax : 1e9;
    double dt_y = vmax>0? geom_.mesh().dy()/vmax : 1e9;
    return safety * std::min(dt_x, dt_y);
}

// ----------------------------- diff dt -----------------------------------
double VelocityPressureSolver::compute_diff_dt() const 
{
    if (nu_ < 1e-12) {
        return std::numeric_limits<double>::max();
    }
    double dx = geom_.mesh().dx();
    double dy = geom_.mesh().dy();
    return 0.5 / (nu_ * (1.0 / (dx * dx) + 1.0 / (dy * dy)));
}

//-------------------------------- adv / diff -----------------------------
void VelocityPressureSolver::advect(Field2D<double>& f,const Field2D<double>& u,const Field2D<double>& v,double dt)
{
    const double dx=geom_.mesh().dx(), dy=geom_.mesh().dy();
    const std::size_t nx=geom_.mesh().nx(), ny=geom_.mesh().ny();
    Field2D<double> f_old = f;
// #ifdef USE_OPENMP
//     #pragma omp parallel for
// #endif
    for(std::size_t j=1;j<ny-1;++j)
        for(std::size_t i=1;i<nx-1;++i){
            if(tag_(i,j)==CellTag::SOLID){ f(i,j)=0; continue; }
            double dfdx = (u(i,j)>0)? (f_old(i,j)-f_old(i-1,j))/dx : (f_old(i+1,j)-f_old(i,j))/dx;
            double dfdy = (v(i,j)>0)? (f_old(i,j)-f_old(i,j-1))/dy : (f_old(i,j+1)-f_old(i,j))/dy;
            f(i,j) = f_old(i,j) - dt*(u(i,j)*dfdx + v(i,j)*dfdy);
        }
}

void VelocityPressureSolver::diffuse(Field2D<double>& f, double dt)
{
    // Коэффициент диффузии за шаг dt
    double alpha = nu_ * dt;

    // Если вязкость или шаг нулевые, диффузии нет
    if (alpha == 0.0) {
        return;
    }

    const std::size_t nx = geom_.mesh().nx();
    const std::size_t ny = geom_.mesh().ny();
    const double dx = geom_.mesh().dx();
    const double dy = geom_.mesh().dy();
    const double dx2 = dx * dx;
    const double dy2 = dy * dy;

    // Создаем копию поля f на начало шага (время n)
    // Это ВАЖНО для корректного явного метода Эйлера
    Field2D<double> f_old = f;

    // Применяем явное обновление ОДИН раз
#ifdef USE_OPENMP
    #pragma omp parallel for collapse(2) // Можно склеить циклы для OpenMP
#endif
    for (std::size_t j = 1; j < ny - 1; ++j) {
        for (std::size_t i = 1; i < nx - 1; ++i) {
            // Пропускаем твердые тела (или обрабатываем их ГУ)
            if (tag_(i, j) == CellTag::SOLID) {
                // Можно явно задать значение в твердых телах, если нужно
                // f(i, j) = 0.0;
                continue;
            }

            // Вычисляем Лапласиан, используя значения из f_old (время n)
            double lap = (f_old(i + 1, j) + f_old(i - 1, j) - 2.0 * f_old(i, j)) / dx2
                       + (f_old(i, j + 1) + f_old(i, j - 1) - 2.0 * f_old(i, j)) / dy2;

            // Обновляем значение f на новом шаге (время n+1)
            // f^{n+1} = f^n + alpha * lap(f^n)
            f(i, j) = f_old(i, j) + alpha * lap;
        }
    }
}


//-------------------------------- project -------------------------------
void VelocityPressureSolver::project(double dt)
{
    const std::size_t nx = geom_.mesh().nx();
    const std::size_t ny = geom_.mesh().ny();
    const double dx = geom_.mesh().dx();
    const double dy = geom_.mesh().dy();

    
    // 1. Вычисляем правую часть для уравнения Пуассона (rhs_) (Rhie-chow like)
    calculatePoissonRHS_RhieChow(rhs_, dt);

    // 3. Установка ГРАНИЧНЫХ УСЛОВИЙ для ДАВЛЕНИЯ ПЕРЕД решением Пуассона
    // (Этот блок остается без изменений)
#ifdef USE_OPENMP
    #pragma omp parallel for
#endif
    for (std::size_t i = 0; i < nx; ++i) {
        p_(i, 0)    = p_(i, 1);
        p_(i, ny-1) = p_(i, ny-2);
    }
#ifdef USE_OPENMP
    #pragma omp parallel for
#endif
     for (std::size_t j = 0; j < ny; ++j) {
         if (tag_(1, j) != CellTag::SOLID) {
             p_(0, j) = p_(1, j);
         }
          if (tag_(nx-1, j) != CellTag::SOLID && tag_(nx-2, j) != CellTag::SOLID ) {
             p_(nx-1, j) = 0.0;
         } else if (tag_(nx-1, j) != CellTag::SOLID) {
             p_(nx-1, j) = p_(nx-2, j);
         }
    }

    // 4. Решаем Пуассона
    ConvergenceInfo p_info = poisson_.solve(p_, rhs_, geom_, max_pressure_iter_, pressure_tol_);
    // --- Обновляем информацию о сходимости НАПРЯМУЮ в структуре ---
    last_monitor_info_.pressureIterations = p_info.iterations;
    last_monitor_info_.pressureResidual = p_info.residual;

    // --- Опциональная отладка ---
    if(p_info.iterations >= max_pressure_iter_ || p_info.residual > pressure_tol_*10) {
        std::cout << "Warning: Pressure solver convergence issues. Iter: " << p_info.iterations
                  << ", Res: " << p_info.residual << std::endl;
    }
    // --- Конец отладки ---

    // 5. Корректируем скорости
// #ifdef USE_OPENMP
//     #pragma omp parallel for collapse(2)
// #endif
//     for (std::size_t j = 1; j < ny - 1; ++j) {
//         for (std::size_t i = 1; i < nx - 1; ++i) {
//             if (tag_(i, j) == CellTag::SOLID) {
//                 continue;
//             }

//              double p_ip1 = (tag_(i+1,j) == CellTag::SOLID) ? p_(i,j) : p_(i+1,j);
//              double p_im1 = (tag_(i-1,j) == CellTag::SOLID) ? p_(i,j) : p_(i-1,j);
//              double p_jp1 = (tag_(i,j+1) == CellTag::SOLID) ? p_(i,j) : p_(i,j+1);
//              double p_jm1 = (tag_(i,j-1) == CellTag::SOLID) ? p_(i,j) : p_(i,j-1);

//              double dpdx_center = (p_ip1 - p_im1) / (2.0 * dx);
//              double dpdy_center = (p_jp1 - p_jm1) / (2.0 * dy);

//             u_(i, j) = u_star_(i, j) - dt / rho_ * dpdx_center;
//             v_(i, j) = v_star_(i, j) - dt / rho_ * dpdy_center;
//         }
//     }
    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 1; i < nx - 1; ++i) {
            if (tag_(i, j) == CellTag::SOLID) {
                u_(i,j) = 0.0;
                v_(i,j) = 0.0;
                continue;
            }

            // --- Расчет градиентов давления в центре (i,j) ---
            double p_L, p_R, p_B, p_T;

            // Соседи по X для dpdx
            if (i == 1) p_L = p_(1, j); // ГУ Неймана p(0)=p(1) -> используем p(1) как эффективное значение в i=0 для градиента в i=1
            else if (tag_(i - 1, j) == CellTag::SOLID) p_L = p_(i, j); // Нулевой градиент у SOLID
            else p_L = p_(i - 1, j);

            if (i == nx - 2) { // Узел перед выходом
                if (tag_(nx - 1, j) == CellTag::SOLID) p_R = p_(i,j); // Нулевой градиент у SOLID
                else p_R = p_(nx-1, j); // Используем давление на выходе (p=0 или Нейман p(N-1)=p(N-2))
            }
            else if (tag_(i + 1, j) == CellTag::SOLID) p_R = p_(i, j); // Нулевой градиент у SOLID
            else p_R = p_(i + 1, j);

            double dpdx_center = (std::abs(dx) > 1e-12) ? (p_R - p_L) / (2.0 * dx) : 0.0;

            // Соседи по Y для dpdy
            // !!! ОБЪЯВЛЯЕМ dpdy_center ЗДЕСЬ !!!
            double dpdy_center = 0.0;

            if (j == 0) { // Нижняя стенка: ГУ Неймана dp/dy=0
                // dpdy_center остается 0.0
            }
            else if (j == ny - 1) { // Верхняя стенка: ГУ Неймана dp/dy=0
                // dpdy_center остается 0.0
            }
            else { // Внутренняя ячейка (j от 1 до ny-2)
                p_T = (tag_(i, j + 1) == CellTag::SOLID) ? p_(i, j) : p_(i, j + 1); // Нулевой градиент у SOLID
                p_B = (tag_(i, j - 1) == CellTag::SOLID) ? p_(i, j) : p_(i, j - 1); // Нулевой градиент у SOLID
                // Вычисляем градиент для внутренних ячеек
                dpdy_center = (std::abs(dy) > 1e-12) ? (p_T - p_B) / (2.0 * dy) : 0.0;
            }

            // --- Коррекция скорости ---
            u_(i, j) = u_star_(i, j) - dt / rho_ * dpdx_center;
            v_(i, j) = v_star_(i, j) - dt / rho_ * dpdy_center; // Теперь dpdy_center определен всегда

            // --- Принудительно v=0 у горизонтальных стенок ---
            // Это важно, так как коррекция v основана на dpdy=0, но v_star мог быть ненулевым
            if (j == 0 || j == ny-1) {
                v_(i, j) = 0.0;
            }
        }
    }
}

//-------------------------------- main step -----------------------------
void VelocityPressureSolver::step(double dt_user)
{
    //double dt = std::min(dt_user, compute_cfl_dt(cfl_));

    double dt_cfl = compute_cfl_dt(cfl_);
    double dt_diff = compute_diff_dt();
    double dt = dt_user;
    // Выбираем минимальный шаг из трех: пользовательский, CFL, диффузионный
    dt = std::min(dt, dt_cfl);
    dt = std::min(dt, dt_diff); // Учитываем лимит диффузии

    apply_bc();
    u_star_ = u_; v_star_ = v_;
    advect(u_star_, u_, v_, dt);
    advect(v_star_, u_, v_, dt);
    diffuse(u_star_, dt);
    diffuse(v_star_, dt);
    project(dt);

    last_monitor_info_.actualDt = dt;
}

} // namespace cfd