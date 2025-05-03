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
    TurbulenceModelType turb_type,
    PoissonType ptype,
    double cfl, 
    double omega,
    unsigned max_p_iter,
    double p_tol)
    : 
    Solver(geom, rho, nu, turb_type),
    u_(geom.mesh().nx(), geom.mesh().ny(), 0.0),
    v_(u_), p_(u_),
    u_star_(u_), v_star_(u_), rhs_(u_),
    tag_(geom.tags()),
    cfl_(cfl),
    max_pressure_iter_(max_p_iter), pressure_tol_(p_tol),
    poisson_(ptype, omega) {}

//---------------------------------------------------------------- inlet --
void VelocityPressureSolver::set_inlet_parabola(double umax)
{
    const std::size_t ny = geom_.mesh().ny();
    const double H = geom_.mesh().dy() * (ny - 1);
    for (std::size_t j = 0; j < ny; ++j) {
        double y = j * geom_.mesh().dy();
        u_(0, j) = 4.0 * umax * y * (H - y) / (H * H);
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
    double max_nu_eff = nu_molecular_;
    const Field2D<double>* nu_t_field = turbulence_model_ ? turbulence_model_->nu_t() : nullptr;

    if (nu_t_field) {
        double current_max_nu_t = 0.0;
#ifdef USE_OPENMP
    #pragma omp parallel for collapse(2)
#endif
        for(std::size_t j=0; j<geom_.mesh().ny(); ++j) {
            for(std::size_t i=0; i<geom_.mesh().nx(); ++i) {
                if(tag_(i,j) == CellTag::FLUID) {
                    current_max_nu_t = std::max(current_max_nu_t, (*nu_t_field)(i,j));
                }
            }
        }
        max_nu_eff += current_max_nu_t;
    }

    if (max_nu_eff <= 1e-12) {
         return std::numeric_limits<double>::max();
    }
    double dx = geom_.mesh().dx();
    double dy = geom_.mesh().dy();
    return 0.5 / (max_nu_eff * (1.0 / (dx * dx) + 1.0 / (dy * dy)));
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

            if (tag_(i, j) == CellTag::SOLID) continue;
            
            double nu_eff = get_effective_viscosity(i, j);

            if (nu_eff < 1e-12) continue;
            double alpha = nu_eff * dt;

            // Вычисляем Лапласиан, используя значения из f_old (время n)
            double lap = (f_old(i + 1, j) + f_old(i - 1, j) - 2.0 * f_old(i, j)) / dx2
                       + (f_old(i, j + 1) + f_old(i, j - 1) - 2.0 * f_old(i, j)) / dy2;

            // Обновляем значение f на новом шаге (время n+1)
            // f^{n+1} = f^n + alpha * lap(f^n)
            f(i, j) = f_old(i, j) + alpha * lap;
        }
    }
}

void VelocityPressureSolver::calculatePoissonRHS_RhieChow(Field2D<double>& rhs, double dt)
{
    const std::size_t nx = geom_.mesh().nx();
    const std::size_t ny = geom_.mesh().ny();
    const double dx = geom_.mesh().dx();
    const double dy = geom_.mesh().dy();

#ifdef USE_OPENMP
    #pragma omp parallel for collapse(2) // Потенциально нужны доп. меры для OMP здесь
#endif
    for (std::size_t j = 1; j < ny - 1; ++j) { // Цикл по внутренним ячейкам
        for (std::size_t i = 1; i < nx - 1; ++i) {
            if (tag_(i, j) == CellTag::SOLID) {
                rhs(i, j) = 0.0; // Используем параметр rhs напрямую
                continue;
            }

            // Скорость u* на восточной грани (i+1/2, j)
            double u_e_left  = u_star_(i, j);
            double u_e_right = (tag_(i + 1, j) == CellTag::SOLID) ? 0.0 : u_star_(i + 1, j);
            double u_e = 0.5 * (u_e_left + u_e_right);

            // Скорость u* на западной грани (i-1/2, j)
            double u_w_left  = (tag_(i - 1, j) == CellTag::SOLID) ? 0.0 : u_star_(i - 1, j);
            double u_w_right = u_star_(i, j);
            double u_w = 0.5 * (u_w_left + u_w_right);

            // Скорость v* на северной грани (i, j+1/2)
            double v_n_bottom = v_star_(i, j);
            double v_n_top    = (tag_(i, j + 1) == CellTag::SOLID) ? 0.0 : v_star_(i, j + 1);
            double v_n = 0.5 * (v_n_bottom + v_n_top);

            // Скорость v* на южной грани (i, j-1/2)
            double v_s_bottom = (tag_(i, j - 1) == CellTag::SOLID) ? 0.0 : v_star_(i, j - 1);
            double v_s_top    = v_star_(i, j);
            double v_s = 0.5 * (v_s_bottom + v_s_top);

            // Дивергенция через скорости на гранях
            double div = (u_e - u_w) / dx + (v_n - v_s) / dy;

            // Правая часть уравнения Пуассона
            rhs(i, j) = rho_ * div / dt; // Используем параметр rhs напрямую
        }
    }
    // Здесь можно добавить обработку RHS на границах, если это необходимо
    // по схеме (хотя обычно решатель Пуассона обрабатывает ГУ для p)
}


//-------------------------------- project -------------------------------
void VelocityPressureSolver::project(double dt)
{
    const std::size_t nx = geom_.mesh().nx();
    const std::size_t ny = geom_.mesh().ny();
    const double dx = geom_.mesh().dx();
    const double dy = geom_.mesh().dy();

    
    // 1. Вычисляем правую часть для уравнения Пуассона (rhs_) (Rhie-chow like)
// #ifdef USE_OPENMP
//     #pragma omp parallel for collapse(2)
// #endif
//     for (std::size_t j = 1; j < ny - 1; ++j) {
//         for (std::size_t i = 1; i < nx - 1; ++i) {
//             if (tag_(i, j) == CellTag::SOLID) {
//                 rhs_(i, j) = 0.0;
//                 continue;
//             }
//             double div=(u_star_(i+1,j)-u_star_(i-1,j))/(2*dx) + (v_star_(i,j+1)-v_star_(i,j-1))/(2*dy);
//             rhs_(i,j)=rho_*div/dt;
//         }
//     }
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
#ifdef USE_OPENMP
    #pragma omp parallel for collapse(2)
#endif
    for (std::size_t j = 1; j < ny - 1; ++j) {
        for (std::size_t i = 1; i < nx - 1; ++i) {
            if (tag_(i, j) == CellTag::SOLID) {
                continue;
            }

             double p_ip1 = (tag_(i+1,j) == CellTag::SOLID) ? p_(i,j) : p_(i+1,j);
             double p_im1 = (tag_(i-1,j) == CellTag::SOLID) ? p_(i,j) : p_(i-1,j);
             double p_jp1 = (tag_(i,j+1) == CellTag::SOLID) ? p_(i,j) : p_(i,j+1);
             double p_jm1 = (tag_(i,j-1) == CellTag::SOLID) ? p_(i,j) : p_(i,j-1);

             double dpdx_center = (p_ip1 - p_im1) / (2.0 * dx);
             double dpdy_center = (p_jp1 - p_jm1) / (2.0 * dy);

            u_(i, j) = u_star_(i, j) - dt / rho_ * dpdx_center;
            v_(i, j) = v_star_(i, j) - dt / rho_ * dpdy_center;
        }
    }
}

//-------------------------------- main step -----------------------------
void VelocityPressureSolver::step(double dt_user)
{
    double dt_cfl = compute_cfl_dt(cfl_);
    double dt_diff = compute_diff_dt();
    double dt = dt_user;
    // Выбираем минимальный шаг из трех: пользовательский, CFL, диффузионный
    dt = std::min(dt, dt_cfl);
    dt = std::min(dt, dt_diff); 

    apply_bc();
    u_star_ = u_; v_star_ = v_;
    advect(u_star_, u_, v_, dt);
    advect(v_star_, u_, v_, dt);
    diffuse(u_star_, dt);
    diffuse(v_star_, dt);

    if (turbulence_model_) {
       turbulence_model_->solve_step(u_, v_, dt); // Передаем текущие u_, v_
    }

    project(dt);

    last_monitor_info_.actualDt = dt;
}

} // namespace cfd