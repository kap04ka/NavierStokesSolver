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
    double u_max_inlet,                 
    double inlet_turb_intensity,      
    double inlet_length_scale_factor, 
    PoissonType ptype,
    double cfl, 
    double omega,
    unsigned max_p_iter,
    double p_tol)
    : 
    Solver(geom, rho, nu, turb_type, u_max_inlet, inlet_turb_intensity, inlet_length_scale_factor),
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

    // for (std::size_t i = 0; i < nx; ++i) {
    //     u_(i, 0) = v_(i, 0) = 0.0;
    //     u_(i, ny - 1) = v_(i, ny - 1) = 0.0;
    // }

    if (nx >= 2) {
        for (std::size_t j = 0; j < ny; ++j) {
            if (tag_(nx - 1, j) == CellTag::FLUID && tag_(nx - 2, j) != CellTag::SOLID) {
                u_(nx - 1, j) = u_(nx - 2, j);
                v_(nx - 1, j) = v_(nx - 2, j);
            }
        }
    }
    
    // и для ячеек FLUID рядом с ними
    for (std::size_t j = 0; j < ny; ++j) { // Идем по всем j, включая 0 и ny-1
        for (std::size_t i = 0; i < nx; ++i) { // Идем по всем i, включая 0 и nx-1
            if (tag_(i,j)==CellTag::SOLID) { // Внутри SOLID всегда 0
               u_(i,j)=v_(i,j)=0.0;
            } else if (tag_(i,j)==CellTag::FLUID) { // Для FLUID проверяем соседей
               // Проверка с защитой от выхода за границы массива
                bool solid_neighbor = false;
                if (i > 0 && tag_(i-1, j) == CellTag::SOLID) solid_neighbor = true;
                if (!solid_neighbor && i < nx - 1 && tag_(i+1, j) == CellTag::SOLID) solid_neighbor = true;
                if (!solid_neighbor && j > 0 && tag_(i, j-1) == CellTag::SOLID) solid_neighbor = true;
                if (!solid_neighbor && j < ny - 1 && tag_(i, j+1) == CellTag::SOLID) solid_neighbor = true;

                if (solid_neighbor) {
                    u_(i,j)=v_(i,j)=0.0; // Задаем прилипание у SOLID препятствий
                }
           }
        }
    }

        
}

//------------------------------ вспом. -----------------------------------
double VelocityPressureSolver::compute_cfl_dt(double safety) const
{
    double umax=0.0, vmax=0.0;
    #ifdef USE_OPENMP
       #pragma omp parallel for reduction(max: umax) reduction(max: vmax) collapse(2)
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
#ifdef USE_OPENMP
    #pragma omp parallel for collapse(2)
#endif
    for(std::size_t j=1;j<ny-1;++j)
        for(std::size_t i=1;i<nx-1;++i){
            if(tag_(i,j)==CellTag::SOLID){ f(i,j)=0; continue; }
            double dfdx = (u(i,j)>0)? (f_old(i,j)-f_old(i-1,j))/dx : (f_old(i+1,j)-f_old(i,j))/dx;
            double dfdy = (v(i,j)>0)? (f_old(i,j)-f_old(i,j-1))/dy : (f_old(i,j+1)-f_old(i,j))/dy;
            f(i,j) = f_old(i,j) - dt*(u(i,j)*dfdx + v(i,j)*dfdy);
        }
}

void VelocityPressureSolver::diffuse_u(double dt) {
    const auto nx = geom_.mesh().nx();
    const auto ny = geom_.mesh().ny();
    const double dx = geom_.mesh().dx();
    const double dy = geom_.mesh().dy();
    const double dx2 = dx * dx;
    const double dy2 = dy * dy;
    const auto& tags = geom_.tags();

    Field2D<double> u_old = u_star_;

    bool turbulent = (turbulence_model_ != nullptr);
    TurbulenceModel* turb_model_ptr = getTurbulenceModel();

#ifdef USE_OPENMP
    // OMP пока убран для ясности логики ГУ
    // #pragma omp parallel for collapse(2)
#endif
    // Идем по внутренним ячейкам, где решаем уравнение
    for (std::size_t j = 1; j < ny - 1; ++j) {
        for (std::size_t i = 1; i < nx - 1; ++i) {
            if (tags(i, j) == CellTag::SOLID) {
                u_star_(i, j) = 0.0;
                continue;
            }

            double nu_eff_ij = get_effective_viscosity(i, j);
            if (nu_eff_ij < 1e-12) continue;

            double lap_ij = 0.0; // Лапласиан для текущей ячейки

            // --- X-компонента Лапласиана (с учетом SOLID) ---
            double u_ip1 = (tags(i + 1, j) == CellTag::SOLID) ? -u_old(i, j) // Ghost cell u=0 -> ug = -ui
                           : u_old(i + 1, j);
            double u_im1 = (tags(i - 1, j) == CellTag::SOLID) ? -u_old(i, j) // Ghost cell u=0 -> ug = -ui
                           : u_old(i - 1, j);
            // Проверяем границы входа/выхода - там ГУ другие
            if (i == 1) u_im1 = u_old(0, j); // Используем значение Дирихле на входе
            if (i == nx - 2) u_ip1 = u_old(i, j); // Используем Неймана на выходе du/dx=0 => u(N)=u(N-1)

            double d2udx2 = (u_ip1 + u_im1 - 2.0 * u_old(i, j)) / dx2;

            // --- Y-компонента Лапласиана (с учетом стенок и режима) ---
            double d2udy2 = 0.0;
            double u_jp1, u_jm1;

            // Верхний сосед j+1
            if (j + 1 >= ny - 1 || tags(i, j + 1) == CellTag::SOLID) { // Верхняя стенка (граница или SOLID)
                if (turbulent && turb_model_ptr) {
                    // Турбулентный: Поток через стенку равен tau_wx/rho
                    std::pair<double, double> tau_w_pair = turb_model_ptr->get_wall_shear_stress(i, j, BoundarySide::TOP);
                    double flux_n = tau_w_pair.first / rho_;
                    // Поток через южную грань (j-1/2) - стандартный
                    double nu_eff_s = 0.5 * (nu_eff_ij + get_effective_viscosity(i, j - 1)); // Сосед j-1 точно FLUID
                    double flux_s = -nu_eff_s * (u_old(i, j) - u_old(i, j - 1)) / dy;
                    // Вся Y-диффузия = (ПотокСевер - ПотокЮг) / dy
                    d2udy2 = (flux_n - flux_s) / dy; // Делим на dy, а не dy2!
                } else { // Ламинарный: u=0 на стенке
                    u_jp1 = 0.0; // Значение НА стенке
                    u_jm1 = u_old(i, j - 1); // Сосед снизу точно FLUID
                    d2udy2 = (u_jp1 + u_jm1 - 2.0 * u_old(i, j)) / dy2; // Стандартный Лапласиан с u_jp1=0
                }
            }
            // Нижний сосед j-1
            else if (j - 1 <= 0 || tags(i, j - 1) == CellTag::SOLID) { // Нижняя стенка (граница или SOLID)
                if (turbulent && turb_model_ptr) {
                    std::pair<double, double> tau_w_pair = turb_model_ptr->get_wall_shear_stress(i, j, BoundarySide::BOTTOM);
                    double flux_s = tau_w_pair.first / rho_;
                    // Поток через северную грань (j+1/2) - стандартный
                    double nu_eff_n = 0.5 * (nu_eff_ij + get_effective_viscosity(i, j + 1)); // Сосед j+1 точно FLUID
                    double flux_n = -nu_eff_n * (u_old(i, j + 1) - u_old(i, j)) / dy;
                    d2udy2 = (flux_n - flux_s) / dy;
                } else { // Ламинарный: u=0 на стенке
                    u_jp1 = u_old(i, j + 1); // Сосед сверху точно FLUID
                    u_jm1 = 0.0; // Значение НА стенке
                    d2udy2 = (u_jp1 + u_jm1 - 2.0 * u_old(i, j)) / dy2;
                }
            }
            // Полностью внутренняя ячейка по Y (и не рядом с SOLID по Y)
            else {
                u_jp1 = u_old(i, j + 1);
                u_jm1 = u_old(i, j - 1);
                d2udy2 = (u_jp1 + u_jm1 - 2.0 * u_old(i, j)) / dy2;
            }

            // Собираем Лапласиан
            // Для турбулентного режима у стенки d2udy2 уже содержит nu_eff (через поток/tau_w)
            // Для ламинарного - нет. Нужно умножить на nu_eff = nu_molecular_
            if (turbulent && turb_model_ptr && (j==1 || j==ny-2 || tags(i,j-1)==CellTag::SOLID || tags(i,j+1)==CellTag::SOLID)) {
                 // В турбулентном режиме у стенки Y-диффузия уже посчитана как (Fn-Fs)/dy
                 // Нужно добавить X-диффузию, умноженную на nu_eff
                 lap_ij = nu_eff_ij * d2udx2 + d2udy2;
            } else {
                 // В ламинарном режиме ИЛИ вдали от стенок в турбулентном - стандартный Лапласиан * nu_eff
                 lap_ij = nu_eff_ij * (d2udx2 + d2udy2);
            }


            // Обновление u_star_ явным Эйлером
            // f_new = f_old + dt * Div(nu_eff * Grad(f_old)) = f_old + dt * lap_ij
            u_star_(i, j) = u_old(i, j) + dt * lap_ij;

        } // end for i
    } // end for j

    // Принудительное применение ГУ Дирихле к u_star_ ПОСЛЕ диффузии
    apply_bc(); // Вызовем apply_bc еще раз, он обнулит u на стенках и у SOLID
                // и применит ГУ Неймана на выходе.
                // Это перезапишет некоторые значения u_star_ после диффузии,
                // обеспечивая точное выполнение ГУ перед проекцией.
}

void VelocityPressureSolver::diffuse_v(double dt) {
    const auto nx = geom_.mesh().nx();
    const auto ny = geom_.mesh().ny();
    const double dx = geom_.mesh().dx();
    const double dy = geom_.mesh().dy();
    const double dx2 = dx * dx;
    const double dy2 = dy * dy;
    const auto& tags = geom_.tags();

    Field2D<double> v_old = v_star_;

    bool turbulent = (turbulence_model_ != nullptr);
    TurbulenceModel* turb_model_ptr = getTurbulenceModel();

    for (std::size_t j = 1; j < ny - 1; ++j) {
         for (std::size_t i = 1; i < nx - 1; ++i) {
             if (tags(i, j) == CellTag::SOLID) {
                 v_star_(i, j) = 0.0;
                 continue;
             }

             double nu_eff_ij = get_effective_viscosity(i, j);
             if (nu_eff_ij < 1e-12) continue;

             double lap_ij = 0.0;
             double d2vdx2 = 0.0;
             double d2vdy2 = 0.0;

             // X-производные (вертикальные стенки - SOLID)
             double v_ip1, v_im1;
             if (i + 1 >= nx - 1 || tags(i + 1, j) == CellTag::SOLID) { // SOLID справа или выход
                  if (tags(i + 1, j) == CellTag::SOLID) { // SOLID
                     if (turbulent && turb_model_ptr) { // Турбулентный SOLID
                         std::pair<double, double> tau_w_pair = turb_model_ptr->get_wall_shear_stress(i, j, BoundarySide::RIGHT);
                         double flux_e = tau_w_pair.second / rho_; // tau_wy
                         // Западный поток
                         double nu_eff_w = 0.5 * (nu_eff_ij + get_effective_viscosity(i - 1, j)); // i-1 точно FLUID
                         double flux_w = -nu_eff_w * (v_old(i, j) - v_old(i - 1, j)) / dx;
                         d2vdx2 = (flux_e - flux_w) / dx;
                     } else { // Ламинарный SOLID
                         v_ip1 = 0.0; // v=0 на стенке
                         v_im1 = v_old(i - 1, j);
                         d2vdx2 = (v_ip1 + v_im1 - 2.0 * v_old(i, j)) / dx2;
                     }
                 } else { // Выход (i = nx-1) - Нейман для v => dv/dx=0
                     v_ip1 = v_old(i, j);
                     v_im1 = v_old(i - 1, j);
                     d2vdx2 = (v_ip1 + v_im1 - 2.0 * v_old(i, j)) / dx2; // d2vdx2=0 ? Нет.
                     // Правильнее поток flux_e = 0
                     double nu_eff_w = 0.5 * (nu_eff_ij + get_effective_viscosity(i - 1, j));
                     double flux_w = -nu_eff_w * (v_old(i, j) - v_old(i - 1, j)) / dx;
                     d2vdx2 = (0.0 - flux_w) / dx;
                 }
             } else if (i - 1 <= 0 || tags(i - 1, j) == CellTag::SOLID) { // SOLID слева или вход
                  if (tags(i - 1, j) == CellTag::SOLID) { // SOLID слева
                      if (turbulent && turb_model_ptr) {
                          std::pair<double, double> tau_w_pair = turb_model_ptr->get_wall_shear_stress(i, j, BoundarySide::LEFT);
                          double flux_w = tau_w_pair.second / rho_; // tau_wy
                          // Восточный поток
                          double nu_eff_e = 0.5 * (nu_eff_ij + get_effective_viscosity(i + 1, j)); // i+1 точно FLUID
                          double flux_e = -nu_eff_e * (v_old(i + 1, j) - v_old(i, j)) / dx;
                          d2vdx2 = (flux_e - flux_w) / dx;
                      } else { // Ламинарный SOLID
                          v_im1 = 0.0; // v=0 на стенке
                          v_ip1 = v_old(i + 1, j);
                          d2vdx2 = (v_ip1 + v_im1 - 2.0 * v_old(i, j)) / dx2;
                      }
                  } else { // Вход (i = 0) - v=0 (Дирихле)
                       v_im1 = 0.0;
                       v_ip1 = v_old(i + 1, j);
                       d2vdx2 = (v_ip1 + v_im1 - 2.0 * v_old(i, j)) / dx2; // Используем v=0 на входе
                  }
             } else { // Полностью внутренняя ячейка по X
                 v_ip1 = v_old(i + 1, j);
                 v_im1 = v_old(i - 1, j);
                 d2vdx2 = (v_ip1 + v_im1 - 2.0 * v_old(i, j)) / dx2;
             }

             // Y-производные (горизонтальные стенки трубы - v=0 всегда)
             double v_jp1 = (j + 1 >= ny - 1 || tags(i, j + 1) == CellTag::SOLID) ? -v_old(i, j) // v=0 ghost cell
                            : v_old(i, j + 1);
             double v_jm1 = (j - 1 <= 0    || tags(i, j - 1) == CellTag::SOLID) ? -v_old(i, j) // v=0 ghost cell
                            : v_old(i, j - 1);
             d2vdy2 = (v_jp1 + v_jm1 - 2.0 * v_old(i, j)) / dy2;

             // Собираем Лапласиан
             if (turbulent && turb_model_ptr && (tags(i+1,j)==CellTag::SOLID || tags(i-1,j)==CellTag::SOLID) ) {
                 // Если у верт. стенки в турб. режиме, d2vdx2 уже посчитан через поток/tau_w
                 lap_ij = d2vdx2 + nu_eff_ij * d2vdy2;
             } else {
                 // Иначе стандартный Лапласиан * nu_eff
                  lap_ij = nu_eff_ij * (d2vdx2 + d2vdy2);
             }


             // Обновление v_star_ явным Эйлером
             v_star_(i, j) = v_old(i, j) + dt * lap_ij;
         } // end for i
    } // end for j

    // Принудительное применение ГУ к v_star_ ПОСЛЕ диффузии
    apply_bc(); // Повторный вызов применит v=0 на стенках/SOLID и Нейман на выходе
}

void VelocityPressureSolver::calculatePoissonRHS_RhieChow(Field2D<double>& rhs, double dt)
{
    const std::size_t nx = geom_.mesh().nx();
    const std::size_t ny = geom_.mesh().ny();
    const double dx = geom_.mesh().dx();
    const double dy = geom_.mesh().dy();

#ifdef USE_OPENMP
    #pragma omp parallel for collapse(2)
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
    ConvergenceInfo p_info = poisson_.solve(p_, rhs_, geom_, PoissonMode::VelocityPressure, max_pressure_iter_, pressure_tol_);
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
    diffuse_u(dt);
    diffuse_v(dt);

    if (turbulence_model_) {
       turbulence_model_->solve_step(u_, v_, dt); // Передаем текущие u_, v_
    }

    project(dt);

    last_monitor_info_.actualDt = dt;
}

} // namespace cfd