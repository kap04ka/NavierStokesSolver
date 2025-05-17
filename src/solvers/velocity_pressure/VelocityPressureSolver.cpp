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
    ParallelizationMode parallel_choice, 
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
    poisson_(ptype, parallel_choice, omega) {}

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


// В src/solvers/velocity_pressure/VelocityPressureSolver.cpp

void VelocityPressureSolver::diffuse_u(double dt) {
    const auto nx = geom_.mesh().nx();
    const auto ny = geom_.mesh().ny();
    const double dx = geom_.mesh().dx();
    const double dy = geom_.mesh().dy();
    const double dx2 = dx * dx;
    const double dy2 = dy * dy;
    const auto& tags = geom_.tags();

    Field2D<double> u_old = u_star_;

    bool turbulent = (this->turbulence_model_ != nullptr);
    TurbulenceModel* turb_model_ptr = this->getTurbulenceModel();

    for (std::size_t j = 1; j < ny - 1; ++j) {
        for (std::size_t i = 1; i < nx - 1; ++i) {
            if (tags(i, j) == CellTag::SOLID) {
                u_star_(i, j) = 0.0;
                continue;
            }

            double nu_eff_ij = get_effective_viscosity(i, j);
            
            // Определяем, примыкает ли ячейка к ВНУТРЕННЕМУ SOLID препятствию СВЕРХУ или СНИЗУ
            // (не к внешним стенкам канала j=0 или j=ny-1)
            bool is_top_channel_boundary_for_logic    = (j == ny - 2);
            bool is_bottom_channel_boundary_for_logic = (j == 1);
            bool solid_obstacle_above = !is_top_channel_boundary_for_logic && (j + 1 < ny -1) && (tags(i,j+1) == CellTag::SOLID);
            bool solid_obstacle_below = !is_bottom_channel_boundary_for_logic && (j - 1 > 0) && (tags(i,j-1) == CellTag::SOLID);

            bool is_potentially_turbulent_wall_for_u = turbulent && turb_model_ptr && 
                                                       (solid_obstacle_above || solid_obstacle_below);

            if (nu_eff_ij < 1e-12 && !is_potentially_turbulent_wall_for_u) {
                 continue;
            }

            double lap_ij = 0.0; // <--- ОБЪЯВЛЕНИЕ lap_ij

            // --- X-компонента Лапласиана (d2udx2_simple) ---
            double u_ip1_val = (tags(i + 1, j) == CellTag::SOLID) ? -u_old(i, j) 
                                                                 : u_old(i + 1, j);
            double u_im1_val = (tags(i - 1, j) == CellTag::SOLID) ? -u_old(i, j) 
                                                                 : u_old(i - 1, j);
            if (i == 1) u_im1_val = u_old(0, j); 
            if (i == nx - 2) u_ip1_val = u_old(i, j); 
            double d2udx2_simple = (u_ip1_val + u_im1_val - 2.0 * u_old(i, j)) / dx2;

            // --- Y-компонента Лапласиана (d2udy2_final_term) ---
            double d2udy2_final_term; 
            bool y_term_includes_nu_eff = false; // <--- ОБЪЯВЛЕНИЕ y_term_includes_nu_eff
            
            if (is_top_channel_boundary_for_logic) { 
                double u_jp1_boundary = 0.0; 
                double u_jm1_neighbor = u_old(i, j - 1);
                d2udy2_final_term = (u_jp1_boundary + u_jm1_neighbor - 2.0 * u_old(i, j)) / dy2;
            } 
            else if (solid_obstacle_above) { 
                if (turbulent && turb_model_ptr) {
                    std::pair<double, double> tau_w_pair = turb_model_ptr->get_wall_shear_stress(i, j, BoundarySide::TOP);
                    double flux_n_wall = tau_w_pair.first / rho_;
                    double nu_eff_s_face = 0.5 * (nu_eff_ij + get_effective_viscosity(i, j - 1));
                    double flux_s_neighbor = nu_eff_s_face * (u_old(i, j) - u_old(i, j - 1)) / dy; 
                    d2udy2_final_term = (flux_n_wall - flux_s_neighbor) / dy; 
                    y_term_includes_nu_eff = true;
                } else { 
                    double u_jp1_boundary = 0.0; // Для ламинарного SOLID - как стенка
                    double u_jm1_neighbor = u_old(i, j - 1);
                    d2udy2_final_term = (u_jp1_boundary + u_jm1_neighbor - 2.0 * u_old(i, j)) / dy2;
                }
            }
            else if (is_bottom_channel_boundary_for_logic) { 
                double u_jp1_neighbor = u_old(i, j + 1);
                double u_jm1_boundary = 0.0; 
                d2udy2_final_term = (u_jp1_neighbor + u_jm1_boundary - 2.0 * u_old(i, j)) / dy2;
            }
            else if (solid_obstacle_below) { 
                if (turbulent && turb_model_ptr) {
                    std::pair<double, double> tau_w_pair = turb_model_ptr->get_wall_shear_stress(i, j, BoundarySide::BOTTOM);
                    double flux_s_wall = tau_w_pair.first / rho_;
                    double nu_eff_n_face = 0.5 * (nu_eff_ij + get_effective_viscosity(i, j + 1));
                    double flux_n_neighbor = nu_eff_n_face * (u_old(i, j + 1) - u_old(i, j)) / dy; 
                    d2udy2_final_term = (flux_n_neighbor - flux_s_wall) / dy;
                    y_term_includes_nu_eff = true;
                } else { 
                    double u_jp1_neighbor = u_old(i, j + 1);
                    double u_jm1_boundary = 0.0; 
                    d2udy2_final_term = (u_jp1_neighbor + u_jm1_boundary - 2.0 * u_old(i, j)) / dy2;
                }
            }
            else { 
                d2udy2_final_term = (u_old(i, j + 1) + u_old(i, j - 1) - 2.0 * u_old(i, j)) / dy2;
            }

            // Собираем Лапласиан
            if (y_term_includes_nu_eff) { 
                 lap_ij = nu_eff_ij * d2udx2_simple + d2udy2_final_term;
            } else { 
                 lap_ij = nu_eff_ij * (d2udx2_simple + d2udy2_final_term);
            }
            
            u_star_(i, j) = u_old(i, j) + dt * lap_ij;
        }
    } 
    apply_bc(); 
}

// В src/solvers/velocity_pressure/VelocityPressureSolver.cpp

void VelocityPressureSolver::diffuse_v(double dt) {
    const auto nx = geom_.mesh().nx();
    const auto ny = geom_.mesh().ny();
    const double dx = geom_.mesh().dx();
    const double dy = geom_.mesh().dy();
    const double dx2 = dx * dx;
    const double dy2 = dy * dy;
    const auto& tags = geom_.tags();

    Field2D<double> v_old = v_star_;

    bool turbulent = (this->turbulence_model_ != nullptr);
    TurbulenceModel* turb_model_ptr = this->getTurbulenceModel();

    for (std::size_t j = 1; j < ny - 1; ++j) {
        for (std::size_t i = 1; i < nx - 1; ++i) {
            if (tags(i, j) == CellTag::SOLID) {
                v_star_(i, j) = 0.0;
                continue;
            }

            double nu_eff_ij = get_effective_viscosity(i, j);
            if (nu_eff_ij < 1e-12) continue;

            double lap_ij = 0.0;

            // --- X-компонента Лапласиана (d2vdx2_term) ---
            double d2vdx2_term;
            bool x_term_is_flux_based = false;

            // Определяем, примыкает ли ячейка (i,j) к ВНУТРЕННЕМУ SOLID препятствию слева или справа
            // (не к границам области i=0 или i=nx-1)
            bool solid_right = (i + 1 < nx - 1) && (tags(i + 1, j) == CellTag::SOLID);
            bool solid_left  = (i - 1 > 0)      && (tags(i - 1, j) == CellTag::SOLID);

            if (solid_right) { // Внутреннее SOLID препятствие справа
                if (turbulent && turb_model_ptr) {
                    std::pair<double, double> tau_w_pair = turb_model_ptr->get_wall_shear_stress(i, j, BoundarySide::RIGHT);
                    double flux_e_wall = tau_w_pair.second / rho_; // Поток v через правую грань (стенка)
                    
                    // Поток через левую грань (w, i-1/2) - стандартный (сосед (i-1,j) не может быть SOLID, т.к. solid_left было бы true)
                    // Если i-1 == 0 (вход), то v_old(0,j) = 0.
                    double nu_eff_w_face = 0.5 * (nu_eff_ij + get_effective_viscosity(i - 1, j));
                    double flux_w_neighbor = nu_eff_w_face * (v_old(i, j) - v_old(i - 1, j)) / dx;
                    
                    d2vdx2_term = (flux_e_wall - flux_w_neighbor) / dx;
                    x_term_is_flux_based = true;
                } else { // Ламинарный SOLID справа
                    double v_ip1_eff = -v_old(i,j); // v=0 на SOLID
                    double v_im1_eff = v_old(i-1,j); // Сосед слева (может быть вход i-1=0, где v=0)
                    d2vdx2_term = (v_ip1_eff + v_im1_eff - 2.0 * v_old(i,j)) / dx2;
                }
            } else if (solid_left) { // Внутреннее SOLID препятствие слева
                if (turbulent && turb_model_ptr) {
                    std::pair<double, double> tau_w_pair = turb_model_ptr->get_wall_shear_stress(i, j, BoundarySide::LEFT);
                    double flux_w_wall = tau_w_pair.second / rho_; // Поток v через левую грань (стенка)

                    // Поток через правую грань (e, i+1/2) - стандартный (сосед (i+1,j) не может быть SOLID)
                    // Если i+1 == nx-1 (выход), то Нейман.
                    double v_ip1_eff_for_flux_e;
                    if (i + 1 == nx -1) { // Сосед справа - выходная граница
                        v_ip1_eff_for_flux_e = v_old(i,j); // Для Неймана dv/dx=0 => v(N)=v(N-1)
                    } else {
                        v_ip1_eff_for_flux_e = v_old(i+1,j);
                    }
                    double nu_eff_e_face = 0.5 * (nu_eff_ij + get_effective_viscosity(i + 1, j));
                    double flux_e_neighbor = nu_eff_e_face * (v_ip1_eff_for_flux_e - v_old(i, j)) / dx;

                    d2vdx2_term = (flux_e_neighbor - flux_w_wall) / dx;
                    x_term_is_flux_based = true;
                } else { // Ламинарный SOLID слева
                    double v_im1_eff = -v_old(i,j); // v=0 на SOLID
                    double v_ip1_eff = v_old(i+1,j); // Сосед справа (может быть выход i+1=nx-1)
                    if (i + 1 == nx - 1) v_ip1_eff = v_old(i,j); // Учет Неймана на выходе
                    d2vdx2_term = (v_ip1_eff + v_im1_eff - 2.0 * v_old(i,j)) / dx2;
                }
            } 
            // Границы области i=0 (вход) и i=nx-1 (выход), если нет SOLID препятствий рядом
            else if (i == 1) { // Ячейка (1,j) примыкает к входу (0,j)
                double v_im1_eff = v_old(0,j); // v=0 на входе (Дирихле)
                double v_ip1_eff = v_old(i+1,j); // Сосед (2,j)
                d2vdx2_term = (v_ip1_eff + v_im1_eff - 2.0 * v_old(i,j)) / dx2;
            } else if (i == nx - 2) { // Ячейка (nx-2,j) примыкает к выходу (nx-1,j)
                // Нейман dv/dx=0 => v(nx-1,j) = v(nx-2,j)
                double v_ip1_eff = v_old(i,j); // v(nx-1) = v(nx-2)
                double v_im1_eff = v_old(i-1,j); // Сосед (nx-3,j)
                d2vdx2_term = (v_ip1_eff + v_im1_eff - 2.0 * v_old(i,j)) / dx2;
            }
            // Полностью внутренняя жидкая ячейка по X
            else {
                d2vdx2_term = (v_old(i + 1, j) + v_old(i - 1, j) - 2.0 * v_old(i, j)) / dx2;
            }


            // --- Y-компонента Лапласиана (d2vdy2_term) ---
            // Для v-компоненты на горизонтальных стенках (внешних или SOLID) всегда v=0 (непротекание).
            // get_wall_shear_stress() НЕ используется.
            double v_jp1_eff, v_jm1_eff;

            // Сосед сверху (j+1)
            if (j + 1 == ny - 1 || tags(i, j + 1) == CellTag::SOLID) { // Верхняя стенка (граница или SOLID)
                v_jp1_eff = -v_old(i,j); // v=0 на стенке -> фиктивная ячейка
            } else {
                v_jp1_eff = v_old(i,j+1);
            }
            // Сосед снизу (j-1)
            if (j - 1 == 0 || tags(i, j - 1) == CellTag::SOLID) { // Нижняя стенка (граница или SOLID)
                v_jm1_eff = -v_old(i,j); // v=0 на стенке -> фиктивная ячейка
            } else {
                v_jm1_eff = v_old(i,j-1);
            }
            double d2vdy2_simple = (v_jp1_eff + v_jm1_eff - 2.0 * v_old(i,j)) / dy2;


            // Собираем Лапласиан
            if (x_term_is_flux_based) { // Если X-компонента была от турбулентной SOLID стенки
                 lap_ij = d2vdx2_term + nu_eff_ij * d2vdy2_simple; 
            } else { // Все остальное (ламинарный, или турбулентный без вертикальных SOLID стенок рядом)
                 lap_ij = nu_eff_ij * (d2vdx2_term + d2vdy2_simple);
            }
            
            v_star_(i, j) = v_old(i, j) + dt * lap_ij;
        } 
    } 
    apply_bc(); 
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