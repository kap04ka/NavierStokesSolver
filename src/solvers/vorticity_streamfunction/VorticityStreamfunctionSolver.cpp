#include "solvers/vorticity_streamfunction/VorticityStreamfunctionSolver.hpp"
#include "utils/export_utils.hpp"
#include <cmath>
#include <algorithm>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <limits>
#include <queue> // Для BFS при идентификации препятствий

namespace cfd {

VorticityStreamfunctionSolver::VorticityStreamfunctionSolver(
    const Geometry& geom,
    double rho,
    double nu_molecular,
    TurbulenceModelType turb_type,
    double u_max_inlet,
    double inlet_turb_intensity,
    double inlet_length_scale_factor,
    PoissonType ptype,
    double cfl_omega,
    double omega_sor_psi,
    unsigned max_psi_iter,
    double psi_tol)
    : Solver(geom, rho, nu_molecular, turb_type, u_max_inlet, inlet_turb_intensity, inlet_length_scale_factor),
      omega_(geom.mesh().nx(), geom.mesh().ny(), 0.0),
      psi_(geom.mesh().nx(), geom.mesh().ny(), 0.0),
      omega_old_(geom.mesh().nx(), geom.mesh().ny(), 0.0),
      u_from_psi_(geom.mesh().nx(), geom.mesh().ny(), 0.0),
      v_from_psi_(geom.mesh().nx(), geom.mesh().ny(), 0.0),
      original_tags_(geom.tags()), // Сохраняем оригинальные теги
      cfl_omega_(cfl_omega),
      max_psi_iter_(max_psi_iter),
      psi_tol_(psi_tol),
      poisson_solver_psi_(ptype, omega_sor_psi) {

    const std::size_t ny_mesh = geom_.mesh().ny();
    const double dy_mesh = geom_.mesh().dy();
    const double H_actual_channel = dy_mesh * (ny_mesh - 1);

    if (ny_mesh <= 1 || H_actual_channel < 1e-9) {
        initial_Q_ = 0.0;
        std::cerr << "Warning (VS_Solver): Invalid channel height for Q calculation. Q set to 0." << std::endl;
    } else {
        initial_Q_ = (2.0 / 3.0) * this->u_max_inlet_bc_ * H_actual_channel;
    }

    initialize_fields();

    std::cout << "Vorticity-Streamfunction Solver initialized." << std::endl;
    std::cout << "  Using fixed psi on obstacles: " << (use_fixed_psi_on_obstacles_ ? "Yes" : "No (Floating - requires modified Poisson)") << std::endl;
    std::cout << "  Channel height H (for Q calc): " << H_actual_channel << std::endl;
    std::cout << "  Umax for inlet profile: " << this->u_max_inlet_bc_ << std::endl;
    std::cout << "  Initial Q calculated: " << initial_Q_ << std::endl;
}

void VorticityStreamfunctionSolver::identify_and_setup_internal_obstacles() {
    internal_obstacles_psi_info_.clear();
    const std::size_t nx = geom_.mesh().nx();
    const std::size_t ny = geom_.mesh().ny();
    const double dy = geom_.mesh().dy();
    const double H_channel = dy * (ny - 1);

    Field2D<char> visited(nx, ny, false); 

    for (std::size_t j_start_obs = 1; j_start_obs < ny - 1; ++j_start_obs) {
        for (std::size_t i_start_obs = 1; i_start_obs < nx - 1; ++i_start_obs) {
            if (original_tags_(i_start_obs, j_start_obs) == CellTag::SOLID && !visited(i_start_obs, j_start_obs)) {
                ObstaclePsiInfo current_obstacle_info;
                std::queue<std::pair<std::size_t, std::size_t>> q_bfs;
                
                q_bfs.push({i_start_obs, j_start_obs});
                visited(i_start_obs, j_start_obs) = true;
                
                double sum_y_indices = 0;
                
                while(!q_bfs.empty()){
                    std::pair<std::size_t, std::size_t> cell = q_bfs.front();
                    q_bfs.pop();

                    current_obstacle_info.cells.push_back(cell);
                    sum_y_indices += static_cast<double>(cell.second);

                    int d_i_neighbors[] = {0, 0, 1, -1};
                    int d_j_neighbors[] = {1, -1, 0, 0};

                    for(int k=0; k<4; ++k){
                        std::size_t ni = cell.first + d_i_neighbors[k];
                        std::size_t nj = cell.second + d_j_neighbors[k];

                        if (ni > 0 && ni < nx - 1 && nj > 0 && nj < ny - 1 &&
                            original_tags_(ni,nj) == CellTag::SOLID && !visited(ni,nj)) 
                        {
                            visited(ni,nj) = true;
                            q_bfs.push({ni,nj});
                        }
                    }
                }

                if (!current_obstacle_info.cells.empty()) {
                    current_obstacle_info.characteristic_y_coord_idx = sum_y_indices / current_obstacle_info.cells.size();
                    
                    double y_char_obs_phys = current_obstacle_info.characteristic_y_coord_idx * dy;
                    if (H_channel > 1e-9) {
                        current_obstacle_info.target_psi_value = 4.0 * this->u_max_inlet_bc_ / (H_channel * H_channel) *
                                           (H_channel * y_char_obs_phys * y_char_obs_phys / 2.0 -
                                            y_char_obs_phys * y_char_obs_phys * y_char_obs_phys / 3.0);
                        current_obstacle_info.target_psi_value = std::max(0.0, std::min(initial_Q_, current_obstacle_info.target_psi_value));
                    } else {
                        current_obstacle_info.target_psi_value = initial_Q_ / 2.0; 
                    }
                    internal_obstacles_psi_info_.push_back(current_obstacle_info);
                     std::cout << "VS_Solver: Identified obstacle (approx y_idx_avg=" << current_obstacle_info.characteristic_y_coord_idx 
                               << "), target_psi=" << current_obstacle_info.target_psi_value << std::endl;
                }
            }
        }
    }
}


void VorticityStreamfunctionSolver::initialize_fields() {
    omega_.fill(0.0);
    psi_.fill(0.0); 

    identify_and_setup_internal_obstacles(); // Находим препятствия и вычисляем для них target_psi

    apply_boundary_conditions_psi();   
    update_velocities_from_psi();      
    apply_boundary_conditions_omega(); 
    
    omega_old_ = omega_;
}

void VorticityStreamfunctionSolver::apply_boundary_conditions_psi() {
    const std::size_t nx = geom_.mesh().nx();
    const std::size_t ny = geom_.mesh().ny();
    const double dy = geom_.mesh().dy();
    const double H_channel = dy * (ny - 1);

    for (std::size_t i = 0; i < nx; ++i) {
        psi_(i, 0) = 0.0;
        psi_(i, ny - 1) = initial_Q_;
    }

    if (H_channel > 1e-9 && ny > 1) {
        for (std::size_t j = 1; j < ny - 1; ++j) {
            const double y = j * dy;
            const double psi_val = 4.0 * u_max_inlet_bc_ / (H_channel * H_channel) *
                                   (H_channel * y * y / 2.0 - y * y * y / 3.0);
            psi_(0, j) = std::clamp(psi_val, 0.0, initial_Q_);
        }
    }

    if (nx >= 3) {
        for (std::size_t j = 1; j < ny - 1; ++j) { 
            psi_(nx - 1, j) = 2.0 * psi_(nx - 2, j) - psi_(nx - 3, j);
        }
    } else if (nx == 2) { 
         for (std::size_t j = 1; j < ny - 1; ++j) {
            psi_(nx - 1, j) = psi_(nx - 2, j);
         }
    }

    if (use_fixed_psi_on_obstacles_) {
        for (const auto& obs_info : internal_obstacles_psi_info_) {
            for (const auto& cell_coord : obs_info.cells) {
                psi_(cell_coord.first, cell_coord.second) = obs_info.target_psi_value;
            }
        }
    }
}

void VorticityStreamfunctionSolver::apply_boundary_conditions_omega() {
    const std::size_t nx = geom_.mesh().nx();
    const std::size_t ny = geom_.mesh().ny();
    const double dx = geom_.mesh().dx();
    const double dy = geom_.mesh().dy();
    const double dx2 = dx * dx;
    const double dy2 = dy * dy;
    const double H_channel = dy * (ny - 1);

    Field2D<double> omega_temp_for_relaxation = omega_old_; // Используем omega_old_ для релаксации

    // 1. Стенки области (j=0, j=ny-1)
    if (ny > 1) {
        for (std::size_t i = 0; i < nx; ++i) { 
            omega_(i, 0)     = -2.0 * (psi_(i, 1) - psi_(i, 0)) / dy2;
            omega_(i, ny - 1) = -2.0 * (psi_(i, ny - 2) - psi_(i, ny - 1)) / dy2;
        }
    }

    // 2. Внутренние SOLID препятствия (используя фиксированное target_psi_value для них)
    if (use_fixed_psi_on_obstacles_) {
        double alpha_omega_wall = 0.05; 
        for (const auto& obs_info : internal_obstacles_psi_info_) {
            for (const auto& solid_cell : obs_info.cells) {
                std::size_t i_solid = solid_cell.first;
                std::size_t j_solid = solid_cell.second;
                
                double sum_omega_thom_contrib = 0.0;
                int active_faces = 0;
                double psi_solid_val = obs_info.target_psi_value; // Используем фиксированное psi

                // Сосед сверху жидкий?
                if (j_solid + 1 < ny - 1 && original_tags_(i_solid, j_solid + 1) == CellTag::FLUID) {
                    sum_omega_thom_contrib += -2.0 * (psi_(i_solid, j_solid + 1) - psi_solid_val) / dy2;
                    active_faces++;
                }
                // Сосед снизу жидкий?
                if (j_solid > 0 && original_tags_(i_solid, j_solid - 1) == CellTag::FLUID) { // j_solid-1 может быть 0 (внешняя стенка)
                    sum_omega_thom_contrib += -2.0 * (psi_(i_solid, j_solid - 1) - psi_solid_val) / dy2;
                    active_faces++;
                }
                // Сосед слева жидкий?
                if (i_solid > 0 && original_tags_(i_solid - 1, j_solid) == CellTag::FLUID) {
                    sum_omega_thom_contrib += -2.0 * (psi_(i_solid - 1, j_solid) - psi_solid_val) / dx2;
                    active_faces++;
                }
                // Сосед справа жидкий?
                if (i_solid + 1 < nx - 1 && original_tags_(i_solid + 1, j_solid) == CellTag::FLUID) {
                    sum_omega_thom_contrib += -2.0 * (psi_(i_solid + 1, j_solid) - psi_solid_val) / dx2;
                    active_faces++;
                }

                if (active_faces > 0) {
                    double omega_thom_calculated = sum_omega_thom_contrib / static_cast<double>(active_faces);
                    omega_(i_solid, j_solid) = alpha_omega_wall * omega_thom_calculated +
                                               (1.0 - alpha_omega_wall) * omega_temp_for_relaxation(i_solid, j_solid);
                } else {
                     omega_(i_solid, j_solid) = omega_temp_for_relaxation(i_solid, j_solid); 
                }
            }
        }
    } // конец if (use_fixed_psi_on_obstacles_)

    if (H_channel > 1e-12 && ny > 1) {
        for (std::size_t j = 1; j < ny - 1; ++j) {
            const double y = j * dy;
            omega_(0, j) = -4.0 * u_max_inlet_bc_ / (H_channel * H_channel) * (H_channel - 2.0 * y);
        }
    }

    for (std::size_t j = 1; j < ny - 1; ++j) omega_(nx - 1, j) = omega_(nx - 2, j);
}

void VorticityStreamfunctionSolver::update_velocities_from_psi() {
    const std::size_t nx = geom_.mesh().nx();
    const std::size_t ny = geom_.mesh().ny();
    const double two_dx = 2.0 * geom_.mesh().dx();
    const double two_dy = 2.0 * geom_.mesh().dy();
    const double dy = geom_.mesh().dy(); 
    const double H_channel = dy * (ny - 1);

    if (nx < 2 || ny < 2) { 
        u_from_psi_.fill(0.0);
        v_from_psi_.fill(0.0);
        if (ny > 0) { for(std::size_t i_idx=0; i_idx<nx; ++i_idx) { u_from_psi_(i_idx,0)=0.0; v_from_psi_(i_idx,0)=0.0; if(ny>1){u_from_psi_(i_idx,ny-1)=0.0; v_from_psi_(i_idx,ny-1)=0.0;}}}
        if (nx > 0) { for(std::size_t j_idx=0; j_idx<ny; ++j_idx) { u_from_psi_(0,j_idx)=0.0; v_from_psi_(0,j_idx)=0.0; if(nx>1){u_from_psi_(nx-1,j_idx)=0.0; v_from_psi_(nx-1,j_idx)=0.0;}}}
        return;
    }

    for (std::size_t j = 0; j < ny; ++j) { // Идем по всем ячейкам для установки ГУ
        for (std::size_t i = 0; i < nx; ++i) {
            if (original_tags_(i,j) == CellTag::SOLID || i == 0 || i == nx-1 || j == 0 || j == ny-1) { // Если SOLID или внешняя граница
                // Сначала устанавливаем ГУ для скорости
                if (j == 0 || j == ny - 1) { // Нижняя и верхняя стенки
                    u_from_psi_(i,j) = 0.0;
                    v_from_psi_(i,j) = 0.0;
                } else if (i == 0) { // Вход (j от 1 до ny-2)
                    double y_coord_cell_center_j = j * dy;
                    u_from_psi_(0, j) = (H_channel > 1e-9) ? 
                                         4.0 * this->u_max_inlet_bc_ / (H_channel * H_channel) *
                                         y_coord_cell_center_j * (H_channel - y_coord_cell_center_j) : 0.0;
                    v_from_psi_(0,j) = 0.0;
                } else if (i == nx - 1 && nx >=2) { // Выход (j от 1 до ny-2)
                    u_from_psi_(nx-1,j) = u_from_psi_(nx-2,j); // Скопируем после основного цикла
                    v_from_psi_(nx-1,j) = v_from_psi_(nx-2,j); // Скопируем после основного цикла
                } else if (original_tags_(i,j) == CellTag::SOLID) { // Внутренние SOLID
                     u_from_psi_(i,j) = 0.0;
                     v_from_psi_(i,j) = 0.0;
                }
            } else { // Внутренние жидкие ячейки
                u_from_psi_(i, j) = (psi_(i, j + 1) - psi_(i, j - 1)) / two_dy;
                v_from_psi_(i, j) = -(psi_(i + 1, j) - psi_(i - 1, j)) / two_dx;
            }
        }
    }
    // Отдельный проход для выходной границы, чтобы использовать уже вычисленные u_from_psi_(nx-2,j)
    if (nx >= 2) {
        for (std::size_t j = 1; j < ny - 1; ++j) {
            u_from_psi_(nx-1,j) = u_from_psi_(nx-2,j);
            v_from_psi_(nx-1,j) = v_from_psi_(nx-2,j);
        }
    }
}


void VorticityStreamfunctionSolver::advect_vorticity(double dt) {
    const std::size_t nx = geom_.mesh().nx();
    const std::size_t ny = geom_.mesh().ny();
    const double dx = geom_.mesh().dx();
    const double dy = geom_.mesh().dy();
    
    Field2D<double> advection_term_values(nx, ny, 0.0);

    for (std::size_t j = 1; j < ny - 1; ++j) {
        for (std::size_t i = 1; i < nx - 1; ++i) {
            if (original_tags_(i,j) == CellTag::SOLID) {
                continue; 
            }

            double u_ij = u_from_psi_(i,j);
            double v_ij = v_from_psi_(i,j);
            double d_omega_dx, d_omega_dy;

            if (u_ij >= 0.0) { 
                d_omega_dx = (omega_old_(i,j) - omega_old_(i-1,j)) / dx;
            } else { 
                d_omega_dx = (omega_old_(i+1,j) - omega_old_(i,j)) / dx;
            }

            if (v_ij >= 0.0) { 
                d_omega_dy = (omega_old_(i,j) - omega_old_(i,j-1)) / dy;
            } else { 
                d_omega_dy = (omega_old_(i,j+1) - omega_old_(i,j)) / dy;
            }
            
            advection_term_values(i,j) = u_ij * d_omega_dx + v_ij * d_omega_dy;
        }
    }
    
    for (std::size_t j = 1; j < ny - 1; ++j) {
        for (std::size_t i = 1; i < nx - 1; ++i) {
             if (original_tags_(i,j) != CellTag::SOLID) {
                omega_(i,j) = omega_old_(i,j) - dt * advection_term_values(i,j);
             }
        }
    }
}

void VorticityStreamfunctionSolver::diffuse_vorticity(double dt) {
    const std::size_t nx = geom_.mesh().nx();
    const std::size_t ny = geom_.mesh().ny();
    const double dx = geom_.mesh().dx();
    const double dy = geom_.mesh().dy();
    // const auto& tags = original_tags_; // Используем original_tags_ для SOLID препятствий

    // Временное поле для хранения полного диффузионного члена D(omega^n)
    Field2D<double> diffusion_term_laplacian(nx, ny, 0.0);

    for (std::size_t j = 1; j < ny - 1; ++j) {
        for (std::size_t i = 1; i < nx - 1; ++i) {
            if (original_tags_(i, j) == CellTag::SOLID) { // Пропускаем внутренние SOLID препятствия
                continue;
            }

            // Эффективная вязкость в центре текущей ячейки (i,j) и соседних
            // nu_eff^n (т.к. используем omega_old_ для производных)
            double nu_eff_ij = get_effective_viscosity(i, j);
            double nu_eff_ip1j = get_effective_viscosity(i + 1, j);
            double nu_eff_im1j = get_effective_viscosity(i - 1, j);
            double nu_eff_ijp1 = get_effective_viscosity(i, j + 1);
            double nu_eff_ijm1 = get_effective_viscosity(i, j - 1);

            // Вязкости на гранях ячейки (i,j)
            double nu_eff_e_face = 0.5 * (nu_eff_ij + nu_eff_ip1j); // Восточная грань (i+1/2, j)
            double nu_eff_w_face = 0.5 * (nu_eff_ij + nu_eff_im1j); // Западная грань (i-1/2, j)
            double nu_eff_n_face = 0.5 * (nu_eff_ij + nu_eff_ijp1); // Северная грань (i, j+1/2)
            double nu_eff_s_face = 0.5 * (nu_eff_ij + nu_eff_ijm1); // Южная грань  (i, j-1/2)

            // Градиенты omega_old_ на гранях
            double domega_dx_e = (omega_old_(i + 1, j) - omega_old_(i, j)) / dx;
            double domega_dx_w = (omega_old_(i, j) - omega_old_(i - 1, j)) / dx;
            double domega_dy_n = (omega_old_(i, j + 1) - omega_old_(i, j)) / dy;
            double domega_dy_s = (omega_old_(i, j) - omega_old_(i, j - 1)) / dy;

            // Потоки вихря через грани
            double flux_omega_e = nu_eff_e_face * domega_dx_e;
            double flux_omega_w = nu_eff_w_face * domega_dx_w;
            double flux_omega_n = nu_eff_n_face * domega_dy_n;
            double flux_omega_s = nu_eff_s_face * domega_dy_s;

            // Диффузионный член для ячейки (i,j)
            diffusion_term_laplacian(i, j) = (flux_omega_e - flux_omega_w) / dx +
                                             (flux_omega_n - flux_omega_s) / dy;
        }
    }
    
    for (std::size_t j = 1; j < ny - 1; ++j) {
        for (std::size_t i = 1; i < nx - 1; ++i) {
             if (original_tags_(i,j) != CellTag::SOLID) { // Только для жидких ячеек
                omega_(i,j) += dt * diffusion_term_laplacian(i,j);
             }
        }
    }
}


double VorticityStreamfunctionSolver::compute_cfl_dt_omega(double safety_factor) const {
    double max_abs_u = 0.0;
    double max_abs_v = 0.0;
    const std::size_t nx = geom_.mesh().nx();
    const std::size_t ny = geom_.mesh().ny();

    for (std::size_t j = 0; j < ny; ++j) { 
        for (std::size_t i = 0; i < nx; ++i) {
            if (original_tags_(i,j) != CellTag::SOLID) { 
                max_abs_u = std::max(max_abs_u, std::abs(u_from_psi_(i,j)));
                max_abs_v = std::max(max_abs_v, std::abs(v_from_psi_(i,j)));
            }
        }
    }

    double dx = geom_.mesh().dx();
    double dy = geom_.mesh().dy();
    double dt_cfl = std::numeric_limits<double>::max();

    if (max_abs_u > 1e-12 && dx > 1e-12) {
        dt_cfl = std::min(dt_cfl, dx / max_abs_u);
    }
    if (max_abs_v > 1e-12 && dy > 1e-12) {
        dt_cfl = std::min(dt_cfl, dy / max_abs_v);
    }
    
    return safety_factor * dt_cfl;
}

double VorticityStreamfunctionSolver::compute_diff_dt_omega() const {
    double dx = geom_.mesh().dx();
    double dy = geom_.mesh().dy();
    if (dx < 1e-9 || dy < 1e-9) return std::numeric_limits<double>::max();

    double max_nu_eff = nu_molecular_; 

    if (this->turbulence_model_ && this->turbulence_model_->nu_t()) {
        const auto* nu_t_field_ptr = this->turbulence_model_->nu_t();
        if (nu_t_field_ptr) { 
            const auto& nu_t_field = *nu_t_field_ptr;
            double current_max_nu_t = 0.0;
            for(std::size_t j=0; j < geom_.mesh().ny(); ++j) {
                for(std::size_t i=0; i < geom_.mesh().nx(); ++i) {
                    if(original_tags_(i,j) == CellTag::FLUID) { 
                        current_max_nu_t = std::max(current_max_nu_t, nu_t_field(i,j));
                    }
                }
            }
            max_nu_eff += current_max_nu_t;
        }
    }
    max_nu_eff = std::max(max_nu_eff, 1e-12); 

    double denom = 2.0 * max_nu_eff * (1.0/(dx*dx) + 1.0/(dy*dy));
    if (denom < 1e-12) return std::numeric_limits<double>::max();
    
    return 1.0 / denom; 
}


void VorticityStreamfunctionSolver::step(double dt_user) {
    omega_old_ = omega_; // Сохраняем omega^n

    // ... (расчет dt) ...
    double dt_cfl  = compute_cfl_dt_omega(cfl_omega_); 
    double dt_diff = compute_diff_dt_omega();        
    double dt = std::min({dt_user, dt_cfl, dt_diff});
    if (dt < 1e-12) { dt = 1e-9; /* ... warning ... */ }
    last_monitor_info_.actualDt = dt;
    
    // --- Решаем уравнение переноса вихря ---
    omega_ = omega_old_; // Начинаем с omega^n для обновления до omega^{n+1}
    advect_vorticity(dt);  // omega_ теперь (omega^n - dt*Adv(omega^n))
    diffuse_vorticity(dt); // omega_ теперь (omega^n - dt*Adv(omega^n) + dt*Diff(omega^n))
                           
    // --- Устанавливаем ГУ для omega^{n+1} ---
    // psi_ здесь все еще psi^n. omega_ на выходе из этого метода будет omega^{n+1}
    apply_boundary_conditions_omega(); 
                                       
    Field2D<double> rhs_psi(geom_.mesh().nx(), geom_.mesh().ny());
    for(std::size_t j=0; j<geom_.mesh().ny(); ++j) {
        for(std::size_t i=0; i<geom_.mesh().nx(); ++i) {
            rhs_psi(i,j) = -omega_(i,j); // Используем omega^{n+1}
        }
    }
    
    apply_boundary_conditions_psi(); 

    ConvergenceInfo psi_conv_info;
    if (use_fixed_psi_on_obstacles_) {

        psi_conv_info = poisson_solver_psi_.solve(psi_, rhs_psi, geom_,
                                                PoissonMode::VorticityStreamfunction, max_psi_iter_, psi_tol_);
    } else {
        throw std::runtime_error("Floating psi on obstacles is not fully robustly implemented yet.");
    }

    last_monitor_info_.streamfunctionIterations = psi_conv_info.iterations;
    last_monitor_info_.streamfunctionResidual = psi_conv_info.residual;
 
    update_velocities_from_psi(); 

    if (this->turbulence_model_) { this->turbulence_model_->solve_step(u_from_psi_, v_from_psi_, dt); }

    static int s_vs_step_count_for_export = 0; // Статический счетчик только для этой функции
    s_vs_step_count_for_export++;

    // Выводим на шагах 1, 2, 3, 4, 5, и потом, например, каждый 10-й или 100-й
    bool condition_to_export = (s_vs_step_count_for_export <= 5) || 
                               (s_vs_step_count_for_export == 10) ||
                               (s_vs_step_count_for_export == 20) ||
                               (s_vs_step_count_for_export == 50) ||
                               (s_vs_step_count_for_export == 100) ||
                               (s_vs_step_count_for_export == 1000); 
                               // (s_vs_step_count_for_export % 100 == 0 && s_vs_step_count_for_export > 5); // для более редкого вывода позже

    if (condition_to_export) {
        std::cout << "Exporting VS fields at step: " << s_vs_step_count_for_export << " (SimTime approx: " << last_monitor_info_.actualDt * s_vs_step_count_for_export << ")" << std::endl;
        std::string step_str = std::to_string(s_vs_step_count_for_export);

        cfd::exportFieldToCSV(psi_, "vs_psi_step_" + step_str + ".csv");
        cfd::exportFieldToCSV(omega_, "vs_omega_step_" + step_str + ".csv");
        cfd::exportFieldToCSV(u_from_psi_, "vs_u_step_" + step_str + ".csv");
        cfd::exportFieldToCSV(v_from_psi_, "vs_v_step_" + step_str + ".csv");

        int nx = geom_.mesh().nx();
        int ny = geom_.mesh().ny();
        Field2D<double> eff_nu(nx, ny);

        for(size_t i = 0; i < nx; ++i)
            for(size_t j = 0; j < ny; ++j)
                eff_nu(i,j) = get_effective_viscosity(i,j);
        cfd::exportFieldToCSV(eff_nu, "eff_nu" + step_str + ".csv");

    }

}

} // namespace cfd