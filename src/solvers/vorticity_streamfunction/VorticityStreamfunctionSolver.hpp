#pragma once

#include "solvers/base/Solver.hpp"
#include "core/Field2D.hpp"
#include "core/CellTag.hpp"
#include "math/PoissonSolverDyn.hpp" 
#include "turbulence/base/TurbulenceModel.hpp"
#include <vector> 
#include <utility> // Для std::pair

namespace cfd {

struct VS_SolverMonitorInfo {
    unsigned streamfunctionIterations = 0;
    double streamfunctionResidual = 0.0;
    double actualDt = 0.0;
};

// Структура для хранения информации об одном препятствии
struct ObstaclePsiInfo {
    std::vector<std::pair<std::size_t, std::size_t>> cells; // Координаты (i,j) ячеек препятствия
    double characteristic_y_coord_idx; // Характерный y-ИНДЕКС для этого препятствия
    double target_psi_value;           // Целевое значение psi для этого препятствия
};

class VorticityStreamfunctionSolver : public Solver {
public:
    VorticityStreamfunctionSolver(
        const Geometry& geom,
        double rho,
        double nu_molecular,
        TurbulenceModelType turb_type,
        double u_max_inlet,
        double inlet_turb_intensity,
        double inlet_length_scale_factor,
        PoissonType ptype,    // Тип решателя Пуассона (Якоби/SOR) для psi
        ParallelizationMode parallel_choice, // Если параллелизм
        double cfl_omega,     // Коэффициент CFL для уравнения вихря
        double omega_sor_psi, // Параметр SOR для решателя psi
        unsigned max_psi_iter,  // Макс. итераций для psi
        double psi_tol          // Точность для psi
    );

    void step(double dt_user) override;
    void initialize_fields();

    [[nodiscard]] const Field2D<double>& vorticity() const noexcept { return omega_; }
    [[nodiscard]] const Field2D<double>& streamfunction() const noexcept { return psi_; }
    [[nodiscard]] const Field2D<double>& u_velocity_from_psi() const noexcept { return u_from_psi_; }
    [[nodiscard]] const Field2D<double>& v_velocity_from_psi() const noexcept { return v_from_psi_; }
    [[nodiscard]] VS_SolverMonitorInfo getMonitorInfo() const noexcept { return last_monitor_info_; }

private:
    Field2D<double> omega_;        
    Field2D<double> psi_;          
    Field2D<double> omega_old_;    
    Field2D<double> u_from_psi_;
    Field2D<double> v_from_psi_;
    const Field2D<CellTag>& original_tags_; // Оригинальные теги из Geometry
    
    double cfl_omega_;             
    unsigned max_psi_iter_;      
    double psi_tol_;               
    double initial_Q_ = 0.0;       
    PoissonSolverDyn poisson_solver_psi_; 
    VS_SolverMonitorInfo last_monitor_info_;

    std::vector<ObstaclePsiInfo> internal_obstacles_psi_info_; // Хранилище для информации о psi на препятствиях
    bool use_fixed_psi_on_obstacles_ = true; // Флаг для выбора стратегии psi на препятствиях

    // Вспомогательные методы
    void apply_boundary_conditions_psi();
    void apply_boundary_conditions_omega();
    void update_velocities_from_psi();
    void advect_vorticity(double dt);
    void diffuse_vorticity(double dt);
    double compute_cfl_dt_omega(double safety_factor) const;
    double compute_diff_dt_omega() const;

    void identify_and_setup_internal_obstacles(); // Объединенный метод
};

} // namespace cfd