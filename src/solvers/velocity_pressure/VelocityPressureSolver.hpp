#pragma once
#include "solvers/base/Solver.hpp"
#include "core/Field2D.hpp"
#include "core/CellTag.hpp"
#include "math/PoissonSolverDyn.hpp" 

namespace cfd {

struct SolverMonitorInfo {
    unsigned pressureIterations = 0;
    double   pressureResidual = 0.0;
    double   actualDt = 0.0;
    // Сюда можно будет добавить другие параметры для мониторинга позже
};

class VelocityPressureSolver : public Solver {
public:
    VelocityPressureSolver(
        const Geometry& geom,
        double rho,
        double nu,
        TurbulenceModelType turb_type,
        double u_max_inlet,
        double inlet_turb_intensity = 0.05,
        double inlet_length_scale_factor = 0.07,
        PoissonType ptype = PoissonType::Jacobi,
        ParallelizationMode parallel_choice = ParallelizationMode::Sequential, 
        double cfl = 0.4,
        double omega = 1.7,
        unsigned max_p_iter = 400,
        double p_tol = 1e-5);

    /** шаг интегрирования; фактический dt ограничивается условием CFL */
    void step(double dt_user) override;

    void set_inlet_parabola(double umax);

    [[nodiscard]] const Field2D<double>& u() const noexcept { return u_; }
    [[nodiscard]] const Field2D<double>& v() const noexcept { return v_; }
    [[nodiscard]] const Field2D<double>& p() const noexcept { return p_; }

    // Метод для получения информации о сходимости 
    [[nodiscard]] SolverMonitorInfo getMonitorInfo() const noexcept { return last_monitor_info_; }

private:
    // вспомогательные процедуры
    void apply_bc();
    void advect  (Field2D<double>& f,const Field2D<double>& u,const Field2D<double>& v,double dt);
    void diffuse_u(double dt); // <<< Диффузия для u
    void diffuse_v(double dt); // <<< Диффузия для v
    void project (double dt);
    double compute_cfl_dt(double safety) const;
    double compute_diff_dt() const; // Добавляем для диффузионного лимита
    void calculatePoissonRHS_RhieChow(Field2D<double>& rhs, double dt);

    Field2D<double>      u_, v_, p_, u_star_, v_star_, rhs_;
    const Field2D<CellTag>& tag_;

    double      cfl_;   // коэффициент безопасности CFL
    unsigned    max_pressure_iter_; // Макс. итераций для Пуассона
    double      pressure_tol_;      // Точность для Пуассона

    PoissonSolverDyn poisson_;

    // Структура мониторинга
    SolverMonitorInfo last_monitor_info_;
};

} // namespace cfd