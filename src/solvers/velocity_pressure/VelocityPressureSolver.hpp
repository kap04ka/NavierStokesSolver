#pragma once
#include "solvers/base/Solver.hpp"
#include "core/Field2D.hpp"
#include "core/CellTag.hpp"
#include "math/PoissonSolver.hpp"

namespace cfd {

class VelocityPressureSolver : public Solver {
public:
    VelocityPressureSolver(const Geometry& geom,
                           double rho,
                           double nu,
                           PoissonType ptype = PoissonType::Jacobi,
                           double cfl = 0.4,
                           double omega = 1.7);

    /** шаг интегрирования; фактический dt ограничивается условием CFL */
    void step(double dt_user) override;

    void set_inlet_parabola(double umax);

    [[nodiscard]] const Field2D<double>& u() const noexcept { return u_; }
    [[nodiscard]] const Field2D<double>& v() const noexcept { return v_; }
    [[nodiscard]] const Field2D<double>& p() const noexcept { return p_; }

private:
    // вспомогательные процедуры
    void apply_bc();
    void advect  (Field2D<double>& f,const Field2D<double>& u,const Field2D<double>& v,double dt);
    void diffuse (Field2D<double>& f,double dt);
    //void rhie_chow_face_flux(double dt);
    //void calculate_face_fluxes();
    void project (double dt);
    double compute_cfl_dt(double safety) const;

    Field2D<double>      u_, v_, p_, u_star_, v_star_, rhs_, uf_, vf_;
    const Field2D<CellTag>& tag_;

    double rho_;   // плотность
    double cfl_;   // коэффициент безопасности CFL
    PoissonSolverDyn poisson_;
};

} // namespace cfd