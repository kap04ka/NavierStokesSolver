#pragma once
#include "math/PoissonSolver.hpp"

namespace cfd {

class SORSolverSeq : public PoissonSolver {
public:
    explicit SORSolverSeq(double omega = 1.7) : omega_(omega) {}

    [[nodiscard]] ConvergenceInfo solve(
        Field2D<double>&       phi,
        const Field2D<double>& rhs,
        const Geometry&        geom,
        PoissonMode            mode,
        unsigned               maxIter = 400,
        double                 tol = 1e-5) override;
private:
    double omega_; // Коэффициент релаксации
};

} // namespace cfd