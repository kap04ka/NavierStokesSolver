#pragma once
#include "math/PoissonSolver.hpp" // Базовый интерфейс

namespace cfd {

class JacobiSolverSeq : public PoissonSolver {
public:
    JacobiSolverSeq() = default;

    [[nodiscard]] ConvergenceInfo solve(
        Field2D<double>&       phi,
        const Field2D<double>& rhs,
        const Geometry&        geom,
        PoissonMode            mode,
        unsigned               maxIter = 400,
        double                 tol = 1e-5) override;
};

} // namespace cfd