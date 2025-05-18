#pragma once

#include "math/PoissonSolver.hpp" // Базовый интерфейс

namespace cfd {

class JacobiSolverCUDA : public PoissonSolver {
public:
    JacobiSolverCUDA() = default;

    [[nodiscard]] ConvergenceInfo solve(
        Field2D<double>&       phi,      // Данные на CPU (Host)
        const Field2D<double>& rhs,      // Данные на CPU (Host)
        const Geometry&        geom,
        PoissonMode            mode,
        unsigned               maxIter = 400,
        double                 tol = 1e-5) override;
};

} // namespace cfd