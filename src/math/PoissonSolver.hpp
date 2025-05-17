#pragma once
#include "core/Field2D.hpp"
#include "core/Geometry.hpp"
#include <utility>

namespace cfd {

/// Доступные схемы решения Пуассона
enum class PoissonType { Jacobi, SOR };

/// Структура для возврата информации о сходимости
struct ConvergenceInfo {
    unsigned iterations = 0;
    double residual = 0;
};

enum class PoissonMode { 
    VelocityPressure, 
    VorticityStreamfunction 
};

/** Абстрактный базовый класс */
class PoissonSolver {
public:
    virtual ~PoissonSolver() = default;
    [[nodiscard]] virtual ConvergenceInfo solve(
        Field2D<double>&       phi,
        const Field2D<double>& rhs,
        const Geometry&        geom,
        PoissonMode            mode,
        unsigned               maxIter = 400,
        double                 tol     = 1e-5) = 0;
};
}
