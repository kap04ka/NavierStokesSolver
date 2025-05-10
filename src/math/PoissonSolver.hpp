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

/* ---------------- Якоби ---------------- */
class JacobiSolver : public PoissonSolver {
public:
    [[nodiscard]] ConvergenceInfo solve(
        Field2D<double>&       phi,
        const Field2D<double>& rhs,
        const Geometry&        geom,
        PoissonMode            mode,
        unsigned               maxIter = 400,
        double                 tol     = 1e-5) override;
};

/* ---------------- ω‑SOR ---------------- */
class SORSolver : public PoissonSolver {
public:
    explicit SORSolver(double omega = 1.7) : omega_(omega) {}

    [[nodiscard]] ConvergenceInfo solve(
        Field2D<double>&       phi,
        const Field2D<double>& rhs,
        const Geometry&        geom,
        PoissonMode            mode,
        unsigned               maxIter = 400,
        double                 tol     = 1e-5) override;
private:
    double omega_;
};

/* --------- динамический выбор --------- */
class PoissonSolverDyn : public PoissonSolver {
public:
    explicit PoissonSolverDyn(PoissonType type, double omega) : type_(type), jac_(), sor_(omega) {}

    [[nodiscard]] ConvergenceInfo solve(
        Field2D<double>&       phi,
        const Field2D<double>& rhs,
        const Geometry&        geom,
        PoissonMode            mode,
        unsigned               maxIter = 400,
        double                 tol     = 1e-5) override;
private:
    PoissonType  type_;
    JacobiSolver jac_;
    SORSolver    sor_;
};

} // namespace cfd