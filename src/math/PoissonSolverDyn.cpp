#include "math/PoissonSolverDyn.hpp"

namespace cfd {

PoissonSolverDyn::PoissonSolverDyn(PoissonType type, 
                                   ParallelizationMode parallel_mode, 
                                   double omega_sor) 
    : type_(type), 
    parallel_mode_(parallel_mode), 
    jac_seq_solver_(),      // Использует конструктор по умолчанию
    jac_omp_solver_(),      // Использует конструктор по умолчанию
    sor_seq_solver_(omega_sor), // Передаем omega для SOR
    sor_omp_solver_(omega_sor) // Если будет
{}

ConvergenceInfo PoissonSolverDyn::solve(
    Field2D<double>&       phi,
    const Field2D<double>& rhs,
    const Geometry&        geom,
    PoissonMode            mode,
    unsigned               maxIter,
    double                 tol)
{
    if (type_ == PoissonType::Jacobi) {
        if (parallel_mode_ == ParallelizationMode::OpenMP) {
            return jac_omp_solver_.solve(phi, rhs, geom, mode, maxIter, tol);
        } else if (parallel_mode_ == ParallelizationMode::CUDA) {
            return jac_cuda_solver_.solve(phi, rhs, geom, mode, maxIter, tol);
        } else { // Sequential
            return jac_seq_solver_.solve(phi, rhs, geom, mode, maxIter, tol);
        }
    } else {
        if (parallel_mode_ == ParallelizationMode::OpenMP) {
            return sor_omp_solver_.solve(phi, rhs, geom, mode, maxIter, tol);
        } else if (parallel_mode_ == ParallelizationMode::CUDA) {
            return sor_cuda_solver_.solve(phi, rhs, geom, mode, maxIter, tol);
        } else { // Sequential
             return sor_seq_solver_.solve(phi, rhs, geom, mode, maxIter, tol);
        }
    }
}

} // namespace cfd