#pragma once
#include "math/PoissonSolver.hpp"
#include "math/sequential/JacobiSolverSeq.hpp"
#include "math/sequential/SORSolverSeq.hpp"
#include "math/omp/JacobiSolverOMP.hpp"
#include "math/omp/SORSolverRedBlackOMP.hpp" // В будущем

namespace cfd {

// Режим параллелизации для PoissonSolverDyn
enum class ParallelizationMode { 
    Sequential, 
    OpenMP 
    // CUDA, // В будущем
    // MPI   // В будущем
};

class PoissonSolverDyn : public PoissonSolver {
public:
    explicit PoissonSolverDyn(PoissonType type, 
                              ParallelizationMode parallel_mode,
                              double omega_sor = 1.7);

    [[nodiscard]] ConvergenceInfo solve(
        Field2D<double>&       phi,
        const Field2D<double>& rhs,
        const Geometry&        geom,
        PoissonMode            mode, // Этот mode (VP/VS) передается дальше
        unsigned               maxIter = 400,
        double                 tol = 1e-5) override;
    
    // Метод для получения типа решателя (если нужен в main.cpp для вывода)
    [[nodiscard]] PoissonType get_solver_type_for_dyn() const { return type_; }


private:
    PoissonType type_;
    ParallelizationMode parallel_mode_;
    
    // Экземпляры всех возможных решателей
    JacobiSolverSeq jac_seq_solver_;
    JacobiSolverOMP jac_omp_solver_;
    SORSolverSeq    sor_seq_solver_; 
    SORSolverRedBlackOMP sor_omp_solver_; // В будущем
};

} // namespace cfd