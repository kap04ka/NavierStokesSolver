#include "math/PoissonSolver.hpp"
#include <cmath>
#include <limits>
#ifdef USE_OPENMP
 #include <omp.h>
#endif

namespace cfd {

//========================= Jacobi =========================================
[[nodiscard]] ConvergenceInfo JacobiSolver::solve(
    Field2D<double>&       phi,
    const Field2D<double>& rhs,
    const Geometry&        geom,
    unsigned               maxIter,
    double                 tol)
{
    const auto& mesh = geom.mesh();
    const auto& tag  = geom.tags();
    std::size_t nx = mesh.nx(), ny = mesh.ny();
    double dx2 = mesh.dx()*mesh.dx();
    double dy2 = mesh.dy()*mesh.dy();
    const double coef = 1.0 /(2.0*(1.0/dx2 + 1.0/dy2));

    Field2D<double> pn(nx, ny, 0.0);
    ConvergenceInfo info;
    info.residual = std::numeric_limits<double>::max();

    for (unsigned it = 0; it < maxIter; ++it) {
        double current_max_err = 0.0;

#ifdef USE_OPENMP
        #pragma omp parallel for reduction(max:current_max_err)
#endif
        for (std::size_t j = 1; j < ny-1; ++j) {
            for (std::size_t i = 1; i < nx-1; ++i) {
                if (tag(i,j)==CellTag::SOLID) { pn(i,j)=0.0; continue; }
                auto nb=[&](std::size_t ii,std::size_t jj){return tag(ii,jj)==CellTag::SOLID?phi(i,j):phi(ii,jj);} ;
                pn(i,j)=coef*((nb(i+1,j)+nb(i-1,j))/dx2 
                            + (nb(i,j+1)+nb(i,j-1))/dy2 - rhs(i,j));
                            current_max_err = std::max(current_max_err, std::fabs(pn(i,j)-phi(i,j)));
            }
        }
        phi.swap(pn);
        info.residual = current_max_err;
        info.iterations = it + 1;
        if (info.residual < tol) break;
    }
    return info;
}

//========================= ω‑SOR ==========================================
[[nodiscard]] ConvergenceInfo SORSolver::solve(
    Field2D<double>&       phi,
    const Field2D<double>& rhs,
    const Geometry&        geom,
    unsigned               maxIter,
    double                 tol)
{
    const auto& mesh = geom.mesh();
    const auto& tag  = geom.tags();
    std::size_t nx = mesh.nx(), ny = mesh.ny();
    double dx2 = mesh.dx()*mesh.dx();
    double dy2 = mesh.dy()*mesh.dy();
    const double coef = 1.0 /(2.0*(1.0/dx2 + 1.0/dy2));

    ConvergenceInfo info;
    info.residual = std::numeric_limits<double>::max();

    for (unsigned it = 0; it < maxIter; ++it) {
        double current_max_err = 0.0;
        for (std::size_t j = 1; j < ny-1; ++j) {
            for (std::size_t i = 1; i < nx-1; ++i) {
                if (tag(i,j)==CellTag::SOLID) continue;
                auto nb=[&](std::size_t ii,std::size_t jj){return tag(ii,jj)==CellTag::SOLID?phi(i,j):phi(ii,jj);} ;
                double p_new = coef*((nb(i+1,j)+nb(i-1,j))/dx2 + (nb(i,j+1)+nb(i,j-1))/dy2 - rhs(i,j));
                double diff  = p_new - phi(i,j);
                phi(i,j) += omega_ * diff;
                current_max_err = std::max(current_max_err, std::fabs(diff));
            }
        }
        info.residual = current_max_err;
        info.iterations = it + 1;
        if (info.residual < tol) break;
    }
    return info;
}

//========================= Wrapper ========================================
[[nodiscard]] ConvergenceInfo PoissonSolverDyn::solve(
    Field2D<double>&       phi,
    const Field2D<double>& rhs,
    const Geometry&        geom,
    unsigned               maxIter,
    double                 tol)
{
    if (type_ == PoissonType::Jacobi)
        return jac_.solve(phi, rhs, geom, maxIter, tol);
    else
        return sor_.solve(phi, rhs, geom, maxIter, tol);
}

} // namespace cfd