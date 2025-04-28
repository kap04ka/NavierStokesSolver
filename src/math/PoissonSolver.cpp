#include "math/PoissonSolver.hpp"
#include <cmath>
#ifdef USE_OPENMP
 #include <omp.h>
#endif

namespace cfd {

//========================= Jacobi =========================================
void JacobiSolver::solve(Field2D<double>&       phi,
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
    for (unsigned it = 0; it < maxIter; ++it) {
        double err = 0.0;
#ifdef USE_OPENMP
        #pragma omp parallel for reduction(max:err)
#endif
        for (std::size_t j = 1; j < ny-1; ++j) {
            for (std::size_t i = 1; i < nx-1; ++i) {
                if (tag(i,j)==CellTag::SOLID) { pn(i,j)=0.0; continue; }
                auto nb=[&](std::size_t ii,std::size_t jj){return tag(ii,jj)==CellTag::SOLID?phi(i,j):phi(ii,jj);} ;
                pn(i,j)=coef*((nb(i+1,j)+nb(i-1,j))/dx2 + (nb(i,j+1)+nb(i,j-1))/dy2 - rhs(i,j));
                err = std::max(err, std::fabs(pn(i,j)-phi(i,j)));
            }
        }
        phi.swap(pn);
        if (err < tol) break;
    }
}

//========================= ω‑SOR ==========================================
void SORSolver::solve(Field2D<double>&       phi,
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

    for (unsigned it = 0; it < maxIter; ++it) {
        double err = 0.0;
#ifdef USE_OPENMP
        #pragma omp parallel for reduction(max:err)
#endif
        for (std::size_t j = 1; j < ny-1; ++j) {
            for (std::size_t i = 1; i < nx-1; ++i) {
                if (tag(i,j)==CellTag::SOLID) continue;
                auto nb=[&](std::size_t ii,std::size_t jj){return tag(ii,jj)==CellTag::SOLID?phi(i,j):phi(ii,jj);} ;
                double p_new = coef*((nb(i+1,j)+nb(i-1,j))/dx2 + (nb(i,j+1)+nb(i,j-1))/dy2 - rhs(i,j));
                double diff  = p_new - phi(i,j);
                phi(i,j) += omega_ * diff;
                err = std::max(err, std::fabs(diff));
            }
        }
        if (err < tol) break;
    }
}

//========================= Wrapper ========================================
void PoissonSolverDyn::solve(Field2D<double>&       phi,
                             const Field2D<double>& rhs,
                             const Geometry&        geom,
                             unsigned               maxIter,
                             double                 tol)
{
    if (type_ == PoissonType::Jacobi)
        jac_.solve(phi, rhs, geom, maxIter, tol);
    else
        sor_.solve(phi, rhs, geom, maxIter, tol);
}

} // namespace cfd