#include "math/omp/JacobiSolverOMP.hpp"
#include <cmath>
#include <limits>
#include <vector>
#include <algorithm>
#include <iostream>

#ifdef USE_OPENMP 
 #include <omp.h>
#endif

namespace cfd {

ConvergenceInfo JacobiSolverOMP::solve(
    Field2D<double>&       phi,
    const Field2D<double>& rhs,
    const Geometry&        geom,
    PoissonMode            mode,
    unsigned               maxIter,
    double                 tol)
{
    const auto& mesh = geom.mesh();
    const auto& tag  = geom.tags();
    std::size_t nx = mesh.nx();
    std::size_t ny = mesh.ny();

    if (nx <= 2 || ny <= 2) return {0, 0.0};

    double dx2 = mesh.dx() * mesh.dx();
    double dy2 = mesh.dy() * mesh.dy();
    if (dx2 < 1e-12 || dy2 < 1e-12) return {0, std::numeric_limits<double>::max()};
    const double coef = 1.0 / (2.0 * (1.0/dx2 + 1.0/dy2));

    Field2D<double> pn = phi;
    ConvergenceInfo info;
    info.residual = std::numeric_limits<double>::max();

    for (unsigned it = 0; it < maxIter; ++it) {
        double iter_global_max_err = 0.0; 
        #ifdef USE_OPENMP
        #pragma omp parallel for collapse(2) reduction(max:iter_global_max_err) schedule(static)
        #endif
        for (std::size_t j = 1; j < ny - 1; ++j) { 
            for (std::size_t i = 1; i < nx - 1; ++i) {
                if (tag(i,j) == CellTag::SOLID) continue; 

                double new_val_at_ij; 
                if (mode == PoissonMode::VelocityPressure) {
                    auto nb_vp = [&](std::size_t ni, std::size_t nj) {
                        // Читаем из phi (которое phi^k)
                        return tag(ni,nj)==CellTag::SOLID ? phi(i,j) : phi(ni,nj);
                    };
                    new_val_at_ij = coef * ( (nb_vp(i+1,j) + nb_vp(i-1,j))/dx2 +
                                             (nb_vp(i,j+1) + nb_vp(i,j-1))/dy2 - rhs(i,j) );
                } else { // PoissonMode::VorticityStreamfunction
                    auto nb_vs = [&](std::size_t ni, std::size_t nj) { 
                        // Читаем из phi (которое phi^k)
                        return phi(ni,nj); 
                    };
                    new_val_at_ij = coef * ( (nb_vs(i+1,j) + nb_vs(i-1,j))/dx2 +
                                             (nb_vs(i,j+1) + nb_vs(i,j-1))/dy2 - rhs(i,j) );
                }
                pn(i,j) = new_val_at_ij; 
                iter_global_max_err = std::max(iter_global_max_err, std::abs(pn(i,j) - phi(i,j)));
            }
        }
        
        phi.swap(pn); 
        info.residual = iter_global_max_err;
        info.iterations = it + 1;
        if (info.residual < tol && it > 0) { break; }
    }
    
    return info;
}

} // namespace cfd