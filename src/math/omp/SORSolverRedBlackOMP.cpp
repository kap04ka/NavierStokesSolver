#include "math/omp/SORSolverRedBlackOMP.hpp"
#include <cmath>
#include <limits>
#include <algorithm> // Для std::max

#ifdef USE_OPENMP
 #include <omp.h>
#endif

namespace cfd {

ConvergenceInfo SORSolverRedBlackOMP::solve(
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

    double dx2 = mesh.dx()*mesh.dx();
    double dy2 = mesh.dy()*mesh.dy();
    if (dx2 < 1e-12 || dy2 < 1e-12) return {0, std::numeric_limits<double>::max()};
    const double coef = 1.0 /(2.0*(1.0/dx2 + 1.0/dy2));

    ConvergenceInfo info;
    info.residual = std::numeric_limits<double>::max();

    Field2D<double> phi_at_iteration_start(nx,ny);

    for (unsigned it = 0; it < maxIter; ++it) {
        double current_iter_max_diff = 0.0;
        phi_at_iteration_start = phi; 

        // --- КРАСНЫЙ ПОЛУШАГ ---
        #ifdef USE_OPENMP
        #pragma omp parallel for collapse(2) schedule(static)
        #endif
        for (std::size_t j = 1; j < ny - 1; ++j) {
            for (std::size_t i = 1; i < nx - 1; ++i) {
                if ((i + j) % 2 != 0) continue; 

                if (tag(i,j) == CellTag::SOLID) continue; 

                double old_phi_ij = phi_at_iteration_start(i,j); 

                auto get_neighbor_val_rb = [&](std::size_t ni, std::size_t nj) {
                    if (mode == PoissonMode::VelocityPressure) {
                        return tag(ni,nj) == CellTag::SOLID ? old_phi_ij : phi_at_iteration_start(ni,nj);
                    } else { 
                        return phi_at_iteration_start(ni,nj);
                    }
                };
                
                double phi_star_val = coef * ( 
                    (get_neighbor_val_rb(i+1,j) + get_neighbor_val_rb(i-1,j))/dx2 +
                    (get_neighbor_val_rb(i,j+1) + get_neighbor_val_rb(i,j-1))/dy2 - rhs(i,j) 
                );
                
                phi(i,j) = old_phi_ij + omega_ * (phi_star_val - old_phi_ij);
            }
        }

        // --- ЧЕРНЫЙ ПОЛУШАГ ---
        #ifdef USE_OPENMP
        #pragma omp parallel for collapse(2) schedule(static)
        #endif
        for (std::size_t j = 1; j < ny - 1; ++j) {
            for (std::size_t i = 1; i < nx - 1; ++i) {
                if ((i + j) % 2 == 0) continue; 

                if (tag(i,j) == CellTag::SOLID) continue;

                double old_phi_ij_for_black_step = phi_at_iteration_start(i,j); 

                auto get_neighbor_val_rb = [&](std::size_t ni, std::size_t nj) {
                    if (mode == PoissonMode::VelocityPressure) {
                        return tag(ni,nj) == CellTag::SOLID ? old_phi_ij_for_black_step : phi(ni,nj);
                    } else { 
                        return phi(ni,nj); 
                    }
                };
                
                double phi_star_val = coef * ( 
                    (get_neighbor_val_rb(i+1,j) + get_neighbor_val_rb(i-1,j))/dx2 +
                    (get_neighbor_val_rb(i,j+1) + get_neighbor_val_rb(i,j-1))/dy2 - rhs(i,j) 
                );
                
                phi(i,j) = old_phi_ij_for_black_step + omega_ * (phi_star_val - old_phi_ij_for_black_step);
            }
        }

        // --- ВЫЧИСЛЕНИЕ ОШИБКИ ПОСЛЕ ОБОИХ ПОЛУШАГОВ ---
        current_iter_max_diff = 0.0;
        #ifdef USE_OPENMP
        #pragma omp parallel for collapse(2) reduction(max:current_iter_max_diff) schedule(static)
        #endif
        for (std::size_t j = 1; j < ny - 1; ++j) {
            for (std::size_t i = 1; i < nx - 1; ++i) {
                if (tag(i,j) == CellTag::SOLID) continue;
                current_iter_max_diff = std::max(current_iter_max_diff, std::abs(phi(i,j) - phi_at_iteration_start(i,j)));
            }
        }
        
        info.residual = current_iter_max_diff;
        info.iterations = it + 1;
        if (info.residual < tol && it > 0) break;
    }
    return info;
}

} // namespace cfd