#include "math/sequential/SORSolverSeq.hpp"
#include <cmath>
#include <limits>
#include <algorithm> // Для std::max

namespace cfd {

ConvergenceInfo SORSolverSeq::solve(
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

    for (unsigned it = 0; it < maxIter; ++it) {
        double current_max_err = 0.0;
        // В SOR обновляем phi на месте, используя уже обновленные значения соседей на текущей итерации 'it'
        for (std::size_t j = 1; j < ny - 1; ++j) { // Только внутренние ячейки
            for (std::size_t i = 1; i < nx - 1; ++i) {
                if (tag(i,j) == CellTag::SOLID) continue; 
                double old_phi_ij = phi(i,j); // Значение до SOR-обновления на этой итерации

                auto get_neighbor_val_sor = [&](std::size_t ni, std::size_t nj) {
                    if (mode == PoissonMode::VelocityPressure) {
                        return tag(ni,nj) == CellTag::SOLID ? old_phi_ij : phi(ni,nj);
                    } else { 
                        return phi(ni,nj);
                    }
                };
                   
                double phi_star_gauss_seidel = coef * ( (get_neighbor_val_sor(i+1,j) + get_neighbor_val_sor(i-1,j))/dx2 +
                                                            (get_neighbor_val_sor(i,j+1) + get_neighbor_val_sor(i,j-1))/dy2 - rhs(i,j) );
                phi(i,j) = old_phi_ij + omega_ * (phi_star_gauss_seidel - old_phi_ij);
                current_max_err = std::max(current_max_err, std::abs(phi(i,j) - old_phi_ij));
            }
        }
        info.residual = current_max_err;
        info.iterations = it + 1;
        if (info.residual < tol && it > 0) break;
    }
    return info;
}

} // namespace cfd