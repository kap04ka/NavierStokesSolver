#include "math/sequential/JacobiSolverSeq.hpp"
#include <cmath>
#include <limits>
#include <algorithm>
#include <iostream>


namespace cfd {

ConvergenceInfo JacobiSolverSeq::solve(
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

    if (nx <= 2 || ny <= 2) return {0, 0.0}; // Нет внутренних точек для обновления

    double dx2 = mesh.dx() * mesh.dx();
    double dy2 = mesh.dy() * mesh.dy();
    if (dx2 < 1e-12 || dy2 < 1e-12) return {0, std::numeric_limits<double>::max()};
    const double coef = 1.0 / (2.0 * (1.0/dx2 + 1.0/dy2));

    Field2D<double> pn = phi;
    ConvergenceInfo info;
    info.residual = std::numeric_limits<double>::max();

    for (unsigned it = 0; it < maxIter; ++it) {
        double current_max_err = 0.0;
        
        for (std::size_t j = 1; j < ny - 1; ++j) {
            for (std::size_t i = 1; i < nx - 1; ++i) {
                if (tag(i,j) == CellTag::SOLID) continue; 

                auto get_neighbor_val = [&](std::size_t ni, std::size_t nj) {
                    if (mode == PoissonMode::VelocityPressure) {
                        return tag(ni,nj) == CellTag::SOLID ? phi(i,j) : phi(ni,nj);
                    } else { 
                        return phi(ni,nj);
                    }
                };
                
                pn(i,j) = coef * ( (get_neighbor_val(i+1,j) + get_neighbor_val(i-1,j))/dx2 +
                                 (get_neighbor_val(i,j+1) + get_neighbor_val(i,j-1))/dy2 - rhs(i,j) );

                current_max_err = std::max(current_max_err, std::abs(pn(i,j) - phi(i,j)));
            }
        }
        phi.swap(pn);                     
        info.residual = current_max_err;
        info.iterations = it + 1;
        if (info.residual < tol && it > 0) { // it > 0 чтобы не выйти на первой итерации, если поля случайно совпали
            break;
        }
    }
    return info;
}

} // namespace cfd