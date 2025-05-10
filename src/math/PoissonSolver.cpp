#include "math/PoissonSolver.hpp"
#include <cmath>     // Для std::abs, std::fabs
#include <limits>    // Для std::numeric_limits
#include <vector>    // Используется Field2D
#include <algorithm> // Для std::max

#ifdef USE_OPENMP // Если используется OpenMP
 #include <omp.h>
#endif

namespace cfd {

//========================= JacobiSolver ===================================
ConvergenceInfo JacobiSolver::solve(
    Field2D<double>&       phi,      // Обновляемое поле (p или psi)
    const Field2D<double>& rhs,      // Правая часть
    const Geometry&        geom,     // Геометрия (сетка + теги)
    PoissonMode            mode,     // Режим работы
    unsigned               maxIter,  // Макс. итераций
    double                 tol)      // Точность
{
    const auto& mesh = geom.mesh();
    const auto& tag  = geom.tags();
    std::size_t nx = mesh.nx();
    std::size_t ny = mesh.ny();

    // Уравнение Пуассона решается только для внутренних точек
    // Граничные ячейки (0 и N-1) предполагаются заданными как Дирихле (снаружи)
    // и не должны изменяться этим итерационным процессом.
    if (nx <= 2 || ny <= 2) { // Нужно хотя бы 3 точки по каждому направлению для одной внутренней точки
        return {0, 0.0}; 
    }

    double dx2 = mesh.dx() * mesh.dx();
    double dy2 = mesh.dy() * mesh.dy();
    if (dx2 < 1e-12 || dy2 < 1e-12) { // Защита от деления на очень малые или нулевые шаги
        return {0, std::numeric_limits<double>::max()};
    }
    const double coef = 1.0 / (2.0 * (1.0/dx2 + 1.0/dy2));

    Field2D<double> pn = phi; // Копируем phi, включая ГУ Дирихле на внешних границах.
                              // pn будет полем на новой итерации (k+1). phi - на старой (k).
    ConvergenceInfo info;
    info.residual = std::numeric_limits<double>::max();

    for (unsigned it = 0; it < maxIter; ++it) {
        double current_max_err = 0.0;
        
        // Обновляем поле pn, используя значения из phi (с предыдущей итерации Якоби)
        // Циклы по внутренним точкам, где применяется 5-точечный шаблон.
        // Внешние границы (i=0, nx-1, j=0, ny-1) не обновляются здесь, их значения в pn остаются из phi.
        for (std::size_t j = 1; j < ny - 1; ++j) {
            for (std::size_t i = 1; i < nx - 1; ++i) {
                if (tag(i,j) == CellTag::SOLID) continue; 
                auto nb_vp = [&](std::size_t ni, std::size_t nj) {
                    return tag(ni,nj) == CellTag::SOLID ? phi(i,j) : phi(ni,nj);
                };
                pn(i,j) = coef * ( (nb_vp(i+1,j) + nb_vp(i-1,j))/dx2 +
                                 (nb_vp(i,j+1) + nb_vp(i,j-1))/dy2 - rhs(i,j) );

                current_max_err = std::max(current_max_err, std::abs(pn(i,j) - phi(i,j)));
            }
        }
        phi.swap(pn); // Обновляем phi для следующей итерации или как результат.
                      // pn теперь можно переиспользовать (или она будет перезаписана phi на след. итерации).
                      // Теперь phi содержит значения (k+1), а pn содержит старые (k)
                      // Нет, Якоби использует значения только с предыдущего слоя.
                      // phi - это (k), pn - это (k+1). После swap, phi становится (k+1) для след. итерации.
        
        info.residual = current_max_err;
        info.iterations = it + 1;
        if (info.residual < tol && it > 0) { // it > 0 чтобы не выйти на первой итерации, если поля случайно совпали
            break;
        }
    }
    return info;
}


//========================= SORSolver ======================================
ConvergenceInfo SORSolver::solve(
    Field2D<double>&       phi,      // Обновляемое поле (p или psi)
    const Field2D<double>& rhs,      // Правая часть
    const Geometry&        geom,     // Геометрия
    PoissonMode            mode,     // Режим работы
    unsigned               maxIter,  // Макс. итераций
    double                 tol)      // Точность
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
                double old_phi_ij = phi(i,j); // Значение до SOR-обновления на этой итерации
                if (tag(i,j) == CellTag::SOLID) continue; 

                auto nb_vp = [&](std::size_t ni, std::size_t nj) {
                    return tag(ni,nj)==CellTag::SOLID ? old_phi_ij : phi(ni,nj);
                };
                   
                double phi_star_gauss_seidel = coef * ( (nb_vp(i+1,j) + nb_vp(i-1,j))/dx2 +
                                                            (nb_vp(i,j+1) + nb_vp(i,j-1))/dy2 - rhs(i,j) );
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

// PoissonSolverDyn::solve остается без изменений, просто пробрасывает 'mode'
ConvergenceInfo PoissonSolverDyn::solve(
    Field2D<double>&       phi,
    const Field2D<double>& rhs,
    const Geometry&        geom,
    PoissonMode            mode,
    unsigned               maxIter,
    double                 tol)
{
    if (type_ == PoissonType::Jacobi)
        return jac_.solve(phi, rhs, geom, mode, maxIter, tol);
    else // SOR
        return sor_.solve(phi, rhs, geom, mode, maxIter, tol);
}

} // namespace cfd