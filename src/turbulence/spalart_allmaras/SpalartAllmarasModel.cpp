#include "turbulence/spalart_allmaras/SpalartAllmarasModel.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>
#include <limits>

namespace cfd {

SpalartAllmarasModel::SpalartAllmarasModel(const Geometry& geom, double rho, double nu_molecular)
    : TurbulenceModel(geom, rho, nu_molecular), // Вызываем конструктор базового класса
      nu_tilde_(geom.mesh().nx(), geom.mesh().ny(), 0.0),
      nu_tilde_old_(geom.mesh().nx(), geom.mesh().ny(), 0.0),
      nu_t_(geom.mesh().nx(), geom.mesh().ny(), 0.0),
      wall_dist_(geom.mesh().nx(), geom.mesh().ny(), 0.0),
      tags_(geom.tags())
{
    std::cout << "Spalart-Allmaras model: Calculating wall distances..." << std::endl;
    calculate_wall_distance();
    initialize_fields();
    std::cout << "Spalart-Allmaras model initialized." << std::endl;
}

void SpalartAllmarasModel::calculate_wall_distance() {
    const std::size_t nx = geom_.mesh().nx();
    const std::size_t ny = geom_.mesh().ny();
    const double dx = geom_.mesh().dx();
    const double dy = geom_.mesh().dy();

    std::vector<std::pair<std::size_t, std::size_t>> wall_cells;
    // Находим все стенки: SOLID ячейки и внешние границы j=0, j=ny-1
    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            if (tags_(i, j) == CellTag::SOLID || j == 0 || j == ny - 1 || i == 0 || i == nx-1) { // Вход/выход тоже считаем стенками для d
                wall_cells.push_back({i, j});
            }
        }
    }

    // Брут-форс вычисление расстояния до ближайшей стенки (может быть медленным для больших сеток)
    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            if (tags_(i, j) == CellTag::FLUID) {
                auto [x_fluid, y_fluid] = geom_.mesh().centre(i,j);
                double min_dist_sq = std::numeric_limits<double>::max();
                for (const auto& wall_cell_idx : wall_cells) {
                    auto [x_wall, y_wall] = geom_.mesh().centre(wall_cell_idx.first, wall_cell_idx.second);
                    double dist_sq = std::pow(x_fluid - x_wall, 2) + std::pow(y_fluid - y_wall, 2);
                    min_dist_sq = std::min(min_dist_sq, dist_sq);
                }
                wall_dist_(i, j) = std::sqrt(min_dist_sq);
            }
        }
    }
}

void SpalartAllmarasModel::initialize_fields() {
    // Начальное значение nu_tilde для внутренних ячеек - небольшое
    // (например, 0.1 * nu_molecular).
    nu_tilde_.fill(0.1 * this->nu_molecular_);
    Field2D<double> dummy_u(geom_.mesh().nx(), geom_.mesh().ny(), 0.0);
    Field2D<double> dummy_v(geom_.mesh().nx(), geom_.mesh().ny(), 0.0);
    apply_bc(dummy_u, dummy_v); 
}

void SpalartAllmarasModel::apply_bc(const Field2D<double>& u, const Field2D<double>& v) {
    (void)u; (void)v; 

    const std::size_t nx = geom_.mesh().nx();
    const std::size_t ny = geom_.mesh().ny();

    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            if (tags_(i, j) == CellTag::SOLID || j == 0 || j == ny - 1) {
                nu_tilde_(i, j) = 0.0;
            }
        }
    }

    for (std::size_t j = 1; j < ny - 1; ++j) {
        nu_tilde_(0, j) = 100.0 * this->nu_molecular_; 
    }

    for (std::size_t j = 1; j < ny - 1; ++j) {
        nu_tilde_(nx - 1, j) = nu_tilde_(nx - 2, j);
    }
}


void SpalartAllmarasModel::solve_step(const Field2D<double>& u, const Field2D<double>& v, double dt) {
    nu_tilde_old_ = nu_tilde_; // Сохраняем значение с предыдущего шага

    const std::size_t nx = geom_.mesh().nx();
    const std::size_t ny = geom_.mesh().ny();
    const double dx = geom_.mesh().dx();
    const double dy = geom_.mesh().dy();
    
    // Временные поля для хранения членов уравнения
    Field2D<double> production_term(nx, ny, 0.0);
    Field2D<double> destruction_term(nx, ny, 0.0);
    Field2D<double> diffusion_term(nx, ny, 0.0);
    Field2D<double> advection_term(nx, ny, 0.0);

    // --- ШАГ I-IV: Вычисление всех членов уравнения переноса для внутренних жидких ячеек ---
    for (std::size_t j = 1; j < ny - 1; ++j) {
        for (std::size_t i = 1; i < nx - 1; ++i) {
            if (tags_(i,j) != CellTag::FLUID) continue;

            // --- Вычисляем вспомогательные величины на основе nu_tilde_old_
            double nu_tilde_val = std::max(0.0, nu_tilde_old_(i,j)); // Защита от отрицательных значений
            double chi = nu_tilde_val / this->nu_molecular_;
            double chi3 = chi * chi * chi;
            double fv1 = chi3 / (chi3 + std::pow(constants_.Cv1, 3));

            // --- Вычисляем член генерации (Production) ---
            // 1. Модуль вихря S
            double dudy = (u(i,j+1) - u(i,j-1)) / (2.0*dy);
            double dvdx = (v(i+1,j) - v(i-1,j)) / (2.0*dx);
            double S = std::abs(dvdx - dudy);

            // 2. Модифицированный вихрь S_tilde
            double d_sq = wall_dist_(i,j) * wall_dist_(i,j);
            double fv2 = 1.0 - chi / (1.0 + chi * fv1);
            double S_tilde = S + (nu_tilde_val / (constants_.Kappa * constants_.Kappa * d_sq)) * fv2;
            
            production_term(i,j) = constants_.Cb1 * S_tilde * nu_tilde_val;

            // --- Вычисляем член разрушения (Destruction) ---
            double r = std::min(nu_tilde_val / (S_tilde * constants_.Kappa * constants_.Kappa * d_sq + 1e-10), 10.0);
            double g = r + constants_.Cw2 * (std::pow(r,6) - r);
            double g6 = std::pow(g,6);
            double Cw3_6 = std::pow(constants_.Cw3, 6);
            double fw = g * std::pow((1.0 + Cw3_6) / (g6 + Cw3_6), 1.0/6.0);

            destruction_term(i,j) = constants_.Cw1 * fw * std::pow(nu_tilde_val / wall_dist_(i,j), 2);

            // --- Вычисляем член диффузии (включая кросс-диффузию) ---
            // Используем потоковую форму, как мы обсуждали для nu_eff
            // Коэффициент диффузии D = (nu_mol + nu_tilde) / sigma
            double sigma_inv = 1.0 / constants_.SigmaNu;
            
            // Вязкости на гранях (nu_mol + nu_tilde_old)
            double diff_coef_ij = this->nu_molecular_ + nu_tilde_val;
            double diff_coef_ip1j = this->nu_molecular_ + nu_tilde_old_(i+1,j);
            double diff_coef_im1j = this->nu_molecular_ + nu_tilde_old_(i-1,j);
            double diff_coef_ijp1 = this->nu_molecular_ + nu_tilde_old_(i,j+1);
            double diff_coef_ijm1 = this->nu_molecular_ + nu_tilde_old_(i,j-1);

            double D_e = 0.5 * (diff_coef_ij + diff_coef_ip1j);
            double D_w = 0.5 * (diff_coef_ij + diff_coef_im1j);
            double D_n = 0.5 * (diff_coef_ij + diff_coef_ijp1);
            double D_s = 0.5 * (diff_coef_ij + diff_coef_ijm1);
            
            // Потоки
            double flux_e = D_e * (nu_tilde_old_(i+1,j) - nu_tilde_val) / dx;
            double flux_w = D_w * (nu_tilde_val - nu_tilde_old_(i-1,j)) / dx;
            double flux_n = D_n * (nu_tilde_old_(i,j+1) - nu_tilde_val) / dy;
            double flux_s = D_s * (nu_tilde_val - nu_tilde_old_(i,j-1)) / dy;

            double standard_diffusion = sigma_inv * ((flux_e - flux_w)/dx + (flux_n - flux_s)/dy);

            // Кросс-диффузионный член
            double d_nut_dx = (nu_tilde_old_(i+1,j) - nu_tilde_old_(i-1,j)) / (2.0*dx);
            double d_nut_dy = (nu_tilde_old_(i,j+1) - nu_tilde_old_(i,j-1)) / (2.0*dy);
            double cross_diffusion = constants_.Cb2 * sigma_inv * (d_nut_dx*d_nut_dx + d_nut_dy*d_nut_dy);
            
            diffusion_term(i,j) = standard_diffusion + cross_diffusion;

            // --- Вычисляем член адвекции ---
            double u_ij = u(i,j);
            double v_ij = v(i,j);
            double d_nut_dx_adv, d_nut_dy_adv;
            if (u_ij >= 0.0) { d_nut_dx_adv = (nu_tilde_val - nu_tilde_old_(i-1,j)) / dx; }
            else { d_nut_dx_adv = (nu_tilde_old_(i+1,j) - nu_tilde_val) / dx; }
            if (v_ij >= 0.0) { d_nut_dy_adv = (nu_tilde_val - nu_tilde_old_(i,j-1)) / dy; }
            else { d_nut_dy_adv = (nu_tilde_old_(i,j+1) - nu_tilde_val) / dy; }
            advection_term(i,j) = u_ij * d_nut_dx_adv + v_ij * d_nut_dy_adv;
        }
    }

    // --- ШАГ V: Обновляем поле nu_tilde_ явным Эйлером ---
    for (std::size_t j = 1; j < ny - 1; ++j) {
        for (std::size_t i = 1; i < nx - 1; ++i) {
            if (tags_(i,j) == CellTag::FLUID) {
                nu_tilde_(i,j) = nu_tilde_old_(i,j) + dt * (
                    production_term(i,j) - 
                    destruction_term(i,j) + 
                    diffusion_term(i,j) - 
                    advection_term(i,j) // Адвекция с минусом, т.к. она слева в уравнении
                );
            }
        }
    }

    // --- ШАГ VI: Применяем ограничение положительности и ГУ ---
    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            nu_tilde_(i,j) = std::max(0.0, nu_tilde_(i,j));
        }
    }
    apply_bc(u, v);

    // --- ШАГ VII: Обновляем поле турбулентной вязкости nu_t_ ---
    for (std::size_t j = 0; j < ny; ++j) {
        for (std::size_t i = 0; i < nx; ++i) {
            double chi = nu_tilde_(i,j) / this->nu_molecular_;
            double chi3 = chi * chi * chi;
            double fv1 = chi3 / (chi3 + std::pow(constants_.Cv1, 3));
            nu_t_(i,j) = nu_tilde_(i,j) * fv1;
        }
    }
}

} // namespace cfd