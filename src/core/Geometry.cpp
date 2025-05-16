#include "Geometry.hpp"
#include <fstream>
#include <cmath>     
#include <iostream>  

namespace cfd {

Geometry::Geometry(std::size_t nx, std::size_t ny, double Lx, double Ly)
    : mesh_(nx, ny, Lx, Ly), tags_(nx, ny, CellTag::FLUID) {}

void Geometry::add_rectangle(std::size_t i0, std::size_t j0,
                             std::size_t i1, std::size_t j1)
{
    for (std::size_t j = j0; j <= j1 && j < tags_.ny(); ++j)
        for (std::size_t i = i0; i <= i1 && i < tags_.nx(); ++i)
            tags_(i, j) = CellTag::SOLID;
}

void Geometry::add_circle(int center_i_idx, int center_j_idx, double radius_phys) {
    if (radius_phys <= 0.0) {
        std::cerr << "Warning (Geometry::add_circle): Radius must be positive. Circle not added." << std::endl;
        return;
    }

    // Проверка корректности индексов центра
    if (center_i_idx < 0 || center_i_idx >= static_cast<int>(mesh_.nx()) ||
        center_j_idx < 0 || center_j_idx >= static_cast<int>(mesh_.ny())) {
        std::cerr << "Warning (Geometry::add_circle): Center cell indices (i=" << center_i_idx 
                  << ", j=" << center_j_idx << ") are out of mesh bounds [" 
                  << mesh_.nx() << "," << mesh_.ny() << "]. Circle not added." << std::endl;
        return;
    }

    // Получаем физические координаты центра круга, используя центр указанной ячейки
    auto [phys_center_x, phys_center_y] = mesh_.centre(
                                                static_cast<std::size_t>(center_i_idx), 
                                                static_cast<std::size_t>(center_j_idx)
                                            );

    const double radius_sq = radius_phys * radius_phys;

    std::cout << "Geometry: Adding circle centered near cell (" << center_i_idx << "," << center_j_idx 
              << ") at physical coords (" << phys_center_x << "," << phys_center_y 
              << ") with R_phys=" << radius_phys << std::endl;

    for (std::size_t j = 0; j < tags_.ny(); ++j) {
        for (std::size_t i = 0; i < tags_.nx(); ++i) {
            auto [cell_cx, cell_cy] = mesh_.centre(i,j);
            double dx_half = mesh_.dx() / 2.0;
            double dy_half = mesh_.dy() / 2.0;

            // Координаты углов ячейки (i,j)
            double x_coords[] = {cell_cx - dx_half, cell_cx + dx_half, cell_cx + dx_half, cell_cx - dx_half};
            double y_coords[] = {cell_cy - dy_half, cell_cy - dy_half, cell_cy + dy_half, cell_cy + dy_half};

            bool cell_is_solid = false;
            // Проверка центра ячейки
            if (std::pow(cell_cx - phys_center_x, 2) + std::pow(cell_cy - phys_center_y, 2) <= radius_sq) {
                cell_is_solid = true;
            }
            // Проверка углов ячейки (если центр не попал, может попасть угол)
            if (!cell_is_solid) {
                for (int k=0; k<4; ++k) {
                    if (std::pow(x_coords[k] - phys_center_x, 2) + std::pow(y_coords[k] - phys_center_y, 2) <= radius_sq) {
                        cell_is_solid = true;
                        break;
                    }
                }
            }

            if (cell_is_solid) {
                tags_(i,j) = CellTag::SOLID;
            }
        }
    }
}


void Geometry::export_tags_csv(const std::string& file) const
{
    std::ofstream f(file);
    for (std::size_t j = 0; j < tags_.ny(); ++j) {
        for (std::size_t i = 0; i < tags_.nx(); ++i)
            f << static_cast<int>(tags_(i, j)) << (i + 1 == tags_.nx() ? '\n' : ',');
    }
}

} // namespace cfd