#include "Geometry.hpp"
#include <fstream>

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

void Geometry::export_tags_csv(const std::string& file) const
{
    std::ofstream f(file);
    for (std::size_t j = 0; j < tags_.ny(); ++j) {
        for (std::size_t i = 0; i < tags_.nx(); ++i)
            f << static_cast<int>(tags_(i, j)) << (i + 1 == tags_.nx() ? '\n' : ',');
    }
}

} // namespace cfd