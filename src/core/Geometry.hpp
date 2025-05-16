#pragma once
#include "Mesh2D.hpp"
#include "Field2D.hpp"
#include "CellTag.hpp"
#include <string>

namespace cfd {

class Geometry {
public:
    Geometry(std::size_t nx, std::size_t ny, double Lx, double Ly);

    void add_rectangle(std::size_t i0, std::size_t j0,
                       std::size_t i1, std::size_t j1);
    void add_circle(int center_i_idx, int center_j_idx, double radius_phys);
    const Mesh2D&           mesh() const noexcept { return mesh_; }
    const Field2D<CellTag>& tags() const noexcept { return tags_; }

    void export_tags_csv(const std::string& path) const;

private:
    Mesh2D              mesh_;
    Field2D<CellTag>    tags_;
};

} // namespace cfd