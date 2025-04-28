#pragma once
#include <cstddef>
#include <utility>

namespace cfd {

class Mesh2D {
public:
    Mesh2D() = default;
    Mesh2D(std::size_t nx, std::size_t ny, double Lx, double Ly)
        : nx_(nx), ny_(ny), dx_(Lx / (nx - 1)), dy_(Ly / (ny - 1)) {}

    [[nodiscard]] std::size_t nx() const noexcept { return nx_; }
    [[nodiscard]] std::size_t ny() const noexcept { return ny_; }
    [[nodiscard]] double dx() const noexcept { return dx_; }
    [[nodiscard]] double dy() const noexcept { return dy_; }

    [[nodiscard]] std::pair<double,double> centre(std::size_t i, std::size_t j) const noexcept
    { return { (i + 0.5) * dx_, (j + 0.5) * dy_ }; }

private:
    std::size_t nx_{0}, ny_{0};
    double      dx_{0.0}, dy_{0.0};
};

} // namespace cfd