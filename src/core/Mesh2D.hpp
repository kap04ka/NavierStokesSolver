#pragma once
#include <cstddef>
#include <utility>
#include <string>
#include <stdexcept>

namespace cfd {

class Mesh2D {
public:
    Mesh2D() = default;
    Mesh2D(std::size_t nx, std::size_t ny, double Lx, double Ly)
        : nx_(nx), ny_(ny), Lx_(Lx), Ly_(Ly), dx_(0.0), dy_(0.0)
        {
            if (nx == 0 || ny == 0 || Lx <= 0.0 || Ly <= 0.0) {
                std::string error_msg = "Mesh dimensions must be positive: nx=" +
                                        std::to_string(nx) + ", ny=" + std::to_string(ny) +
                                        ", Lx=" + std::to_string(Lx) +
                                        ", Ly=" + std::to_string(Ly);
                throw std::invalid_argument(error_msg);
            }
            dx_ = Lx / nx; 
            dy_ = Ly / ny;
        }
    

    [[nodiscard]] std::size_t nx() const noexcept { return nx_; }
    [[nodiscard]] std::size_t ny() const noexcept { return ny_; }
    [[nodiscard]] double dx() const noexcept { return dx_; }
    [[nodiscard]] double dy() const noexcept { return dy_; }
    [[nodiscard]] double Lx() const noexcept { return Lx_; }
    [[nodiscard]] double Ly() const noexcept { return Ly_; }

    [[nodiscard]] std::pair<double,double> centre(std::size_t i, std::size_t j) const noexcept
    { return { (i + 0.5) * dx_, (j + 0.5) * dy_ }; }

private:
    std::size_t nx_{0}, ny_{0};
    double      Lx_{0.0}, Ly_{0.0};
    double      dx_{0.0}, dy_{0.0};
};

} // namespace cfd