#pragma once
#include "core/Field2D.hpp"
#include "core/Geometry.hpp"
#include <memory> 

namespace cfd {

enum class TurbulenceModelType {
    None,     // Ламинарный режим
    KEpsilon,  // k-epsilon
    SpalartAllmaras // SA модель
};

enum class BoundarySide { BOTTOM, TOP, LEFT, RIGHT };

class TurbulenceModel {
public:
    TurbulenceModel(const Geometry& geom, double rho, double nu_molecular)
        : geom_(geom), rho_(rho), nu_molecular_(nu_molecular) {}

    virtual ~TurbulenceModel() = default;

    virtual void solve_step(const Field2D<double>& u, const Field2D<double>& v, double dt) = 0;

    virtual void apply_bc(const Field2D<double>& u, const Field2D<double>& v) = 0;

    [[nodiscard]] virtual const Field2D<double>* nu_t() const = 0;


    [[nodiscard]] virtual const Field2D<double>* k() const { return nullptr; }

    [[nodiscard]] virtual const Field2D<double>* epsilon() const { return nullptr; }

    [[nodiscard]] virtual std::pair<double, double> get_wall_shear_stress(std::size_t i, std::size_t j, BoundarySide side) const {
        (void)i; (void)j; (void)side; 
        return {0.0, 0.0};
    }

protected:
    const Geometry& geom_;
    double rho_;
    double nu_molecular_; 
};
} 