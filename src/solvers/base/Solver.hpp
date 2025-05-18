#pragma once
#include "core/Geometry.hpp"
#include "turbulence/base/TurbulenceModel.hpp"
#include "turbulence/k_epsilon/KEpsilonModel.hpp"
#include <memory>
#include <iostream>

namespace cfd {
class Solver {
public:
    Solver(const Geometry& geom, double rho, double nu_molecular,
        TurbulenceModelType turb_type,
        double u_max_inlet,
        double inlet_turb_intensity = 0.05,
        double inlet_length_scale_factor = 0.07)
        : geom_(geom),
          rho_(rho),
          nu_molecular_(nu_molecular),
          u_max_inlet_bc_(u_max_inlet),
          turbulence_model_(nullptr)
        {
            switch (turb_type) {
                case TurbulenceModelType::None:
                    std::cout << "Solver created with Turbulence Model: None (Laminar)" << std::endl;
                    break;
                case TurbulenceModelType::KEpsilon:
                    std::cout << "Solver created with Turbulence Model: k-epsilon" << std::endl;
                    // Создаем KEpsilonModel, передавая нужные параметры
                    turbulence_model_ = std::make_unique<KEpsilonModel>(
                        geom, rho, nu_molecular,
                        u_max_inlet_bc_,
                        inlet_turb_intensity,
                        inlet_length_scale_factor
                    );
                    break;
                default:
                    throw std::runtime_error("Unknown turbulence model type requested in Solver constructor.");
            } 
        }
    virtual ~Solver() = default;
    virtual void step(double dt) = 0;
    [[nodiscard]] double get_effective_viscosity(std::size_t i, std::size_t j) const {
            const Field2D<double>* nu_t_field = turbulence_model_ ? turbulence_model_->nu_t() : nullptr;
            if (nu_t_field) {
                return nu_molecular_ + std::max(0.0, (*nu_t_field)(i, j));
            } else {
                return nu_molecular_;
            }
        }
    [[nodiscard]] const TurbulenceModel* getTurbulenceModel() const { return turbulence_model_.get(); }
    [[nodiscard]] TurbulenceModel* getTurbulenceModel() { return turbulence_model_.get(); }

protected:
    const Geometry& geom_;
    double          rho_;
    double          nu_molecular_;
    double          u_max_inlet_bc_;
    std::unique_ptr<TurbulenceModel> turbulence_model_;
};
} // namespace cfd