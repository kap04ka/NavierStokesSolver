#pragma once
#include "core/Geometry.hpp"
#include "turbulence/base/TurbulenceModel.hpp"
#include <memory>

namespace cfd {
class Solver {
public:
    Solver(const Geometry& geom, double rho, double nu_molecular, TurbulenceModelType turb_type)
        :   geom_(geom),
            rho_(rho),
            nu_molecular_(nu_molecular),
            turbulence_model_(createTurbulenceModel(turb_type, geom, rho, nu_molecular))
            {}
    virtual ~Solver() = default;
    virtual void step(double dt) = 0;
    // Метод для получения эффективной вязкости
    [[nodiscard]] double get_effective_viscosity(std::size_t i, std::size_t j) const {
            const Field2D<double>* nu_t_field = turbulence_model_ ? turbulence_model_->nu_t() : nullptr;
            if (nu_t_field) {
                // Возвращаем сумму молекулярной и турбулентной (с ограничением >= 0)
                return nu_molecular_ + std::max(0.0, (*nu_t_field)(i, j));
            } else {
                // Возвращаем только молекулярную для ламинарного режима
                return nu_molecular_;
            }
        }
protected:
    const Geometry& geom_;
    double          rho_;
    double          nu_molecular_;
    
    std::unique_ptr<TurbulenceModel> turbulence_model_;
};
} // namespace cfd