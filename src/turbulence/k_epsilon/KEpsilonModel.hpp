#pragma once
#include "turbulence/base/TurbulenceModel.hpp"
#include <iostream> // Для заглушки

namespace cfd {

class KEpsilonModel : public TurbulenceModel {
public:
    KEpsilonModel(const Geometry& geom, double rho, double nu_molecular);

    void solve_step(const Field2D<double>& u, const Field2D<double>& v, double dt) override;

    const Field2D<double>* nu_t() const override;

    // Добавить поля k_, epsilon_, nu_t_ и т.д.
private:
     std::unique_ptr<Field2D<double>> k_;
     std::unique_ptr<Field2D<double>> epsilon_;
     std::unique_ptr<Field2D<double>> nu_t_;
     // Поля для градиентов, Pk и т.д. по мере необходимости
};

} // namespace cfd
