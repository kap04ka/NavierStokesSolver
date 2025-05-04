#pragma once
#include "turbulence/base/TurbulenceModel.hpp"
#include "core/Field2D.hpp"
#include <vector>
#include <string>
#include <memory> // Для unique_ptr

namespace cfd {

struct KEpsilonConstants {
    // Стандартные значения для Standard k-epsilon
    static constexpr double Cmu  = 0.09;
    static constexpr double C1e  = 1.44;
    static constexpr double C2e  = 1.92;
    static constexpr double sigk = 1.0;
    static constexpr double sige = 1.3;
    // Пристенные функции
    static constexpr double kappa = 0.41; // Постоянная Кармана
    static constexpr double E_log = 9.8;  // Константа для лог. закона (гладкие стенки)
    // Минимальные значения для стабильности
    static constexpr double k_min = 1e-9;
    static constexpr double epsilon_min = 1e-9;
    // Можно добавить другие константы, если модель будет расширяться
};


class KEpsilonModel : public TurbulenceModel {
    public:
    KEpsilonModel(const Geometry& geom, double rho, double nu_molecular,
                  double u_max_for_init,
                  double inlet_turb_intensity = 0.05,
                  double inlet_length_scale_factor = 0.07, 
                  const KEpsilonConstants& constants = KEpsilonConstants{});

    void solve_step(const Field2D<double>& u, const Field2D<double>& v, double dt) override;
    const Field2D<double>* nu_t() const override;
    void apply_bc(const Field2D<double>& u, const Field2D<double>& v) override; // Объявление было добавлено ранее
    std::pair<double, double> get_wall_shear_stress(std::size_t i, std::size_t j, BoundarySide side) const override; // Объявление было добавлено ранее
    const Field2D<double>* k() const { return k_.get(); }
    const Field2D<double>* epsilon() const { return epsilon_.get(); }

    // Добавить поля k_, epsilon_, nu_t_ и т.д.
private:

    void calculate_Pk(const Field2D<double>& u, const Field2D<double>& v);
    void solve_k_equation(const Field2D<double>& u, const Field2D<double>& v, double dt);
    void solve_epsilon_equation(const Field2D<double>& u, const Field2D<double>& v, double dt);
    void update_nu_t();
    void apply_inlet_outlet_bc(const Field2D<double>& u, const Field2D<double>& v);
    void apply_wall_functions(const Field2D<double>& u, const Field2D<double>& v);
    double calculate_u_tau(double Up, double yp) const;

    std::unique_ptr<Field2D<double>> k_;
    std::unique_ptr<Field2D<double>> epsilon_;
    std::unique_ptr<Field2D<double>> nu_t_;
    std::unique_ptr<Field2D<double>> Pk_;
    std::unique_ptr<Field2D<double>> k_old_;
    std::unique_ptr<Field2D<double>> epsilon_old_;
    std::unique_ptr<Field2D<double>> tau_wx_;
    std::unique_ptr<Field2D<double>> tau_wy_;

    const KEpsilonConstants constants_;

    const double inlet_turbulence_intensity_;
    const double inlet_length_scale_;
};

} // namespace cfd
