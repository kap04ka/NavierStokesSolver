#pragma once

#include "turbulence/base/TurbulenceModel.hpp"
#include "core/Field2D.hpp"

namespace cfd {

// Константы для стандартной модели Спаларта-Аллмареса
struct SA_Constants {
    // Основные константы
    static constexpr double Cb1 = 0.1355;
    static constexpr double Cb2 = 0.622;
    static constexpr double SigmaNu = 2.0 / 3.0;
    static constexpr double Kappa = 0.41;

    // Константы для функции разрушения
    static constexpr double Cw1 = Cb1 / (Kappa * Kappa) + (1.0 + Cb2) / SigmaNu;
    static constexpr double Cw2 = 0.3;
    static constexpr double Cw3 = 2.0;

    // Константа для функции демпфирования
    static constexpr double Cv1 = 7.1;
};


class SpalartAllmarasModel : public TurbulenceModel {
public:
    SpalartAllmarasModel(const Geometry& geom, double rho, double nu_molecular);

    void solve_step(const Field2D<double>& u, const Field2D<double>& v, double dt) override;

    [[nodiscard]] const Field2D<double>* nu_t() const override { return &nu_t_; }

    [[nodiscard]] const Field2D<double>* k() const override { return nullptr; }
    [[nodiscard]] const Field2D<double>* epsilon() const override { return nullptr; }
    void apply_bc(const Field2D<double>& u, const Field2D<double>& v) override;
    [[nodiscard]] std::pair<double, double> get_wall_shear_stress(std::size_t i, std::size_t j, BoundarySide side) const override {
        return {0.0, 0.0};
    }
    
private:
    void initialize_fields();
    void calculate_wall_distance(); // Вычисление расстояния до стенки

    // Основные поля
    Field2D<double> nu_tilde_;     // Рабочая переменная модели (ню с тильдой)
    Field2D<double> nu_tilde_old_; // Значение на предыдущем временном шаге
    Field2D<double> nu_t_;         // Итоговая турбулентная вязкость
    
    // Вспомогательные поля
    Field2D<double> wall_dist_;    // Расстояние до ближайшей стенки для каждой ячейки
    
    const Field2D<CellTag>& tags_; // Ссылка на теги геометрии
    SA_Constants constants_;       // Константы модели
    
};

} // namespace cfd