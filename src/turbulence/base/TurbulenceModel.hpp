#pragma once
#include "core/Field2D.hpp"
#include "core/Geometry.hpp"
#include <memory> // Для std::unique_ptr

namespace cfd {

// Перечисление типов моделей (можно вынести в отдельный файл)
enum class TurbulenceModelType {
    None,     // Ламинарный режим
    KEpsilon  // k-epsilon
    // Добавить другие позже
};

enum class BoundarySide { BOTTOM, TOP, LEFT, RIGHT };

class TurbulenceModel {
public:
    // Конструктор может принимать общие параметры
    TurbulenceModel(const Geometry& geom, double rho, double nu_molecular)
        : geom_(geom), rho_(rho), nu_molecular_(nu_molecular) {}

    virtual ~TurbulenceModel() = default;

    // Основной метод для выполнения шага модели турбулентности
    // Принимает осредненные скорости и dt
    virtual void solve_step(const Field2D<double>& u, const Field2D<double>& v, double dt) = 0;

    // Методы для граничных условий (пока пустые)
    virtual void apply_bc(const Field2D<double>& u, const Field2D<double>& v) = 0;

    // Метод для получения поля турбулентной кин. вязкости nu_t
    // Возвращает ссылку на поле (или nullptr, если модель выключена/не инициализирована)
    [[nodiscard]] virtual const Field2D<double>* nu_t() const = 0;

    // --- Опциональные Виртуальные Методы (для доступа к состоянию) ---
    // Предоставляют реализацию по умолчанию

    // Доступ к полю k (может вернуть nullptr)
    [[nodiscard]] virtual const Field2D<double>* k() const { return nullptr; }

    // Доступ к полю epsilon (может вернуть nullptr)
    [[nodiscard]] virtual const Field2D<double>* epsilon() const { return nullptr; }

    // Получение касательного напряжения на стенке для fluid-ячейки (i,j)
    // у грани 'side' (относительно ячейки i,j)
    // Возвращает пару {tau_x, tau_y}
    [[nodiscard]] virtual std::pair<double, double> get_wall_shear_stress(std::size_t i, std::size_t j, BoundarySide side) const {
        // Реализация по умолчанию возвращает нули
        // Производный класс (KEpsilonModel) должен переопределить это
        (void)i; (void)j; (void)side; // Отметить как неиспользуемые
        return {0.0, 0.0};
    }

protected:
    const Geometry& geom_;
    double rho_;
    double nu_molecular_; // Молекулярная вязкость
};
} // namespace cfd