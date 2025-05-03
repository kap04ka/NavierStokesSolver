#pragma once
#include "core/Field2D.hpp"
#include "core/Geometry.hpp"
#include <memory> // Для std::unique_ptr

namespace cfd {

class TurbulenceModel {
public:
    // Конструктор может принимать общие параметры
    TurbulenceModel(const Geometry& geom, double rho, double nu_molecular)
        : geom_(geom), rho_(rho), nu_molecular_(nu_molecular) {}

    virtual ~TurbulenceModel() = default;

    // Основной метод для выполнения шага модели турбулентности
    // Принимает осредненные скорости и dt
    virtual void solve_step(const Field2D<double>& u, const Field2D<double>& v, double dt) = 0;

    // Метод для получения поля турбулентной кин. вязкости nu_t
    // Возвращает ссылку на поле (или nullptr, если модель выключена/не инициализирована)
    virtual const Field2D<double>* nu_t() const = 0;

    // Методы для граничных условий (пока пустые)
    virtual void apply_bc() {} // Можно сделать чисто виртуальными, если нужно

    // Можно добавить методы для доступа к k, epsilon, если нужно для визуализации/отладки
    // virtual const Field2D<double>* k() const { return nullptr; }
    // virtual const Field2D<double>* epsilon() const { return nullptr; }

protected:
    const Geometry& geom_;
    double rho_;
    double nu_molecular_; // Молекулярная вязкость
};

// Перечисление типов моделей (можно вынести в отдельный файл)
enum class TurbulenceModelType {
    None,     // Ламинарный режим
    KEpsilon  // k-epsilon
    // Добавить другие позже
};

// Фабричная функция для создания модели (опционально, но удобно)
std::unique_ptr<TurbulenceModel> createTurbulenceModel(
    TurbulenceModelType type,
    const Geometry& geom,
    double rho,
    double nu_molecular);

} // namespace cfd