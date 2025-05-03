#include "turbulence/base/TurbulenceModel.hpp"
#include "turbulence/k_epsilon/KEpsilonModel.hpp" // Подключим реализацию k-eps
#include <stdexcept>

namespace cfd {

// Реализация фабрики
std::unique_ptr<TurbulenceModel> createTurbulenceModel(
    TurbulenceModelType type,
    const Geometry& geom,
    double rho,
    double nu_molecular)
{
    switch (type) {
        case TurbulenceModelType::None:
            return nullptr; // Возвращаем пустой указатель для ламинарного режима
        case TurbulenceModelType::KEpsilon:
            return std::make_unique<KEpsilonModel>(geom, rho, nu_molecular);
        // Добавить другие модели позже
        default:
            throw std::runtime_error("Unknown turbulence model type requested.");
    }
}

} // namespace cfd