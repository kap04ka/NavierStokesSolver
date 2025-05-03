#include "turbulence/k_epsilon/KEpsilonModel.hpp"
#include <stdexcept> // Для runtime_error
#include <cmath>      // Для std::max, std::fabs (могут понадобиться позже)
#include <algorithm>  // Для std::max

namespace cfd {

// ... (Конструктор KEpsilonModel как был) ...
KEpsilonModel::KEpsilonModel(const Geometry& geom, double rho, double nu_molecular)
    : TurbulenceModel(geom, rho, nu_molecular)
{
    std::cout << "Initializing k-epsilon turbulence model framework..." << std::endl;
    k_       = std::make_unique<Field2D<double>>(geom.mesh().nx(), geom.mesh().ny(), 1e-6);
    epsilon_ = std::make_unique<Field2D<double>>(geom.mesh().nx(), geom.mesh().ny(), 1e-6);
    nu_t_    = std::make_unique<Field2D<double>>(geom.mesh().nx(), geom.mesh().ny(), 0.0);
}


void KEpsilonModel::solve_step(const Field2D<double>& /*u*/, const Field2D<double>& /*v*/, double /*dt*/)
{
    if (!k_ || !epsilon_ || !nu_t_) {
         throw std::runtime_error("k-epsilon fields not initialized!");
    }

    // --- ЗАГЛУШКА ---
    std::cout << "Warning: KEpsilonModel::solve_step() is not implemented yet!" << std::endl;
    // Здесь будет реальное решение уравнений k и epsilon

    // --- Обновление nu_t (примерно, пока k и eps не решаются) ---
    const double Cmu = 0.09;
    // Получаем ссылку на поле тегов для удобства
    const auto& cell_tags = geom_.tags(); 

    for(std::size_t j=0; j<geom_.mesh().ny(); ++j) {
        for(std::size_t i=0; i<geom_.mesh().nx(); ++i) {
            // Используем cell_tags вместо несуществующей tag_
            if(cell_tags(i,j) == CellTag::FLUID) {
                // Обеспечиваем положительность k и epsilon перед использованием
                // Используем -> для доступа к методам объекта через unique_ptr
                if (k_->operator()(i,j) < 1e-9) k_->operator()(i,j) = 1e-9;
                if (epsilon_->operator()(i,j) < 1e-9) epsilon_->operator()(i,j) = 1e-9;

                // Обновляем nu_t, используя operator() через unique_ptr
                nu_t_->operator()(i,j) = Cmu * (k_->operator()(i,j) * k_->operator()(i,j)) / epsilon_->operator()(i,j);

            } else { // Предполагаем, что остальное - SOLID
                k_->operator()(i,j) = 0.0;
                epsilon_->operator()(i,j) = 0.0; // Или другое ГУ для эпсилон в стене? Зависит от модели.
                nu_t_->operator()(i,j) = 0.0;
            }
        }
    }
    // --- Конец Заглушки ---
}

const Field2D<double>* KEpsilonModel::nu_t() const
{
    // Возвращаем сырой указатель из unique_ptr
    return nu_t_.get();
}

} // namespace cfd