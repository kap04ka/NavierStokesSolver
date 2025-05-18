#pragma once
#include <vector>
#include <algorithm>
#include <cstddef>
#include <stdexcept>

namespace cfd {

template<class T>
class Field2D {
public:
    Field2D() = default;
    Field2D(std::size_t nx, std::size_t ny, const T& init = T{})
        : nx_(nx), ny_(ny), data_(nx * ny, init) {}

    [[nodiscard]] std::size_t nx() const noexcept { return nx_; }
    [[nodiscard]] std::size_t ny() const noexcept { return ny_; }

    T& operator()(std::size_t i, std::size_t j)               { return data_[idx(i, j)]; }
    const T& operator()(std::size_t i, std::size_t j) const   { return data_[idx(i, j)]; }

    T* raw_data_ptr_for_write() { // Для неконстантного доступа (запись)
        if (data_.empty()) return nullptr; // Защита от пустого вектора
        return data_.data(); 
    }
    const T* raw_data_ptr() const { 
        if (data_.empty()) return nullptr;
        return data_.data(); 
    }

    void fill(const T& val) { std::fill(data_.begin(), data_.end(), val); }
    void swap(Field2D& other) noexcept
    {
        std::swap(nx_,   other.nx_);
        std::swap(ny_,   other.ny_);
        data_.swap(other.data_);
    }

private:
    std::size_t idx(std::size_t i, std::size_t j) const
    {
        if (i >= nx_ || j >= ny_) throw std::out_of_range("Field2D index");
        return j * nx_ + i;
    }

    std::size_t   nx_ {0}, ny_ {0};
    std::vector<T> data_;
};

} // namespace cfd