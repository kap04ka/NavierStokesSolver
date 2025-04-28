#pragma once
#include "core/Geometry.hpp"

namespace cfd {
class Solver {
public:
    Solver(const Geometry& geom, double nu)
        : geom_(geom), nu_(nu) {}
    virtual ~Solver() = default;

    virtual void step(double dt) = 0;
protected:
    const Geometry& geom_;
    double          nu_;
};
} // namespace cfd