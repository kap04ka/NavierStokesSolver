#include "core/Geometry.hpp"
#include "solvers/velocity_pressure/VelocityPressureSolver.hpp"
#include "math/PoissonSolver.hpp"
#include <iostream>
#include <cmath>

int main()
{
    constexpr std::size_t NX = 60, NY = 30;
    cfd::Geometry geom(NX, NY, 1.0, 0.5);
    geom.add_rectangle(20, 10, 40, 20);

    cfd::PoissonType ptype = cfd::PoissonType::Jacobi;
    double cfl = 0.1;

    cfd::VelocityPressureSolver solver(geom, 1000.0, 1e-3, ptype, cfl);
    solver.set_inlet_parabola(1.0);
    geom.export_tags_csv("mask_new.csv");


    const double dt = 0.001;
    for (int n = 0; n < 3000; ++n) {
        solver.step(dt);
        if (n % 1000 == 0) {
            std::cout<<"Pressure drop Δp = "<<solver.p()(1,NY/2)-solver.p()(NX-2,NY/2)<< std::endl;
            std::cout << "near obstacle u = " << solver.u()(NX / 2, NY - 5) << std::endl;
            std::cout << "near left u = " << solver.u()(5, NY / 2) << std::endl;
            std::cout << "near right u = " << solver.u()(NX - 5, NY / 2) << std::endl;
            std::cout << "centre u = " << solver.u()(NX / 2, NY / 2) << std::endl;
            std::cout << "right u = " << solver.u()(NX - 2, NY / 2) << std::endl;
            std::cout << "left u = " << solver.u()(0, NY / 2) << std::endl;
        }
    }

    double uc = solver.u()(NX / 2, NY / 2);
    double ur = solver.u()(NX - 2, NY / 2);
    double vr = solver.v()(NX - 2, NY / 2);
    double ul = solver.u()(0, NY / 2);
    double vl = solver.v()(0, NY / 2);
    std::cout << "\nRESULTS\n" << std::endl;
    std::cout << "centre u = " << uc << std::endl;
    std::cout << "right u = " << ur << std::endl;
    std::cout << "right v = " << vr << std::endl;
    std::cout << "left u = " << ul << std::endl;
    std::cout << "left v = " << vl << std::endl;
    return std::fabs(uc) < 1e-3 ? 0 : 1;
}