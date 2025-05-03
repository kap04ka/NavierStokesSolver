// #include "solvers/velocity_pressure/VelocityPressureSolver.hpp"
// #include <iostream>
// int main()
// {
//     constexpr std::size_t NX=60, NY=30;
//     cfd::Geometry geom(NX,NY,1.0,0.5);
//     // внутренний квадрат‑блок для проверки
//     geom.add_rectangle(25,10,34,19);
//     geom.export_tags_csv("mask_new.csv");

//     cfd::VelocityPressureSolver solver(geom,1000.0,1e-3);
//     solver.set_inlet_parabola(1.0);

//     double dt=0.001;
//     for(int n=0;n<3000;++n) solver.step(dt);

//     std::cout<<"Centreline u ≈ "<<solver.u()(NX - 3,NY/2)<<"\n";
//     std::cout<<"Pressure drop Δp = "<<solver.p()(1,NY/2)-solver.p()(NX-2,NY/2)<<"\n";
//     return 0;
// }