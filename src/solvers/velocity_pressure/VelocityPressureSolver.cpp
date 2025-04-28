#include "solvers/velocity_pressure/VelocityPressureSolver.hpp"
#include <algorithm>
#include <cmath>
#ifdef USE_OPENMP
 #include <omp.h>
#endif

namespace cfd {

//------------------------------------------------------------------------
VelocityPressureSolver::VelocityPressureSolver(const Geometry& geom,
                                               double rho,
                                               double nu,
                                               PoissonType ptype,
                                               double cfl, 
                                               double omega)
    : Solver(geom, nu),
      u_(geom.mesh().nx(), geom.mesh().ny(), 0.0),
      v_(u_), p_(u_),
      u_star_(u_), v_star_(u_), rhs_(u_),
      tag_(geom.tags()),
      rho_(rho), cfl_(cfl), poisson_(ptype, omega) {}

//---------------------------------------------------------------- inlet --
void VelocityPressureSolver::set_inlet_parabola(double umax)
{
    const std::size_t ny = geom_.mesh().ny();
    const double H = geom_.mesh().dy() * (ny - 1);
    for (std::size_t j = 0; j < ny; ++j) {
        double y = j * geom_.mesh().dy();
        u_(0, j) = 4.0 * umax * y * (H - y) / (H * H);
        v_(0, j) = 0.0;
    }
}

//---------------------------------------------------------------- bc -----
void VelocityPressureSolver::apply_bc()
{
    const std::size_t nx = geom_.mesh().nx();
    const std::size_t ny = geom_.mesh().ny();

    for (std::size_t i = 0; i < nx; ++i) {
        u_(i, 0) = v_(i, 0) = 0.0;
        u_(i, ny - 1) = v_(i, ny - 1) = 0.0;
    }
    for (std::size_t j = 0; j < ny; ++j) {
        u_(nx - 1, j) = u_(nx - 2, j);
        v_(nx - 1, j) = v_(nx - 2, j);
        p_(nx - 1, j) = 0.0;
    }
    for (std::size_t j = 1; j < ny - 1; ++j)
        for (std::size_t i = 1; i < nx - 1; ++i) {
            if (tag_(i,j)==CellTag::SOLID) { u_(i,j)=v_(i,j)=0.0; continue; }
            if (tag_(i-1,j)==CellTag::SOLID||tag_(i+1,j)==CellTag::SOLID||
                tag_(i,j-1)==CellTag::SOLID||tag_(i,j+1)==CellTag::SOLID)
                u_(i,j)=v_(i,j)=0.0;
        }
}

//------------------------------ вспом. -----------------------------------
double VelocityPressureSolver::compute_cfl_dt(double safety) const
{
    double umax=0.0, vmax=0.0;
    for(std::size_t j=0;j<geom_.mesh().ny();++j)
        for(std::size_t i=0;i<geom_.mesh().nx();++i){
            umax = std::max(umax, std::fabs(u_(i,j)));
            vmax = std::max(vmax, std::fabs(v_(i,j)));
        }
    double dt_x = umax>0? geom_.mesh().dx()/umax : 1e9;
    double dt_y = vmax>0? geom_.mesh().dy()/vmax : 1e9;
    return safety * std::min(dt_x, dt_y);
}

//-------------------------------- adv / diff -----------------------------
void VelocityPressureSolver::advect(Field2D<double>& f,const Field2D<double>& u,const Field2D<double>& v,double dt)
{
    const double dx=geom_.mesh().dx(), dy=geom_.mesh().dy();
    const std::size_t nx=geom_.mesh().nx(), ny=geom_.mesh().ny();
#ifdef USE_OPENMP
    #pragma omp parallel for
#endif
    for(std::size_t j=1;j<ny-1;++j)
        for(std::size_t i=1;i<nx-1;++i){
            if(tag_(i,j)==CellTag::SOLID){ f(i,j)=0; continue; }
            double dfdx = (u(i,j)>0)? (f(i,j)-f(i-1,j))/dx : (f(i+1,j)-f(i,j))/dx;
            double dfdy = (v(i,j)>0)? (f(i,j)-f(i,j-1))/dy : (f(i,j+1)-f(i,j))/dy;
            f(i,j) -= dt*(u(i,j)*dfdx + v(i,j)*dfdy);
        }
}

void VelocityPressureSolver::diffuse(Field2D<double>& f,double dt)
{
    double alpha = nu_*dt;
    if(alpha==0.0) return;
    const std::size_t nx=geom_.mesh().nx(), ny=geom_.mesh().ny();
    const double dx2=geom_.mesh().dx()*geom_.mesh().dx();
    const double dy2=geom_.mesh().dy()*geom_.mesh().dy();
    for(int it=0;it<8;++it){
#ifdef USE_OPENMP
        #pragma omp parallel for
#endif
        for(std::size_t j=1;j<ny-1;++j)
            for(std::size_t i=1;i<nx-1;++i){
                if(tag_(i,j)==CellTag::SOLID) continue;
                double lap=(f(i+1,j)+f(i-1,j)-2*f(i,j))/dx2 + (f(i,j+1)+f(i,j-1)-2*f(i,j))/dy2;
                f(i,j)+=alpha*lap;
            }
    }
}

//-------------------------------- project -------------------------------
void VelocityPressureSolver::project(double dt)
{
    const double dx=geom_.mesh().dx(), dy=geom_.mesh().dy();
    const std::size_t nx=geom_.mesh().nx(), ny=geom_.mesh().ny();

    for(std::size_t j=1;j<ny-1;++j)
        for(std::size_t i=1;i<nx-1;++i){
            if(tag_(i,j)==CellTag::SOLID){ rhs_(i,j)=0; continue; }
            double div=(u_star_(i+1,j)-u_star_(i-1,j))/(2*dx) + (v_star_(i,j+1)-v_star_(i,j-1))/(2*dy);
            rhs_(i,j)=rho_*div/dt;
        }

    poisson_.solve(p_, rhs_, geom_);

    for(std::size_t j=1;j<ny-1;++j)
        for(std::size_t i=1;i<nx-1;++i){
            if(tag_(i,j)==CellTag::SOLID){ u_(i,j)=v_(i,j)=0; continue; }
            u_(i,j)=u_star_(i,j)-dt/rho_*(p_(i+1,j)-p_(i-1,j))/(2*dx);
            v_(i,j)=v_star_(i,j)-dt/rho_*(p_(i,j+1)-p_(i,j-1))/(2*dy);
        }
}

//-------------------------------- main step -----------------------------
void VelocityPressureSolver::step(double dt_user)
{
    double dt = std::min(dt_user, compute_cfl_dt(cfl_));

    apply_bc();
    u_star_ = u_; v_star_ = v_;
    advect(u_star_, u_, v_, dt);
    advect(v_star_, u_, v_, dt);
    diffuse(u_star_, dt);
    diffuse(v_star_, dt);
    project(dt);
}

} // namespace cfd