#include "math/cuda/SORSolverRedBlackCUDA.hpp"
#include "core/CellTag.hpp"
#include "math/PoissonSolver.hpp"

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <stdexcept>
#include <string>
#include <limits>
#include <iostream>

// -----------------------------------------------------------------------------
// Helper macro for CUDA error checking
#ifndef CUDA_CHECK_THROW_SOR
#define CUDA_CHECK_THROW_SOR(call)                                                   \
    do {                                                                             \
        cudaError_t err = call;                                                      \
        if (err != cudaSuccess) {                                                    \
            throw std::runtime_error(std::string("CUDA Error in SOR: ") +          \
                                     cudaGetErrorString(err) + " (" + __FILE__ +   \
                                     ":" + std::to_string(__LINE__) + ")");      \
        }                                                                            \
    } while (0)
#endif

namespace cfd {
__device__ inline double atomicMaxDouble(double* address, double val)
{
    unsigned long long* addr_ull = reinterpret_cast<unsigned long long*>(address);
    unsigned long long old = *addr_ull, assumed;
    do {
        if (__longlong_as_double(old) >= val) break;
        assumed = old;
        old = atomicCAS(addr_ull, assumed, __double_as_longlong(val));
    } while (assumed != old);
    return __longlong_as_double(old);
}

__device__ inline double get_neighbor_gpu(
    int ci, int cj, int ni, int nj,
    const double* phi, const CellTag* tags,
    int nx, PoissonMode mode)
{
    int nIdx = nj * nx + ni;
    int cIdx = cj * nx + ci;
    if (mode == PoissonMode::VelocityPressure) {
        return (tags[nIdx] == CellTag::SOLID) ? phi[cIdx] : phi[nIdx];
    }
    return phi[nIdx];
}

__global__ void sor_red_black_kernel(
    double* __restrict__ phi,
    const double* __restrict__ rhs,
    const CellTag* __restrict__ tags,
    int nx, int ny,
    double dx2, double dy2, double coef, double omega,
    PoissonMode mode,
    int parity,                 // 0 – red, 1 – black
    double* residual_g)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    extern __shared__ double sdata[];                 
    int tid = threadIdx.y * blockDim.x + threadIdx.x;
    sdata[tid] = 0.0;

    if (i >= 1 && i < nx - 1 && j >= 1 && j < ny - 1) {
        if (((i + j) & 1) == parity) {
            int idx = j * nx + i;
            if (tags[idx] != CellTag::SOLID) {
                double old_phi = phi[idx];
                double v_ip1 = get_neighbor_gpu(i, j, i + 1, j, phi, tags, nx, mode);
                double v_im1 = get_neighbor_gpu(i, j, i - 1, j, phi, tags, nx, mode);
                double v_jp1 = get_neighbor_gpu(i, j, i, j + 1, phi, tags, nx, mode);
                double v_jm1 = get_neighbor_gpu(i, j, i, j - 1, phi, tags, nx, mode);

                double phi_star = coef * ((v_ip1 + v_im1) / dx2 + (v_jp1 + v_jm1) / dy2 - rhs[idx]);
                double new_phi  = old_phi + omega * (phi_star - old_phi);
                phi[idx] = new_phi;
                sdata[tid] = fabs(new_phi - old_phi);
            }
        }
    }

    __syncthreads();
    for (int stride = blockDim.x * blockDim.y / 2; stride > 0; stride >>= 1) {
        if (tid < stride) sdata[tid] = fmax(sdata[tid], sdata[tid + stride]);
        __syncthreads();
    }

    if (tid == 0) atomicMaxDouble(residual_g, sdata[0]);
}

ConvergenceInfo SORSolverRedBlackCUDA::solve(
    Field2D<double>&       phi_h,
    const Field2D<double>& rhs_h,
    const Geometry&        geom,
    PoissonMode            mode,
    unsigned               maxIter,
    double                 tol)
{
    const auto& mesh  = geom.mesh();
    const auto& tag_h = geom.tags();

    int nx = static_cast<int>(mesh.nx());
    int ny = static_cast<int>(mesh.ny());
    if (nx <= 2 || ny <= 2) return {0, 0.0};

    double dx = mesh.dx();
    double dy = mesh.dy();
    if (dx < 1e-12 || dy < 1e-12) return {0, std::numeric_limits<double>::max()};

    double dx2 = dx * dx;
    double dy2 = dy * dy;
    double coef = 1.0 / (2.0 * (1.0 / dx2 + 1.0 / dy2));

    std::size_t N = static_cast<std::size_t>(nx) * ny;
    std::size_t bytes     = N * sizeof(double);
    std::size_t tag_bytes = N * sizeof(CellTag);

    // Allocate device arrays
    double *d_phi = nullptr, *d_rhs = nullptr, *d_res = nullptr;
    CellTag *d_tags = nullptr;

    CUDA_CHECK_THROW_SOR(cudaMalloc(&d_phi, bytes));
    CUDA_CHECK_THROW_SOR(cudaMalloc(&d_rhs, bytes));
    CUDA_CHECK_THROW_SOR(cudaMalloc(&d_tags, tag_bytes));
    CUDA_CHECK_THROW_SOR(cudaMalloc(&d_res, sizeof(double)));

    CUDA_CHECK_THROW_SOR(cudaMemcpy(d_phi,  phi_h.raw_data_ptr(), bytes,     cudaMemcpyHostToDevice));
    CUDA_CHECK_THROW_SOR(cudaMemcpy(d_rhs,  rhs_h.raw_data_ptr(), bytes,     cudaMemcpyHostToDevice));
    CUDA_CHECK_THROW_SOR(cudaMemcpy(d_tags, tag_h.raw_data_ptr(), tag_bytes, cudaMemcpyHostToDevice));

    dim3 block(16, 16);
    dim3 grid((nx + block.x - 1) / block.x, (ny + block.y - 1) / block.y);
    size_t shmem = block.x * block.y * sizeof(double);

    ConvergenceInfo info{0, std::numeric_limits<double>::max()};

    for (unsigned it = 0; it < maxIter; ++it) {
        CUDA_CHECK_THROW_SOR(cudaMemset(d_res, 0, sizeof(double)));

        sor_red_black_kernel<<<grid, block, shmem>>>(d_phi, d_rhs, d_tags, nx, ny,
                                                     dx2, dy2, coef, omega_, mode, 0, d_res);
        sor_red_black_kernel<<<grid, block, shmem>>>(d_phi, d_rhs, d_tags, nx, ny,
                                                     dx2, dy2, coef, omega_, mode, 1, d_res);
        CUDA_CHECK_THROW_SOR(cudaGetLastError());

        double h_residual;
        CUDA_CHECK_THROW_SOR(cudaMemcpy(&h_residual, d_res, sizeof(double), cudaMemcpyDeviceToHost));

        info.residual   = h_residual;
        info.iterations = it + 1;

        if (h_residual < tol && it > 0) break;
    }

    CUDA_CHECK_THROW_SOR(cudaMemcpy(phi_h.raw_data_ptr_for_write(), d_phi, bytes, cudaMemcpyDeviceToHost));

    cudaFree(d_phi);
    cudaFree(d_rhs);
    cudaFree(d_tags);
    cudaFree(d_res);

    return info;
}

} // namespace cfd
