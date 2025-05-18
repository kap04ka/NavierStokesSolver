#include "math/cuda/JacobiSolverCUDA.hpp"
#include "core/CellTag.hpp"
#include "math/PoissonSolver.hpp"

#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include <stdexcept>
#include <string>
#include <limits>
#include <iostream>

#ifndef CUDA_CHECK_THROW_JACOBI
#define CUDA_CHECK_THROW_JACOBI(call)                                                \
    do {                                                                             \
        cudaError_t err = call;                                                      \
        if (err != cudaSuccess) {                                                    \
            throw std::runtime_error(std::string("CUDA Error in Jacobi: ") +        \
                                     cudaGetErrorString(err) + " (" + __FILE__ +    \
                                     ":" + std::to_string(__LINE__) + ")");      \
        }                                                                            \
    } while (0)
#endif

namespace cfd {

// -------- helpers -------------------------------------------------------------
__device__ inline double atomicMaxDouble(double* address, double val) {
    unsigned long long int* addr_as_ull = reinterpret_cast<unsigned long long int*>(address);
    unsigned long long int old = *addr_as_ull, assumed;
    do {
        if (__longlong_as_double(old) >= val) break; 
        assumed = old;
        old = atomicCAS(addr_as_ull, assumed, __double_as_longlong(val));
    } while (assumed != old);
    return __longlong_as_double(old);
}

__device__ inline double get_neighbor_val_for_kernel(
    int cur_i, int cur_j,
    int nei_i, int nei_j,
    const double* phi, const CellTag* tags,
    int nx, PoissonMode mode)
{
    int idx_neighbor = nei_j * nx + nei_i;
    int idx_cur      = cur_j * nx + cur_i;

    if (mode == PoissonMode::VelocityPressure) {
        return (tags[idx_neighbor] == CellTag::SOLID) ? phi[idx_cur] : phi[idx_neighbor];
    }
    return phi[idx_neighbor];
}

__global__ void jacobi_update_kernel(
    const double* __restrict__ phi_k,
    double* __restrict__ phi_next,
    const double* __restrict__ rhs,
    const CellTag* __restrict__ tags,
    int nx, int ny,
    double dx2, double dy2, double coef,
    PoissonMode mode,
    double* residual_global)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    extern __shared__ double sdata[];   
    unsigned int tid = threadIdx.y * blockDim.x + threadIdx.x;
    sdata[tid] = 0.0;

    if (i < nx && j < ny) {
        int idx = j * nx + i;

        if (i == 0 || i == nx - 1 || j == 0 || j == ny - 1 || tags[idx] == CellTag::SOLID) {
            phi_next[idx] = phi_k[idx];
        } else {
            double val_ip1 = get_neighbor_val_for_kernel(i, j, i + 1, j, phi_k, tags, nx, mode);
            double val_im1 = get_neighbor_val_for_kernel(i, j, i - 1, j, phi_k, tags, nx, mode);
            double val_jp1 = get_neighbor_val_for_kernel(i, j, i, j + 1, phi_k, tags, nx, mode);
            double val_jm1 = get_neighbor_val_for_kernel(i, j, i, j - 1, phi_k, tags, nx, mode);

            double new_phi = coef * ((val_ip1 + val_im1) / dx2 +
                                      (val_jp1 + val_jm1) / dy2 - rhs[idx]);

            phi_next[idx] = new_phi;
            sdata[tid] = fabs(new_phi - phi_k[idx]);
        }
    }

    __syncthreads();
    for (unsigned int stride = blockDim.x * blockDim.y / 2; stride > 0; stride >>= 1) {
        if (tid < stride) {
            sdata[tid] = fmax(sdata[tid], sdata[tid + stride]);
        }
        __syncthreads();
    }

    if (tid == 0) {
        atomicMaxDouble(residual_global, sdata[0]);
    }
}

// -----------------------------------------------------------------------------
ConvergenceInfo JacobiSolverCUDA::solve(
    Field2D<double>&       phi_h,
    const Field2D<double>& rhs_h,
    const Geometry&        geom,
    PoissonMode            mode,
    unsigned               maxIter,
    double                 tol)
{
    const auto& mesh   = geom.mesh();
    const auto& tag_h  = geom.tags();

    int nx = static_cast<int>(mesh.nx());
    int ny = static_cast<int>(mesh.ny());
    if (nx <= 2 || ny <= 2) return {0, 0.0};

    double dx = mesh.dx();
    double dy = mesh.dy();
    if (dx < 1e-9 || dy < 1e-9) return {0, std::numeric_limits<double>::max()};

    double dx2  = dx * dx;
    double dy2  = dy * dy;
    double coef = 1.0 / (2.0 * (1.0 / dx2 + 1.0 / dy2));

    std::size_t N_total = static_cast<std::size_t>(nx) * ny;
    std::size_t field_bytes = N_total * sizeof(double);
    std::size_t tag_bytes   = N_total * sizeof(CellTag);

    double *d_phi_A = nullptr, *d_phi_B = nullptr, *d_rhs = nullptr;
    CellTag* d_tags = nullptr;
    double* d_residual = nullptr; 

    CUDA_CHECK_THROW_JACOBI(cudaMalloc(&d_phi_A, field_bytes));
    CUDA_CHECK_THROW_JACOBI(cudaMalloc(&d_phi_B, field_bytes));
    CUDA_CHECK_THROW_JACOBI(cudaMalloc(&d_rhs,   field_bytes));
    CUDA_CHECK_THROW_JACOBI(cudaMalloc(&d_tags,  tag_bytes));
    CUDA_CHECK_THROW_JACOBI(cudaMalloc(&d_residual, sizeof(double)));

    CUDA_CHECK_THROW_JACOBI(cudaMemcpy(d_phi_A, phi_h.raw_data_ptr(), field_bytes, cudaMemcpyHostToDevice));
    CUDA_CHECK_THROW_JACOBI(cudaMemcpy(d_rhs,   rhs_h.raw_data_ptr(), field_bytes, cudaMemcpyHostToDevice));
    CUDA_CHECK_THROW_JACOBI(cudaMemcpy(d_tags,  tag_h.raw_data_ptr(), tag_bytes,   cudaMemcpyHostToDevice));

    dim3 block(16, 16);
    dim3 grid((nx + block.x - 1) / block.x, (ny + block.y - 1) / block.y);
    size_t shmem = block.x * block.y * sizeof(double); 

    ConvergenceInfo info{0, std::numeric_limits<double>::max()};

    double* d_cur = d_phi_A;
    double* d_next = d_phi_B;

    for (unsigned it = 0; it < maxIter; ++it) {
        CUDA_CHECK_THROW_JACOBI(cudaMemset(d_residual, 0, sizeof(double)));

        jacobi_update_kernel<<<grid, block, shmem>>>(
            d_cur, d_next, d_rhs, d_tags, nx, ny,
            dx2, dy2, coef, mode, d_residual);
        CUDA_CHECK_THROW_JACOBI(cudaGetLastError());

        double h_residual;
        CUDA_CHECK_THROW_JACOBI(cudaMemcpy(&h_residual, d_residual, sizeof(double), cudaMemcpyDeviceToHost));

        info.iterations = it + 1;
        info.residual   = h_residual;

        if (h_residual < tol && it > 0) {
            std::swap(d_cur, d_next); 
            break;
        }

        std::swap(d_cur, d_next);
    }

    CUDA_CHECK_THROW_JACOBI(cudaMemcpy(phi_h.raw_data_ptr_for_write(), d_cur, field_bytes, cudaMemcpyDeviceToHost));

    cudaFree(d_phi_A);
    cudaFree(d_phi_B);
    cudaFree(d_rhs);
    cudaFree(d_tags);
    cudaFree(d_residual);

    return info;
}

} // namespace cfd
