/*
 * Copyright (C) 2026 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "ColorGradientKernelsCUDA.cuh"

#include <core/DataTypes.h>
#include <core/math/Vector3.h>
#include <gpu/FieldAccessor.h>
#include <gpu/FieldIndexing.h>
#include <gpu/GPUField.h>
#include <gpu/Kernel.h>

#include <utils/math/bspline.hpp>

#include <thrust/copy.h>
#include <thrust/device_vector.h>

#include <vector>

namespace walberla {
namespace lbm {
namespace color_gradient {

namespace {
// LCOV_EXCL_START
template <typename FloatType>
__global__ void kernel_count_negative_densities(
    gpu::FieldAccessor<FloatType> rho_a, gpu::FieldAccessor<FloatType> rho_b,
    unsigned int *count) {
  rho_a.set(blockIdx, threadIdx);
  rho_b.set(blockIdx, threadIdx);
  if (rho_a.isValidPosition() and rho_b.isValidPosition() and
      (rho_a.get() < FloatType{0} or rho_b.get() < FloatType{0})) {
    atomicAdd(count, 1u);
  }
}

template <typename FloatType>
__global__ void kernel_add_density_weighted_force(
    gpu::FieldAccessor<FloatType> force_a,
    gpu::FieldAccessor<FloatType> force_b,
    gpu::FieldAccessor<FloatType> rho_a, gpu::FieldAccessor<FloatType> rho_b,
    FloatType force_x, FloatType force_y, FloatType force_z) {
  force_a.set(blockIdx, threadIdx);
  force_b.set(blockIdx, threadIdx);
  rho_a.set(blockIdx, threadIdx);
  rho_b.set(blockIdx, threadIdx);
  if (force_a.isValidPosition() and force_b.isValidPosition() and
      rho_a.isValidPosition() and rho_b.isValidPosition()) {
    // same arithmetic as the CPU path in
    // LBWalberlaImplColorGradient::reset_force_and_apply_external_force()
    auto const density_a = rho_a.get();
    auto const density_b = rho_b.get();
    auto const total_density = density_a + density_b;
    if (total_density <= FloatType{0}) {
      return;
    }
    auto const inv_density = FloatType{1} / total_density;
    auto const fraction_a = density_a * inv_density;
    auto const fraction_b = density_b * inv_density;
    force_a.get(0u) += fraction_a * force_x;
    force_a.get(1u) += fraction_a * force_y;
    force_a.get(2u) += fraction_a * force_z;
    force_b.get(0u) += fraction_b * force_x;
    force_b.get(1u) += fraction_b * force_y;
    force_b.get(2u) += fraction_b * force_z;
  }
}

/**
 * @brief Interpolate both densities, and optionally their gradients, at
 * one position per thread with B-splines of order @p order. Same corner,
 * distance and weights as @c Utils::Interpolation::bspline_3d and
 * @c bspline_3d_gradient (grid spacing 1, offset 0.5), and the same summation
 * order. Stencil nodes outside the field and its ghost layers are skipped, as
 * on the CPU; this check is needed because the field accessors have none.
 */
template <int order, bool with_gradient>
__global__ void kernel_interpolate_densities(
    gpu::FieldAccessor<double> rho_a, gpu::FieldAccessor<double> rho_b,
    double const *__restrict__ pos, double *__restrict__ out, uint n_pos,
    int size_x, int size_y, int size_z, int gl) {
  constexpr uint n_values = with_gradient ? 8u : 2u;
  auto const index = blockIdx.x * blockDim.x + threadIdx.x;
  // the accessors point to the first ghost cell (-gl, -gl, -gl)
  rho_a.set({0u, 0u, 0u}, {0u, 0u, 0u});
  rho_b.set({0u, 0u, 0u}, {0u, 0u, 0u});
  if (index >= n_pos) {
    return;
  }
  int const size[3] = {size_x, size_y, size_z};
  int corner[3];
  double weight[3][order];
  double weight_d[3][order];
#pragma unroll
  for (int dim = 0; dim < 3; ++dim) {
    auto const fractional_index = pos[3u * index + dim] - 0.5;
    int nmp;
    double distance;
    if constexpr (order % 2 == 0) {
      nmp = static_cast<int>(floor(fractional_index));
      distance = fractional_index - nmp - 0.5;
    } else {
      nmp = static_cast<int>(floor(fractional_index + 0.5));
      distance = fractional_index - nmp;
    }
    corner[dim] = nmp - (order - 1) / 2;
#pragma unroll
    for (int i = 0; i < order; ++i) {
      weight[dim][i] = Utils::bspline<order>(i, distance);
      weight_d[dim][i] = Utils::bspline_d<order>(i, distance);
    }
  }
  auto const in_field = [&](int dim, int cell) {
    return cell >= -gl and cell < size[dim] + gl;
  };
  double density_a = 0., density_b = 0.;
  double gradient_a[3] = {0., 0., 0.};
  double gradient_b[3] = {0., 0., 0.};
  for (int i = 0; i < order; ++i) {
    auto const x = corner[0] + i;
    if (not in_field(0, x)) {
      continue;
    }
    for (int j = 0; j < order; ++j) {
      auto const y = corner[1] + j;
      if (not in_field(1, y)) {
        continue;
      }
      auto const weight_xy = weight[0][i] * weight[1][j];
      for (int k = 0; k < order; ++k) {
        auto const z = corner[2] + k;
        if (not in_field(2, z)) {
          continue;
        }
        auto const value_a = rho_a.getNeighbor(x + gl, y + gl, z + gl);
        auto const value_b = rho_b.getNeighbor(x + gl, y + gl, z + gl);
        auto const w = weight_xy * weight[2][k];
        density_a += value_a * w;
        density_b += value_b * w;
        if constexpr (with_gradient) {
          double const w_d[3] = {
              weight_d[0][i] * weight[1][j] * weight[2][k],
              weight[0][i] * weight_d[1][j] * weight[2][k],
              weight[0][i] * weight[1][j] * weight_d[2][k]};
#pragma unroll
          for (int dim = 0; dim < 3; ++dim) {
            gradient_a[dim] += value_a * w_d[dim];
            gradient_b[dim] += value_b * w_d[dim];
          }
        }
      }
    }
  }
  auto *const result = out + n_values * index;
  result[0] = density_a;
  result[1] = density_b;
  if constexpr (with_gradient) {
#pragma unroll
    for (int dim = 0; dim < 3; ++dim) {
      result[2 + dim] = gradient_a[dim];
      result[5 + dim] = gradient_b[dim];
    }
  }
}

template <int order, bool with_gradient>
std::vector<double> interpolate(gpu::GPUField<double> const *rho_a,
                                gpu::GPUField<double> const *rho_b,
                                std::vector<double> const &pos,
                                uint_t n_ghost_layers) {
  constexpr std::size_t n_values = with_gradient ? 8u : 2u;
  auto const n_pos = static_cast<uint>(pos.size() / 3u);
  if (n_pos == 0u) {
    return {};
  }
  thrust::device_vector<double> dev_pos(pos.begin(), pos.end());
  thrust::device_vector<double> dev_out(n_values * n_pos, 0.);
  auto const threads_per_block = 64u;
  auto const n_blocks = (n_pos + threads_per_block - 1u) / threads_per_block;
  using Indexing = gpu::FieldIndexing<double>;
  kernel_interpolate_densities<order, with_gradient>
      <<<n_blocks, threads_per_block, 0u, nullptr>>>(
          Indexing::withGhostLayerXYZ(*rho_a, n_ghost_layers).gpuAccess(),
          Indexing::withGhostLayerXYZ(*rho_b, n_ghost_layers).gpuAccess(),
          thrust::raw_pointer_cast(dev_pos.data()),
          thrust::raw_pointer_cast(dev_out.data()), n_pos,
          static_cast<int>(rho_a->xSize()), static_cast<int>(rho_a->ySize()),
          static_cast<int>(rho_a->zSize()), static_cast<int>(n_ghost_layers));
  std::vector<double> out(dev_out.size());
  thrust::copy(dev_out.begin(), dev_out.end(), out.begin());
  return out;
}
// LCOV_EXCL_STOP
} // namespace

std::vector<double> interpolate_densities(gpu::GPUField<double> const *rho_a,
                                          gpu::GPUField<double> const *rho_b,
                                          std::vector<double> const &pos,
                                          uint_t n_ghost_layers) {
  return interpolate<2, false>(rho_a, rho_b, pos, n_ghost_layers);
}

std::vector<double>
interpolate_densities_and_gradients(gpu::GPUField<double> const *rho_a,
                                    gpu::GPUField<double> const *rho_b,
                                    std::vector<double> const &pos,
                                    uint_t n_ghost_layers) {
  return interpolate<3, true>(rho_a, rho_b, pos, n_ghost_layers);
}

uint_t count_negative_densities(gpu::GPUField<double> const *rho_a,
                                gpu::GPUField<double> const *rho_b) {
  thrust::device_vector<unsigned int> count(1u, 0u);
  auto kernel = gpu::make_kernel(kernel_count_negative_densities<double>);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<double>::xyz(*rho_a));
  kernel.addFieldIndexingParam(gpu::FieldIndexing<double>::xyz(*rho_b));
  kernel.addParam(thrust::raw_pointer_cast(count.data()));
  kernel();
  return static_cast<uint_t>(count[0]);
}

void add_density_weighted_force(gpu::GPUField<double> *force_a,
                                gpu::GPUField<double> *force_b,
                                gpu::GPUField<double> const *rho_a,
                                gpu::GPUField<double> const *rho_b,
                                Vector3<double> const &force,
                                uint_t n_ghost_layers) {
  using Indexing = gpu::FieldIndexing<double>;
  auto kernel = gpu::make_kernel(kernel_add_density_weighted_force<double>);
  kernel.addFieldIndexingParam(
      Indexing::withGhostLayerXYZ(*force_a, n_ghost_layers));
  kernel.addFieldIndexingParam(
      Indexing::withGhostLayerXYZ(*force_b, n_ghost_layers));
  kernel.addFieldIndexingParam(
      Indexing::withGhostLayerXYZ(*rho_a, n_ghost_layers));
  kernel.addFieldIndexingParam(
      Indexing::withGhostLayerXYZ(*rho_b, n_ghost_layers));
  kernel.addParam(force[0]);
  kernel.addParam(force[1]);
  kernel.addParam(force[2]);
  kernel();
}

} // namespace color_gradient
} // namespace lbm
} // namespace walberla
