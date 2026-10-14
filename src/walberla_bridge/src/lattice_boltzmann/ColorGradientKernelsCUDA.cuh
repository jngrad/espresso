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

#pragma once

/**
 * @file
 * Hand-written CUDA kernels of the two-component color-gradient LB that the
 * lbmpy-generated kernels and field accessors don't provide. The GPU model
 * only runs in double precision.
 */

#include <core/DataTypes.h>
#include <core/math/Vector3.h>
#include <gpu/GPUField.h>

#include <vector>

namespace walberla {
namespace lbm {
namespace color_gradient {

/**
 * @brief Count the interior cells where the density of at least one
 * component is negative.
 */
uint_t count_negative_densities(gpu::GPUField<double> const *rho_a,
                                gpu::GPUField<double> const *rho_b);

/**
 * @brief Split a force density between the components in proportion to
 * their densities, @f$ f_c = \rho_c / (\rho_a + \rho_b)\, f @f$, and add it
 * to both force fields in every cell, including @p n_ghost_layers ghost
 * layers. Cells without fluid (@f$ \rho_a + \rho_b \le 0 @f$) get nothing.
 */
void add_density_weighted_force(gpu::GPUField<double> *force_a,
                                gpu::GPUField<double> *force_b,
                                gpu::GPUField<double> const *rho_a,
                                gpu::GPUField<double> const *rho_b,
                                Vector3<double> const &force,
                                uint_t n_ghost_layers);

/**
 * @brief Interpolate both component densities at block-local positions
 * with order-2 B-splines, as the CPU friction coupling does.
 * @param rho_a, rho_b  density fields of the components
 * @param pos           positions, 3 values each, relative to the block corner
 * @param n_ghost_layers number of ghost layers of the fields
 * @return @f$ (\rho_a, \rho_b) @f$ per position
 */
std::vector<double> interpolate_densities(gpu::GPUField<double> const *rho_a,
                                          gpu::GPUField<double> const *rho_b,
                                          std::vector<double> const &pos,
                                          uint_t n_ghost_layers);

/**
 * @brief Interpolate both component densities and their gradients at
 * block-local positions with order-3 B-splines, as the CPU solvation
 * coupling does. Stencil nodes outside the field (including its ghost
 * layers) are skipped.
 * @return @f$ (\rho_a, \rho_b, \nabla\rho_a, \nabla\rho_b) @f$, 8 values
 * per position
 */
std::vector<double>
interpolate_densities_and_gradients(gpu::GPUField<double> const *rho_a,
                                    gpu::GPUField<double> const *rho_b,
                                    std::vector<double> const &pos,
                                    uint_t n_ghost_layers);

} // namespace color_gradient
} // namespace lbm
} // namespace walberla
