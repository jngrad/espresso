/*
 * Copyright (C) 2019-2026 The ESPResSo project
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
 * Out-of-class position-based interpolation definitions for
 * @ref walberla::LBWalberlaImplColorGradient.
 */

#include <utils/Vector.hpp>
#include <utils/interpolation/bspline_3d.hpp>
#include <utils/interpolation/bspline_3d_gradient.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <iterator>
#include <memory>
#include <optional>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace walberla {

#if defined(__CUDACC__) and defined(WALBERLA_BUILD_WITH_CUDA)
namespace detail {
/**
 * @brief Positions relative to the corner of the rank's only block (as on
 * the GPU there is one block per rank), and whether each position lies in
 * the block extended by one cell. The CPU kernels skip the other positions.
 */
inline auto block_local_positions(::LatticeWalberla const &lattice,
                                  std::vector<Utils::Vector3d> const &pos) {
  assert(lattice.get_blocks()->getNumberOfBlocks() == 1u);
  auto const &block = *(lattice.get_blocks()->begin());
  auto const &aabb = block.getAABB();
  auto const extended = aabb.getExtended(real_c(1u));
  std::vector<double> local;
  std::vector<bool> inside;
  local.reserve(3u * pos.size());
  inside.reserve(pos.size());
  for (auto const &p : pos) {
    inside.emplace_back(
        extended.contains(real_c(p[0]), real_c(p[1]), real_c(p[2])));
    for (std::size_t i : {0u, 1u, 2u}) {
      local.emplace_back(p[i] - aabb.min()[i]);
    }
  }
  return std::make_pair(std::move(local), std::move(inside));
}
} // namespace detail
#endif

/**
 * @brief Distribute forces density-weighted to the two component lattices
 * at given positions. Uses B-spline interpolation to spread each force
 * over the surrounding lattice nodes, split between components by local
 * rho_a/rho_b. Writes into force_to_be_applied (double-buffered). On the GPU,
 * the densities are interpolated on the device, the split is computed on the
 * host, and each component's share is spread on the device.
 */
template <typename FloatType, lbmpy::Arch Architecture>
void LBWalberlaImplColorGradient<FloatType, Architecture>::add_forces_at_pos(
    std::vector<Utils::Vector3d> const &pos,
    std::vector<Utils::Vector3d> const &forces) {
  assert(pos.size() == forces.size());
  if (pos.empty()) {
    return;
  }
  if constexpr (Architecture == lbmpy::Arch::CPU) {
    auto const kernel = make_density_weighted_force_interpolation_kernel();
    for (std::size_t i = 0ul; i < pos.size(); ++i) {
      kernel(pos[i], forces[i]);
    }
  }
#if defined(__CUDACC__) and defined(WALBERLA_BUILD_WITH_CUDA)
  if constexpr (Architecture == lbmpy::Arch::GPU) {
    auto const &lattice = get_lattice();
    auto const gl = lattice.get_ghost_layers();
    auto &block = *(lattice.get_blocks()->begin());
    auto const [local, inside] = detail::block_local_positions(lattice, pos);
    auto const densities = lbm::color_gradient::interpolate_densities(
        block.template getData<ScalarField>(m_rho_field_id[0]),
        block.template getData<ScalarField>(m_rho_field_id[1]), local, gl);
    std::vector<double> positions, forces_a, forces_b;
    for (std::size_t i = 0ul; i < pos.size(); ++i) {
      auto const rho_a = densities[2u * i + 0u];
      auto const rho_b = densities[2u * i + 1u];
      auto const total_rho = rho_a + rho_b;
      if (not inside[i] or total_rho == 0.) {
        continue;
      }
      auto const inv_rho = 1. / total_rho;
      for (std::size_t j : {0u, 1u, 2u}) {
        positions.emplace_back(local[3u * i + j]);
        forces_a.emplace_back((rho_a * inv_rho) * forces[i][j]);
        forces_b.emplace_back((rho_b * inv_rho) * forces[i][j]);
      }
    }
    if (positions.empty()) {
      return;
    }
    lbm::accessor::Interpolation::add_force(
        block.template getData<VectorField>(m_force_to_be_applied_id[0]),
        positions, forces_a, gl);
    lbm::accessor::Interpolation::add_force(
        block.template getData<VectorField>(m_force_to_be_applied_id[1]),
        positions, forces_b, gl);
  }
#endif
}

/**
 * @brief Distribute solvation force on the two lattices
 * for the two components at given positions. On the GPU, densities and
 * their gradients are interpolated on the device, the forces are computed
 * on the host, and spread on the device.
 */
template <typename FloatType, lbmpy::Arch Architecture>
void LBWalberlaImplColorGradient<FloatType, Architecture>::
    add_solvation_forces_at_pos(std::vector<Utils::Vector3d> const &pos,
                                std::vector<double> const &delta_mus) {
  assert(pos.size() == delta_mus.size());
  if (pos.empty())
    return;
  if constexpr (Architecture == lbmpy::Arch::CPU) {
    auto const kernel = make_solvation_force_interpolation_kernel();
    for (std::size_t i = 0ul; i < pos.size(); ++i) {
      kernel(pos[i], delta_mus[i]);
    }
  }
#if defined(__CUDACC__) and defined(WALBERLA_BUILD_WITH_CUDA)
  if constexpr (Architecture == lbmpy::Arch::GPU) {
    auto const &lattice = get_lattice();
    auto const gl = lattice.get_ghost_layers();
    auto &block = *(lattice.get_blocks()->begin());
    auto const [local, inside] = detail::block_local_positions(lattice, pos);
    auto const values = lbm::color_gradient::interpolate_densities_and_gradients(
        block.template getData<ScalarField>(m_rho_field_id[0]),
        block.template getData<ScalarField>(m_rho_field_id[1]), local, gl);
    std::vector<double> positions, forces_a, forces_b;
    for (std::size_t i = 0ul; i < pos.size(); ++i) {
      auto const *const v = values.data() + 8u * i;
      auto const rho_a = v[0];
      auto const rho_b = v[1];
      auto const total_rho = rho_a + rho_b;
      if (not inside[i] or total_rho == 0.) {
        continue;
      }
      auto const inv_rho_sq = 1. / (total_rho * total_rho);
      // same expressions as the CPU kernel; the device multiplies by the
      // spreading weight
      auto const f_a =
          +delta_mus[i] * rho_b * inv_rho_sq * Utils::Vector3d{v[2], v[3], v[4]};
      auto const f_b =
          -delta_mus[i] * rho_a * inv_rho_sq * Utils::Vector3d{v[5], v[6], v[7]};
      for (std::size_t j : {0u, 1u, 2u}) {
        positions.emplace_back(local[3u * i + j]);
        forces_a.emplace_back(f_a[j]);
        forces_b.emplace_back(f_b[j]);
      }
    }
    if (positions.empty()) {
      return;
    }
    lbm::accessor::Interpolation::add_force(
        block.template getData<VectorField>(m_force_to_be_applied_id[0]),
        positions, forces_a, gl);
    lbm::accessor::Interpolation::add_force(
        block.template getData<VectorField>(m_force_to_be_applied_id[1]),
        positions, forces_b, gl);
  }
#endif
}

template <typename FloatType, lbmpy::Arch Architecture>
auto LBWalberlaImplColorGradient<
    FloatType, Architecture>::make_solvation_force_interpolation_kernel()
    /* Only for two-component LB, adds solvation force to fluid*/
    const {
  auto const &lattice = *m_lattice;
  auto const &blocks = *lattice.get_blocks();
  assert(lattice.get_ghost_layers() == 1u);
  return [&](Utils::Vector3d const &pos, double delta_mu) {
    if (not get_block_extended(lattice, pos, 1u)) {
      return;
    }
    auto const grid_spacing = Utils::Vector3d::broadcast(1.);
    auto const offset = Utils::Vector3d::broadcast(.5);
    /* Accumulate grad(rho_a) and grad(rho_b) at particle position using
     * bspline_3d_gradient*/
    Utils::Vector3d grad_rho_a{}, grad_rho_b{};
    Utils::Interpolation::bspline_3d_gradient<3>( // SOLVATION_BSPLINE_ORDER
        pos,
        [&](std::array<int, 3> const node, Utils::Vector3d const &grad_weight) {
          auto block = get_block_extended(lattice, node, 1u);
          if (!block)
            return;
          auto cell = to_cell(node);
          blocks.transformGlobalToBlockLocalCell(cell, *block);
          auto const cell_rho_a = static_cast<double>(
              block
                  ->template uncheckedFastGetData<ScalarField>(
                      m_rho_field_id[0])
                  ->get(cell));
          auto const cell_rho_b = static_cast<double>(
              block
                  ->template uncheckedFastGetData<ScalarField>(
                      m_rho_field_id[1])
                  ->get(cell));
          grad_rho_a += cell_rho_a * grad_weight;
          grad_rho_b += cell_rho_b * grad_weight;
        },
        grid_spacing, offset);

    /* Accumulate rho_a and rho_b at particle position using bspline_3d*/
    double rho_a{}, rho_b{};
    Utils::Interpolation::bspline_3d<3>( // SOLVATION_BSPLINE_ORDER
        pos,
        [&](std::array<int, 3> const node, double weight) {
          auto block = get_block_extended(lattice, node, 1u);
          if (!block)
            return;
          auto cell = to_cell(node);
          blocks.transformGlobalToBlockLocalCell(cell, *block);
          auto const cell_rho_a = static_cast<double>(
              block
                  ->template uncheckedFastGetData<ScalarField>(
                      m_rho_field_id[0])
                  ->get(cell));
          auto const cell_rho_b = static_cast<double>(
              block
                  ->template uncheckedFastGetData<ScalarField>(
                      m_rho_field_id[1])
                  ->get(cell));
          rho_a += cell_rho_a * weight;
          rho_b += cell_rho_b * weight;
        },
        grid_spacing, offset);

    /*spread force back onto fluid using bspline_3d*/
    interpolate_bspline_at_pos(pos, [&](std::array<int, 3> const node,
                                        double weight) {
      auto block = get_block_extended(lattice, node, 0u);
      if (!block)
        block = get_block_extended(lattice, node, 1u);
      if (!block)
        return;
      auto cell = to_cell(node);
      blocks.transformGlobalToBlockLocalCell(cell, *block);
      auto const total_rho = rho_a + rho_b;
      if (total_rho == 0.)
        /* No fluid to couple to */
        return;
      auto const inv_rho_sq = 1. / (total_rho * total_rho);
      // f^a = +delta_mu * rho_b * inv_rho_sq * grad_rho_a * weight
      // f^b = -delta_mu * rho_a * inv_rho_sq * grad_rho_b * weight
      auto const f_a = to_vector3<FloatType>(+delta_mu * rho_b * inv_rho_sq *
                                             grad_rho_a * weight);
      auto const f_b = to_vector3<FloatType>(-delta_mu * rho_a * inv_rho_sq *
                                             grad_rho_b * weight);
      auto field_a = block->template uncheckedFastGetData<VectorField>(
          m_force_to_be_applied_id[0]);
      auto field_b = block->template uncheckedFastGetData<VectorField>(
          m_force_to_be_applied_id[1]);
      lbm::accessor::Vector::add(field_a, f_a, cell);
      lbm::accessor::Vector::add(field_b, f_b, cell);
    });
  };
}

template <typename FloatType, lbmpy::Arch Architecture>
auto LBWalberlaImplColorGradient<
    FloatType, Architecture>::make_density_weighted_force_interpolation_kernel()
    /* Only for two-component LB, e.g. friction coupling*/
    const {
  auto const &lattice = *m_lattice;
  auto const &blocks = *lattice.get_blocks();
  assert(lattice.get_ghost_layers() == 1u);
  return [&](Utils::Vector3d const &pos, Utils::Vector3d const &force) {
    if (not get_block_extended(lattice, pos, 1u)) {
      return;
    }
    auto const grid_spacing = Utils::Vector3d::broadcast(1.);
    auto const offset = Utils::Vector3d::broadcast(.5);
    /* Accumulate rho_a and rho_b at particle position using bspline_3d*/
    double rho_a{}, rho_b{};
    Utils::Interpolation::bspline_3d<2>(
        pos,
        [&](std::array<int, 3> const node, double weight) {
          auto block = get_block_extended(lattice, node, 1u);
          if (!block)
            return;
          auto cell = to_cell(node);
          blocks.transformGlobalToBlockLocalCell(cell, *block);
          auto const cell_rho_a = static_cast<double>(
              block
                  ->template uncheckedFastGetData<ScalarField>(
                      m_rho_field_id[0])
                  ->get(cell));
          auto const cell_rho_b = static_cast<double>(
              block
                  ->template uncheckedFastGetData<ScalarField>(
                      m_rho_field_id[1])
                  ->get(cell));
          rho_a += cell_rho_a * weight;
          rho_b += cell_rho_b * weight;
        },
        grid_spacing, offset);

    /*spread force back onto fluid using bspline_3d*/
    interpolate_bspline_at_pos(
        pos, [&](std::array<int, 3> const node, double weight) {
          /* no conversion factor m_zc_to_lb since densities for two-component
           * lbs are not saved zero_centered*/
          auto block = get_block_extended(lattice, node, 0u);
          if (!block)
            block = get_block_extended(lattice, node, 1u);
          if (block) {
            auto cell = to_cell(node);
            blocks.transformGlobalToBlockLocalCell(cell, *block);
            auto const total_rho = rho_a + rho_b;
            if (total_rho == 0.)
              /* No fluid to couple to */
              return;
            auto const inv_rho = 1. / total_rho;
            auto const weighted_force_a =
                to_vector3<FloatType>(weight * (rho_a * inv_rho) * force);
            auto const weighted_force_b =
                to_vector3<FloatType>(weight * (rho_b * inv_rho) * force);
            auto field_a = block->template uncheckedFastGetData<VectorField>(
                m_force_to_be_applied_id[0]);
            auto field_b = block->template uncheckedFastGetData<VectorField>(
                m_force_to_be_applied_id[1]);
            lbm::accessor::Vector::add(field_a, weighted_force_a, cell);
            lbm::accessor::Vector::add(field_b, weighted_force_b, cell);
          }
        });
  };
}

template <typename FloatType, lbmpy::Arch Architecture>
auto LBWalberlaImplColorGradient<FloatType, Architecture>::
    make_solvation_particle_force_kernel() const {
  auto const &lattice = *m_lattice;
  auto const &blocks = *lattice.get_blocks();
  assert(lattice.get_ghost_layers() == 1u);
  return [&](Utils::Vector3d const &pos, double delta_mu) -> Utils::Vector3d {
    if (not get_block_extended(lattice, pos, 1u))
      return {};
    auto const grid_spacing = Utils::Vector3d::broadcast(1.);
    auto const offset = Utils::Vector3d::broadcast(.5);
    Utils::Vector3d grad_rho_a{}, grad_rho_b{};
    Utils::Interpolation::bspline_3d_gradient<3>( // SOLVATION_BSPLINE_ORDER
        pos,
        [&](std::array<int, 3> const node, Utils::Vector3d const &grad_weight) {
          auto block = get_block_extended(lattice, node, 1u);
          if (!block)
            return;
          auto cell = to_cell(node);
          blocks.transformGlobalToBlockLocalCell(cell, *block);
          auto const cell_rho_a = static_cast<double>(
              block->template uncheckedFastGetData<ScalarField>(m_rho_field_id[0])->get(cell));
          auto const cell_rho_b = static_cast<double>(
              block->template uncheckedFastGetData<ScalarField>(m_rho_field_id[1])->get(cell));
          grad_rho_a += cell_rho_a * grad_weight;
          grad_rho_b += cell_rho_b * grad_weight;
        },
        grid_spacing, offset);
    double rho_a{}, rho_b{};
    Utils::Interpolation::bspline_3d<3>( // SOLVATION_BSPLINE_ORDER
        pos,
        [&](std::array<int, 3> const node, double weight) {
          auto block = get_block_extended(lattice, node, 1u);
          if (!block)
            return;
          auto cell = to_cell(node);
          blocks.transformGlobalToBlockLocalCell(cell, *block);
          auto const cell_rho_a = static_cast<double>(
              block->template uncheckedFastGetData<ScalarField>(m_rho_field_id[0])->get(cell));
          auto const cell_rho_b = static_cast<double>(
              block->template uncheckedFastGetData<ScalarField>(m_rho_field_id[1])->get(cell));
          rho_a += cell_rho_a * weight;
          rho_b += cell_rho_b * weight;
        },
        grid_spacing, offset);
    auto const total_rho = rho_a + rho_b;
    if (total_rho == 0.)
      return {};
    auto const inv_rho_sq = 1. / (total_rho * total_rho);
    return -delta_mu * inv_rho_sq * (rho_b * grad_rho_a - rho_a * grad_rho_b);
  };
}

template <typename FloatType, lbmpy::Arch Architecture>
std::vector<Utils::Vector3d>
LBWalberlaImplColorGradient<FloatType, Architecture>::
    get_solvation_particle_forces_at_pos(
        std::vector<Utils::Vector3d> const &pos,
        std::vector<double> const &delta_mus) {
  assert(pos.size() == delta_mus.size());
  if (pos.empty())
    return {};
  std::vector<Utils::Vector3d> forces{};
  forces.reserve(pos.size());
  if constexpr (Architecture == lbmpy::Arch::CPU) {
    auto const kernel = make_solvation_particle_force_kernel();
    for (std::size_t i = 0ul; i < pos.size(); ++i)
      forces.emplace_back(kernel(pos[i], delta_mus[i]));
  }
#if defined(__CUDACC__) and defined(WALBERLA_BUILD_WITH_CUDA)
  if constexpr (Architecture == lbmpy::Arch::GPU) {
    auto const &lattice = get_lattice();
    auto &block = *(lattice.get_blocks()->begin());
    auto const [local, inside] = detail::block_local_positions(lattice, pos);
    auto const values = lbm::color_gradient::interpolate_densities_and_gradients(
        block.template getData<ScalarField>(m_rho_field_id[0]),
        block.template getData<ScalarField>(m_rho_field_id[1]), local,
        lattice.get_ghost_layers());
    for (std::size_t i = 0ul; i < pos.size(); ++i) {
      auto const *const v = values.data() + 8u * i;
      auto const rho_a = v[0];
      auto const rho_b = v[1];
      auto const total_rho = rho_a + rho_b;
      if (not inside[i] or total_rho == 0.) {
        forces.emplace_back(Utils::Vector3d{});
        continue;
      }
      auto const inv_rho_sq = 1. / (total_rho * total_rho);
      auto const grad_rho_a = Utils::Vector3d{v[2], v[3], v[4]};
      auto const grad_rho_b = Utils::Vector3d{v[5], v[6], v[7]};
      forces.emplace_back(-delta_mus[i] * inv_rho_sq *
                          (rho_b * grad_rho_a - rho_a * grad_rho_b));
    }
  }
#endif
  return forces;
}

template <typename FloatType, lbmpy::Arch Architecture>
bool LBWalberlaImplColorGradient<FloatType, Architecture>::add_force_at_pos(
    Utils::Vector3d const &pos, Utils::Vector3d const &force) {
  if (!m_lattice->pos_in_local_halo(pos))
    return false;
  if constexpr (Architecture == lbmpy::Arch::CPU) {
    auto const kernel = make_density_weighted_force_interpolation_kernel();
    kernel(pos, force);
  } else {
    add_forces_at_pos({pos}, {force});
  }
  return true;
}

} // namespace walberla
