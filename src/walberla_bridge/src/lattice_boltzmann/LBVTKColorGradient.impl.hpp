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
 * Out-of-class VTK writer registration definition for
 * @ref walberla::LBWalberlaImplColorGradient.
 */

#include "LBVTKHelpers.hpp"

#include <memory>
#include <optional>
#include <string>

namespace walberla {

template <typename FloatType, lbmpy::Arch Architecture>
void LBWalberlaImplColorGradient<FloatType, Architecture>::
    register_vtk_field_writers(walberla::vtk::VTKOutput &vtk_obj,
                               LatticeModel::units_map const &units,
                               int flag_observables) {
  // On the GPU, every field gets its own host copy, filled before each write;
  // the writers read the host copy. Copy 0 is the phase field, copies 1 and 2
  // are the densities of components a and b.
  auto const add_scalar_writer = [&](BlockDataID const &field_id,
                                     [[maybe_unused]] std::size_t copy,
                                     std::string const &name,
                                     FloatType unit_conversion) {
    if constexpr (Architecture == lbmpy::Arch::CPU) {
      vtk_obj.addCellDataWriter(
          std::make_shared<ScalarVTKWriter<FloatType, ScalarField, float>>(
              field_id, name, unit_conversion));
    }
#if defined(__CUDACC__) and defined(WALBERLA_BUILD_WITH_CUDA)
    else {
      auto const &blocks = m_lattice->get_blocks();
      auto &cpu_field_id = (copy == 0u) ? m_phasefield_cpu_field_id
                                        : m_rho_cpu_field_id[copy - 1u];
      if (not cpu_field_id) {
        cpu_field_id = field::addToStorage<ScalarFieldCpu>(
            blocks, name + "_cpu", FloatType{0}, field::fzyx,
            m_lattice->get_ghost_layers(), m_host_field_allocator);
      }
      vtk_obj.addBeforeFunction(
          gpu::fieldCpyFunctor<ScalarFieldCpu, ScalarField>(
              blocks, *cpu_field_id, field_id));
      vtk_obj.addCellDataWriter(
          std::make_shared<ScalarVTKWriter<FloatType, ScalarFieldCpu, float>>(
              *cpu_field_id, name, unit_conversion));
    }
#endif
  };

  if (flag_observables & static_cast<int>(OutputVTK::phasefield)) {
    add_scalar_writer(m_phasefield_id, 0u, "phasefield", FloatType{1});
  }
  if (flag_observables & static_cast<int>(OutputVTK::density_a)) {
    add_scalar_writer(m_rho_field_id[0], 1u, "density_a",
                      FloatType_c(zero_centered_to_md(units.at("density"))));
  }
  if (flag_observables & static_cast<int>(OutputVTK::density_b)) {
    add_scalar_writer(m_rho_field_id[1], 2u, "density_b",
                      FloatType_c(zero_centered_to_md(units.at("density"))));
  }
  if (flag_observables & static_cast<int>(OutputVTK::velocity_vector)) {
    auto const unit_conversion = FloatType_c(units.at("velocity"));
    if constexpr (Architecture == lbmpy::Arch::CPU) {
      vtk_obj.addCellDataWriter(
          std::make_shared<VelocityVTKWriter<FloatType, VectorField, float>>(
              m_velocity_field_id, "velocity_vector", unit_conversion));
    }
#if defined(__CUDACC__) and defined(WALBERLA_BUILD_WITH_CUDA)
    else {
      auto const &blocks = m_lattice->get_blocks();
      if (not m_vel_cpu_field_id) {
        m_vel_cpu_field_id = field::addToStorage<VectorFieldCpu>(
            blocks, "vel_cpu", FloatType{0}, field::fzyx,
            m_lattice->get_ghost_layers(), m_host_field_allocator);
      }
      vtk_obj.addBeforeFunction(
          gpu::fieldCpyFunctor<VectorFieldCpu, VectorField>(
              blocks, *m_vel_cpu_field_id, m_velocity_field_id));
      vtk_obj.addCellDataWriter(
          std::make_shared<VelocityVTKWriter<FloatType, VectorFieldCpu, float>>(
              *m_vel_cpu_field_id, "velocity_vector", unit_conversion));
    }
#endif
  }
}

} // namespace walberla
