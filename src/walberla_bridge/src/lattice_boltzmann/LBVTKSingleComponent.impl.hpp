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
 * @ref walberla::LBWalberlaImplSingleComponent.
 */

#include "LBVTKHelpers.hpp"

#include <memory>
#include <optional>
#include <string>

namespace walberla {

template <typename FloatType, lbmpy::Arch Architecture>
void LBWalberlaImplSingleComponent<FloatType, Architecture>::
    register_vtk_field_writers(walberla::vtk::VTKOutput &vtk_obj,
                               LatticeModel::units_map const &units,
                               int flag_observables) {
  if (flag_observables & static_cast<int>(OutputVTK::density)) {
    auto const unit_conversion = FloatType_c(units.at("density"));
    auto const blocks = m_lattice->get_blocks();
    WALBERLA_ASSERT_NOT_NULLPTR(blocks);
    auto density_writer = std::make_shared<DensityVTKWriter<FloatType, float>>(
        m_pdf_field_id[0], "density", unit_conversion);
    auto before_function = [this, blocks, density_writer]() {
      for (auto &block : *blocks) {
        auto *pdf_field = block.template getData<PdfField>(m_pdf_field_id[0]);
        auto const bci = pdf_field->xyzSize();
        density_writer->set_content(
            lbm::accessor::Density::get(pdf_field, m_density, bci));
        density_writer->set_dims(
            Vector3<uint_t>(bci.xSize(), bci.ySize(), bci.zSize()));
      }
    };
    vtk_obj.addBeforeFunction(std::move(before_function));
    vtk_obj.addCellDataWriter(density_writer);
  }
  if (flag_observables & static_cast<int>(OutputVTK::velocity_vector)) {
    auto const unit_conversion = FloatType_c(units.at("velocity"));
    auto const blocks = m_lattice->get_blocks();
    WALBERLA_ASSERT_NOT_NULLPTR(blocks);
    auto velocity_writer =
        std::make_shared<VelocityVTKWriter<FloatType, float>>(
            m_pdf_field_id[0], "velocity_vector", unit_conversion);
    auto before_function = [this, blocks, velocity_writer]() {
      for (auto &block : *blocks) {
        auto *velocity_field =
            block.template getData<VectorField>(m_velocity_field_id);

        auto const offset = m_lattice->get_block_corner(block, true);
        auto const *flag_field =
            block.template getData<FlagField>(m_flag_field_id);
        auto const boundary_flag = flag_field->getFlag(Boundary_flag);
        WALBERLA_FOR_ALL_CELLS_XYZ(flag_field, {
          if (flag_field->isFlagSet(x, y, z, boundary_flag)) {
            Cell const global(offset[0] + x, offset[1] + y, offset[2] + z);
            auto const &vel = m_boundary[0]->get_node_value_at_boundary(global);
            Cell const local(x, y, z);
            lbm::accessor::Vector::set(velocity_field, vel, local);
          }
        }) // WALBERLA_FOR_ALL_CELLS_XYZ

        auto const bci = velocity_field->xyzSize();
        velocity_writer->set_content(
            lbm::accessor::Vector::get(velocity_field, bci));
        velocity_writer->set_dims(
            Vector3<uint_t>(bci.xSize(), bci.ySize(), bci.zSize()));
      }
    };

    vtk_obj.addBeforeFunction(std::move(before_function));
    vtk_obj.addCellDataWriter(velocity_writer);
  }
  if (flag_observables & static_cast<int>(OutputVTK::pressure_tensor)) {
    auto const unit_conversion = FloatType_c(units.at("pressure"));
    auto const blocks = m_lattice->get_blocks();
    WALBERLA_ASSERT_NOT_NULLPTR(blocks);
    auto pressure_writer =
        std::make_shared<PressureTensorVTKWriter<FloatType, float>>(
            m_pdf_field_id[0], "pressure_tensor", unit_conversion);
    auto before_function = [this, blocks, pressure_writer]() {
      for (auto &block : *blocks) {
        auto *pdf_field = block.template getData<PdfField>(m_pdf_field_id[0]);
        auto const bci = pdf_field->xyzSize();
        auto values =
            lbm::accessor::PressureTensor::get(pdf_field, m_density, bci);
        for (std::size_t n = 0u; n < values.size(); n += 9u) {
          pressure_tensor_correction(
              std::span<FloatType, 9ul>(&values[n], 9ul));
        }
        pressure_writer->set_content(std::move(values));
        pressure_writer->set_dims(
            Vector3<uint_t>(bci.xSize(), bci.ySize(), bci.zSize()));
      }
    };
    vtk_obj.addBeforeFunction(std::move(before_function));
    vtk_obj.addCellDataWriter(pressure_writer);
  }
  if (flag_observables & static_cast<int>(OutputVTK::boundary)) {
    vtk_obj.addCellDataWriter(
        std::make_shared<BoundaryVTKWriter<FlagField, float>>(
            m_flag_field_id, "boundary", Boundary_flag));
  }
}

} // namespace walberla
