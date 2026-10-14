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
 * VTK writer helper class templates shared by LBWalberlaImpl and
 * LBWalberlaImplSingleComponent.
 */

#include <field/iterators/IteratorMacros.h>

#include <memory>
#include <string>

namespace walberla {

/**
 * @brief Base class for LB field VTK writers.
 * Provides unit conversion and field access for cell-based VTK output.
 * All field is copied to a Vector container before writing.
 * @tparam FloatType   Internal LB precision (float or double).
 * @tparam VecType     Vector type to copy the field data into.
 * @tparam F_SIZE_ARG  Number of components per cell (1, 3, or 9).
 * @tparam OutputType  VTK output precision (default: float).
 */
template <typename FloatType, typename VecType, uint_t F_SIZE_ARG,
          typename OutputType>
class VTKWriter : public vtk::BlockCellDataWriter<OutputType, F_SIZE_ARG> {
public:
  VTKWriter(ConstBlockDataID const &block_id, std::string const &id,
            FloatType unit_conversion)
      : vtk::BlockCellDataWriter<OutputType, F_SIZE_ARG>(id),
        m_conversion(unit_conversion), m_content{} {}

protected:
  void configure() override { WALBERLA_ASSERT_NOT_NULLPTR(this->block_); }

  std::size_t get_first_index(cell_idx_t const x, cell_idx_t const y,
                              cell_idx_t const z) {
    return (static_cast<std::size_t>(x) * m_dims[2] * m_dims[1] +
            static_cast<std::size_t>(y) * m_dims[2] +
            static_cast<std::size_t>(z)) *
           F_SIZE_ARG;
  }

  FloatType m_conversion;
  VecType m_content;
  Vector3<uint_t> m_dims;

public:
  void set_content(VecType content) { m_content = std::move(content); }

  void set_dims(Vector3<uint_t> dims) { m_dims = dims; }
};

template <typename FloatType, typename OutputType = float>
class DensityVTKWriter
    : public VTKWriter<FloatType, std::vector<FloatType>, 1u, OutputType> {
public:
  using Base = VTKWriter<FloatType, std::vector<FloatType>, 1u, OutputType>;
  using Base::Base;
  using Base::evaluate;

protected:
  OutputType evaluate(cell_idx_t const x, cell_idx_t const y,
                      cell_idx_t const z, cell_idx_t const) override {
    WALBERLA_ASSERT(!this->m_content.empty());
    auto const density = this->m_content[this->get_first_index(x, y, z)];
    return numeric_cast<OutputType>(this->m_conversion * density);
  }
};

template <typename FloatType, typename ScalarField, typename OutputType = float>
class ScalarVTKWriter
    : public VTKWriter<FloatType, std::vector<FloatType>, 1u, OutputType> {
public:
  using Base = VTKWriter<FloatType, std::vector<FloatType>, 1u, OutputType>;
  using Base::Base;
  using Base::evaluate;

protected:
  OutputType evaluate(cell_idx_t const x, cell_idx_t const y,
                      cell_idx_t const z, cell_idx_t const) override {
    WALBERLA_ASSERT(!this->m_content.empty());
    auto const field_value = this->m_content[this->get_first_index(x, y, z)];
    return numeric_cast<OutputType>(field_value);
  }
};

template <typename FloatType, typename OutputType = float>
class VelocityVTKWriter
    : public VTKWriter<FloatType, std::vector<FloatType>, 3u, OutputType> {
public:
  using Base = VTKWriter<FloatType, std::vector<FloatType>, 3u, OutputType>;
  using Base::Base;
  using Base::evaluate;

protected:
  OutputType evaluate(cell_idx_t const x, cell_idx_t const y,
                      cell_idx_t const z, cell_idx_t const f) override {
    WALBERLA_ASSERT(!this->m_content.empty());
    auto velocity = this->m_content[this->get_first_index(x, y, z) + f];
    return numeric_cast<OutputType>(this->m_conversion * velocity);
  }
};

template <typename FloatType, typename OutputType = float>
class PressureTensorVTKWriter
    : public VTKWriter<FloatType, std::vector<FloatType>, 9u, OutputType> {
public:
  using Base = VTKWriter<FloatType, std::vector<FloatType>, 9u, OutputType>;
  using Base::Base;
  using Base::evaluate;

protected:
  OutputType evaluate(cell_idx_t const x, cell_idx_t const y,
                      cell_idx_t const z, cell_idx_t const f) override {
    WALBERLA_ASSERT(!this->m_content.empty());
    auto pressure = this->m_content[this->get_first_index(x, y, z) + f];
    return numeric_cast<OutputType>(this->m_conversion * pressure);
  }
};

template <typename FlagField_T, typename OutputType = float>
class BoundaryVTKWriter : public vtk::BlockCellDataWriter<OutputType, 1u> {
public:
  using Base = vtk::BlockCellDataWriter<OutputType, 1u>;
  using Base::evaluate;
  BoundaryVTKWriter(ConstBlockDataID const &flag_field_id,
                    std::string const &id, FlagUID const &boundary_flag)
      : vtk::BlockCellDataWriter<OutputType, 1u>(id),
        m_flag_field_id(flag_field_id), m_flag_field(nullptr),
        m_boundary_flag(boundary_flag) {}

protected:
  void configure() override {
    WALBERLA_ASSERT_NOT_NULLPTR(this->block_);
    m_flag_field = this->block_->template getData<FlagField_T>(m_flag_field_id);
    m_boundary_flag_value = m_flag_field->getFlag(m_boundary_flag);
  }

  OutputType evaluate(cell_idx_t const x, cell_idx_t const y,
                      cell_idx_t const z, cell_idx_t const) override {
    WALBERLA_ASSERT_NOT_NULLPTR(m_flag_field);
    return m_flag_field->isFlagSet(x, y, z, m_boundary_flag_value)
               ? OutputType{1}
               : OutputType{0};
  }

  ConstBlockDataID const m_flag_field_id;
  FlagField_T const *m_flag_field;
  FlagUID const m_boundary_flag;
  typename FlagField_T::flag_t m_boundary_flag_value;
};

} // namespace walberla
