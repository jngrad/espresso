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
 * Out-of-class slice access definitions for
 * @ref walberla::LBWalberlaImplColorGradient.
 */

#include <utils/Vector.hpp>

#include <span>
#include <vector>

namespace walberla {

template <typename FloatType, lbmpy::Arch Architecture>
std::vector<double>
LBWalberlaImplColorGradient<FloatType, Architecture>::get_slice_velocity(
    Utils::Vector3i const &lower_corner,
    Utils::Vector3i const &upper_corner) const {
  std::vector<double> out;
  for_each_block_in_slice(
      get_lattice(), lower_corner, upper_corner,
      [&](auto &block, auto const &bci, auto const &ci,
          auto const &block_offset) {
        if (out.empty())
          out.resize(3u * ci.numCells());
        auto const field =
            block.template getData<VectorField>(m_velocity_field_id);
        auto values = lbm::accessor::Vector::get(field, bci);

        auto kernel = [&values, &out, this](unsigned const block_index,
                                            unsigned const local_index,
                                            Utils::Vector3i const &node) {
          if (m_boundary[0]->node_is_boundary(node)) {
            auto const &vec = m_boundary[0]->get_node_value_at_boundary(node);
            for (uint_t f = 0u; f < 3u; ++f) {
              out[3u * local_index + f] = vec[f];
            }
          } else {
            for (uint_t f = 0u; f < 3u; ++f) {
              out[3u * local_index + f] =
                  double_c(values[3u * block_index + f]);
            }
          }
        };

        copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);
      });
  return out;
}

template <typename FloatType, lbmpy::Arch Architecture>
void LBWalberlaImplColorGradient<FloatType, Architecture>::set_slice_velocity(
    Utils::Vector3i const &lower_corner, Utils::Vector3i const &upper_corner,
    std::vector<double> const &velocity) {
  throw std::runtime_error(
      "set_slice_velocity is not supported for two-component LB. "
      "Set densities and populations instead to control the barycentric "
      "velocity.");
}

template <typename FloatType, lbmpy::Arch Architecture>
std::vector<double> LBWalberlaImplColorGradient<FloatType, Architecture>::
    get_slice_last_applied_force(Utils::Vector3i const &lower_corner,
                                 Utils::Vector3i const &upper_corner) const {
  std::vector<double> out;
  for_each_block_in_slice(
      get_lattice(), lower_corner, upper_corner,
      [&](auto const &block, auto const &bci, auto const &ci,
          auto const &block_offset) {
        if (out.empty())
          out.resize(6u * ci.numCells());
        auto const field_a =
            block.template getData<VectorField>(m_last_applied_force_field_id[0]);
        auto const field_b =
            block.template getData<VectorField>(m_last_applied_force_field_id[1]);
        auto const values_a = lbm::accessor::Vector::get(field_a, bci);
        auto const values_b = lbm::accessor::Vector::get(field_b, bci);

        auto kernel = [&values_a, &values_b, &out](unsigned const block_index,
                                                    unsigned const local_index,
                                                    Utils::Vector3i const &) {
          for (uint_t f = 0u; f < 3u; ++f) {
            out[6u * local_index + f] = values_a[3u * block_index + f];
            out[6u * local_index + 3u + f] = values_b[3u * block_index + f];
          }
        };

        copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);
      });
  return out;
}

template <typename FloatType, lbmpy::Arch Architecture>
void LBWalberlaImplColorGradient<FloatType, Architecture>::
    set_slice_last_applied_force(Utils::Vector3i const &lower_corner,
                                 Utils::Vector3i const &upper_corner,
                                 std::vector<double> const &force) {
  m_pending_ghost_comm.set(GhostComm::LAF);
  for_each_block_in_slice(
      get_lattice(), lower_corner, upper_corner,
      [&](auto &block, auto const &bci, auto const &ci,
          auto const &block_offset) {
        assert(force.size() == 6u * ci.numCells());
        auto force_field_a =
            block.template getData<VectorField>(m_last_applied_force_field_id[0]);
        auto force_field_b =
            block.template getData<VectorField>(m_last_applied_force_field_id[1]);
        std::vector<FloatType> values_a(3u * bci.numCells());
        std::vector<FloatType> values_b(3u * bci.numCells());

        auto kernel = [&values_a, &values_b, &force](
                          unsigned const block_index, unsigned const local_index,
                          Utils::Vector3i const &) {
          for (uint_t f = 0u; f < 3u; ++f) {
            values_a[3u * block_index + f] =
                numeric_cast<FloatType>(force[6u * local_index + f]);
            values_b[3u * block_index + f] =
                numeric_cast<FloatType>(force[6u * local_index + 3u + f]);
          }
        };

        copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);
        lbm::accessor::Vector::set(force_field_a, values_a, bci);
        lbm::accessor::Vector::set(force_field_b, values_b, bci);
      });
}

template <typename FloatType, lbmpy::Arch Architecture>
std::vector<double>
LBWalberlaImplColorGradient<FloatType, Architecture>::get_slice_population(
    Utils::Vector3i const &lower_corner,
    Utils::Vector3i const &upper_corner) const {
  std::vector<double> out;
  auto const pop_per_node = 2u * this->stencil_size();
  for_each_block_in_slice(
      get_lattice(), lower_corner, upper_corner,
      [&](auto const &block, auto const &bci, auto const &ci,
          auto const &block_offset) {
        if (out.empty())
          out.resize(pop_per_node * ci.numCells());

        auto const pdf_field_a =
            block.template getData<PdfField>(m_pdf_field_id[0]);
        auto const values_a = lbm::accessor::Population::get(pdf_field_a, bci);

        auto const pdf_field_b =
            block.template getData<PdfField>(m_pdf_field_id[1]);
        auto const values_b = lbm::accessor::Population::get(pdf_field_b, bci);

        auto kernel = [&values_a, &values_b, &out, this](
                          unsigned const block_index,
                          unsigned const local_index, Utils::Vector3i const &) {
          auto const ss = this->stencil_size();
          for (uint_t f = 0u; f < ss; ++f) {
            out[2u * ss * local_index + f] = values_a[ss * block_index + f];
            out[2u * ss * local_index + ss + f] =
                values_b[ss * block_index + f];
          }
        };

        copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);
      });
  return out;
}

template <typename FloatType, lbmpy::Arch Architecture>
void LBWalberlaImplColorGradient<FloatType, Architecture>::set_slice_population(
    Utils::Vector3i const &lower_corner, Utils::Vector3i const &upper_corner,
    std::vector<double> const &population) {
  m_pending_ghost_comm.set(GhostComm::PDF);
  for_each_block_in_slice(
      get_lattice(), lower_corner, upper_corner,
      [&](auto &block, auto const &bci, auto const &ci,
          auto const &block_offset) {
        auto const pop_per_node = 2u * this->stencil_size();
        assert(population.size() == pop_per_node * ci.numCells());
        (void)pop_per_node;
        auto pdf_field_a = block.template getData<PdfField>(m_pdf_field_id[0]);
        auto pdf_field_b = block.template getData<PdfField>(m_pdf_field_id[1]);
        std::vector<FloatType> values_a(this->stencil_size() * bci.numCells());
        std::vector<FloatType> values_b(this->stencil_size() * bci.numCells());

        auto kernel = [&values_a, &values_b, &population, this](
                          unsigned const block_index,
                          unsigned const local_index, Utils::Vector3i const &) {
          auto const ss = this->stencil_size();
          for (uint_t f = 0u; f < ss; ++f) {
            values_a[ss * block_index + f] =
                numeric_cast<FloatType>(population[2u * ss * local_index + f]);
            values_b[ss * block_index + f] = numeric_cast<FloatType>(
                population[2u * ss * local_index + ss + f]);
          }
        };

        copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);
        // populations only, see set_node_population()
        lbm::accessor::Population::set(pdf_field_a, values_a, bci);
        lbm::accessor::Population::set(pdf_field_b, values_b, bci);
      });
}

template <typename FloatType, lbmpy::Arch Architecture>
std::vector<double>
LBWalberlaImplColorGradient<FloatType, Architecture>::get_slice_density(
    Utils::Vector3i const &lower_corner,
    Utils::Vector3i const &upper_corner) const {
  std::vector<double> out;
  for_each_block_in_slice(
      get_lattice(), lower_corner, upper_corner,
      [&](auto const &block, auto const &bci, auto const &ci,
          auto const &block_offset) {
        if (out.empty())
          out.resize(2u * ci.numCells());

        auto const rho_a_field =
            block.template getData<ScalarField>(m_rho_field_id[0]);
        auto const rho_b_field =
            block.template getData<ScalarField>(m_rho_field_id[1]);
        // field values in x-major, z-minor order
        auto const values_a = ek::accessor::Scalar::get(rho_a_field, bci);
        auto const values_b = ek::accessor::Scalar::get(rho_b_field, bci);

        auto kernel = [&values_a, &values_b, &out](unsigned const block_index,
                                                   unsigned const local_index,
                                                   Utils::Vector3i const &) {
          out[2u * local_index + 0u] = double_c(values_a[block_index]);
          out[2u * local_index + 1u] = double_c(values_b[block_index]);
        };

        copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);
      });
  return out;
}

template <typename FloatType, lbmpy::Arch Architecture>
std::vector<double>
LBWalberlaImplColorGradient<FloatType, Architecture>::get_slice_phasefield(
    Utils::Vector3i const &lower_corner,
    Utils::Vector3i const &upper_corner) const {
  std::vector<double> out;
  for_each_block_in_slice(
      get_lattice(), lower_corner, upper_corner,
      [&](auto const &block, auto const &bci, auto const &ci,
          auto const &block_offset) {
        if (out.empty())
          out.resize(ci.numCells());

        auto const phasefield_field =
            block.template getData<ScalarField>(m_phasefield_id);
        // field values in x-major, z-minor order
        auto const values = ek::accessor::Scalar::get(phasefield_field, bci);

        auto kernel = [&values, &out](unsigned const block_index,
                                      unsigned const local_index,
                                      Utils::Vector3i const &) {
          out[local_index] = double_c(values[block_index]);
        };

        copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);
      });
  return out;
}

template <typename FloatType, lbmpy::Arch Architecture>
void LBWalberlaImplColorGradient<FloatType, Architecture>::set_slice_density(
    Utils::Vector3i const &lower_corner, Utils::Vector3i const &upper_corner,
    std::vector<double> const &density) {
  m_pending_ghost_comm.set(GhostComm::RHO);
  for_each_block_in_slice(
      get_lattice(), lower_corner, upper_corner,
      [&](auto &block, auto const &bci, auto const &ci,
          auto const &block_offset) {
        assert(density.size() == 2u * ci.numCells());

        auto rho_a_field =
            block.template getData<ScalarField>(m_rho_field_id[0]);
        auto rho_b_field =
            block.template getData<ScalarField>(m_rho_field_id[1]);
        // Collect values from input into block-local arrays
        std::vector<FloatType> values_a(bci.numCells());
        std::vector<FloatType> values_b(bci.numCells());

        auto kernel = [&values_a, &values_b, &density](
                          unsigned const block_index,
                          unsigned const local_index, Utils::Vector3i const &) {
          values_a[block_index] =
              numeric_cast<FloatType>(density[2u * local_index + 0u]);
          values_b[block_index] =
              numeric_cast<FloatType>(density[2u * local_index + 1u]);
        };

        copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);

        // field values in x-major, z-minor order
        ek::accessor::Scalar::set(rho_a_field, values_a, bci);
        ek::accessor::Scalar::set(rho_b_field, values_b, bci);
      });
}

template <typename FloatType, lbmpy::Arch Architecture>
std::vector<double>
LBWalberlaImplColorGradient<FloatType, Architecture>::get_slice_pressure_tensor(
    Utils::Vector3i const &lower_corner,
    Utils::Vector3i const &upper_corner) const {
  std::vector<double> out;
  for_each_block_in_slice(
      get_lattice(), lower_corner, upper_corner,
      [&](auto const &block, auto const &bci, auto const &ci,
          auto const &block_offset) {
        if (out.empty())
          out.resize(9u * ci.numCells());

        if constexpr (Architecture == lbmpy::Arch::CPU) {
          auto const pdf_field_a =
              block.template getData<PdfField>(m_pdf_field_id[0]);
          auto const pdf_field_b =
              block.template getData<PdfField>(m_pdf_field_id[1]);

          // Reads pre-collision populations via a neighbor-offset pull (see
          // raw_pressure_tensor_moment()), so the third copy_block_buffer
          // argument -- the node in block-local-offset space -- is needed
          // here (unlike other slice accessors, which only need the
          // pre-fetched bulk buffer indices).
          auto kernel = [&out, &block_offset, pdf_field_a, pdf_field_b, this](
                            unsigned const /* block_index */,
                            unsigned const local_index,
                            Utils::Vector3i const &node) {
            auto const local = node - block_offset;
            auto const cell = Cell{local[0], local[1], local[2]};
            auto const tensor = pressure_tensor_from_populations(
                pdf_field_a, pdf_field_b, cell);
            for (uint_t f = 0u; f < 9u; ++f) {
              out[9u * local_index + f] = double_c(tensor[f]);
            }
          };

          copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);
        }
#if defined(__CUDACC__) and defined(WALBERLA_BUILD_WITH_CUDA)
        else {
          // same neighbor-offset pull in the generated accessor, which adds
          // the zero-centered background (see corrected_pressure_tensor())
          std::array<std::vector<FloatType>, 2> values;
          for (std::size_t const component : {0u, 1u}) {
            values[component] = lbm::accessor::PressureTensor::get(
                block.template getData<PdfField>(m_pdf_field_id[component]),
                FloatType{1}, bci);
          }
          auto kernel = [&out, &values, this](unsigned const block_index,
                                              unsigned const local_index,
                                              Utils::Vector3i const &) {
            Matrix3<FloatType> tensor(FloatType{0});
            for (std::size_t const component : {0u, 1u}) {
              Matrix3<FloatType> moment;
              for (uint_t f = 0u; f < 9u; ++f) {
                moment[f] = values[component][9u * block_index + f];
              }
              tensor += corrected_pressure_tensor(moment, FloatType{1},
                                                  component);
            }
            for (uint_t f = 0u; f < 9u; ++f) {
              out[9u * local_index + f] = double_c(tensor[f]);
            }
          };

          copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);
        }
#endif
      });
  return out;
}

} // namespace walberla
