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
 * Out-of-class node access definitions for
 * @ref walberla::LBWalberlaImplColorGradient.
 */

#include <utils/Vector.hpp>

#include <array>
#include <optional>
#include <utility>
#include <vector>

namespace walberla {

template <typename FloatType, lbmpy::Arch Architecture>
std::optional<Utils::Vector3d>
LBWalberlaImplColorGradient<FloatType, Architecture>::get_node_velocity(
    Utils::Vector3i const &node, bool consider_ghosts) const {
  assert(not(consider_ghosts and m_pending_ghost_comm.test(GhostComm::VEL)));
  assert(not(consider_ghosts and m_pending_ghost_comm.test(GhostComm::UBB)));
  if (m_has_boundaries) {
    auto const is_boundary = this->get_node_is_boundary(node, consider_ghosts);
    if (is_boundary and *is_boundary) {
      return this->get_node_velocity_at_boundary(node, consider_ghosts);
    }
  }
  auto const bc = get_block_and_cell(get_lattice(), node, consider_ghosts);
  if (!bc)
    return std::nullopt;

  auto field = bc->block->template uncheckedFastGetData<VectorField>(
      m_velocity_field_id);
  auto const vec = lbm::accessor::Vector::get(field, bc->cell);
  return to_vector3d(vec);
}

template <typename FloatType, lbmpy::Arch Architecture>
bool LBWalberlaImplColorGradient<FloatType, Architecture>::set_node_velocity(
    Utils::Vector3i const & /* node */, Utils::Vector3d const & /* v */) {
  throw std::runtime_error(
      "set_node_velocity is not supported for two-component LB. "
      "Set densities and populations instead to control the barycentric "
      "velocity.");
}

template <typename FloatType, lbmpy::Arch Architecture>
std::optional<Utils::Vector3d>
LBWalberlaImplColorGradient<FloatType, Architecture>::get_node_velocity_raw(
    Utils::Vector3i const &node, bool consider_ghosts) const {
  assert(not(consider_ghosts and m_pending_ghost_comm.test(GhostComm::VEL)));
  auto const bc = get_block_and_cell(get_lattice(), node, consider_ghosts);
  if (!bc)
    return std::nullopt;

  auto field = bc->block->template uncheckedFastGetData<VectorField>(
      m_velocity_field_id);
  return to_vector3d(lbm::accessor::Vector::get(field, bc->cell));
}

template <typename FloatType, lbmpy::Arch Architecture>
bool LBWalberlaImplColorGradient<FloatType, Architecture>::
    set_node_velocity_raw(Utils::Vector3i const &node,
                          Utils::Vector3d const &v) {
  m_pending_ghost_comm.set(GhostComm::VEL);
  auto bc = get_block_and_cell(get_lattice(), node, false);
  if (!bc)
    return false;

  auto vel_field =
      bc->block->template getData<VectorField>(m_velocity_field_id);
  lbm::accessor::Vector::set(vel_field, to_vector3<FloatType>(v), bc->cell);

  return true;
}

template <typename FloatType, lbmpy::Arch Architecture>
std::optional<std::vector<double>>
LBWalberlaImplColorGradient<FloatType, Architecture>::get_node_density(
    Utils::Vector3i const &node, bool consider_ghosts) const {
  assert(not(consider_ghosts and m_pending_ghost_comm.test(GhostComm::RHO)));
  auto bc = get_block_and_cell(get_lattice(), node, consider_ghosts);
  if (!bc)
    return std::nullopt;

  auto const rho_a_field =
      bc->block->template uncheckedFastGetData<ScalarField>(m_rho_field_id[0]);
  auto const rho_b_field =
      bc->block->template uncheckedFastGetData<ScalarField>(m_rho_field_id[1]);
  return std::vector<double>{
      double_c(ek::accessor::Scalar::get(rho_a_field, bc->cell)),
      double_c(ek::accessor::Scalar::get(rho_b_field, bc->cell))};
}

template <typename FloatType, lbmpy::Arch Architecture>
bool LBWalberlaImplColorGradient<FloatType, Architecture>::set_node_density(
    Utils::Vector3i const &node, std::vector<double> const &rho) {
  m_pending_ghost_comm.set(GhostComm::RHO);
  auto bc = get_block_and_cell(get_lattice(), node, false);
  if (!bc)
    return false;

  auto rho_a_field =
      bc->block->template getData<ScalarField>(m_rho_field_id[0]);
  auto rho_b_field =
      bc->block->template getData<ScalarField>(m_rho_field_id[1]);
  ek::accessor::Scalar::set(rho_a_field, FloatType_c(rho[0]), bc->cell);
  ek::accessor::Scalar::set(rho_b_field, FloatType_c(rho[1]), bc->cell);

  return true;
}

template <typename FloatType, lbmpy::Arch Architecture>
std::optional<double>
LBWalberlaImplColorGradient<FloatType, Architecture>::get_node_phasefield(
    Utils::Vector3i const &node, bool consider_ghosts) const {
  assert(not(consider_ghosts and m_pending_ghost_comm.test(GhostComm::PHI)));
  auto bc = get_block_and_cell(get_lattice(), node, consider_ghosts);
  if (!bc)
    return std::nullopt;

  auto const phasefield_field =
      bc->block->template uncheckedFastGetData<ScalarField>(m_phasefield_id);
  return double_c(ek::accessor::Scalar::get(phasefield_field, bc->cell));
}

template <typename FloatType, lbmpy::Arch Architecture>
bool LBWalberlaImplColorGradient<FloatType, Architecture>::set_node_phasefield(
    Utils::Vector3i const &node, double phasefield) {
  m_pending_ghost_comm.set(GhostComm::PHI);
  auto bc = get_block_and_cell(get_lattice(), node, false);
  if (!bc)
    return false;

  auto phasefield_field =
      bc->block->template getData<ScalarField>(m_phasefield_id);
  ek::accessor::Scalar::set(phasefield_field, FloatType_c(phasefield),
                            bc->cell);

  return true;
}

template <typename FloatType, lbmpy::Arch Architecture>
std::optional<std::vector<double>>
LBWalberlaImplColorGradient<FloatType, Architecture>::get_node_population(
    Utils::Vector3i const &node, bool consider_ghosts) const {
  assert(not(consider_ghosts and m_pending_ghost_comm.test(GhostComm::PDF)));
  auto bc = get_block_and_cell(get_lattice(), node, consider_ghosts);
  if (!bc)
    return std::nullopt;

  auto pdf_field_a = bc->block->template getData<PdfField>(m_pdf_field_id[0]);
  auto pdf_field_b = bc->block->template getData<PdfField>(m_pdf_field_id[1]);
  auto const pop_a = lbm::accessor::Population::get(pdf_field_a, bc->cell);
  auto const pop_b = lbm::accessor::Population::get(pdf_field_b, bc->cell);
  std::vector<double> population(2u * Stencil::Size);
  for (uint_t f = 0u; f < Stencil::Size; ++f) {
    population[f] = double_c(pop_a[f]);
    population[Stencil::Size + f] = double_c(pop_b[f]);
  }
  return {std::move(population)};
}

template <typename FloatType, lbmpy::Arch Architecture>
bool LBWalberlaImplColorGradient<FloatType, Architecture>::set_node_population(
    Utils::Vector3i const &node, std::vector<double> const &population) {
  m_pending_ghost_comm.set(GhostComm::PDF);
  auto bc = get_block_and_cell(get_lattice(), node, false);
  if (!bc)
    return false;

  // Only the populations are written: the SC overload that also updates the
  // velocity field assumes zero-centered single-component populations and
  // would store a wrong velocity. rho, phi and u re-derive on next integrate.
  auto pdf_field_a = bc->block->template getData<PdfField>(m_pdf_field_id[0]);
  auto pdf_field_b = bc->block->template getData<PdfField>(m_pdf_field_id[1]);
  std::array<FloatType, Stencil::Size> pop_a;
  std::array<FloatType, Stencil::Size> pop_b;
  for (uint_t f = 0u; f < Stencil::Size; ++f) {
    pop_a[f] = FloatType_c(population[f]);
    pop_b[f] = FloatType_c(population[Stencil::Size + f]);
  }
  lbm::accessor::Population::set(pdf_field_a, pop_a, bc->cell);
  lbm::accessor::Population::set(pdf_field_b, pop_b, bc->cell);

  return true;
}

template <typename FloatType, lbmpy::Arch Architecture>
std::optional<std::vector<Utils::Vector3d>>
LBWalberlaImplColorGradient<FloatType, Architecture>::
    get_node_force_to_be_applied(Utils::Vector3i const &node) const {
  auto const bc = get_block_and_cell(get_lattice(), node, true);
  if (!bc)
    return std::nullopt;
  auto const field_a =
      bc->block->template getData<VectorField>(m_force_to_be_applied_id[0]);
  auto const field_b =
      bc->block->template getData<VectorField>(m_force_to_be_applied_id[1]);
  return std::vector<Utils::Vector3d>{
      to_vector3d(lbm::accessor::Vector::get(field_a, bc->cell)),
      to_vector3d(lbm::accessor::Vector::get(field_b, bc->cell))};
}

template <typename FloatType, lbmpy::Arch Architecture>
std::optional<std::vector<Utils::Vector3d>>
LBWalberlaImplColorGradient<FloatType, Architecture>::
    get_node_last_applied_force(Utils::Vector3i const &node,
                                bool consider_ghosts) const {
  assert(not(consider_ghosts and m_pending_ghost_comm.test(GhostComm::LAF)));
  auto const bc = get_block_and_cell(get_lattice(), node, consider_ghosts);
  if (!bc)
    return std::nullopt;
  auto const field_a =
      bc->block->template getData<VectorField>(m_last_applied_force_field_id[0]);
  auto const field_b =
      bc->block->template getData<VectorField>(m_last_applied_force_field_id[1]);
  return std::vector<Utils::Vector3d>{
      to_vector3d(lbm::accessor::Vector::get(field_a, bc->cell)),
      to_vector3d(lbm::accessor::Vector::get(field_b, bc->cell))};
}

template <typename FloatType, lbmpy::Arch Architecture>
bool LBWalberlaImplColorGradient<FloatType, Architecture>::
    set_node_last_applied_force(Utils::Vector3i const &node,
                                std::vector<Utils::Vector3d> const &force) {
  m_pending_ghost_comm.set(GhostComm::LAF);
  auto bc = get_block_and_cell(get_lattice(), node, false);
  if (!bc)
    return false;

  // CG has no standalone velocity-set accessor (velocity is produced inside
  // the fused collide/stream kernels), so unlike SC we cannot use Force::set;
  // write the field directly. The stored velocity re-derives on next integrate.
  auto force_field_a =
      bc->block->template getData<VectorField>(m_last_applied_force_field_id[0]);
  auto force_field_b =
      bc->block->template getData<VectorField>(m_last_applied_force_field_id[1]);
  lbm::accessor::Vector::set(force_field_a, to_vector3<FloatType>(force[0]),
                             bc->cell);
  lbm::accessor::Vector::set(force_field_b, to_vector3<FloatType>(force[1]),
                             bc->cell);

  return true;
}

template <typename FloatType, lbmpy::Arch Architecture>
std::optional<Utils::VectorXd<9>>
LBWalberlaImplColorGradient<FloatType, Architecture>::get_node_pressure_tensor(
    Utils::Vector3i const &node) const {
  auto bc = get_block_and_cell(get_lattice(), node, false);
  if (!bc)
    return std::nullopt;

  if constexpr (Architecture == lbmpy::Arch::CPU) {
    auto const pdf_field_a =
        bc->block->template getData<PdfField>(m_pdf_field_id[0]);
    auto const pdf_field_b =
        bc->block->template getData<PdfField>(m_pdf_field_id[1]);
    return to_vector9d(
        pressure_tensor_from_populations(pdf_field_a, pdf_field_b, bc->cell));
  }
#if defined(__CUDACC__) and defined(WALBERLA_BUILD_WITH_CUDA)
  else {
    Matrix3<FloatType> tensor(FloatType{0});
    for (std::size_t const component : {0u, 1u}) {
      auto const pdf_field =
          bc->block->template getData<PdfField>(m_pdf_field_id[component]);
      tensor += corrected_pressure_tensor(
          lbm::accessor::PressureTensor::get(pdf_field, FloatType{1}, bc->cell),
          FloatType{1}, component);
    }
    return to_vector9d(tensor);
  }
#endif
}

} // namespace walberla
