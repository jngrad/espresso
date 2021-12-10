/*
 * Copyright (C) 2021 The ESPResSo project
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
#ifndef SCRIPT_INTERFACE_WALBERLA_FLUID_ROUTINES_WALBERLA_HPP
#define SCRIPT_INTERFACE_WALBERLA_FLUID_ROUTINES_WALBERLA_HPP

#include "config.hpp"

#ifdef LB_WALBERLA

#include <walberla_bridge/LBWalberlaBase.hpp>

#include "core/communication.hpp"
#include "core/errorhandling.hpp"
#include "core/grid_based_algorithms/lb_interface.hpp"
#include "core/grid_based_algorithms/lb_walberla_instance.hpp"
#include "core/grid_based_algorithms/lb_walberla_interface.hpp"
#include "core/integrate.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"

#include <utils/Vector.hpp>
#include <utils/math/int_pow.hpp>

#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace ScriptInterface::walberla {

class LBFluidRoutinesWalberla : public AutoParameters<LBFluidRoutinesWalberla> {
  std::shared_ptr<::LBWalberlaBase> m_lb_fluid;
  Utils::Vector3i m_index;
  double m_conv_dens;
  double m_conv_press;
  double m_conv_force;
  double m_conv_velocity;

public:
  LBFluidRoutinesWalberla() {
    add_parameters(
        {{"_index", AutoParameter::read_only, [this]() { return (m_index); }},
         {"velocity",
          [this](const Variant &v) {
            auto const velocity =
                get_value<Utils::Vector3d>(v) * m_conv_velocity;
            Walberla::mpi_set_node_velocity(m_index, velocity);
          },
          [this]() {
            return Walberla::mpi_get_node_velocity(m_index) / m_conv_velocity;
          }},
         {"velocity_at_boundary",
          [this](const Variant &v) {
            if (is_none(v)) {
              Walberla::mpi_remove_node_from_boundary(m_index);
            } else {
              auto const velocity =
                  get_value<Utils::Vector3d>(v) * m_conv_velocity;
              Walberla::mpi_set_node_velocity_at_boundary(m_index, velocity);
            }
          },
          [this]() {
            if (Walberla::mpi_get_node_is_boundary(m_index)) {
              return Variant{
                  Walberla::mpi_get_node_velocity_at_boundary(m_index) /
                  m_conv_velocity};
            } else {
              return Variant{None{}};
            }
          }},
         {"density",
          [this](const Variant &v) {
            auto const density = get_value<double>(v) * m_conv_dens;
            Walberla::mpi_set_node_density(m_index, density);
          },
          [this]() {
            return Walberla::mpi_get_node_density(m_index) / m_conv_dens;
          }},
         {"_population",
          [this](const Variant &v) {
            auto const population = get_value<std::vector<double>>(v);
            Walberla::mpi_set_node_pop(m_index, population);
          },
          [this]() { return Walberla::mpi_get_node_pop(m_index); }},
         {"is_boundary", AutoParameter::read_only,
          [this]() { return Walberla::mpi_get_node_is_boundary(m_index); }},
         {"boundary_force", AutoParameter::read_only,
          [this]() {
            return Walberla::mpi_get_node_boundary_force(m_index) /
                   m_conv_force;
          }},
         {"_pressure_tensor", AutoParameter::read_only,
          [this]() {
            return Walberla::mpi_get_node_pressure_tensor(m_index) /
                   m_conv_press;
          }},
         {"last_applied_force",
          [this](const Variant &v) {
            auto const last_applied_force =
                get_value<Utils::Vector3d>(v) * m_conv_force;
            Walberla::mpi_set_node_last_applied_force(m_index,
                                                      last_applied_force);
          },
          [this]() {
            return Walberla::mpi_get_node_last_applied_force(m_index) /
                   m_conv_force;
          }}});
  }

  void do_construct(VariantMap const &params) override {
    try {
      m_lb_fluid = lb_walberla();
      auto const &lb_params = lb_walberla_params();
      auto const tau = lb_params->get_tau();
      auto const agrid = lb_params->get_agrid();
      auto const grid_size = m_lb_fluid->lattice().get_grid_dimensions();
      m_index = get_value<Utils::Vector3i>(params, "index");
      if (not(m_index < grid_size and m_index >= Utils::Vector3i{})) {
        throw std::out_of_range("Index error");
      }
      m_conv_dens = Utils::int_pow<3>(agrid);
      m_conv_press = Utils::int_pow<1>(agrid) * Utils::int_pow<2>(tau);
      m_conv_force = Utils::int_pow<2>(tau) / Utils::int_pow<1>(agrid);
      m_conv_velocity = Utils::int_pow<1>(tau) / Utils::int_pow<1>(agrid);
    } catch (const std::exception &e) {
      runtimeErrorMsg() << "LatticeWalberla failed: " << e.what();
      m_lb_fluid.reset();
    }
  }
};
} // namespace ScriptInterface::walberla

#endif // LB_WALBERLA
#endif
