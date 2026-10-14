/*
 * Copyright (C) 2021-2026 The ESPResSo project
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
#include <config/config.hpp>

#ifdef ESPRESSO_WALBERLA

#include "LBFluid.hpp"
#include "LBWalberlaNodeState.hpp"
#include "WalberlaCheckpoint.hpp"

#include "core/BoxGeometry.hpp"
#include "core/lb/LBWalberla.hpp"
#include "core/lees_edwards/lees_edwards.hpp"
#include "core/lees_edwards/protocols.hpp"
#include "core/system/System.hpp"

#include <script_interface/code_info/CodeInfo.hpp>
#include <script_interface/communication.hpp>

#include <walberla_bridge/LatticeWalberla.hpp>
#include <walberla_bridge/lattice_boltzmann/LBWalberlaColorGradientBase.hpp>
#include <walberla_bridge/lattice_boltzmann/LeesEdwardsPack.hpp>
#include <walberla_bridge/lattice_boltzmann/lb_walberla_init.hpp>
#include <walberla_bridge/utils/ResourceManager.hpp>
#include <walberla_bridge/walberla_init.hpp>

#include <utils/Vector.hpp>
#include <utils/matrix.hpp>
#include <utils/mpi/reduce_optional.hpp>

#include <boost/mpi.hpp>
#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/collectives/broadcast.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace ScriptInterface::walberla {
namespace {
std::vector<double> get_kinematic_viscosity(VariantMap const &params) {
  auto const &visc_variant = params.at("kinematic_viscosity");
  if (is_type<double>(visc_variant)) {
    return {get_value<double>(params, "kinematic_viscosity")};
  }
  return get_value<std::vector<double>>(params, "kinematic_viscosity");
}
} // namespace

std::unordered_map<std::string, int> const LBVTKHandle::obs_map = {
    {"density", static_cast<int>(OutputVTK::density)},
    {"velocity_vector", static_cast<int>(OutputVTK::velocity_vector)},
    {"pressure_tensor", static_cast<int>(OutputVTK::pressure_tensor)},
    {"boundary", static_cast<int>(OutputVTK::boundary)},
    {"phasefield", static_cast<int>(OutputVTK::phasefield)},
    {"density_a", static_cast<int>(OutputVTK::density_a)},
    {"density_b", static_cast<int>(OutputVTK::density_b)},
};

Variant LBFluid::do_call_method(std::string const &name,
                                VariantMap const &params) {
  if (name == "activate") {
    auto &system = get_system();
    context()->parallel_try_catch(
        [&]() { lb_throw_if_expired(m_mpi_cart_comm_observer); });
    system.lb.set<::LB::LBWalberla>(m_instance, m_lb_params);
    system.lb.update_collision_model();
    m_is_active = true;
    return {};
  }
  if (name == "deactivate") {
    get_system().lb.reset();
    m_is_active = false;
    return {};
  }
  if (not name.starts_with("get_")) {
    context()->parallel_try_catch(
        [&]() { lb_throw_if_expired(m_mpi_cart_comm_observer); });
  }
  if (name == "add_force_at_pos") {
    auto const &box_geo = *get_system().box_geo;
    auto const pos = get_value<Utils::Vector3d>(params, "pos");
    auto const f = get_value<Utils::Vector3d>(params, "force");
    auto const folded_pos = box_geo.folded_position(pos);
    m_instance->add_force_at_pos(folded_pos * m_conv_dist, f * m_conv_force);
    return {};
  }
  if (name == "get_interpolated_velocity") {
    auto const pos = get_value<Utils::Vector3d>(params, "pos");
    return get_interpolated_velocity(pos);
  }
  if (name == "get_boundary_force_from_shape") {
    return get_boundary_force_from_shape(
        get_value<std::vector<int>>(params, "raster"));
  }
  if (name == "get_boundary_force") {
    return get_boundary_force();
  }
  if (name == "get_pressure_tensor") {
    return get_average_pressure_tensor();
  }
  if (name == "load_checkpoint") {
    auto const path = get_value<std::filesystem::path>(params, "path");
    auto const mode = get_value<int>(params, "mode");
    load_checkpoint(path, mode);
    return {};
  }
  if (name == "save_checkpoint") {
    auto const path = get_value<std::filesystem::path>(params, "path");
    auto const mode = get_value<int>(params, "mode");
    save_checkpoint(path, mode);
    return {};
  }
  if (name == "clear_boundaries") {
    m_instance->clear_boundaries();
    get_system().on_lb_boundary_conditions_change();
    return {};
  }
  if (name == "add_boundary_from_shape") {
    m_instance->update_boundary_from_shape(
        get_value<std::vector<int>>(params, "raster"),
        get_value<std::vector<double>>(params, "values"));
    return {};
  }
  if (name == "get_lattice_speed") {
    return 1. / m_conv_speed;
  }
  if (name == "init_two_component") {
    if (m_color_gradient == nullptr) {
      throw std::runtime_error(
          "this LB instance is single-component; 'init_two_component' is only "
          "valid for two-component (color-gradient) LB");
    }
    m_color_gradient->init_pdfs_from_components();
    m_instance->ghost_communication();
    return {};
  }

  return Base::do_call_method(name, params);
}

void LBFluid::make_instance(VariantMap const &params) {
  auto lb_visc = get_kinematic_viscosity(params);
  auto const dens = get_value<double>(params, "density");
  auto const gpu = get_value_or(params, "gpu", false);
  // the two-component model only runs in double precision on the GPU
  auto const precision =
      get_value_or(params, "single_precision", gpu and lb_visc.size() == 1u);
  auto const lb_lattice = m_lattice->lattice();
  for (auto &vi : lb_visc) {
    vi *= m_conv_visc;
  }
  auto const lb_dens = m_conv_dens * dens;
  using LBFactory = std::shared_ptr<LBWalberlaBase>(
      std::shared_ptr<::LatticeWalberla> const &, std::vector<double>, double,
      bool);
  LBFactory *make_new_instance = &new_lb_walberla_cpu;
  if (gpu) {
    std::vector<std::string> required_features;
    required_features.emplace_back("CUDA");
    CodeInfo::check_features(required_features);
#ifdef ESPRESSO_CUDA
    auto const blocks_per_mpi_rank = get_value<Utils::Vector3i>(
        m_lattice->get_parameter("blocks_per_mpi_rank"));
    if (blocks_per_mpi_rank != Utils::Vector3i{{1, 1, 1}}) {
      throw std::runtime_error(
          "Using more than one block per MPI rank is not supported for GPU LB");
    }
    make_new_instance = &new_lb_walberla_gpu;
#endif
  }
  m_instance = make_new_instance(lb_lattice, lb_visc, lb_dens, precision);
  m_color_gradient =
      dynamic_cast<LBWalberlaColorGradientBase *>(m_instance.get());
}

void LBFluid::do_construct(VariantMap const &params) {
  m_lattice = get_value<std::shared_ptr<LatticeWalberla>>(params, "lattice");
  m_vtk_writers =
      get_value_or<decltype(m_vtk_writers)>(params, "vtk_writers", {});
  auto const tau = get_value<double>(params, "tau");
  auto const agrid = get_value<double>(m_lattice->get_parameter("agrid"));
  auto const visc = get_kinematic_viscosity(params);
  auto const dens = get_value<double>(params, "density");
  auto const kT = get_value<double>(params, "kT");
  auto const ext_f = get_value<Utils::Vector3d>(params, "ext_force_density");
  m_lb_params = std::make_shared<::LB::LBWalberlaParams>(agrid, tau);
  m_is_active = false;
  auto const seed = get_value<int>(params, "seed");
  context()->parallel_try_catch([&]() {
    if (tau <= 0.) {
      throw std::domain_error("Parameter 'tau' must be > 0");
    }
    m_conv_dist = 1. / agrid;
    m_conv_visc = Utils::int_pow<1>(tau) / Utils::int_pow<2>(agrid);
    m_conv_energy = Utils::int_pow<2>(tau) / Utils::int_pow<2>(agrid);
    m_conv_dens = Utils::int_pow<3>(agrid);
    m_conv_speed = Utils::int_pow<1>(tau) / Utils::int_pow<1>(agrid);
    m_conv_press = Utils::int_pow<2>(tau) * Utils::int_pow<1>(agrid);
    m_conv_force = Utils::int_pow<2>(tau) / Utils::int_pow<1>(agrid);
    m_conv_force_dens = Utils::int_pow<2>(tau) * Utils::int_pow<2>(agrid);
    // surface tension is energy per area: [M T^-2]
    m_conv_sigma = Utils::int_pow<2>(tau);
    auto lb_visc = visc;
    for (auto &vi : lb_visc) {
      vi *= m_conv_visc;
    }
    auto const lb_dens = m_conv_dens * dens;
    auto const lb_kT = m_conv_energy * kT;
    auto const lb_ext_f = m_conv_force_dens * ext_f;
    if (seed < 0) {
      throw std::domain_error("Parameter 'seed' must be >= 0");
    }
    if (lb_kT < 0.) {
      throw std::domain_error("Parameter 'kT' must be >= 0");
    }
    if (lb_dens <= 0.) {
      throw std::domain_error("Parameter 'density' must be > 0");
    }
    for (auto const &vi : lb_visc) {
      if (vi < 0.) {
        throw std::domain_error("Parameter 'kinematic_viscosity' must be >= 0");
      }
    }
    // single precision doesn't conserve mass and momentum of the
    // two-component model well enough
    if (visc.size() != 1u and get_value_or(params, "gpu", false) and
        get_value_or(params, "single_precision", false)) {
      throw std::domain_error("single precision is not supported for the "
                              "two-component LB on the GPU");
    }
    make_instance(params);
    m_mpi_cart_comm_observer = ::walberla::get_mpi_cart_comm_observer();
    m_instance->set_collision_model(lb_kT, seed);
    m_sigma = get_value_or<double>(params, "sigma", 0.);
    m_beta = get_value_or<double>(params, "beta", 0.7);
    if (m_color_gradient) {
      m_color_gradient->set_collision_model_color_gradient(
          m_sigma * m_conv_sigma, m_beta);
    }
    // restore the RNG counter when deserializing a thermalized fluid;
    // must come after the collision model setup, which resets the counter
    if (lb_kT > 0. and params.contains("rng_state") and
        is_type<int>(params.at("rng_state"))) {
      auto const rng_state = get_value<int>(params, "rng_state");
      if (rng_state < 0) {
        throw std::domain_error("Parameter 'rng_state' must be >= 0");
      }
      m_instance->set_rng_state(static_cast<uint64_t>(rng_state));
    }
    m_instance->set_external_force(lb_ext_f);
    m_instance->ghost_communication();
    for (auto &vtk : m_vtk_writers) {
      vtk->attach_to_lattice(m_instance, get_lattice_to_md_units_conversion());
    }
  });
}

Variant
LBFluid::get_boundary_force_from_shape(std::vector<int> const &raster) const {
  auto const local =
      m_instance->get_boundary_force_from_shape(raster) / m_conv_force;
  return mpi_reduce_sum(context()->get_comm(), local);
}

Variant LBFluid::get_boundary_force() const {
  auto const local = m_instance->get_boundary_force() / m_conv_force;
  return mpi_reduce_sum(context()->get_comm(), local);
}

std::vector<Variant> LBFluid::get_average_pressure_tensor() const {
  auto const local = m_instance->get_pressure_tensor() / m_conv_press;
  auto const tensor_flat = mpi_reduce_sum(context()->get_comm(), local);
  auto tensor = Utils::Matrix<double, 3, 3>{};
  std::ranges::copy(tensor_flat, tensor.m_data.begin());
  return std::vector<Variant>{tensor.row<0>().as_vector(),
                              tensor.row<1>().as_vector(),
                              tensor.row<2>().as_vector()};
}

Variant LBFluid::get_interpolated_velocity(Utils::Vector3d const &pos) const {
  auto const &box_geo = *get_system().box_geo;
  auto const lb_pos = box_geo.folded_position(pos) * m_conv_dist;
  auto const result = m_instance->get_velocity_at_pos(lb_pos);
  return Utils::Mpi::reduce_optional(context()->get_comm(), result) /
         m_conv_speed;
}

/*
 * Checkpoint file layout (ASCII or binary, same field order):
 *   grid_size, pop_size (= n_components * stencil_size)
 *   two-component only: has_rng [, rng_state]
 *   per node in (i, j, k) order:
 *     populations (pop_size values)
 *     last_applied_force (one 3-vector per component)
 *     two-component only: density (2 values), phasefield, velocity
 *     is_boundary [, slip_velocity]
 * Single-component files are unchanged by the two-component extension.
 * The two-component derived fields are stored rather than recomputed, since
 * the particle coupling reads them before the next LB step.
 */

void LBFluid::load_checkpoint(std::filesystem::path const &path, int mode) {
  auto &lb_obj = *m_instance;
  auto *const cg_obj = m_color_gradient;
  auto const n_components = (cg_obj) ? std::size_t{2u} : std::size_t{1u};
  auto const pop_size = n_components * lb_obj.stencil_size();
  std::optional<uint64_t> rng_state;

  auto const read_metadata = [&lb_obj, cg_obj, pop_size,
                              &rng_state](CheckpointFile &cpfile) {
    auto const expected_grid_size = lb_obj.get_lattice().get_grid_dimensions();
    auto const expected_pop_size = pop_size;
    Utils::Vector3i read_grid_size;
    std::size_t read_pop_size;
    cpfile.read(read_grid_size);
    cpfile.read(read_pop_size);
    if (read_grid_size != expected_grid_size) {
      std::stringstream message;
      message << "grid dimensions mismatch, read [" << read_grid_size << "], "
              << "expected [" << expected_grid_size << "].";
      throw std::runtime_error(message.str());
    }
    if (read_pop_size != expected_pop_size) {
      throw std::runtime_error("population size mismatch, read " +
                               std::to_string(read_pop_size) + ", expected " +
                               std::to_string(expected_pop_size) + ".");
    }
    if (cg_obj) {
      int has_rng;
      cpfile.read(has_rng);
      if (has_rng) {
        uint64_t counter;
        cpfile.read(counter);
        rng_state = counter;
      }
    }
  };

  auto const read_data = [&lb_obj, cg_obj, pop_size,
                          n_components](CheckpointFile &cpfile) {
    auto const grid_size = lb_obj.get_lattice().get_grid_dimensions();
    auto const i_max = grid_size[0];
    auto const j_max = grid_size[1];
    auto const k_max = grid_size[2];
    LBWalberlaNodeState cpnode;
    cpnode.populations.resize(pop_size);
    cpnode.last_applied_force.resize(n_components);
    if (cg_obj) {
      cpnode.density.resize(2u);
    }
    for (int i = 0; i < i_max; i++) {
      for (int j = 0; j < j_max; j++) {
        for (int k = 0; k < k_max; k++) {
          auto const ind = Utils::Vector3i{{i, j, k}};
          cpfile.read(cpnode.populations);
          for (auto &force : cpnode.last_applied_force) {
            cpfile.read(force);
          }
          if (cg_obj) {
            cpfile.read(cpnode.density);
            cpfile.read(cpnode.phasefield);
            cpfile.read(cpnode.velocity);
          }
          cpfile.read(cpnode.is_boundary);
          if (cpnode.is_boundary) {
            cpfile.read(cpnode.slip_velocity);
          }
          lb_obj.set_node_population(ind, cpnode.populations);
          lb_obj.set_node_last_applied_force(ind, cpnode.last_applied_force);
          if (cg_obj) {
            lb_obj.set_node_density(ind, cpnode.density);
            cg_obj->set_node_phasefield(ind, cpnode.phasefield);
            cg_obj->set_node_velocity_raw(ind, cpnode.velocity);
          }
          if (cpnode.is_boundary) {
            lb_obj.set_node_velocity_at_boundary(ind, cpnode.slip_velocity);
          }
        }
      }
    }
  };

  auto const on_success = [&lb_obj, &rng_state]() {
    lb_obj.ghost_communication();
    lb_obj.reallocate_ubb_field();
    if (rng_state and lb_obj.get_rng_state()) {
      lb_obj.set_rng_state(*rng_state);
    }
  };

  load_checkpoint_common(*context(), "LB", path, mode, read_metadata, read_data,
                         on_success);
}

void LBFluid::save_checkpoint(std::filesystem::path const &path, int mode) {
  auto &lb_obj = *m_instance;
  auto *const cg_obj = m_color_gradient;
  auto const n_components = (cg_obj) ? std::size_t{2u} : std::size_t{1u};
  auto const pop_size = n_components * lb_obj.stencil_size();

  auto const write_metadata = [&lb_obj, cg_obj, pop_size,
                               mode](std::shared_ptr<CheckpointFile> cpfile_ptr,
                                     Context const &context) {
    auto const grid_size = lb_obj.get_lattice().get_grid_dimensions();
    if (context.is_head_node()) {
      cpfile_ptr->write(grid_size);
      cpfile_ptr->write(pop_size);
      if (cg_obj) {
        auto const rng_state = lb_obj.get_rng_state();
        cpfile_ptr->write(static_cast<int>(rng_state.has_value()));
        if (rng_state) {
          cpfile_ptr->write(*rng_state);
        }
      }
      unit_test_handle(mode);
    }
  };

  auto const on_failure = [](std::shared_ptr<CheckpointFile> const &,
                             Context const &context) {
    if (context.is_head_node()) {
      auto failure = true;
      boost::mpi::broadcast(context.get_comm(), failure, 0);
    }
  };

  auto const write_data = [&lb_obj, cg_obj,
                           mode](std::shared_ptr<CheckpointFile> cpfile_ptr,
                                 Context const &context) {
    auto const get_node_checkpoint =
        [&](Utils::Vector3i const &ind) -> std::optional<LBWalberlaNodeState> {
      auto const pop = lb_obj.get_node_population(ind);
      auto const laf = lb_obj.get_node_last_applied_force(ind);
      auto const lbb = lb_obj.get_node_is_boundary(ind);
      auto const vbb = lb_obj.get_node_velocity_at_boundary(ind);
      if (not(pop and laf and lbb and ((*lbb) ? vbb.has_value() : true))) {
        return std::nullopt;
      }
      LBWalberlaNodeState cpnode;
      cpnode.populations = *pop;
      cpnode.last_applied_force = *laf;
      if (cg_obj) {
        auto const dens = lb_obj.get_node_density(ind);
        auto const phi = cg_obj->get_node_phasefield(ind);
        auto const vel = cg_obj->get_node_velocity_raw(ind);
        if (not(dens and phi and vel)) {
          return std::nullopt;
        }
        cpnode.density = *dens;
        cpnode.phasefield = *phi;
        cpnode.velocity = *vel;
      } else {
        // SC stores a single force
        cpnode.last_applied_force.resize(1u);
      }
      cpnode.is_boundary = *lbb;
      if (*lbb) {
        cpnode.slip_velocity = *vbb;
      }
      return {cpnode};
    };

    auto failure = false;
    auto const &comm = context.get_comm();
    auto const is_head_node = context.is_head_node();
    auto const unit_test_mode = (mode != static_cast<int>(CptMode::ascii)) and
                                (mode != static_cast<int>(CptMode::binary));
    auto const grid_size = lb_obj.get_lattice().get_grid_dimensions();
    auto const i_max = grid_size[0];
    auto const j_max = grid_size[1];
    auto const k_max = grid_size[2];
    LBWalberlaNodeState cpnode;
    for (int i = 0; i < i_max; i++) {
      for (int j = 0; j < j_max; j++) {
        for (int k = 0; k < k_max; k++) {
          auto const ind = Utils::Vector3i{{i, j, k}};
          auto const result = get_node_checkpoint(ind);
          if (!unit_test_mode) {
            assert(1 == boost::mpi::all_reduce(comm, static_cast<int>(!!result),
                                               std::plus<>()) &&
                   "Incorrect number of return values");
          }
          if (is_head_node) {
            if (result) {
              cpnode = *result;
            } else {
              comm.recv(boost::mpi::any_source, 42, cpnode);
            }
            auto &cpfile = *cpfile_ptr;
            cpfile.write(cpnode.populations);
            for (auto const &force : cpnode.last_applied_force) {
              cpfile.write(force);
            }
            if (cg_obj) {
              cpfile.write(cpnode.density);
              cpfile.write(cpnode.phasefield);
              cpfile.write(cpnode.velocity);
            }
            cpfile.write(cpnode.is_boundary);
            if (cpnode.is_boundary) {
              cpfile.write(cpnode.slip_velocity);
            }
            boost::mpi::broadcast(comm, failure, 0);
          } else {
            if (result) {
              comm.send(0, 42, *result);
            }
            boost::mpi::broadcast(comm, failure, 0);
            if (failure) {
              return;
            }
          }
        }
      }
    }
  };

  save_checkpoint_common(*context(), "LB", path, mode, write_metadata,
                         write_data, on_failure);
}

} // namespace ScriptInterface::walberla

#endif // ESPRESSO_WALBERLA
