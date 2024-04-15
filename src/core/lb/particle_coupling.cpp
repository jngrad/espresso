/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include "lb/particle_coupling.hpp"
#include "BoxGeometry.hpp"
#include "LocalBox.hpp"
#include "Particle.hpp"
#include "cell_system/CellStructure.hpp"
#include "communication.hpp"
#include "config/config.hpp"
#include "errorhandling.hpp"
#include "random.hpp"
#include "system/System.hpp"
#include "thermostat.hpp"

#include "Implementation.hpp"
#include "LBWalberla.hpp"

#include <utils/Counter.hpp>
#include <utils/Vector.hpp>

#include <boost/mpi.hpp>
#include <boost/serialization/optional.hpp>

#ifdef CALIPER
#include <caliper/cali.h>
#endif

#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <vector>

static Thermostat::GammaType lb_handle_particle_anisotropy(Particle const &p,
                                                           double lb_gamma) {
#ifdef THERMOSTAT_PER_PARTICLE
  auto const &partcl_gamma = p.gamma();
#ifdef PARTICLE_ANISOTROPY
  auto const default_gamma = Thermostat::GammaType::broadcast(lb_gamma);
#else
  auto const default_gamma = lb_gamma;
#endif // PARTICLE_ANISOTROPY
  return Thermostat::handle_particle_gamma(partcl_gamma, default_gamma);
#else
  return lb_gamma;
#endif // THERMOSTAT_PER_PARTICLE
}

Utils::Vector3d lb_drag_force(LB::Solver const &lb, double lb_gamma,
                              Particle const &p,
                              Utils::Vector3d const &shifted_pos,
                              Utils::Vector3d const &vel_offset) {
  /* calculate fluid velocity at particle's position
     this is done by linear interpolation (eq. (11) @cite ahlrichs99a) */
  auto const v_fluid = lb.get_coupling_interpolated_velocity(shifted_pos);
  auto const v_drift = v_fluid + vel_offset;
  auto const gamma = lb_handle_particle_anisotropy(p, lb_gamma);

  /* calculate viscous force (eq. (9) @cite ahlrichs99a) */
  return Utils::hadamard_product(gamma, v_drift - p.v());
}

Utils::Vector3d lb_drag_force(Particle const &p, double lb_gamma,
                              Utils::Vector3d const &v_fluid,
                              Utils::Vector3d const &vel_offset) {
  auto const v_drift = v_fluid + vel_offset;
  auto const gamma = lb_handle_particle_anisotropy(p, lb_gamma);

  /* calculate viscous force (eq. (9) @cite ahlrichs99a) */
  return Utils::hadamard_product(gamma, v_drift - p.v());
}

/**
 * @brief Check if a position is within the local box + halo.
 *
 * @param local_box Local geometry
 * @param pos Position to check
 * @param halo Halo
 *
 * @return True iff the point is inside of the box up to halo.
 */
static bool in_local_domain(LocalBox const &local_box,
                            Utils::Vector3d const &pos, double halo = 0.) {
  auto const halo_vec = Utils::Vector3d::broadcast(halo);
  auto const lower_corner = local_box.my_left() - halo_vec;
  auto const upper_corner = local_box.my_right() + halo_vec;

  return pos >= lower_corner and pos < upper_corner;
}

static bool in_box(Utils::Vector3d const &pos,
                   Utils::Vector3d const &lower_corner,
                   Utils::Vector3d const &upper_corner) {
  return pos >= lower_corner and pos < upper_corner;
}

bool in_local_halo(LocalBox const &local_box, Utils::Vector3d const &pos,
                   double agrid) {
  auto const halo = 0.5 * agrid;
  return in_local_domain(local_box, pos, halo);
}

static void positions_in_halo_impl(Utils::Vector3d const &pos_folded,
                                   Utils::Vector3d const &fully_inside_lower,
                                   Utils::Vector3d const &fully_inside_upper,
                                   Utils::Vector3d const &halo_lower_corner,
                                   Utils::Vector3d const &halo_upper_corner,
                                   BoxGeometry const &box_geo,
                                   std::vector<Utils::Vector3d> &res) {
  for (int i : {-1, 0, 1}) {
    for (int j : {-1, 0, 1}) {
      for (int k : {-1, 0, 1}) {
        Utils::Vector3d shift{{double(i), double(j), double(k)}};
        Utils::Vector3d pos_shifted =
            pos_folded + Utils::hadamard_product(box_geo.length(), shift);

        if (box_geo.type() == BoxType::LEES_EDWARDS) {
          auto le = box_geo.lees_edwards_bc();
          auto normal_shift = (pos_shifted - pos_folded)[le.shear_plane_normal];
          if (normal_shift > std::numeric_limits<double>::epsilon())
            pos_shifted[le.shear_direction] += le.pos_offset;
          if (normal_shift < -std::numeric_limits<double>::epsilon())
            pos_shifted[le.shear_direction] -= le.pos_offset;
        }

        if (in_box(pos_shifted, halo_lower_corner, halo_upper_corner)) {
          res.emplace_back(pos_shifted);
        }
      }
    }
  }
}

/**
 * @brief Return a vector of positions shifted by +,- box length in each
 * coordinate
 */
std::vector<Utils::Vector3d> positions_in_halo(Utils::Vector3d const &pos,
                                               BoxGeometry const &box_geo,
                                               LocalBox const &local_box,
                                               double agrid) {
  auto const halo = 0.5 * agrid;
  auto const halo_vec = Utils::Vector3d::broadcast(halo);
  auto const fully_inside_lower = local_box.my_left() + 2. * halo_vec;
  auto const fully_inside_upper = local_box.my_right() - 2. * halo_vec;
  auto const pos_folded = box_geo.folded_position(pos);
  if (in_box(pos_folded, fully_inside_lower, fully_inside_upper)) {
    return {pos_folded};
  }
  auto const halo_lower_corner = local_box.my_left() - halo_vec;
  auto const halo_upper_corner = local_box.my_right() + halo_vec;
  std::vector<Utils::Vector3d> res;
  positions_in_halo_impl(pos_folded, fully_inside_lower, fully_inside_upper,
                         halo_lower_corner, halo_upper_corner, box_geo, res);
  return res;
}

namespace LB {

Utils::Vector3d ParticleCoupling::get_noise_term(Particle const &p) const {
  if (not m_thermalized) {
    return Utils::Vector3d{};
  }
  using std::sqrt;
  using Utils::sqrt;
  auto const gamma = lb_handle_particle_anisotropy(p, m_thermostat.gamma);
  auto const noise = Random::noise_uniform<RNGSalt::PARTICLES>(
      m_thermostat.rng_counter(), m_thermostat.rng_seed(), p.id());
  return m_noise_pref_wo_gamma * Utils::hadamard_product(sqrt(gamma), noise);
}

void ParticleCoupling::commit(System::System const &system,
                              std::vector<Particle *> const &particles) {
  if (particles.empty()) {
    return;
  }
  auto const halo = 0.5 * m_lb.get_agrid();
  auto const halo_vec = Utils::Vector3d::broadcast(halo);
  auto const fully_inside_lower = m_local_box.my_left() + 2. * halo_vec;
  auto const fully_inside_upper = m_local_box.my_right() - 2. * halo_vec;
  auto const halo_lower_corner = m_local_box.my_left() - halo_vec;
  auto const halo_upper_corner = m_local_box.my_right() + halo_vec;
  std::vector<Utils::Vector3d> positions_velocity_coupling;
  std::vector<Utils::Vector3d> positions_force_coupling;
  std::vector<Utils::Vector3d> force_coupling_forces;
  std::vector<unsigned> positions_force_coupling_counter;
  std::vector<Particle *> coupled_particles;
  for (auto ptr : particles) {
    auto &p = *ptr;
    auto span_size = 1u;
    auto const folded_pos = m_box_geo.folded_position(p.pos());
    if (in_box(folded_pos, fully_inside_lower, fully_inside_upper)) {
      positions_force_coupling.emplace_back(folded_pos);
    } else {
      auto const old_size = positions_force_coupling.size();
      positions_in_halo_impl(folded_pos, fully_inside_lower, fully_inside_upper,
                             halo_lower_corner, halo_upper_corner, m_box_geo,
                             positions_force_coupling);
      auto const new_size = positions_force_coupling.size();
      span_size = static_cast<unsigned>(new_size - old_size);
    }
    auto found = false;
#ifdef ENGINE
    if (not p.swimming().is_engine_force_on_fluid)
#endif
      for (auto end = positions_force_coupling.end(), it = end - span_size;
           it != end; ++it) {
        auto const &pos = *it;
        if (pos >= halo_lower_corner and pos < halo_upper_corner) {
          positions_velocity_coupling.emplace_back(pos);
          found = true;
          break;
        }
      }
    if (found) {
      coupled_particles.emplace_back(ptr);
      positions_force_coupling_counter.emplace_back(span_size);
    } else {
      positions_force_coupling.erase(positions_force_coupling.end() - span_size,
                                     positions_force_coupling.end());
    }
  }

  if (positions_velocity_coupling.empty()) {
    return;
  }
  auto interpolated_velocities =
      m_lb.get_coupling_interpolated_velocities(positions_velocity_coupling);

  auto const &domain_lower_corner = m_local_box.my_left();
  auto const &domain_upper_corner = m_local_box.my_right();
  auto it_interpolated_velocities = interpolated_velocities.begin();
  auto it_positions_force_coupling = positions_force_coupling.begin();
  auto it_positions_force_coupling_counter =
      positions_force_coupling_counter.begin();
  for (auto ptr : coupled_particles) {
    auto &p = *ptr;
    Utils::Vector3d force_on_particle = {};
#ifdef ENGINE
    if (not p.swimming().is_engine_force_on_fluid)
#endif
    {
      auto const &v_fluid = *it_interpolated_velocities;
      auto const vel_offset = lb_drift_velocity_offset(p);
      auto const drag_force =
          lb_drag_force(p, m_thermostat.gamma, v_fluid, vel_offset);
      auto const random_force = get_noise_term(p);
      force_on_particle = drag_force + random_force;
      ++it_interpolated_velocities;
    }

    auto force_on_fluid = -force_on_particle;
#ifdef ENGINE
    if (p.swimming().is_engine_force_on_fluid) {
      force_on_fluid = p.calc_director() * p.swimming().f_swim;
    }
#endif

    auto const span_size = *it_positions_force_coupling_counter;
    ++it_positions_force_coupling_counter;
    for (auto i = 0u; i < span_size; ++i) {
      auto &pos = *it_positions_force_coupling;
      if (pos >= domain_lower_corner and pos < domain_upper_corner) {
        /* Particle is in our LB volume, so this node
         * is responsible to adding its force */
        p.force() += force_on_particle;
      }
      force_coupling_forces.emplace_back(force_on_fluid);
      ++it_positions_force_coupling;
    }
  }
  m_lb.add_forces_at_pos(positions_force_coupling, force_coupling_forces);
}

void ParticleCoupling::kernel(Particle &p) {
  auto const agrid = m_lb.get_agrid();

  // Calculate coupling force
  Utils::Vector3d force_on_particle = {};
  auto const halo_pos = positions_in_halo(m_box_geo.folded_position(p.pos()),
                                          m_box_geo, m_local_box, agrid);

#ifdef ENGINE
  if (not p.swimming().is_engine_force_on_fluid)
#endif
    for (auto const &pos : halo_pos) {
      if (in_local_halo(m_local_box, pos, agrid)) {
        auto const vel_offset = lb_drift_velocity_offset(p);
        auto const drag_force =
            lb_drag_force(m_lb, m_thermostat.gamma, p, pos, vel_offset);
        auto const random_force = get_noise_term(p);
        force_on_particle = drag_force + random_force;
        break;
      }
    }

  auto force_on_fluid = -force_on_particle;
#ifdef ENGINE
  if (p.swimming().is_engine_force_on_fluid) {
    force_on_fluid = p.calc_director() * p.swimming().f_swim;
  }
#endif

  // couple positions including shifts by one box length to add
  // forces to ghost layers
  for (auto const &pos : halo_pos) {
    if (in_local_domain(m_local_box, pos)) {
      /* Particle is in our LB volume, so this node
       * is responsible to adding its force */
      p.force() += force_on_particle;
    }
    m_lb.add_force_density(pos, force_on_fluid);
  }
}

#if defined(THERMOSTAT_PER_PARTICLE) and defined(PARTICLE_ANISOTROPY)
static void lb_coupling_sanity_checks(Particle const &p) {
  /*
  lb does (at the moment) not support rotational particle coupling.
  Consequently, anisotropic particles are also not supported.
  */
  auto const &p_gamma = p.gamma();
  if (p_gamma[0] != p_gamma[1] or p_gamma[1] != p_gamma[2]) {
    runtimeErrorMsg() << "anisotropic particle (id " << p.id()
                      << ") coupled to LB.";
  }
}
#endif

} // namespace LB

void System::System::lb_couple_particles(double time_step) {
#ifdef CALIPER
  CALI_CXX_MARK_FUNCTION;
#endif
  assert(thermostat->lb != nullptr);
  if (thermostat->lb->couple_to_md) {
    if (not lb.is_solver_set()) {
      runtimeErrorMsg() << "The LB thermostat requires a LB fluid";
      return;
    }
    auto const real_particles = cell_structure->local_particles();
    auto const ghost_particles = cell_structure->ghost_particles();
    LB::ParticleCoupling coupling{*thermostat->lb, lb, *box_geo, *local_geo,
                                  time_step};
    LB::CouplingBookkeeping bookkeeping{*cell_structure};
    std::vector<Particle *> particles{};
    auto const is_gpu = lb.is_gpu();
    for (auto const *particle_range : {&real_particles, &ghost_particles}) {
      for (auto &p : *particle_range) {
        if (not LB::is_tracer(p) and bookkeeping.should_be_coupled(p)) {
#if defined(THERMOSTAT_PER_PARTICLE) and defined(PARTICLE_ANISOTROPY)
          LB::lb_coupling_sanity_checks(p);
#endif
          if (is_gpu) {
            particles.emplace_back(&p);
          } else {
            coupling.kernel(p);
          }
        }
      }
    }
    if (is_gpu) {
      coupling.commit(*this, particles);
    }
  }
}
