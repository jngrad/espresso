/*
 * Copyright (C) 2010-2022 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

/** \file
 *  Molecular dynamics integrator.
 *
 *  For more information about the integrator
 *  see \ref integrate.hpp "integrate.hpp".
 */

#include "integrate.hpp"
#include "integrators/brownian_inline.hpp"
#include "integrators/steepest_descent.hpp"
#include "integrators/stokesian_dynamics_inline.hpp"
#include "integrators/velocity_verlet_inline.hpp"
#include "integrators/velocity_verlet_npt.hpp"

#include "Integrator.hpp"
#include "ParticleRange.hpp"
#include "accumulators.hpp"
#include "bond_breakage/bond_breakage.hpp"
#include "bonded_interactions/rigid_bond.hpp"
#include "cells.hpp"
#include "collision.hpp"
#include "communication.hpp"
#include "cuda_interface.hpp"
#include "errorhandling.hpp"
#include "event.hpp"
#include "forces.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "grid_based_algorithms/lb_particle_coupling.hpp"
#include "interactions.hpp"
#include "lees_edwards/lees_edwards.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "npt.hpp"
#include "partCfg_global.hpp"
#include "particle_node.hpp"
#include "rattle.hpp"
#include "rotation.hpp"
#include "signalhandling.hpp"
#include "thermostat.hpp"
#include "virtual_sites.hpp"

#include <profiler/profiler.hpp>
#include <utils/mpi/all_compare.hpp>

#include <boost/mpi/collectives/reduce.hpp>
#include <boost/range/algorithm/min_element.hpp>

#include <mpi.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <csignal>
#include <functional>
#include <stdexcept>
#include <utility>

#ifdef VALGRIND_MARKERS
#include <callgrind.h>
#endif

static Integrator integrator{};

namespace {
volatile std::sig_atomic_t ctrl_C = 0;
} // namespace

namespace LeesEdwards {
/** @brief Currently active Lees-Edwards protocol. */
static std::shared_ptr<ActiveProtocol> protocol = nullptr;

/**
 * @brief Update the Lees-Edwards parameters of the box geometry
 * for the current simulation time.
 */
static void update_box_params() {
  if (box_geo.type() == BoxType::LEES_EDWARDS) {
    assert(protocol != nullptr);
    auto const sim_time = integrator.sim_time;
    box_geo.lees_edwards_update(get_pos_offset(sim_time, *protocol),
                                get_shear_velocity(sim_time, *protocol));
  }
}

void set_protocol(std::shared_ptr<ActiveProtocol> new_protocol) {
  box_geo.set_type(BoxType::LEES_EDWARDS);
  protocol = std::move(new_protocol);
  LeesEdwards::update_box_params();
  ::integrator.recalc_forces = true;
  cell_structure.set_resort_particles(Cells::RESORT_LOCAL);
}

void unset_protocol() {
  protocol = nullptr;
  box_geo.set_type(BoxType::CUBOID);
  ::integrator.recalc_forces = true;
  cell_structure.set_resort_particles(Cells::RESORT_LOCAL);
}

template <class Kernel> void run_kernel() {
  if (box_geo.type() == BoxType::LEES_EDWARDS) {
    auto const kernel = Kernel{box_geo};
    auto const particles = cell_structure.local_particles();
    std::for_each(particles.begin(), particles.end(),
                  [&kernel](auto &p) { kernel(p); });
  }
}
} // namespace LeesEdwards

void integrator_sanity_checks() {
  if (integrator.time_step < 0.0) {
    runtimeErrorMsg() << "time_step not set";
  }
  switch (integrator.type) {
  case INTEG_METHOD_STEEPEST_DESCENT:
    if (thermo_switch != THERMO_OFF)
      runtimeErrorMsg()
          << "The steepest descent integrator is incompatible with thermostats";
    break;
  case INTEG_METHOD_NVT:
    if (thermo_switch & (THERMO_NPT_ISO | THERMO_BROWNIAN | THERMO_SD))
      runtimeErrorMsg() << "The VV integrator is incompatible with the "
                           "currently active combination of thermostats";
    break;
#ifdef NPT
  case INTEG_METHOD_NPT_ISO:
    if (thermo_switch != THERMO_OFF and thermo_switch != THERMO_NPT_ISO)
      runtimeErrorMsg() << "The NpT integrator requires the NpT thermostat";
    if (box_geo.type() == BoxType::LEES_EDWARDS)
      runtimeErrorMsg() << "The NpT integrator cannot use Lees-Edwards";
    break;
#endif
  case INTEG_METHOD_BD:
    if (thermo_switch != THERMO_BROWNIAN)
      runtimeErrorMsg() << "The BD integrator requires the BD thermostat";
    break;
#ifdef STOKESIAN_DYNAMICS
  case INTEG_METHOD_SD:
    if (thermo_switch != THERMO_OFF and thermo_switch != THERMO_SD)
      runtimeErrorMsg() << "The SD integrator requires the SD thermostat";
    break;
#endif
  default:
    runtimeErrorMsg() << "Unknown value for integrator.type";
  }
}

static void resort_particles_if_needed(ParticleRange const &particles) {
  auto const offset = LeesEdwards::verlet_list_offset(
      box_geo, cell_structure.get_le_pos_offset_at_last_resort());
  if (cell_structure.check_resort_required(particles, integrator.skin,
                                           offset)) {
    cell_structure.set_resort_particles(Cells::RESORT_LOCAL);
  }
}

void on_thermostat_param_change() { integrator.reinit_thermo = true; }

/** called every time the simulation is continued/started, i.e.
 *  when switching from the script interface to the simulation core.
 */
static void on_integration_start() {
  /********************************************/
  /* sanity checks                            */
  /********************************************/

  integrator_sanity_checks();
#ifdef NPT
  integrator_npt_sanity_checks(integrator);
#endif
  long_range_interactions_sanity_checks();
  lb_lbfluid_sanity_checks(integrator.time_step);

  /********************************************/
  /* end sanity checks                        */
  /********************************************/

  lb_lbfluid_on_integration_start();

#ifdef CUDA
  MPI_Bcast(gpu_get_global_particle_vars_pointer_host(),
            sizeof(CUDA_global_part_vars), MPI_BYTE, 0, comm_cart);
#endif

  /* Prepare the thermostat */
  if (integrator.reinit_thermo) {
    thermo_init(integrator.time_step);
    integrator.reinit_thermo = false;
    integrator.recalc_forces = true;
  }

#ifdef NPT
  npt_ensemble_init(box_geo, integrator);
#endif

  partCfg().invalidate();
  invalidate_fetch_cache();

#ifdef ADDITIONAL_CHECKS
  if (!Utils::Mpi::all_compare(comm_cart, cell_structure.use_verlet_list)) {
    runtimeErrorMsg() << "Nodes disagree about use of verlet lists.";
  }
#ifdef ELECTROSTATICS
  {
    auto const &actor = electrostatics_actor;
    if (!Utils::Mpi::all_compare(comm_cart, static_cast<bool>(actor)) or
        (actor and !Utils::Mpi::all_compare(comm_cart, (*actor).which())))
      runtimeErrorMsg() << "Nodes disagree about Coulomb long-range method";
  }
#endif
#ifdef DIPOLES
  {
    auto const &actor = magnetostatics_actor;
    if (!Utils::Mpi::all_compare(comm_cart, static_cast<bool>(actor)) or
        (actor and !Utils::Mpi::all_compare(comm_cart, (*actor).which())))
      runtimeErrorMsg() << "Nodes disagree about dipolar long-range method";
  }
#endif
#endif /* ADDITIONAL_CHECKS */

  on_observable_calc();
}

/** @brief Calls the hook for propagation kernels before the force calculation
 *  @return whether or not to stop the integration loop early.
 */
static bool integrator_step_1(ParticleRange const &particles) {
  auto early_exit = false;
  switch (integrator.type) {
  case INTEG_METHOD_STEEPEST_DESCENT:
    early_exit = steepest_descent_step(particles);
    break;
  case INTEG_METHOD_NVT:
    velocity_verlet_step_1(particles, integrator.time_step);
    break;
#ifdef NPT
  case INTEG_METHOD_NPT_ISO:
    velocity_verlet_npt_step_1(particles, integrator.time_step);
    break;
#endif
  case INTEG_METHOD_BD:
    // the Ermak-McCammon's Brownian Dynamics requires a single step
    // so, just skip here
    break;
#ifdef STOKESIAN_DYNAMICS
  case INTEG_METHOD_SD:
    stokesian_dynamics_step_1(particles, integrator.time_step);
    break;
#endif // STOKESIAN_DYNAMICS
  default:
    throw std::runtime_error("Unknown value for integrator.type");
  }
  return early_exit;
}

/** Calls the hook of the propagation kernels after force calculation */
static void integrator_step_2(ParticleRange const &particles, double kT) {
  switch (integrator.type) {
  case INTEG_METHOD_STEEPEST_DESCENT:
    // Nothing
    break;
  case INTEG_METHOD_NVT:
    velocity_verlet_step_2(particles, integrator.time_step);
    break;
#ifdef NPT
  case INTEG_METHOD_NPT_ISO:
    velocity_verlet_npt_step_2(particles, integrator.time_step);
    break;
#endif
  case INTEG_METHOD_BD:
    // the Ermak-McCammon's Brownian Dynamics requires a single step
    brownian_dynamics_propagator(brownian, particles, integrator.time_step, kT);
    resort_particles_if_needed(particles);
    break;
#ifdef STOKESIAN_DYNAMICS
  case INTEG_METHOD_SD:
    // Nothing
    break;
#endif // STOKESIAN_DYNAMICS
  default:
    throw std::runtime_error("Unknown value for integrator.type");
  }
}

int integrate(int n_steps, int reuse_forces) {
  ESPRESSO_PROFILER_CXX_MARK_FUNCTION;

  // Prepare particle structure and run sanity checks of all active algorithms
  on_integration_start();

  // If any method vetoes (e.g. P3M not initialized), immediately bail out
  if (check_runtime_errors(comm_cart))
    return INTEG_ERROR_RUNTIME;

  // Additional preparations for the first integration step
  if (reuse_forces == INTEG_REUSE_FORCES_NEVER or
      (integrator.recalc_forces and
       reuse_forces != INTEG_REUSE_FORCES_ALWAYS)) {
    ESPRESSO_PROFILER_MARK_BEGIN("Initial Force Calculation");
    lb_lbcoupling_deactivate();

#ifdef VIRTUAL_SITES
    virtual_sites()->update();
#endif

    // Communication step: distribute ghost positions
    cells_update_ghosts(global_ghost_flags());

    force_calc(cell_structure, integrator, temperature);

    if (integrator.type != INTEG_METHOD_STEEPEST_DESCENT) {
#ifdef ROTATION
      convert_initial_torques(cell_structure.local_particles());
#endif
    }

    ESPRESSO_PROFILER_MARK_END("Initial Force Calculation");
  }

  lb_lbcoupling_activate();

  if (check_runtime_errors(comm_cart))
    return INTEG_ERROR_RUNTIME;

  // Keep track of the number of Verlet updates (i.e. particle resorts)
  int n_verlet_updates = 0;

  // Keep track of whether an interrupt signal was caught (only in singleton
  // mode, since signal handlers are unreliable with more than 1 MPI rank)
  auto const singleton_mode = comm_cart.size() == 1;
  auto caught_sigint = false;
  auto caught_error = false;

#ifdef VALGRIND_MARKERS
  CALLGRIND_START_INSTRUMENTATION;
#endif
  // Integration loop
  ESPRESSO_PROFILER_CXX_MARK_LOOP_BEGIN(integration_loop, "Integration loop");
  int integrated_steps = 0;
  for (int step = 0; step < n_steps; step++) {
    ESPRESSO_PROFILER_CXX_MARK_LOOP_ITERATION(integration_loop, step);

    auto particles = cell_structure.local_particles();

#ifdef BOND_CONSTRAINT
    if (n_rigidbonds)
      save_old_position(particles, cell_structure.ghost_particles());
#endif

    LeesEdwards::update_box_params();
    auto const early_exit = integrator_step_1(particles);
    if (early_exit)
      break;

    LeesEdwards::run_kernel<LeesEdwards::Push>();

#ifdef NPT
    if (integrator.type != INTEG_METHOD_NPT_ISO)
#endif
    {
      resort_particles_if_needed(particles);
    }

    // Propagate philox RNG counters
    philox_counter_increment();

#ifdef BOND_CONSTRAINT
    // Correct particle positions that participate in a rigid/constrained bond
    if (n_rigidbonds) {
      correct_position_shake(cell_structure);
    }
#endif

#ifdef VIRTUAL_SITES
    virtual_sites()->update();
#endif

    if (cell_structure.get_resort_particles() >= Cells::RESORT_LOCAL)
      n_verlet_updates++;

    // Communication step: distribute ghost positions
    cells_update_ghosts(global_ghost_flags());

    particles = cell_structure.local_particles();

    force_calc(cell_structure, integrator, temperature);

#ifdef VIRTUAL_SITES
    virtual_sites()->after_force_calc();
#endif
    integrator_step_2(particles, temperature);
    LeesEdwards::run_kernel<LeesEdwards::UpdateOffset>();
#ifdef BOND_CONSTRAINT
    // SHAKE velocity updates
    if (n_rigidbonds) {
      correct_velocity_shake(cell_structure);
    }
#endif

    // propagate one-step functionalities
    if (integrator.type != INTEG_METHOD_STEEPEST_DESCENT) {
      if (lb_lbfluid_get_lattice_switch() != ActiveLB::NONE) {
        auto const tau = lb_lbfluid_get_tau();
        auto const lb_steps_per_md_step =
            static_cast<int>(std::round(tau / integrator.time_step));
        integrator.fluid_step += 1;
        if (integrator.fluid_step >= lb_steps_per_md_step) {
          integrator.fluid_step = 0;
          lb_lbfluid_propagate();
        }
        lb_lbcoupling_propagate();
      }

#ifdef VIRTUAL_SITES
      virtual_sites()->after_lb_propagation(integrator.time_step);
#endif

#ifdef COLLISION_DETECTION
      handle_collisions();
#endif
      BondBreakage::process_queue();
    }

    integrated_steps++;

    if (check_runtime_errors(comm_cart)) {
      caught_error = true;
      break;
    }

    // Check if SIGINT has been caught.
    if (singleton_mode and ctrl_C == 1) {
      caught_sigint = true;
      break;
    }

  } // for-loop over integration steps
  LeesEdwards::update_box_params();
  ESPRESSO_PROFILER_CXX_MARK_LOOP_END(integration_loop);

#ifdef VALGRIND_MARKERS
  CALLGRIND_STOP_INSTRUMENTATION;
#endif

#ifdef VIRTUAL_SITES
  virtual_sites()->update();
#endif

  // Verlet list statistics
  if (n_verlet_updates > 0)
    integrator.verlet_reuse = n_steps / static_cast<double>(n_verlet_updates);
  else
    integrator.verlet_reuse = 0.;

#ifdef NPT
  if (integrator.type == INTEG_METHOD_NPT_ISO) {
    synchronize_npt_state();
  }
#endif
  if (caught_sigint) {
    ctrl_C = 0;
    return INTEG_ERROR_SIGINT;
  }
  if (caught_error) {
    return INTEG_ERROR_RUNTIME;
  }
  return integrated_steps;
}

int integrate_with_signal_handler(int n_steps, int reuse_forces,
                                  bool update_accumulators) {

  assert(n_steps >= 0);

  // Override the signal handler so that the integrator obeys Ctrl+C
  SignalHandler sa(SIGINT, [](int) { ctrl_C = 1; });

  if (not update_accumulators or n_steps == 0) {
    return integrate(n_steps, reuse_forces);
  }

  auto const is_head_node = comm_cart.rank() == 0;

  /* if skin wasn't set, do an educated guess now */
  if (!integrator.skin_set) {
    auto const max_cut = maximal_cutoff(n_nodes);
    if (max_cut <= 0.0) {
      if (is_head_node) {
        throw std::runtime_error(
            "cannot automatically determine skin, please set it manually");
      }
      return INTEG_ERROR_RUNTIME;
    }
    /* maximal skin that can be used without resorting is the maximal
     * range of the cell system minus what is needed for interactions. */
    auto const max_range = *boost::min_element(::cell_structure.max_cutoff());
    auto const new_skin = std::min(0.4 * max_cut, max_range - max_cut);
    ::set_skin(new_skin);
  }

  // re-acquire MpiCallbacks listener on worker nodes
  if (not is_head_node) {
    return 0;
  }

  using Accumulators::auto_update;
  using Accumulators::auto_update_next_update;

  for (int i = 0; i < n_steps;) {
    /* Integrate to either the next accumulator update, or the
     * end, depending on what comes first. */
    auto const steps = std::min((n_steps - i), auto_update_next_update());
    auto const retval = mpi_call(Communication::Result::main_rank, integrate,
                                 steps, reuse_forces);
    if (retval < 0) {
      return retval; // propagate error code
    }

    reuse_forces = INTEG_REUSE_FORCES_ALWAYS;

    auto_update(steps);

    i += steps;
  }

  return 0;
}

REGISTER_CALLBACK_MAIN_RANK(integrate)

double interaction_range() {
  /* Consider skin only if there are actually interactions */
  auto const max_cut = maximal_cutoff(n_nodes == 1);
  return (max_cut > 0.) ? max_cut + integrator.skin : INACTIVE_CUTOFF;
}

double get_time_step() { return integrator.time_step; }

double get_sim_time() { return integrator.sim_time; }

void increment_sim_time(double amount) { integrator.sim_time += amount; }

void set_time_step(double value) {
  if (value <= 0.)
    throw std::domain_error("time_step must be > 0.");
  if (lb_lbfluid_get_lattice_switch() != ActiveLB::NONE)
    check_tau_time_step_consistency(lb_lbfluid_get_tau(), value);
  ::integrator.time_step = value;
  on_timestep_change();
}

void set_skin(double value) {
  ::integrator.skin = value;
  ::integrator.skin_set = true;
  on_skin_change();
}

void set_time(double value) {
  ::integrator.sim_time = value;
  ::integrator.recalc_forces = true;
  LeesEdwards::update_box_params();
}

void set_integ_switch(int value) {
  ::integrator.type = value;
  ::integrator.recalc_forces = true;
}

void set_recalc_forces(bool value) { ::integrator.recalc_forces = value; }

Integrator const &get_integrator() { return ::integrator; }
