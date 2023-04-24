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
 *  Hook procedures.
 *
 *  Implementation of event.hpp.
 */
#include "event.hpp"

#include "bonded_interactions/thermalized_bond.hpp"
#include "cell_system/CellStructureType.hpp"
#include "cells.hpp"
#include "collision.hpp"
#include "communication.hpp"
#include "config/config.hpp"
#include "cuda_init.hpp"
#include "cuda_interface.hpp"
#include "cuda_utils.hpp"
#include "electrostatics/coulomb.hpp"
#include "electrostatics/icc.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/electrokinetics.hpp"
#include "grid_based_algorithms/lb_boundaries.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "immersed_boundaries.hpp"
#include "integrate.hpp"
#include "interactions.hpp"
#include "magnetostatics/dipoles.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "partCfg_global.hpp"
#include "particle_node.hpp"
#include "thermostat.hpp"
#include "virtual_sites.hpp"

#ifdef ELECTROSTATICS
/** whether electrostatics actor has to be reinitialized on observable calc */
static bool reinit_electrostatics = false;
#endif
#ifdef DIPOLES
/** whether magnetostatics actor has to be reinitialized on observable calc */
static bool reinit_magnetostatics = false;
#endif

void on_program_start() {
#ifdef CUDA
  if (this_node == 0) {
    try {
      cuda_init();
    } catch (cuda_runtime_error const &) {
      // pass
    }
  }
#endif

  init_node_grid();

  /* initially go for regular decomposition */
  cells_re_init(CellStructureType::CELL_STRUCTURE_REGULAR);

  /* make sure interaction 0<->0 always exists */
  make_particle_type_exist(0);
}

void on_observable_calc() {
  /* Prepare particle structure: Communication step: number of ghosts and ghost
   * information */
  cells_update_ghosts(global_ghost_flags());
  update_dependent_particles();
#ifdef ELECTROSTATICS
  if (reinit_electrostatics) {
    Coulomb::on_observable_calc();
    reinit_electrostatics = false;
  }
#endif /* ELECTROSTATICS */

#ifdef DIPOLES
  if (reinit_magnetostatics) {
    Dipoles::on_observable_calc();
    reinit_magnetostatics = false;
  }
#endif /* DIPOLES */

#ifdef ELECTROKINETICS
  if (ek_initialized) {
    ek_integrate_electrostatics();
  }
#endif /* ELECTROKINETICS */

  clear_particle_node();
}

void on_particle_charge_change() {
#ifdef ELECTROSTATICS
  reinit_electrostatics = true;
#endif

  /* the particle information is no longer valid */
  partCfg().invalidate();
}

void on_particle_change() {
  if (cell_structure.decomposition_type() ==
      CellStructureType::CELL_STRUCTURE_HYBRID) {
    cell_structure.set_resort_particles(Cells::RESORT_GLOBAL);
  } else {
    cell_structure.set_resort_particles(Cells::RESORT_LOCAL);
  }
#ifdef ELECTROSTATICS
  reinit_electrostatics = true;
#endif
#ifdef DIPOLES
  reinit_magnetostatics = true;
#endif
  set_recalc_forces(true);

  /* the particle information is no longer valid */
  partCfg().invalidate();

  /* the particle information is no longer valid */
  invalidate_fetch_cache();
}

void on_coulomb_and_dipoles_change() {
#ifdef ELECTROSTATICS
  reinit_electrostatics = true;
  Coulomb::on_coulomb_change();
#endif
#ifdef DIPOLES
  reinit_magnetostatics = true;
  Dipoles::on_dipoles_change();
#endif
  on_short_range_ia_change();
}

void on_coulomb_change() {
#ifdef ELECTROSTATICS
  reinit_electrostatics = true;
  Coulomb::on_coulomb_change();
#endif
  on_short_range_ia_change();
}

void on_dipoles_change() {
#ifdef DIPOLES
  reinit_magnetostatics = true;
  Dipoles::on_dipoles_change();
#endif
  on_short_range_ia_change();
}

void on_non_bonded_ia_change() {
  maximal_cutoff_nonbonded();
  on_short_range_ia_change();
}

void on_short_range_ia_change() {
  cells_re_init(cell_structure.decomposition_type());
  set_recalc_forces(true);
}

void on_constraint_change() { set_recalc_forces(true); }

void on_lbboundary_change() {
#if defined(LB_BOUNDARIES) || defined(LB_BOUNDARIES_GPU)
  LBBoundaries::lb_init_boundaries();

  set_recalc_forces(true);
#endif
}

void on_boxl_change(bool skip_method_adaption) {
  grid_changed_box_l(box_geo);
  /* Electrostatics cutoffs mostly depend on the system size,
   * therefore recalculate them. */
  cells_re_init(cell_structure.decomposition_type());

  if (not skip_method_adaption) {
    /* Now give methods a chance to react to the change in box length */
#ifdef ELECTROSTATICS
    Coulomb::on_boxl_change();
#endif

#ifdef DIPOLES
    Dipoles::on_boxl_change();
#endif

    lb_lbfluid_init();
#ifdef LB_BOUNDARIES
    LBBoundaries::lb_init_boundaries();
#endif
  }
}

void on_cell_structure_change() {
  clear_particle_node();

  /* Now give methods a chance to react to the change in cell structure.
   * Most ES methods need to reinitialize, as they depend on skin,
   * node grid and so on. */
#ifdef ELECTROSTATICS
  Coulomb::on_cell_structure_change();
#endif

#ifdef DIPOLES
  Dipoles::on_cell_structure_change();
#endif
}

void on_temperature_change() { lb_lbfluid_reinit_parameters(); }

void on_periodicity_change() {
#ifdef ELECTROSTATICS
  Coulomb::on_periodicity_change();
#endif

#ifdef DIPOLES
  Dipoles::on_periodicity_change();
#endif

#ifdef STOKESIAN_DYNAMICS
  if (get_integrator().type == INTEG_METHOD_SD) {
    if (box_geo.periodic(0) || box_geo.periodic(1) || box_geo.periodic(2))
      runtimeErrorMsg() << "Stokesian Dynamics requires periodicity "
                        << "(False, False, False)\n";
  }
#endif
  on_skin_change();
}

void on_skin_change() {
  cells_re_init(cell_structure.decomposition_type());
  on_coulomb_and_dipoles_change();
}

void on_timestep_change() {
  lb_lbfluid_reinit_parameters();
  on_thermostat_param_change();
}

void on_forcecap_change() { set_recalc_forces(true); }

void on_node_grid_change() {
  grid_changed_n_nodes();
#ifdef ELECTROSTATICS
  Coulomb::on_node_grid_change();
#endif
#ifdef DIPOLES
  Dipoles::on_node_grid_change();
#endif
  cells_re_init(cell_structure.decomposition_type());
}

/**
 * @brief Returns the ghost flags required for running pair
 *        kernels for the global state, e.g. the force calculation.
 * @return Required data parts;
 */
unsigned global_ghost_flags() {
  /* Position and Properties are always requested. */
  unsigned data_parts = Cells::DATA_PART_POSITION | Cells::DATA_PART_PROPERTIES;

  if (lattice_switch == ActiveLB::CPU)
    data_parts |= Cells::DATA_PART_MOMENTUM;

  if (thermo_switch & THERMO_DPD)
    data_parts |= Cells::DATA_PART_MOMENTUM;

  if (n_thermalized_bonds) {
    data_parts |= Cells::DATA_PART_MOMENTUM;
    data_parts |= Cells::DATA_PART_BONDS;
  }

#ifdef COLLISION_DETECTION
  if (collision_params.mode != CollisionModeType::OFF) {
    data_parts |= Cells::DATA_PART_BONDS;
  }
#endif

  return data_parts;
}

void update_dependent_particles() {
#ifdef VIRTUAL_SITES
  virtual_sites()->update();
  cells_update_ghosts(global_ghost_flags());
#endif

#ifdef ELECTROSTATICS
  update_icc_particles();
#endif

  // Here we initialize volume conservation
  // This function checks if the reference volumes have been set and if
  // necessary calculates them
  immersed_boundaries.init_volume_conservation(cell_structure);
}
