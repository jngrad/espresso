/*
 * Copyright (C) 2010-2023 The ESPResSo project
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
#ifndef ESPRESSO_SRC_CORE_INTEGRATOR_HPP
#define ESPRESSO_SRC_CORE_INTEGRATOR_HPP

/** \name Integrator switches */
/**@{*/
#define INTEG_METHOD_NPT_ISO 0
#define INTEG_METHOD_NVT 1
#define INTEG_METHOD_STEEPEST_DESCENT 2
#define INTEG_METHOD_BD 3
#define INTEG_METHOD_SD 7
/**@}*/

/** \name Integrator error codes */
/**@{*/
#define INTEG_ERROR_RUNTIME -1
#define INTEG_ERROR_SIGINT -2
/**@}*/

/** \name Integrator flags */
/**@{*/
/// recalculate forces unconditionally (mostly used for timing)
#define INTEG_REUSE_FORCES_NEVER -1
/// recalculate forces if @ref recalc_forces is set
#define INTEG_REUSE_FORCES_CONDITIONALLY 0
/// do not recalculate forces (mostly when reading checkpoints with forces)
#define INTEG_REUSE_FORCES_ALWAYS 1
/**@}*/

struct Integrator {
  /** Switch determining which integrator to use. */
  int type = INTEG_METHOD_NVT;
  int fluid_step = 0;
  /** Time step for the integration. */
double time_step = -1.;
  /** Verlet list skin. */
  double skin = 0.;
/** Average number of integration steps the Verlet list has been re-using. */
 double verlet_reuse = 0.;
/** Actual simulation time. */
  double sim_time = 0.0;
  /** If true, the forces will be recalculated before the next integration. */
  bool recalc_forces = true;
/** True iff the user has changed the skin setting. */
 bool skin_set = false;
/** Whether the thermostat has to be reinitialized before integration. */
bool reinit_thermo = true;
};

Integrator const &get_integrator();

#endif
