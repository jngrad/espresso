/*
 * Copyright (C) 2021-2023 The ESPResSo project
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

#include "Particle.hpp"
#include "cell_system/CellStructure.hpp"
#include "cells.hpp"
#include "system/System.hpp"

#include <boost/mpi/communicator.hpp>

#include <optional>

inline auto copy_particle_to_head_node(boost::mpi::communicator const &comm,
                                       System::System &system, int p_id) {
  std::optional<Particle> result{};
  auto p = system.cell_structure->get_local_particle(p_id);
  if (p and not p->is_ghost()) {
    if (comm.rank() == 0) {
      result = *p;
    } else {
      comm.send(0, p_id, *p);
    }
  }
  if (comm.rank() == 0 and not result) {
    Particle p{};
    comm.recv(boost::mpi::any_source, p_id, p);
    result = p;
  }
  return result;
}
