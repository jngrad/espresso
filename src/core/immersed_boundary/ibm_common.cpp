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

#include "ibm_common.hpp"

#include "cell_system/CellStructure.hpp"
#include "communication.hpp"
#include "system/System.hpp"

#include <utils/Vector.hpp>
#include <utils/serialization/optional.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>

#include <optional>
#include <stdexcept>

Utils::Vector3d get_ibm_particle_position(int pid) {
  auto &cell_structure = *System::get_system().cell_structure;
  auto *p = cell_structure.get_local_particle(pid);
  std::optional<Particle> opt_part{std::nullopt};

  if (p and not p->is_ghost()) {
    opt_part = *p;
  }
  opt_part = boost::mpi::all_reduce(comm_cart, opt_part,
                                    [](std::optional<Particle> const &acc,
                                       std::optional<Particle> const &item) {
                                      if (acc) {
                                        return acc;
                                      }
                                      return item;
                                    });
  if (opt_part)
    return opt_part.value().pos();
  throw std::runtime_error("Immersed Boundary: Particle not found");
}
