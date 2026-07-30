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

#include "ReactionAlgorithm.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/cell_system/CellSystem.hpp"

#include "core/Observable_stat.hpp"
#include "core/cell_system/CellStructure.hpp"
#include "core/communication.hpp"
#include "core/system/System.hpp"

#include <boost/mpi/collectives.hpp>
#include <boost/serialization/serialization.hpp>

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace ScriptInterface {
namespace ReactionMethods {

Variant ReactionAlgorithm::do_call_method(std::string const &name,
                                          VariantMap const &params) {
  if (name == "count_number_of_particles_per_type") {
    auto const &cs = get_value<std::shared_ptr<CellSystem::CellSystem>>(params, "cell_system")->get_cell_structure();
    auto const types = get_value<std::vector<int>>(params, "types");
    std::vector<int> local_numbers;
    std::vector<int> global_numbers(types.size());
    for (auto const &type : types) {
      if (type < 0) {
        throw std::runtime_error("Types may not be negative");
      }
      int counter = 0;
      for (auto const &p : cs.local_particles()) {
        if (p.type() == type) {
          counter++;
        }
      }
      local_numbers.emplace_back(counter);
    }
    boost::mpi::reduce(::comm_cart, local_numbers, global_numbers, std::plus(), 0);
    return global_numbers;
  }
  if (context()->is_head_node()) {
    throw std::runtime_error("unknown method '" + name + "()'");
  }
  return {};
}

} /* namespace ReactionMethods */
} /* namespace ScriptInterface */
