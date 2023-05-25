/*
 * Copyright (C) 2023 The ESPResSo project
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

#include "FancyBrew.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"
#include "script_interface/communication.hpp"

#include <memory>
#include <string>

namespace ScriptInterface::Barista {

FancyBrew::FancyBrew() {
  add_parameters({
      {"base", AutoParameter::read_only,
       [this]() { return Variant(m_instance->get_base()); }},
      {"toppings",
       [this](Variant const &v) {
         auto const toppings = get_value<std::vector<std::string>>(v);
         context()->parallel_try_catch(
             [&]() { m_instance->set_toppings(toppings); });
       },
       [this]() {
         return make_vector_of_variants(m_instance->get_toppings());
       }},
  });
}

void FancyBrew::do_construct(VariantMap const &params) {
  auto const base = get_value<std::string>(params, "base");
  context()->parallel_try_catch(
      [&]() { m_instance = std::make_shared<::Barista::FancyBrew>(base); });
  if (params.count("toppings") == 1) {
    do_set_parameter("toppings", params.at("toppings"));
  }
}

Variant FancyBrew::do_call_method(std::string const &name,
                                  VariantMap const &params) {
  if (name == "brew_coffee") {
    return m_instance->brew_coffee();
  }
  if (name == "set_toppings") {
    do_set_parameter("toppings", params.at("toppings"));
    return {};
  }
  if (name == "run_cleaning_cycle") {
    auto const local_water_used = m_instance->run_cleaning_cycle();
    return mpi_reduce_sum(context()->get_comm(), local_water_used);
  }
  return {};
}

} // namespace ScriptInterface::Barista
