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

#include "barista/FancyBrew.hpp"

#include <initializer_list>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace Barista {

FancyBrew::FancyBrew(std::string const &base) : m_base{base}, m_toppings{} {
  std::set<std::string> const allowed_bases = {"espresso", "ristretto",
                                               "lungo"};
  if (allowed_bases.count(base) != 1) {
    throw std::invalid_argument("Base not recognized: " + base);
  }
}

void FancyBrew::set_toppings(std::vector<std::string> const &toppings) {
  std::set<std::string> const allowed_toppings = {"milk", "sugar", "cream"};
  for (auto const &topping : toppings) {
    if (allowed_toppings.count(topping) != 1) {
      throw std::invalid_argument("Topping not recognized: " + topping);
    }
  }
  m_toppings = toppings;
}

double FancyBrew::run_cleaning_cycle() const {
  auto constexpr water_usage_extraction_chamber = 0.1;
  auto constexpr water_usage_milk_steamer = 0.2;
  auto water_used = water_usage_extraction_chamber;
  for (auto const &topping : m_toppings) {
    if (topping == "milk") {
      water_used += water_usage_milk_steamer;
      break;
    }
  }
  return water_used;
}

} // namespace Barista
