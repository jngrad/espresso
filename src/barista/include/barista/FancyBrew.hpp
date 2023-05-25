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

#pragma once

#include <string>
#include <vector>

namespace Barista {

class FancyBrew {
  std::string m_base;
  std::vector<std::string> m_toppings;

public:
  FancyBrew(std::string const &base);
  auto const &get_base() const { return m_base; }
  auto const &get_toppings() const { return m_toppings; }
  void set_toppings(std::vector<std::string> const &toppings);
  double run_cleaning_cycle() const;

  std::string brew_coffee() const {
    auto msg = "Here is your " + m_base;
    if (not m_toppings.empty()) {
      msg += " with ";
      for (auto const &topping : m_toppings) {
        msg += topping + ", ";
      }
      msg = msg.substr(0, msg.length() - 2);
    }
    return msg;
  }
};

} // namespace Barista
