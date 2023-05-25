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

#define BOOST_TEST_MODULE fancy brew
#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include "barista/FancyBrew.hpp"

#include <string>

BOOST_AUTO_TEST_CASE(brewing) {
  auto beverage = Barista::FancyBrew("lungo");
  BOOST_CHECK_EQUAL(beverage.get_base(), "lungo");
  BOOST_REQUIRE_EQUAL(beverage.get_toppings().size(), 0ul);

  beverage.set_toppings({std::string("milk"), std::string("sugar")});
  BOOST_REQUIRE_EQUAL(beverage.get_toppings().size(), 2ul);
  BOOST_CHECK_EQUAL(beverage.get_toppings()[0], "milk");
  BOOST_CHECK_EQUAL(beverage.get_toppings()[1], "sugar");

  beverage.set_toppings({std::string("sugar")});
  BOOST_REQUIRE_EQUAL(beverage.get_toppings().size(), 1ul);
  BOOST_CHECK_EQUAL(beverage.get_toppings()[0], "sugar");

  BOOST_CHECK_EQUAL(beverage.brew_coffee(), "Here is your lungo with sugar");
}
