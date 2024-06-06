/*
 * Copyright (C) 2024 The ESPResSo project
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

#define BOOST_TEST_MODULE cuda test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

void gpu_interface_test();

boost::test_tools::assertion_result has_gpu(boost::unit_test::test_unit_id);

BOOST_AUTO_TEST_CASE(gpu_interface, *boost::unit_test::precondition(has_gpu)) {
  gpu_interface_test();
}
