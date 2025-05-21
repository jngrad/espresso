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

#define BOOST_TEST_MODULE ResourceCleanup test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "config/config.hpp"

#include "communication.hpp"
#include "cuda/utils.hpp"
#include "system/GpuParticleData.hpp"
#include "system/ResourceCleanup.hpp"
#include "system/System.hpp"
#include "system/System.impl.hpp"

#include <boost/mpi.hpp>

#include <memory>
#include <vector>

struct GlobalConfig : public MpiContainerUnitTest {
  GlobalConfig() {
  ::this_node=-2;
  }
  ~GlobalConfig() {
  ::System::reset_system();
  }
};

BOOST_TEST_GLOBAL_CONFIGURATION(GlobalConfig);
BOOST_AUTO_TEST_SUITE(suite)

BOOST_AUTO_TEST_CASE(checks) {
  BOOST_REQUIRE_EQUAL(0, 0);
}

BOOST_AUTO_TEST_SUITE_END()
