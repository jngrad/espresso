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

struct GlobalConfig  {
  std::shared_ptr<boost::mpi::environment> m_mpi_env;
  GlobalConfig() {
  ::this_node=-2;
    m_mpi_env = std::make_shared<boost::mpi::environment>(
        boost::unit_test::framework::master_test_suite().argc,
        boost::unit_test::framework::master_test_suite().argv,
        boost::mpi::threading::multiple);
    Communication::init(m_mpi_env);
  }
  ~GlobalConfig() {
      printf("%i: ~MpiContainerUnitTest() Communication::deinit() use_count=%li\n", this_node, m_mpi_env.use_count());
      Communication::deinit();
      printf("%i: ~MpiContainerUnitTest() m_mpi_env.reset()       use_count=%li\n", this_node, m_mpi_env.use_count());
      m_mpi_env.reset();
      printf("%i: ~MpiContainerUnitTest() ending                  use_count=%li\n", this_node, m_mpi_env.use_count());
  }
};

BOOST_TEST_GLOBAL_CONFIGURATION(GlobalConfig);
BOOST_AUTO_TEST_SUITE(suite)

BOOST_AUTO_TEST_CASE(checks) {
  BOOST_REQUIRE_EQUAL(0, 0);
}

BOOST_AUTO_TEST_SUITE_END()
