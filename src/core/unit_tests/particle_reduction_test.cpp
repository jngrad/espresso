/*
 * Copyright (C) 2025 The ESPResSo project
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

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_MODULE particle_reduction test
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "Particle.hpp"
#include "cell_system/CellStructure.hpp"
#include "communication.hpp"
#include "particle_node.hpp"
#include "particle_reduction.hpp"
#include "system/System.hpp"

#include <utils/Vector.hpp>

#include <boost/mpi.hpp>

#include <functional>

#ifdef SHARED_MEMORY_PARALLELISM
BOOST_AUTO_TEST_CASE(make_kokkos_reduction_) {
  auto &system = System::get_system();
  auto const &cell_structure = *system.cell_structure;
  auto const &cells = cell_structure.decomposition().local_cells();

  // Particle setting
  ::make_new_particle(0, Utils::Vector3d{0., 0., 0.});
  Particle p;
  p.id() = 0;
#ifdef MASS
  p.mass() = 2.;
#endif
  p.v() = {3., 0., 0.};

  auto const kernel = [&cells](int i, Utils::Vector3d &res) {
    for (auto &p : cells[i]->particles()) {
      res += p.mass() * p.v();
    }
  };
  auto const reduce_op = [](Utils::Vector3d &a, Utils::Vector3d const &b) {
    a = a + b;
  };
  auto const reducer =
      Reduction::make_kokkos_reducer<Utils::Vector3d>(kernel, reduce_op);
  for (auto i = 0; i < cells.size(); ++i) {
    auto ref = Utils::Vector3d{0., 0., 0.};
    kernel(i, ref);
    auto res = Utils::Vector3d{0., 0., 0.};
    reducer(i, res);
    // The lambda `kernel(i, res)` is executed inside `reducer(i, res)`,
    // so make sure that both results are equal.
    BOOST_CHECK_EQUAL(ref, res);
  }
}

int main(int argc, char **argv) {
  auto const mpi_handle = MpiContainerUnitTest(argc, argv);

  auto system = System::System::create();
  system->set_box_l(Utils::Vector3d{10., 10., 10.});
  system->set_cell_structure_topology(CellStructureType::REGULAR);
  ::System::set_system(system);

  return boost::unit_test::unit_test_main(init_unit_test, argc, argv);
}
#else // SHARED_MEMORY_PARALLELISM
int main(int, char **) {}
#endif
