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

#define BOOST_TEST_MODULE fft test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "p3m/fft.hpp"

#include <utils/Vector.hpp>

#include <array>
#include <cstddef>
#include <limits>
#include <optional>
#include <span>
#include <stdexcept>
#include <vector>

#if defined(P3M) || defined(DP3M)
std::optional<std::vector<int>> find_comm_groups(Utils::Vector3i const &,
                                                 Utils::Vector3i const &,
                                                 std::span<const int>,
                                                 std::span<int>, std::span<int>,
                                                 std::span<int>, int);

BOOST_AUTO_TEST_CASE(fft_find_comm_groups_mismatch) {
  int my_pos[3] = {0};
  int nodelist[4] = {0};
  int nodepos[12] = {0};
  int rank = 0;
  {
    auto const optional = find_comm_groups({0, 1, 2}, {1, 2, 3}, nodelist,
                                           nodelist, nodepos, my_pos, rank);
    BOOST_CHECK(not optional.has_value());
  }
  {
    auto const optional = find_comm_groups({3, 2, 1}, {2, 3, 1}, nodelist,
                                           nodelist, nodepos, my_pos, rank);
    BOOST_CHECK(not optional.has_value());
  }
  {
    auto const optional = find_comm_groups({2, 3, 1}, {3, 2, 1}, nodelist,
                                           nodelist, nodepos, my_pos, rank);
    BOOST_CHECK(not optional.has_value());
  }
}

BOOST_AUTO_TEST_CASE(fft_exceptions) {
  fft_allocator<int> allocator{};
  BOOST_CHECK_EQUAL(allocator.allocate(0ul), nullptr);
  BOOST_CHECK_THROW(
      allocator.allocate(std::numeric_limits<std::size_t>::max() / sizeof(int) +
                         1ul),
      std::bad_array_new_length);
}
#endif // defined(P3M) || defined(DP3M)
