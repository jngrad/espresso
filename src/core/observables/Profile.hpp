/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#ifndef OBSERVABLES_PROFILE_HPP
#define OBSERVABLES_PROFILE_HPP

#include <array>
#include <boost/range/algorithm.hpp>
#include <utility>
#include <vector>

#include <utils/math/make_lin_space.hpp>

namespace Observables {

class Profile {
protected:
  /** Range of the profile. */
  std::array<std::pair<double, double>, 3> m_limits;
  /** Number of bins for each coordinate. */
  std::array<size_t, 3> m_bins;

public:
  Profile(double min_0, double max_0, double min_1, double max_1, double min_2,
          double max_2, int n_0_bins, int n_1_bins, int n_2_bins)
      : m_limits({{{min_0, max_0}, {min_1, max_1}, {min_2, max_2}}}),
        m_bins({{static_cast<size_t>(n_0_bins), static_cast<size_t>(n_1_bins),
                 static_cast<size_t>(n_2_bins)}}) {}

  /** Calculate the bin edges for each dimension */
  std::array<std::vector<double>, 3> edges() {
    std::array<std::vector<double>, 3> ret;
    for (int i = 0; i < 3; ++i) {
      boost::copy(Utils::make_lin_space(m_limits[i].first, m_limits[i].second,
                                        m_bins[i] + 1),
                  std::back_inserter(ret[i]));
    }
    return ret;
  }
};

} // Namespace Observables
#endif
