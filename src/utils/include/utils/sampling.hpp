/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include "utils/Vector.hpp"
#include "utils/math/make_lin_space.hpp"
#include "utils/math/sqr.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <utility>
#include <vector>

namespace Utils {

/**
 * @brief Generate sampling positions for a cylindrical histogram.
 * @param r_limits Range for radial coordinate as std::pair.
 * @param phi_limits Range for azimuthal angle as std::pair.
 * @param z_limits Range for z coordinate as std::pair.
 * @param n_r_bins Number of bins in radial direction.
 * @param n_phi_bins Number of bins in azimuthal direction.
 * @param n_z_bins Number of bins in z direction.
 * @param sampling_density The number of samples per unit volume.
 * @retval Cartesian sampling coordinates.
 */
std::vector<Vector3d> get_cylindrical_sampling_positions(
    std::pair<double, double> const &r_limits,
    std::pair<double, double> const &phi_limits,
    std::pair<double, double> const &z_limits, std::size_t n_r_bins,
    std::size_t n_phi_bins, std::size_t n_z_bins, double sampling_density) {
  auto constexpr endpoint = false;
  auto const delta_r =
      (r_limits.second - r_limits.first) / static_cast<double>(n_r_bins);
  auto const delta_phi =
      (phi_limits.second - phi_limits.first) / static_cast<double>(n_phi_bins);

  // For the smallest bin we chose samples along the z-direction for a single
  // azimuthal angle per bin such that we fulfill the sampling density
  // requirement.
  auto const smallest_bin_volume =
      Utils::sqr(r_limits.first + delta_r) * delta_phi / 2.;
  auto const min_n_samples =
      std::max(n_z_bins, static_cast<std::size_t>(std::round(
                             smallest_bin_volume * sampling_density)));
  auto const delta_z =
      (z_limits.second - z_limits.first) / static_cast<double>(min_n_samples);

  std::vector<double> r_range;
  std::ranges::copy(make_lin_space(r_limits.first + .5 * delta_r,
                                   r_limits.second, n_r_bins, endpoint),
                    std::back_inserter(r_range));

  auto const phi_range =
      make_lin_space(phi_limits.first + .5 * delta_phi, phi_limits.second,
                     n_phi_bins, endpoint);
  auto const z_range = make_lin_space(z_limits.first + .5 * delta_z,
                                      z_limits.second, min_n_samples, endpoint);

  // Create the sampling positions for the innermost bin.
  std::vector<Vector3d> sampling_positions;
  for (auto const z : z_range) {
    for (auto const phi : phi_range) {
      sampling_positions.push_back(Vector3d{{r_range.front(), phi, z}});
    }
  }

  // Scale the number of samples for larger bins
  auto phis = [n_phi_bins, phi_limits](long r_bin) {
    auto const phis_range = make_lin_space(
        phi_limits.first, phi_limits.second,
        n_phi_bins * (static_cast<std::size_t>(r_bin) + 1u), endpoint);
    return phis_range;
  };
  // Calculate the sampling positions
  // Along z
  for (auto const z : z_range) {
    // Along r
    for (auto r_it = ++r_range.begin(); r_it != r_range.end(); ++r_it) {
      // Along phi
      for (auto const phi : phis(std::distance(r_range.begin(), r_it))) {
        sampling_positions.push_back(Vector3d{{*r_it, phi, z}});
      }
    }
  }

  return sampling_positions;
}

} // namespace Utils
