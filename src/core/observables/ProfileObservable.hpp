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
#ifndef OBSERVABLES_PROFILEOBSERVABLE_HPP
#define OBSERVABLES_PROFILEOBSERVABLE_HPP

#include "Observable.hpp"

#include <algorithm>
#include <array>
#include <iterator>
#include <vector>

namespace Observables {

namespace CoordSystem {
struct Cartesian {};
struct Cylindrical {};
struct Spherical {};
} // namespace CoordSystem

namespace {
/** Calculate the dimensionality of a profile observable.
 *  @param x, y, z  Dimensionality of the profile, use 0 for deleted dimensions
 *  @return A 3D, 2D or 1D vector with the number of bins.
 */
inline std::vector<size_t> profile_shape(size_t x, size_t y, size_t z) {
  std::array<size_t, 3> const dimensions = {x, y, z};
  std::vector<size_t> non_zero_dimensions;
  std::copy_if(dimensions.begin(), dimensions.end(),
               std::back_inserter(non_zero_dimensions),
               [](auto value) { return value != 0; });
  return non_zero_dimensions;
}
} // namespace

/** Profile observable
 *
 *  Measures a scalar field or a vector field in various coordinate systems
 *  (defaults to Cartesian). For measurements in 1D or 2D (e.g. an average
 *  along an axis or in a plane), set the number of bins on the degenerate
 *  axes to 0.
 *
 *  @tparam CoordinateSystem  Coordinate system of the profile, can be any
 *                            of the structs defined in @ref CoordSystem
 */
template <typename CoordinateSystem = CoordSystem::Cartesian>
class ProfileObservable : virtual public Observable {
public:
  ProfileObservable(double min_x, double max_x, double min_y, double max_y,
                    double min_z, double max_z, int n_x_bins, int n_y_bins,
                    int n_z_bins)
      : min_x(min_x), max_x(max_x), min_y(min_y), max_y(max_y), min_z(min_z),
        max_z(max_z), n_x_bins(static_cast<size_t>(n_x_bins)),
        n_y_bins(static_cast<size_t>(n_y_bins)),
        n_z_bins(static_cast<size_t>(n_z_bins)) {}
  double min_x, max_x;
  double min_y, max_y;
  double min_z, max_z;
  size_t n_x_bins, n_y_bins, n_z_bins;
  std::vector<size_t> shape() const override {
    assert(n_x_bins + n_y_bins + n_z_bins != 0);
    return profile_shape(n_x_bins, n_y_bins, n_z_bins);
  }
};

/** Cylindrical profile observable */
template <>
class ProfileObservable<CoordSystem::Cylindrical> : virtual public Observable {
public:
  ProfileObservable(Utils::Vector3d const &center, Utils::Vector3d const &axis,
                    double min_r, double max_r, double min_phi, double max_phi,
                    double min_z, double max_z, int n_r_bins, int n_phi_bins,
                    int n_z_bins)
      : center(center), axis(axis), min_r(min_r), max_r(max_r),
        min_phi(min_phi), max_phi(max_phi), min_z(min_z), max_z(max_z),
        n_r_bins(static_cast<size_t>(n_r_bins)),
        n_phi_bins(static_cast<size_t>(n_phi_bins)),
        n_z_bins(static_cast<size_t>(n_z_bins)) {}
  Utils::Vector3d center;
  Utils::Vector3d axis;
  double min_r, max_r;
  double min_phi, max_phi;
  double min_z, max_z;
  // Number of bins for each coordinate.
  size_t n_r_bins, n_phi_bins, n_z_bins;
  std::vector<size_t> shape() const override {
    assert(n_r_bins + n_phi_bins + n_z_bins != 0);
    return profile_shape(n_r_bins, n_phi_bins, n_z_bins);
  }
};

/** Spherical profile observable */
template <>
class ProfileObservable<CoordSystem::Spherical> : virtual public Observable {
public:
  ProfileObservable(double min_r, double max_r, double min_theta,
                    double max_theta, double min_phi, double max_phi,
                    int n_r_bins, int n_theta_bins, int n_phi_bins)
      : min_r(min_r), max_r(max_r), min_theta(min_theta), max_theta(max_theta),
        min_phi(min_phi), max_phi(max_phi),
        n_r_bins(static_cast<size_t>(n_r_bins)),
        n_theta_bins(static_cast<size_t>(n_theta_bins)),
        n_phi_bins(static_cast<size_t>(n_phi_bins)) {}
  double min_r, max_r;
  double min_theta, max_theta;
  double min_phi, max_phi;
  // Number of bins for each coordinate.
  size_t n_r_bins, n_theta_bins, n_phi_bins;
  std::vector<size_t> shape() const override {
    assert(n_r_bins + n_theta_bins + n_phi_bins != 0);
    return profile_shape(n_r_bins, n_theta_bins, n_phi_bins);
  }
};

} // Namespace Observables
#endif
