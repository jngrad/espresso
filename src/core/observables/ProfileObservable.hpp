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

namespace Observables {

/** 3D distribution function */
class ProfileObservable3d : virtual public Observable {
public:
  ProfileObservable3d(double min_x, double max_x, double min_y, double max_y,
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
    return {n_x_bins, n_y_bins, n_z_bins};
  }
};

using ProfileObservable = ProfileObservable3d;

/** 2D distribution function */
class ProfileObservable2d : virtual public Observable {
public:
  ProfileObservable2d(double min_x, double max_x, double min_y, double max_y,
                      int n_x_bins, int n_y_bins)
      : min_x(min_x), max_x(max_x), min_y(min_y), max_y(max_y),
        n_x_bins(static_cast<size_t>(n_x_bins)),
        n_y_bins(static_cast<size_t>(n_y_bins)) {}
  double min_x, max_x;
  double min_y, max_y;
  size_t n_x_bins, n_y_bins;
  std::vector<size_t> shape() const override {
    return {n_x_bins, n_y_bins};
  }
};

/** 1D distribution function */
class ProfileObservable1d : virtual public Observable {
public:
  ProfileObservable1d(double min_x, double max_x, int n_x_bins)
      : min_x(min_x), max_x(max_x), n_x_bins(static_cast<size_t>(n_x_bins)) {}
  double min_x, max_x;
  size_t n_x_bins;
  std::vector<size_t> shape() const override {
    return {n_x_bins};
  }
};

} // Namespace Observables
#endif
