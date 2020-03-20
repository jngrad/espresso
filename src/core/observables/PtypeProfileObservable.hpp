/*
 * Copyright (C) 2010-2020 The ESPResSo project
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
#ifndef OBSERVABLES_PTYPEPROFILEOBSERVABLE_HPP
#define OBSERVABLES_PTYPEPROFILEOBSERVABLE_HPP

#include "Observable.hpp"
#include "Particle.hpp"
#include "ProfileObservable.hpp"
#include "PtypeObservable.hpp"
#include "integrate.hpp"
#include <vector>

namespace Observables {

/** Distribution function which acts on a given list of particle types */
template <typename CoordinateSystem = CoordSystem::Cartesian>
class PtypeProfileObservable : public PtypeObservable,
                               public ProfileObservable<CoordinateSystem> {
public:
  template <class... Args>
  PtypeProfileObservable(std::vector<std::vector<int>> const &types,
                         int n_x_bins, int n_y_bins, int n_z_bins, double min_x,
                         double max_x, double min_y, double max_y, double min_z,
                         double max_z)
      : PtypeObservable(types), ProfileObservable<CoordinateSystem>(
                                    min_x, max_x, min_y, max_y, min_z, max_z,
                                    n_x_bins, n_y_bins, n_z_bins) {}
};

template <>
class PtypeProfileObservable<CoordSystem::Cylindrical>
    : public PtypeObservable,
      public ProfileObservable<CoordSystem::Cylindrical> {
public:
  PtypeProfileObservable(std::vector<std::vector<int>> const &types,
                         Utils::Vector3d const &center,
                         Utils::Vector3d const &axis, int n_x_bins,
                         int n_y_bins, int n_z_bins, double min_x, double max_x,
                         double min_y, double max_y, double min_z, double max_z)
      : PtypeObservable(types), ProfileObservable<CoordSystem::Cylindrical>(
                                    center, axis, min_x, max_x, min_y, max_y,
                                    min_z, max_z, n_x_bins, n_y_bins,
                                    n_z_bins) {}
};

} // Namespace Observables
#endif
