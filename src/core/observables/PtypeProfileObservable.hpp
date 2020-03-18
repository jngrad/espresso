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
#ifndef OBSERVABLES_PTYPEPROFILEOBSERVABLE_HPP
#define OBSERVABLES_PTYPEPROFILEOBSERVABLE_HPP

#include "Observable.hpp"
#include "Particle.hpp"
#include "PtypeObservable.hpp"
#include "ProfileObservable.hpp"
#include "integrate.hpp"
#include <vector>

namespace Observables {

/** 1D distribution function which acts on a given list of particle types */
class PtypeProfileObservable1d : public PtypeObservable, public ProfileObservable1d {
public:
  PtypeProfileObservable1d(std::vector<std::vector<int>> const &types,
                           int n_x_bins, double min_x, double max_x)
      : PtypeObservable(types),
        ProfileObservable1d(min_x, max_x, n_x_bins) {}
};

/** 2D distribution function which acts on a given list of particle types */
class PtypeProfileObservable2d : public PtypeObservable, public ProfileObservable2d {
public:
  PtypeProfileObservable2d(std::vector<std::vector<int>> const &types,
                           int n_x_bins, int n_y_bins, double min_x,
                           double min_y, double max_x, double max_y)
      : PtypeObservable(types),
        ProfileObservable2d(min_x, max_x, min_y, max_y, n_x_bins, n_y_bins) {}
};

/** 3D distribution function which acts on a given list of particle types */
class PtypeProfileObservable3d : public PtypeObservable, public ProfileObservable3d {
public:
  PtypeProfileObservable3d(std::vector<std::vector<int>> const &types, int n_x_bins,
                           int n_y_bins, int n_z_bins, double min_x, double min_y,
                           double min_z, double max_x, double max_y, double max_z)
      : PtypeObservable(types),
        ProfileObservable3d(min_x, max_x, min_y, max_y, min_z, max_z, n_x_bins,
                            n_y_bins, n_z_bins) {}
};

} // Namespace Observables
#endif
