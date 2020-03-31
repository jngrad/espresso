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

#include "ProfileObservableBase.hpp"

namespace Observables {

/** Cartesian profile observable */
class ProfileObservable : virtual public ProfileObservableBase {
public:
  ProfileObservable(double min_x, double max_x, double min_y, double max_y,
                    double min_z, double max_z, int n_x_bins, int n_y_bins,
                    int n_z_bins)
      : ProfileObservableBase(min_x, max_x, min_y, max_y, min_z, max_z,
                              n_x_bins, n_y_bins, n_z_bins) {}
  double min_x() const { return m_limits[0].first; }
  double max_x() const { return m_limits[0].second; }
  double min_y() const { return m_limits[1].first; }
  double max_y() const { return m_limits[1].second; }
  double min_z() const { return m_limits[2].first; }
  double max_z() const { return m_limits[2].second; }
  size_t n_x_bins() const { return m_bins[0]; }
  size_t n_y_bins() const { return m_bins[1]; }
  size_t n_z_bins() const { return m_bins[2]; }
};

} // Namespace Observables
#endif
