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
#ifndef OBSERVABLES_CYLINDRICALPROFILEOBSERVABLE_HPP
#define OBSERVABLES_CYLINDRICALPROFILEOBSERVABLE_HPP

#include "ProfileObservableBase.hpp"

#include <utils/Vector.hpp>

namespace Observables {

/** Cylindrical profile observable */
class CylindricalProfileObservable : virtual public ProfileObservableBase {
public:
  CylindricalProfileObservable(Utils::Vector3d const &center,
                               Utils::Vector3d const &axis, double min_r,
                               double max_r, double min_phi, double max_phi,
                               double min_z, double max_z, int n_r_bins,
                               int n_phi_bins, int n_z_bins)
      : ProfileObservableBase(min_r, max_r, min_phi, max_phi, min_z, max_z,
                              n_r_bins, n_phi_bins, n_z_bins),
        center(center), axis(axis) {}
  double min_r() const { return m_limits[0].first; }
  double max_r() const { return m_limits[0].second; }
  double min_phi() const { return m_limits[1].first; }
  double max_phi() const { return m_limits[1].second; }
  double min_z() const { return m_limits[2].first; }
  double max_z() const { return m_limits[2].second; }
  size_t n_r_bins() const { return m_bins[0]; }
  size_t n_phi_bins() const { return m_bins[1]; }
  size_t n_z_bins() const { return m_bins[2]; }
  Utils::Vector3d center;
  Utils::Vector3d axis;
};

} // Namespace Observables
#endif
