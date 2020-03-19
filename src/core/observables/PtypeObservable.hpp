/*
 * Copyright (C) 2016-2020 The ESPResSo project
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

#ifndef OBSERVABLES_PTYPEOBSERVABLE_HPP
#define OBSERVABLES_PTYPEOBSERVABLE_HPP

#include "Observable.hpp"
#include "Particle.hpp"

#include <set>
#include <vector>

namespace Observables {

/** %Particle-based observable.
 *
 *  Base class for observables extracting raw data from particle subsets and
 *  returning either the data or a statistic derived from it.
 */
class PtypeObservable : virtual public Observable {
  /** Types of particles measured by this observable */
  std::vector<int> m_types;

  virtual std::vector<double>
  evaluate(std::vector<Utils::Span<const Particle *const>> particles) const = 0;

public:
  explicit PtypeObservable(std::vector<std::vector<int>> types_lists) {
    std::set<int> types;
    for (auto const &type_list : types_lists) {
      types.insert(type_list.begin(), type_list.end());
    }
    m_types = std::vector<int>(types.begin(), types.end());
  }
  std::vector<double> operator()() const final;

  std::vector<int> &types() { return m_types; }
  std::vector<int> const &types() const { return m_types; }
};

} // Namespace Observables
#endif
