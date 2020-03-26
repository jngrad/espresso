/*
 * Copyright (C) 2010-2020 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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
#ifndef OBSERVABLES_RDF_HPP
#define OBSERVABLES_RDF_HPP

#include "Observable.hpp"

#include "partCfg_global.hpp"

#include <utils/math/int_pow.hpp>

#include <algorithm>
#include <vector>

namespace Observables {

class TypesObservable : public Observable {
public:
  TypesObservable(std::vector<int> const &p1_types, std::vector<int> const &p2_types,
      double min_r, double max_r, int n_r_bins)
      : p1_types(p1_types), p2_types(p2_types), min_r(min_r), max_r(max_r),
        n_r_bins(static_cast<size_t>(n_r_bins)) {}
  /// Types of particles to find the distribution for.
  std::vector<int> p1_types;
  /// Types of particles the others are distributed around.
  std::vector<int> p2_types;

  double min_r, max_r;
  size_t n_r_bins;

  virtual std::vector<double> evaluate() const = 0;

  std::vector<double> operator()() const final {
    return this->evaluate();
  }

  std::vector<size_t> shape() const override {
    return {n_r_bins};
  }
};

/** Calculate the radial distribution function.
 *
 *  Calculates the radial distribution function of particles with
 *  types given in the @p p1_types list around particles with types given
 *  in the @p p2_types list. The range is given by @p r_min and @p r_max and
 *  the distribution function is binned into @p r_bin bins, which are
 *  equidistant. The result is stored in the array @p rdf.
 *
 *  @param p1_types list with types of particles to find the distribution for.
 *  @param p2_types list with types of particles the others are distributed
 *                  around.
 */
class RDF : public TypesObservable {
public:
  RDF(std::vector<int> const &types_a, std::vector<int> const &types_b,
      double min_r, double max_r, int n_r_bins)
      : TypesObservable(order_unique(types_a), order_unique(types_b),
                        min_r, max_r, n_r_bins) {}

  std::vector<int> order_unique(std::vector<int> const &types_a) const {
    std::set<int> const p1_set(types_a.begin(), types_a.end());
    return std::vector<int>{p1_set.begin(), p1_set.end()};
  }

  std::vector<double> evaluate() const {
    bool const mixed = p1_types != p2_types;
    std::vector<const Particle *> p1_range;
    std::vector<const Particle *> p2_range;
    for (Particle const &p : partCfg()) {
      if(std::find(p1_types.begin(), p1_types.end(), p.p.type) != p1_types.end()) {
        p1_range.push_back(&p);
      }
      if (mixed) {
        if(std::find(p2_types.begin(), p2_types.end(), p.p.type) != p2_types.end()) {
          p2_range.push_back(&p);
        }
      }
    }
    if (!mixed) {
      p2_range = p1_range;
    }
    auto const bin_width = (max_r - min_r) / static_cast<double>(n_r_bins);
    auto const inv_bin_width = 1.0 / bin_width;
    std::vector<double> res(n_values(), 0.0);
    long int cnt = 0;
    //for (auto const *p1 : p1_range) {
    //  for (auto const *p2 : p2_range) {
    for (auto it = p1_range.begin(); it != p1_range.end(); ++it) {
      auto jt = (mixed)? p2_range.begin() : std::next(it);
      auto jt_end = (mixed)? p2_range.end() : p1_range.end();
      for (; jt != jt_end; ++jt) {
        //auto const dist = get_mi_vector(p1->r.p, p2->r.p, box_geo).norm();
        auto const dist = get_mi_vector((*it)->r.p, (*jt)->r.p, box_geo).norm();
        if (dist > min_r && dist < max_r) {
          auto const ind = static_cast<int>((dist - min_r) * inv_bin_width);
          res[ind]++;
        }
       cnt++;
      }
    }
    if (cnt == 0)
      return {};
    // normalization
    auto const volume = box_geo.volume();
    for (int i = 0; i < n_r_bins; ++i) {
      auto const r_in = i * bin_width + min_r;
      auto const r_out = r_in + bin_width;
      auto const bin_volume = (4.0 / 3.0) * Utils::pi() *
                              (Utils::int_pow<3>(r_out) - Utils::int_pow<3>(r_in));
      res[i] *= volume / (bin_volume * cnt);
    }

    return res;
  }
};

} // Namespace Observables

#endif
