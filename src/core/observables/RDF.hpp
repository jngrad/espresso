/*
 * Copyright (C) 2020 The ESPResSo project
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

#include "PtypeProfileObservable.hpp"

#include <map>
#include <vector>

namespace Observables {

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
 *  @param rdf      Array to store the result (size: @p r_bins).
 */
class RDF : public PtypeProfileObservable1d {
public:
  RDF(std::vector<std::vector<int>> const &p_types,
      double min_r, double max_r, int n_r_bins)
      : PtypeProfileObservable1d(p_types, n_r_bins, min_r, max_r),
        p1_types(p_types[0]), p2_types(p_types[1]) {};
  /// Types of particles to find the distribution for.
  std::vector<int> p1_types;
  /// Types of particles the others are distributed around.
  std::vector<int> p2_types;
  
  
  std::vector<size_t> shape() const override { return {n_x_bins}; }

  std::vector<double>
  evaluate(std::vector<Utils::Span<const Particle *const>> particles) const override {
    std::map<int, int> types_map;
    for (int i = 0; i < types().size(); ++i) {
      types_map[types()[i]] = i;
    }
    auto const bin_width = (max_x - min_x) / static_cast<double>(n_x_bins);
    auto const inv_bin_width = 1.0 / bin_width;
    std::vector<double> res(n_values(), 0.0);
    long int cnt = 0;
    for (auto const p1_type: p1_types) {
      auto const &p1_range = particles[types_map[p1_type]];
      for (auto const p2_type: p2_types) {
        auto const &p2_range = particles[types_map[p2_type]];
        for (auto const *p1: p1_range) {
          for (auto const *p2: p2_range) {
            auto const dist = get_mi_vector(p1->r.p, p2->r.p, box_geo).norm();
            if (dist > min_x && dist < max_x) {
              auto const ind = static_cast<int>((dist - min_x) * inv_bin_width);
              res[ind]++;
            }
            cnt++;
          }
        }
      }
    }
    if (cnt == 0)
      return {};
    /* normalization */
    auto const volume = box_geo.volume();
    for (int i = 0; i < n_x_bins; ++i) {
      auto const r_in = i * bin_width + min_x;
      auto const r_out = r_in + bin_width;
      auto const bin_volume = (4.0 / 3.0) * Utils::pi() *
                              ((r_out * r_out * r_out) - (r_in * r_in * r_in));
      res[i] *= volume / (bin_volume * cnt);
    }
      
/*
  long int cnt = 0;
  int ind;
  bool mixed_flag = false;
  if (n_p1 == n_p2) {
    for (int i = 0; i < n_p1; i++)
      if (p1_types[i] != p2_types[i])
        mixed_flag = true;
  } else {
    mixed_flag = true;
  }

  auto const bin_width = (r_max - r_min) / (double)r_bins;
  auto const inv_bin_width = 1.0 / bin_width;
  for (int i = 0; i < r_bins; i++)
    rdf[i] = 0.0;
  // particle loop: p1_types
  for (auto it = partCfg.begin(); it != partCfg.end(); ++it) {
    for (int t1 = 0; t1 < n_p1; t1++) {
      if (it->p.type == p1_types[t1]) {
        // distinguish mixed and identical rdf's
        auto jt = mixed_flag ? partCfg.begin() : std::next(it);

        // particle loop: p2_types
        for (; jt != partCfg.end(); ++jt) {
          for (int t2 = 0; t2 < n_p2; t2++) {
            if (jt->p.type == p2_types[t2]) {
              auto const dist = get_mi_vector(it->r.p, jt->r.p, box_geo).norm();
              if (dist > r_min && dist < r_max) {
                ind = (int)((dist - r_min) * inv_bin_width);
                rdf[ind]++;
              }
              cnt++;
            }
          }
        }
      }
    }
  }*/
    
    
    return res;
  }
};

} // Namespace Observables

#endif
