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
#include "RDF.hpp"

#include "fetch_particles.hpp"
#include "particle_data.hpp"

#include <utils/math/int_pow.hpp>

#include <boost/range/algorithm/transform.hpp>
#include <cmath>
#include <iostream>
#include <chrono>
using namespace std::chrono;

namespace Observables {
std::vector<double> RDF::operator()() const {
auto t1 = std::chrono::high_resolution_clock::now();
  std::vector<Particle> particles1 = fetch_particles(ids1());
  std::vector<Particle> particles2;
  if (!ids2().empty()) {
    particles2 = fetch_particles(ids2());
  }
  std::vector<std::pair<const Particle *const,const Particle *const>> particle_pairs_ptrs;
  bool const mixed_flag = !ids2().empty();
  for (auto it = particles1.begin(); it != particles1.end(); ++it) {
    for (auto jt = mixed_flag ? particles2.begin() : std::next(it),
              jend = mixed_flag ? particles2.end() : particles1.end();
         jt != jend; ++jt) {
      auto const p1 = std::addressof(*it), p2 = std::addressof(*jt);
      if (p1 != p2)
            particle_pairs_ptrs.emplace_back(std::pair<const Particle *const, const Particle *const>{p1 ,p2});
    }
  }
  auto ret = this->evaluate(particle_pairs_ptrs);
  auto t2 = std::chrono::high_resolution_clock::now();
  std::cout << "obs: " << std::chrono::duration_cast<std::chrono::microseconds>( t2 - t1 ).count() << " " << ids1().size() << " " << ids2().size() << std::endl;

  return ret;
}

std::vector<double>
RDF::evaluate(Utils::Span<const std::pair<const Particle *const, const Particle *const>> particle_pairs) const {
  auto const bin_width = (max_r - min_r) / static_cast<double>(n_r_bins);
  auto const inv_bin_width = 1.0 / bin_width;
  std::vector<double> res(n_values(), 0.0);
  long int cnt = 0;
  for (auto pair: particle_pairs) {
    auto const p1 = pair.first, p2 = pair.second;
    auto const dist = get_mi_vector(p1->r.p, p2->r.p, box_geo).norm();
    if (dist > min_r && dist < max_r) {
        auto const ind =
            static_cast<int>(std::floor((dist - min_r) * inv_bin_width));
        res[ind]++;
      }
      cnt++;
  }
  if (cnt == 0)
    return res;
  // normalization
  auto const volume = box_geo.volume();
  for (int i = 0; i < n_r_bins; ++i) {
    auto const r_in = i * bin_width + min_r;
    auto const r_out = r_in + bin_width;
    auto const bin_volume =
        (4.0 / 3.0) * Utils::pi() *
        (Utils::int_pow<3>(r_out) - Utils::int_pow<3>(r_in));
    res[i] *= volume / (bin_volume * static_cast<double>(cnt));
  }

  return res;
}
} // namespace Observables
