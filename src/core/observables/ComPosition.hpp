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
#ifndef OBSERVABLES_COMPOSITION_HPP
#define OBSERVABLES_COMPOSITION_HPP

#include "PidObservable.hpp"

#include <numeric>
#include <vector>

namespace Observables {
template <typename ParticleTrait, typename ParticleTraitType>
class ComReduce : public PidObservable {
  ParticleTrait particle_trait;
public:
  using PidObservable::PidObservable;
  int n_values() const override { return 3; }

  explicit ComReduce(std::vector<int> ids, ParticleTrait &&particle_trait):
    PidObservable(std::move(ids)), particle_trait(particle_trait) {}

  std::vector<double>
  evaluate(Utils::Span<const Particle *const> particles) const override {
    double total_mass;
    ParticleTraitType acc;
    std::tie(total_mass, acc) =
        std::accumulate(particles.begin(), particles.end(),
                        std::make_pair(0.0, ParticleTraitType{}),
                        particle_trait);
    return acc / total_mass;
  };
};

auto particle_trait = [](std::pair<double, Utils::Vector3d> const &acc,
                         Particle const *const p) {
  if (!p->p.is_virtual) {
    auto const mass = p->p.mass;
    return std::make_pair(acc.first + mass, acc.second + mass * p->r.p);
  }
  return acc;
}

using ComPosition = ComReduce<particle_trait, Utils::Vector3d>;

} // Namespace Observables


/*
#include <numeric>
#include <vector>
#include <array>
#include <utility>
#include <tuple>


struct P {
  std::array<double,3> pos;
  double m;
  bool v;
};

inline std::array<double,3> add_arrays(std::array<double,3>  A, std::array<double,3>  B) {
  return {A[0] + B[0],A[1] + B[1],A[2] + B[2]};
}

25 lines
std::array<double,3> foo(std::vector<P const * > &  v) {
  auto kernel = [](std::pair<double, std::array<double,3>> const &acc, P const *const p) {
    if (!p->v) {
      auto const m = p->m;
      return std::make_pair(acc.first + m, add_arrays(acc.second, p->pos));
    }
    return acc;
  };
    double a;
    std::array<double,3> b;
    std::tie(a, b) =
        std::accumulate(v.begin(), v.end(),
                        std::make_pair(0.0, std::array<double,3>{}),
                        kernel);
    return b;
}

32 lines
std::array<double,3> foo2(std::vector<P const * > &  v) {
    double a = 0;
    std::array<double,3> b{};
    for (auto p: v) {
        if (!p->v) {
          a += p->m;
          b = add_arrays(b, p->pos);
        }
    }
    return b;
}

*/
#endif
