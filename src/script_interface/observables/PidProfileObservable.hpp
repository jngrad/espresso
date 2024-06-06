/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#ifndef SCRIPT_INTERFACE_OBSERVABLES_PIDPROFILEOBSERVABLE_HPP
#define SCRIPT_INTERFACE_OBSERVABLES_PIDPROFILEOBSERVABLE_HPP

#include "script_interface/auto_parameters/AutoParameters.hpp"
#include "script_interface/observables/Observable.hpp"

#include "core/observables/DensityProfile.hpp"
#include "core/observables/FluxDensityProfile.hpp"
#include "core/observables/ForceDensityProfile.hpp"
#include "core/observables/PidProfileObservable.hpp"

#include <algorithm>
#include <cstddef>
#include <iterator>
#include <memory>
#include <string>
#include <type_traits>
#include <vector>

namespace ScriptInterface {
namespace Observables {

template <typename CoreObs>
class PidProfileObservable
    : public AutoParameters<PidProfileObservable<CoreObs>, Observable> {
  using Base = AutoParameters<PidProfileObservable<CoreObs>, Observable>;

public:
  using Base::Base;
  PidProfileObservable() {
    this->add_parameters(
        {{"ids", AutoParameter::read_only,
          [this]() { return pid_profile_observable()->ids(); }},
         {"n_x_bins", AutoParameter::read_only,
          [this]() {
            return static_cast<int>(pid_profile_observable()->n_bins()[0]);
          }},
         {"n_y_bins", AutoParameter::read_only,
          [this]() {
            return static_cast<int>(pid_profile_observable()->n_bins()[1]);
          }},
         {"n_z_bins", AutoParameter::read_only,
          [this]() {
            return static_cast<int>(pid_profile_observable()->n_bins()[2]);
          }},
         {"min_x", AutoParameter::read_only,
          [this]() { return pid_profile_observable()->limits()[0].first; }},
         {"min_y", AutoParameter::read_only,
          [this]() { return pid_profile_observable()->limits()[1].first; }},
         {"min_z", AutoParameter::read_only,
          [this]() { return pid_profile_observable()->limits()[2].first; }},
         {"max_x", AutoParameter::read_only,
          [this]() { return pid_profile_observable()->limits()[0].second; }},
         {"max_y", AutoParameter::read_only,
          [this]() { return pid_profile_observable()->limits()[1].second; }},
         {"max_z", AutoParameter::read_only,
          [this]() { return pid_profile_observable()->limits()[2].second; }}});
  }

  void do_construct(VariantMap const &params) override {
    ObjectHandle::context()->parallel_try_catch([&]() {
      m_observable =
          make_shared_from_args<CoreObs, std::vector<int>, int, int, int,
                                double, double, double, double, double, double>(
              params, "ids", "n_x_bins", "n_y_bins", "n_z_bins", "min_x",
              "max_x", "min_y", "max_y", "min_z", "max_z");
    });
  }

  Variant do_call_method(std::string const &method,
                         VariantMap const &parameters) override {
    if (method == "edges") {
      std::vector<Variant> variant_edges;
      std::ranges::copy(pid_profile_observable()->edges(),
                        std::back_inserter(variant_edges));
      return variant_edges;
    }
    return Base::do_call_method(method, parameters);
  }

  std::shared_ptr<::Observables::PidProfileObservable>
  pid_profile_observable() const {
    return m_observable;
  }

  std::shared_ptr<::Observables::Observable> observable() const override {
    return m_observable;
  }

private:
  std::shared_ptr<CoreObs> m_observable;
};

} /* namespace Observables */
} /* namespace ScriptInterface */

#endif
