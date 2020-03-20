/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#include <memory>

#include "core/observables/PidProfileObservable.hpp"

namespace ScriptInterface {
namespace Observables {

template <typename CoreObs,
          typename CoordinateSystem = ::Observables::CoordSystem::Cartesian>
class PidProfileObservable
    : public AutoParameters<PidProfileObservable<CoreObs, CoordinateSystem>,
                            Observable> {
public:
  PidProfileObservable() {
    this->add_parameters(
        {{"ids",
          [this](const Variant &v) {
            pid_profile_observable()->ids() = get_value<std::vector<int>>(v);
          },
          [this]() { return pid_profile_observable()->ids(); }},
         {"n_x_bins",
          [this](const Variant &v) {
            pid_profile_observable()->n_x_bins =
                static_cast<size_t>(get_value<int>(v));
          },
          [this]() {
            return static_cast<int>(pid_profile_observable()->n_x_bins);
          }},
         {"n_y_bins",
          [this](const Variant &v) {
            pid_profile_observable()->n_y_bins =
                static_cast<size_t>(get_value<int>(v));
          },
          [this]() {
            return static_cast<int>(pid_profile_observable()->n_y_bins);
          }},
         {"n_z_bins",
          [this](const Variant &v) {
            pid_profile_observable()->n_z_bins =
                static_cast<size_t>(get_value<int>(v));
          },
          [this]() {
            return static_cast<int>(pid_profile_observable()->n_z_bins);
          }},
         {"min_x",
          [this](const Variant &v) {
            pid_profile_observable()->min_x = get_value<double>(v);
          },
          [this]() { return pid_profile_observable()->min_x; }},
         {"min_y",
          [this](const Variant &v) {
            pid_profile_observable()->min_y = get_value<double>(v);
          },
          [this]() { return pid_profile_observable()->min_y; }},
         {"min_z",
          [this](const Variant &v) {
            pid_profile_observable()->min_z = get_value<double>(v);
          },
          [this]() { return pid_profile_observable()->min_z; }},
         {"max_x",
          [this](const Variant &v) {
            pid_profile_observable()->max_x = get_value<double>(v);
          },
          [this]() { return pid_profile_observable()->max_x; }},
         {"max_y",
          [this](const Variant &v) {
            pid_profile_observable()->max_y = get_value<double>(v);
          },
          [this]() { return pid_profile_observable()->max_y; }},
         {"max_z",
          [this](const Variant &v) {
            pid_profile_observable()->max_z = get_value<double>(v);
          },
          [this]() { return pid_profile_observable()->max_z; }}});
  }

  void construct(VariantMap const &params) override {
    m_observable =
        make_shared_from_args<CoreObs, std::vector<int>, int, int, int, double,
                              double, double, double, double, double>(
            params, "ids", "n_x_bins", "n_y_bins", "n_z_bins", "min_x", "min_y",
            "min_z", "max_x", "max_y", "max_z");
  }

  std::shared_ptr<::Observables::PidProfileObservable<CoordinateSystem>>
  pid_profile_observable() const {
    return m_observable;
  }

  std::shared_ptr<::Observables::Observable> observable() const override {
    return m_observable;
  }

private:
  std::shared_ptr<CoreObs> m_observable;
};

template <typename CoreObs>
class PidProfileObservable<CoreObs, ::Observables::CoordSystem::Cylindrical>
    : public AutoParameters<
          PidProfileObservable<CoreObs,
                               ::Observables::CoordSystem::Cylindrical>,
          Observable> {
public:
  PidProfileObservable() {
    this->add_parameters(
        {{"ids",
          [this](const Variant &v) {
            pid_profile_observable()->ids() = get_value<std::vector<int>>(v);
          },
          [this]() { return pid_profile_observable()->ids(); }},
         {"center",
          [this](const Variant &v) {
            pid_profile_observable()->center = get_value<::Utils::Vector3d>(v);
          },
          [this]() { return pid_profile_observable()->center; }},
         {"axis",
          [this](const Variant &v) {
            pid_profile_observable()->axis = get_value<Utils::Vector3d>(v);
          },
          [this]() { return pid_profile_observable()->axis; }},
         {"n_r_bins",
          [this](const Variant &v) {
            pid_profile_observable()->n_r_bins =
                static_cast<size_t>(get_value<int>(v));
          },
          [this]() {
            return static_cast<int>(pid_profile_observable()->n_r_bins);
          }},
         {"n_phi_bins",
          [this](const Variant &v) {
            pid_profile_observable()->n_phi_bins =
                static_cast<size_t>(get_value<int>(v));
          },
          [this]() {
            return static_cast<int>(pid_profile_observable()->n_phi_bins);
          }},
         {"n_z_bins",
          [this](const Variant &v) {
            pid_profile_observable()->n_z_bins =
                static_cast<size_t>(get_value<int>(v));
          },
          [this]() {
            return static_cast<int>(pid_profile_observable()->n_z_bins);
          }},
         {"min_r",
          [this](const Variant &v) {
            pid_profile_observable()->min_r = get_value<double>(v);
          },
          [this]() { return pid_profile_observable()->min_r; }},
         {"min_phi",
          [this](const Variant &v) {
            pid_profile_observable()->min_phi = get_value<double>(v);
          },
          [this]() { return pid_profile_observable()->min_phi; }},
         {"min_z",
          [this](const Variant &v) {
            pid_profile_observable()->min_z = get_value<double>(v);
          },
          [this]() { return pid_profile_observable()->min_z; }},
         {"max_r",
          [this](const Variant &v) {
            pid_profile_observable()->max_r = get_value<double>(v);
          },
          [this]() { return pid_profile_observable()->max_r; }},
         {"max_phi",
          [this](const Variant &v) {
            pid_profile_observable()->max_phi = get_value<double>(v);
          },
          [this]() { return pid_profile_observable()->max_phi; }},
         {"max_z",
          [this](const Variant &v) {
            pid_profile_observable()->max_z = get_value<double>(v);
          },
          [this]() { return pid_profile_observable()->max_z; }}});
  }

  void construct(VariantMap const &params) override {
    m_observable =
        make_shared_from_args<CoreObs, std::vector<int>, Utils::Vector3d,
                              Utils::Vector3d, int, int, int, double, double,
                              double, double, double, double>(
            params, "ids", "center", "axis", "n_r_bins", "n_phi_bins",
            "n_z_bins", "min_r", "min_phi", "min_z", "max_r", "max_phi",
            "max_z");
  }

  std::shared_ptr<::Observables::PidProfileObservable<
      ::Observables::CoordSystem::Cylindrical>>
  pid_profile_observable() const {
    return m_observable;
  }

  std::shared_ptr<::Observables::Observable> observable() const override {
    return m_observable;
  }

private:
  std::shared_ptr<CoreObs> m_observable;
};

template <typename CoreObs>
class PidProfileObservable<CoreObs, ::Observables::CoordSystem::Spherical>
    : public AutoParameters<
          PidProfileObservable<CoreObs, ::Observables::CoordSystem::Spherical>,
          Observable> {
public:
  PidProfileObservable() {
    this->add_parameters(
        {{"ids",
          [this](const Variant &v) {
            pid_profile_observable()->ids() = get_value<std::vector<int>>(v);
          },
          [this]() { return pid_profile_observable()->ids(); }},
         {"n_r_bins",
          [this](const Variant &v) {
            pid_profile_observable()->n_r_bins =
                static_cast<size_t>(get_value<int>(v));
          },
          [this]() {
            return static_cast<int>(pid_profile_observable()->n_r_bins);
          }},
         {"n_theta_bins",
          [this](const Variant &v) {
            pid_profile_observable()->n_theta_bins =
                static_cast<size_t>(get_value<int>(v));
          },
          [this]() {
            return static_cast<int>(pid_profile_observable()->n_theta_bins);
          }},
         {"n_phi_bins",
          [this](const Variant &v) {
            pid_profile_observable()->n_phi_bins =
                static_cast<size_t>(get_value<int>(v));
          },
          [this]() {
            return static_cast<int>(pid_profile_observable()->n_phi_bins);
          }},
         {"min_r",
          [this](const Variant &v) {
            pid_profile_observable()->min_r = get_value<double>(v);
          },
          [this]() { return pid_profile_observable()->min_r; }},
         {"min_theta",
          [this](const Variant &v) {
            pid_profile_observable()->min_theta = get_value<double>(v);
          },
          [this]() { return pid_profile_observable()->min_theta; }},
         {"min_phi",
          [this](const Variant &v) {
            pid_profile_observable()->min_phi = get_value<double>(v);
          },
          [this]() { return pid_profile_observable()->min_phi; }},
         {"max_r",
          [this](const Variant &v) {
            pid_profile_observable()->max_r = get_value<double>(v);
          },
          [this]() { return pid_profile_observable()->max_r; }},
         {"max_theta",
          [this](const Variant &v) {
            pid_profile_observable()->max_theta = get_value<double>(v);
          },
          [this]() { return pid_profile_observable()->max_theta; }},
         {"max_phi",
          [this](const Variant &v) {
            pid_profile_observable()->max_phi = get_value<double>(v);
          },
          [this]() { return pid_profile_observable()->max_phi; }}});
  }

  void construct(VariantMap const &params) override {
    m_observable =
        make_shared_from_args<CoreObs, std::vector<int>, int, int, int, double,
                              double, double, double, double, double>(
            params, "ids", "n_r_bins", "n_theta_bins", "n_phi_bins", "min_r",
            "min_theta", "min_phi", "max_r", "max_theta", "max_phi");
  }

  std::shared_ptr<::Observables::PidProfileObservable<
      ::Observables::CoordSystem::Spherical>>
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
