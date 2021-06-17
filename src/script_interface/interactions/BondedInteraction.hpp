/*
 * Copyright (C) 2021 The ESPResSo project
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

#ifndef SCRIPT_INTERFACE_INTERACTIONS_BONDED_INTERACTION_HPP
#define SCRIPT_INTERFACE_INTERACTIONS_BONDED_INTERACTION_HPP

#include "core/bonded_interactions/bonded_interaction_data.hpp"
#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"

#include <memory>

namespace ScriptInterface {
namespace Interactions {

template <class T> class BondedInteractionInterface {
protected:
  std::shared_ptr<::Bonded_IA_Parameters> m_bonded_ia;

public:
  std::shared_ptr<::Bonded_IA_Parameters> bonded_ia() { return m_bonded_ia; }
  std::shared_ptr<const ::Bonded_IA_Parameters> bonded_ia() const {
    return m_bonded_ia;
  }
};

class BondedInteraction : public AutoParameters<BondedInteraction>,
                          public BondedInteractionInterface<BondedInteraction> {
  Variant do_call_method(std::string const &method,
                         VariantMap const &parameters) override {
    if (method == "set_sip_from_bond_id") {
      auto const bond_id = get_value<int>(parameters.at("bond_id"));
      m_bonded_ia = ::bonded_ia_params.at(bond_id);
      return none;
    }
    return none;
  }
};

class HarmonicBond : public BondedInteraction {
  using CoreBondedInteraction = ::HarmonicBond;

public:
  HarmonicBond() {
    m_bonded_ia =
        std::make_shared<::Bonded_IA_Parameters>(CoreBondedInteraction());
    add_parameters({
        {"k",
         [this](Variant const &value) {
           get_struct().k = get_value<double>(value);
         },
         [this]() { return get_struct().k; }},
        {"r_0",
         [this](Variant const &value) {
           get_struct().r = get_value<double>(value);
         },
         [this]() { return get_struct().r; }},
        {"r_cut",
         [this](Variant const &value) {
           get_struct().r_cut = get_value<double>(value);
         },
         [this]() { return get_struct().r_cut; }},
    });
  }

public:
  CoreBondedInteraction &get_struct() {
    return boost::get<CoreBondedInteraction>(*bonded_ia());
  }
};

class BondedCoulombSR : public BondedInteraction {
  using CoreBondedInteraction = ::BondedCoulombSR;

public:
  BondedCoulombSR() {
    m_bonded_ia = std::make_shared<::Bonded_IA_Parameters>(::BondedCoulombSR());
    add_parameters({
        {"q1q2",
         [this](Variant const &value) {
           get_struct().q1q2 = get_value<double>(value);
         },
         [this]() { return get_struct().q1q2; }},
    });
  }

public:
  CoreBondedInteraction &get_struct() {
    return boost::get<CoreBondedInteraction>(*bonded_ia());
  }
};

} // namespace Interactions
} // namespace ScriptInterface

#endif
