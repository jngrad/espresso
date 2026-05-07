/*
 * Copyright (C) 2021-2026 The ESPResSo project
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

#pragma once

#include "script_interface/ScriptInterface.hpp"

#include <stdexcept>
#include <string>

namespace ScriptInterface {
namespace ReactionMethods {

class ReactionAlgorithm : public AutoParameters<ReactionAlgorithm> {
public:
  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override;

protected:
private:
  std::string get_internal_state() const override {
    throw std::runtime_error("Reaction methods do not support checkpointing");
  }
};

} /* namespace ReactionMethods */
} /* namespace ScriptInterface */
