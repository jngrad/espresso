/*
 * Copyright (C) 2023 The ESPResSo project
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

#include "barista/FancyBrew.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"

#include <memory>
#include <string>

namespace ScriptInterface::Barista {

class FancyBrew : public AutoParameters<FancyBrew> {
  std::shared_ptr<::Barista::FancyBrew> m_instance;

public:
  FancyBrew();
  void do_construct(VariantMap const &params) override;
  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override;
};

} // namespace ScriptInterface::Barista
