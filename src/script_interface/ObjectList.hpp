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

#ifndef SCRIPT_INTERFACE_OBJECT_LIST_HPP
#define SCRIPT_INTERFACE_OBJECT_LIST_HPP

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/get_value.hpp"

#include <memory>
#include <string>
#include <type_traits>
#include <vector>

namespace ScriptInterface {
/**
 * @brief Owning list of ObjectHandles
 * @tparam ManagedType Type of the managed objects, needs to be
 *         derived from ObjectHandle
 */
template <
    typename ManagedType,
    class = std::enable_if_t<std::is_base_of<ObjectHandle, ManagedType>::value>>
class ObjectList : public ObjectHandle {
public:
  virtual void add_in_core(const std::shared_ptr<ManagedType> &obj_ptr) = 0;
  virtual void remove_in_core(const std::shared_ptr<ManagedType> &obj_ptr) = 0;
  Variant do_call_method(std::string const &method,
                         VariantMap const &parameters) override {

    if (method == "add") {
      auto obj_ptr =
          get_value<std::shared_ptr<ManagedType>>(parameters.at("object"));

      add_in_core(obj_ptr);
      m_elements.push_back(obj_ptr);
      return none;
    }

    if (method == "remove") {
      auto obj_ptr =
          get_value<std::shared_ptr<ManagedType>>(parameters.at("object"));

      remove_in_core(obj_ptr);
      m_elements.erase(
          std::remove(m_elements.begin(), m_elements.end(), obj_ptr),
          m_elements.end());
      return none;
    }

    if (method == "get_elements") {
      std::vector<Variant> ret;
      ret.reserve(m_elements.size());

      for (auto const &e : m_elements)
        ret.emplace_back(e);

      return ret;
    }

    if (method == "clear") {
      for (auto const &e : m_elements) {
        remove_in_core(e);
      }

      m_elements.clear();
      return none;
    }

    if (method == "size") {
      return static_cast<int>(m_elements.size());
    }

    if (method == "empty") {
      return m_elements.empty();
    }

    return none;
  }

private:
  std::vector<std::shared_ptr<ManagedType>> m_elements;
};
} // Namespace ScriptInterface
#endif
