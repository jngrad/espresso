/*
 * Copyright (C) 2010-2025 The ESPResSo project
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

#include "None.hpp"

#include <utils/Vector.hpp>
#include <utils/serialization/filesystem.hpp>
#include <utils/serialization/unordered_map.hpp>
#include <utils/serialization/variant.hpp>

#include <boost/serialization/access.hpp>
#include <boost/serialization/serialization.hpp>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>

#include <cstddef>
#include <filesystem>
#include <memory>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <variant>
#include <vector>

namespace Utils {
using Vector3b = Utils::Vector<bool, 3>;
}

namespace ScriptInterface {

template <class... Ts> struct recursive_variant;

namespace detail {
// CRTP helper class; declare base class of @ref recursive_variant here,
// such that `using BaseClass = detail::recursive_variant_base<Ts...>`
// can be used in the derived class to avoid writing the entire typelist again
template <class... Ts>
using recursive_variant_add_containers =
    std::variant<Ts..., std::vector<recursive_variant<Ts...>>,
                 std::unordered_map<int, recursive_variant<Ts...>>,
                 std::unordered_map<std::string, recursive_variant<Ts...>>>;
} // namespace detail

template <class... Ts>
struct recursive_variant : detail::recursive_variant_add_containers<Ts...> {
  using BaseClass = detail::recursive_variant_add_containers<Ts...>;
  using BaseClass::BaseClass;

private:
  friend class boost::serialization::access;

  template <typename Archive>
  void serialize(Archive &ar, unsigned const /*version*/) {
    BaseClass &self = *this;
    ar & self;
  }
};

/**
 * @brief Helper typedef to generate consistent variant types.
 *
 * This is a custom recursive variant type designed specifically for the script
 * interface. It features all basic types required to interface with the core.
 * Recursive types are added by @c detail::recursive_variant_add_containers.
 * The template parameter @c ObjectType is used to ensure instantiations like
 * @ref Variant and @ref PackedVariant hold the same types in the same order.
 *
 * @tparam ObjectType Type of the script interface object handle or reference.
 */
template <typename ObjectType>
using make_recursive_variant =
    recursive_variant<None, bool, int, std::size_t, double, std::string,
                      std::filesystem::path, ObjectType, Utils::Vector3b,
                      Utils::Vector3i, Utils::Vector2d, Utils::Vector3d,
                      Utils::Vector4d, std::vector<int>, std::vector<double>>;

class ObjectHandle;
using ObjectRef = std::shared_ptr<ObjectHandle>;

/**
 * @brief None-"literal".
 */
constexpr const None none{};

/**
 * @brief Possible types for parameters.
 */
using Variant = make_recursive_variant<ObjectRef>;

using VariantMap = std::unordered_map<std::string, Variant>;

/**
 * @brief Make a Variant from argument.
 *
 * This is a convenience function, so that rather involved constructors from
 * @ref Variant are not needed in the script interface.
 */
template <typename T> Variant make_variant(const T &x) { return Variant(x); }

template <typename K, typename V>
auto make_unordered_map_of_variants(std::unordered_map<K, V> const &v) {
  return std::unordered_map<K, Variant>{v.begin(), v.end()};
}

template <typename T> auto make_vector_of_variants(std::vector<T> const &v) {
  return std::vector<Variant>{v.begin(), v.end()};
}

/**
 * @brief Check is a Variant holds a specific type.
 *
 * @tparam T type to check for
 * @param v Variant to check in
 * @return true, if v holds a T.
 */
template <class T> constexpr bool is_type(Variant const &v) {
  return std::holds_alternative<T>(v);
}

constexpr bool is_none(Variant const &v) { return is_type<None>(v); }
} // namespace ScriptInterface
