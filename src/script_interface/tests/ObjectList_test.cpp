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

#define BOOST_TEST_MODULE ObjectList test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/range/algorithm/find.hpp>

#include "script_interface/LocalContext.hpp"
#include "script_interface/ObjectList.hpp"

#include <algorithm>
#include <memory>
#include <vector>

using namespace ScriptInterface;

struct ObjectListImpl : ObjectList<ObjectHandle> {
  std::vector<ObjectRef> mock_core;
  void add(ObjectRef const &e) { do_call_method("add", {{"object", e}}); }
  void remove(ObjectRef const &e) { do_call_method("remove", {{"object", e}}); }
  int size() { return get_value<int>(do_call_method("size", {})); }
  bool empty() { return get_value<bool>(do_call_method("empty", {})); }

private:
  void add_in_core(const ObjectRef &obj_ptr) override {
    mock_core.push_back(obj_ptr);
  }
  void remove_in_core(const ObjectRef &obj_ptr) override {
    mock_core.erase(std::remove(mock_core.begin(), mock_core.end(), obj_ptr),
                    mock_core.end());
  }
};

BOOST_AUTO_TEST_CASE(default_construction) {
  // A defaulted ObjectList has no elements.
  BOOST_CHECK(ObjectListImpl{}.empty());
}

BOOST_AUTO_TEST_CASE(adding_elements) {
  // Added elements are on the back of the list of elements.
  ObjectRef e = std::make_shared<ObjectHandle>();
  ObjectListImpl list;
  list.add(e);
  // And is added to the core
  BOOST_CHECK(list.mock_core.back() == e);
  BOOST_CHECK(e.use_count() == 3);
}

BOOST_AUTO_TEST_CASE(removing_elements) {
  // An element that is removed from the list is
  // no longer an element of the list.
  ObjectRef e = std::make_shared<ObjectHandle>();
  ObjectListImpl list;
  list.add(e);
  BOOST_CHECK(list.size() == 1);
  list.remove(e);
  BOOST_CHECK(list.size() == 0);
  BOOST_CHECK(e.use_count() == 1);
  // And is removed from the core
  BOOST_CHECK(boost::find(list.mock_core, e) == list.mock_core.end());
}

BOOST_AUTO_TEST_CASE(clearing_elements) {
  // A cleared list is empty.
  ObjectListImpl list;
  list.add(std::make_shared<ObjectHandle>());
  list.add(std::make_shared<ObjectHandle>());
  list.do_call_method("clear", {});
  BOOST_CHECK(list.empty());
  BOOST_CHECK(list.mock_core.empty());
}

BOOST_AUTO_TEST_CASE(serialization) {
  // In a context
  Utils::Factory<ObjectHandle> f;
  f.register_new<ObjectHandle>("ObjectHandle");
  f.register_new<ObjectListImpl>("ObjectList");
  auto ctx = std::make_shared<LocalContext>(f);
  // A list of some elements
  auto list = std::dynamic_pointer_cast<ObjectListImpl>(
      ctx->make_shared("ObjectList", {}));
  // with a bunch of elements
  list->add(std::make_shared<ObjectHandle>());
  list->add(std::make_shared<ObjectHandle>());
  // can be partially serialized to a string
  auto const s = list->serialize();
  // and is restored to an empty list
  auto const list2 = std::dynamic_pointer_cast<ObjectListImpl>(
      ObjectHandle::deserialize(s, *ctx));
  BOOST_CHECK_EQUAL(list2->size(), 0);
  // and the elements are not restored to the core
  BOOST_CHECK_EQUAL(list2->mock_core.size(), 0);
}
