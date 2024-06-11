/*
 * Copyright (C) 2010-2022 The ESPResSo project
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
#define BOOST_TEST_MODULE AutoParameters test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/variant.hpp>

#include "script_interface/auto_parameters/AutoParameters.hpp"

using ScriptInterface::AutoParameters;

struct A : AutoParameters<A> {
  A(int i_, int j_) : AutoParameters({{"i", i}, {"j", j}}), i(i_), j(j_) {}

  int i;
  const int j;
};

BOOST_AUTO_TEST_CASE(basic) {
  A a{0, 42};

  auto const &parameters = a.valid_parameters();

  BOOST_CHECK(parameters.size() == 2u);
  BOOST_CHECK(std::ranges::find(parameters, "i") != parameters.end());
  BOOST_CHECK(std::ranges::find(parameters, "j") != parameters.end());

  BOOST_CHECK(0 == boost::get<int>(a.get_parameter("i")));
  BOOST_CHECK(42 == boost::get<int>(a.get_parameter("j")));

  a.set_parameter("i", 12);

  BOOST_CHECK(12 == boost::get<int>(a.get_parameter("i")));
  BOOST_CHECK(42 == boost::get<int>(a.get_parameter("j")));
}

struct B : public A {
  B(int i_, int j_, int k_, int l_) : A(i_, j_), k(k_), l(l_) {
    add_parameters({{"k", k}, {"i", l} /* override i accessor */});
  }
  int k;
  int l;
};

BOOST_AUTO_TEST_CASE(add_parameters) {
  B b{1, 2, 3, 4};
  A &a = b;

  BOOST_CHECK_EQUAL(a.i, 1);
  BOOST_CHECK_EQUAL(b.i, 1);
  BOOST_CHECK_EQUAL(boost::get<int>(b.get_parameter("j")), 2);
  BOOST_CHECK_EQUAL(boost::get<int>(b.get_parameter("k")), 3);
  BOOST_CHECK_EQUAL(boost::get<int>(b.get_parameter("i")), 4);
  BOOST_CHECK_EQUAL(boost::get<int>(a.get_parameter("i")), 4);
  b.set_parameter("k", 12);
  BOOST_CHECK_EQUAL(boost::get<int>(b.get_parameter("k")), 12);
}

BOOST_AUTO_TEST_CASE(exceptions) {
  A a{0, 42};

  BOOST_CHECK_THROW(a.get_parameter("unknown"),
                    AutoParameters<A>::UnknownParameter);
  BOOST_CHECK_THROW(a.set_parameter("unknown", 12),
                    AutoParameters<A>::UnknownParameter);
  BOOST_CHECK_THROW(a.set_parameter("j", 12), AutoParameters<A>::WriteError);
}
