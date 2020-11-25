/*
 * Copyright (C) 2017-2019 The ESPResSo project
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

#define BOOST_TEST_MODULE coordinate_transformation test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <utils/Vector.hpp>
#include <utils/constants.hpp>
#include <utils/math/coordinate_transformation.hpp>
#include <utils/math/vec_rotate.hpp>

#include <cmath>

using Utils::Vector3d;

BOOST_AUTO_TEST_CASE(cartesian_to_cylinder_test) {
  Vector3d const cart_coord{{1.0, 3.3, 2.0}};
  auto const transformed_x = transform_coordinate_cartesian_to_cylinder(
      cart_coord, Vector3d{{1, 0, 0}});
  auto const transformed_y = transform_coordinate_cartesian_to_cylinder(
      cart_coord, Vector3d{{0, 1, 0}});
  auto const transformed_z = transform_coordinate_cartesian_to_cylinder(
      cart_coord, Vector3d{{0, 0, 1}});
  // For x as the symmetry axis we rotate the cartesian coordinates around the
  // y-axis by -pi/2.
  auto const expected_x = transform_coordinate_cartesian_to_cylinder(
      vec_rotate(Vector3d{{0.0, 1.0, 0.0}}, -Utils::pi() / 2.0, cart_coord),
      Vector3d{{0, 0, 1}});
  // For y as the symmetry axis we rotate the cartesian coordinates around the
  // x-axis by pi/2.
  auto const expected_y = transform_coordinate_cartesian_to_cylinder(
      vec_rotate(Vector3d{{1.0, 0.0, 0.0}}, Utils::pi() / 2.0, cart_coord),
      Vector3d{{0, 0, 1}});
  auto const expected_z = Vector3d{
      {std::sqrt(cart_coord[0] * cart_coord[0] + cart_coord[1] * cart_coord[1]),
       std::atan2(cart_coord[1], cart_coord[0]), cart_coord[2]}};

  for (int i = 0; i < 3; ++i) {
    BOOST_CHECK(transformed_x[i] == expected_x[i]);
    BOOST_CHECK(transformed_y[i] == expected_y[i]);
    BOOST_CHECK(transformed_z[i] == expected_z[i]);
  }
}

BOOST_AUTO_TEST_CASE(cartesian_to_cylinder_orientation_test) {
  constexpr auto eps = 1e-15;
  // tilted orthogonal basis
  auto const x =
      (Vector3d{{1, 0, 0}} - (1. / 3.) * Vector3d{{1, 1, 1}}).normalize();
  auto const y = (Vector3d{{0, 1, -1}}).normalize();
  auto const z = (Vector3d{{1, 1, 1}}).normalize();
  // check simple transformation
  {
    auto const rot_x = transform_coordinate_cartesian_to_cylinder(x, z);
    auto const rot_y = transform_coordinate_cartesian_to_cylinder(y, z);
    auto const rot_z = transform_coordinate_cartesian_to_cylinder(z, z);
    auto const ref_x = Vector3d{{1.0, 0.0, 0.0}};
    auto const ref_y = Vector3d{{1.0, Utils::pi() / 2.0, 0.0}};
    auto const ref_z = Vector3d{{0.0, rot_z[1], 1.0}};
    for (int i = 0; i < 3; ++i) {
      BOOST_CHECK_SMALL(rot_x[i] - ref_x[i], eps);
      BOOST_CHECK_SMALL(rot_y[i] - ref_y[i], eps);
      BOOST_CHECK_SMALL(rot_z[i] - ref_z[i], eps);
    }
  }
  // check transformation with orientation
  {
    auto const rot_x = transform_coordinate_cartesian_to_cylinder(x, z, y);
    auto const rot_y = transform_coordinate_cartesian_to_cylinder(y, z, y);
    auto const rot_z = transform_coordinate_cartesian_to_cylinder(z, z, y);
    auto const ref_x = Vector3d{{1.0, -Utils::pi() / 2.0, 0.0}};
    auto const ref_y = Vector3d{{1.0, 0.0, 0.0}};
    auto const ref_z = Vector3d{{0.0, rot_z[1], 1.0}};
    for (int i = 0; i < 3; ++i) {
      BOOST_CHECK_SMALL(rot_x[i] - ref_x[i], eps);
      BOOST_CHECK_SMALL(rot_y[i] - ref_y[i], eps);
      BOOST_CHECK_SMALL(rot_z[i] - ref_z[i], eps);
    }
  }
  // check transformation with orientation for another angle
  {
    auto const rx = vec_rotate(z, Utils::pi() / 3.0, x);
    auto const ry = vec_rotate(z, Utils::pi() / 3.0, y);
    auto const rz = z;
    auto const rot_rx = transform_coordinate_cartesian_to_cylinder(rx, z, y);
    auto const rot_ry = transform_coordinate_cartesian_to_cylinder(ry, z, y);
    auto const rot_rz = transform_coordinate_cartesian_to_cylinder(rz, z, y);
    auto const ref_rx = Vector3d{{1.0, Utils::pi() * (1. / 3. - 1. / 2.), 0.0}};
    auto const ref_ry = Vector3d{{1.0, Utils::pi() / 3.0, 0.0}};
    auto const ref_rz = Vector3d{{0.0, rot_rz[1], 1.0}};
    for (int i = 0; i < 3; ++i) {
      BOOST_CHECK_SMALL(rot_rx[i] - ref_rx[i], eps);
      BOOST_CHECK_SMALL(rot_ry[i] - ref_ry[i], eps);
      BOOST_CHECK_SMALL(rot_rz[i] - ref_rz[i], eps);
    }
  }
}

BOOST_AUTO_TEST_CASE(cylinder_to_cartesian_test) {
  Vector3d const cylinder_coord{{1.2, 3.123, 42.0}};
  auto const transformed_x = transform_coordinate_cylinder_to_cartesian(
      cylinder_coord, Vector3d{{1, 0, 0}});
  auto const transformed_y = transform_coordinate_cylinder_to_cartesian(
      cylinder_coord, Vector3d{{0, 1, 0}});
  auto const transformed_z = transform_coordinate_cylinder_to_cartesian(
      cylinder_coord, Vector3d{{0, 0, 1}});
  // We transform from cylinder zu cartesian and have to rotate back. See test
  // cartesian_to_cylinder_test.
  auto const expected_x =
      vec_rotate(Vector3d{{0.0, 1.0, 0.0}}, Utils::pi() / 2.0,
                 transform_coordinate_cylinder_to_cartesian(
                     cylinder_coord, Vector3d{{0, 0, 1}}));
  auto const expected_y =
      vec_rotate(Vector3d{{1.0, 0.0, 0.0}}, -Utils::pi() / 2.0,
                 transform_coordinate_cylinder_to_cartesian(
                     cylinder_coord, Vector3d{{0, 0, 1}}));
  // x = r * cos(phi); y = r * sin(phi); z = z
  auto const expected_z = Vector3d{
      {cylinder_coord[0] * std::cos(cylinder_coord[1]),
       cylinder_coord[0] * std::sin(cylinder_coord[1]), cylinder_coord[2]}};
  for (int i = 0; i < 3; ++i) {
    BOOST_CHECK(transformed_x[i] == expected_x[i]);
    BOOST_CHECK(transformed_y[i] == expected_y[i]);
    BOOST_CHECK(transformed_z[i] == expected_z[i]);
  }
}
