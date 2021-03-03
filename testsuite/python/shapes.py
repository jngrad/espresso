#
# Copyright (C) 2010-2019 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
import unittest as ut
import numpy as np

import espressomd.shapes


class ShapeTests(ut.TestCase):
    def test_Union(self):
        union = espressomd.shapes.Union()
        wall1 = espressomd.shapes.Wall(normal=[0, 0, 1], dist=0)
        wall2 = espressomd.shapes.Wall(normal=[0, 0, -1], dist=-10)
        print("self.assertTrue(union.call_method('empty'))")
        self.assertTrue(union.call_method('empty'))
        print("union.add([wall1, wall2])")
        union.add([wall1, wall2])
        print("self.assertFalse(union.call_method('empty'))")
        self.assertFalse(union.call_method('empty'))
        print("self.assertEqual(union.size(), 2)")
        self.assertEqual(union.size(), 2)

        # check object retrieval
        print("pwall1, pwall2 = union.call_method('get_elements')")
        pwall1, pwall2 = union.call_method('get_elements')
        print("self.assertIsInstance(pwall1, espressomd.shapes.Wall)")
        self.assertIsInstance(pwall1, espressomd.shapes.Wall)
        print("self.assertIsInstance(pwall2, espressomd.shapes.Wall)")
        self.assertIsInstance(pwall2, espressomd.shapes.Wall)
        print("np.copy(pwall1.normal), np.copy(wall1.normal))")
        np.testing.assert_almost_equal(
            np.copy(pwall1.normal), np.copy(wall1.normal))
        print("np.copy(pwall2.normal), np.copy(wall2.normal))")
        np.testing.assert_almost_equal(
            np.copy(pwall2.normal), np.copy(wall2.normal))
        print("np.testing.assert_almost_equal(pwall1.dist, wall1.dist)")
        np.testing.assert_almost_equal(pwall1.dist, wall1.dist)
        print("np.testing.assert_almost_equal(pwall2.dist, wall2.dist)")
        np.testing.assert_almost_equal(pwall2.dist, wall2.dist)

        print("position=[1, 2, 4.5])[0], 4.5)")
        self.assertAlmostEqual(union.calc_distance(
            position=[1, 2, 4.5])[0], 4.5)
        print(position=[1, 2, 5.0])[0], 5.0)"")
        self.assertAlmostEqual(union.calc_distance(
            position=[1, 2, 5.0])[0], 5.0)
        print("position=[1, 2, 6.5])[0], 3.5)")
        self.assertAlmostEqual(union.calc_distance(
            position=[1, 2, 6.5])[0], 3.5)

        # negative distances are not well-defined for a union of shapes
        with self.assertRaises(ValueError):
            print("union.calc_distance(position=[1, 2, 11.5])")
            union.calc_distance(position=[1, 2, 11.5])
        print("union.clear()")
        union.clear()
        print("self.assertTrue(union.call_method('empty'))")
        self.assertTrue(union.call_method('empty'))
        print("self.assertEqual(union.size(), 0)")
        self.assertEqual(union.size(), 0)
        print("self.assertEqual(union.calc_distance(position=[1, 2, 6.5])[0], np.inf)")
        self.assertEqual(union.calc_distance(position=[1, 2, 6.5])[0], np.inf)

        print("union.add([wall1, wall2])")
        union.add([wall1, wall2])
        print("union.remove(wall2)")
        union.remove(wall2)
        print("position=[1, 2, 6.5])[0], 6.5)")
        self.assertAlmostEqual(union.calc_distance(
            position=[1, 2, 6.5])[0], 6.5)


if __name__ == "__main__":
    ut.main()
