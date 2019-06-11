#
# Copyright (C) 2017-2018 The ESPResSo project
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
from __future__ import print_function
import unittest as ut
import unittest_system as uts
import numpy as np


class RescaleTest(uts.TestCaseSystem):

    """Test the global box and particle rescaling.

    """

    def setUp(self):
        N = 100
        self.system.cell_system.skin = 0.0
        self.system.time_step = 0.01
        self.system.box_l = 3 * [10]
        self.system.part.add(pos=self.system.box_l * np.random.random((N, 3)))

    def test_iso(self):
        """Test 'isotropic' case (dir="xyz").
        """
        scale = 1.3
        new_box_l = scale * self.system.box_l[0]

        old_pos = self.system.part[:].pos
        self.system.change_volume_and_rescale_particles(new_box_l)
        new_pos = self.system.part[:].pos

        max_diff = np.max(np.abs(new_pos / old_pos - scale))
        self.assertAlmostEqual(0., max_diff, places=10)

    def dir_test(self, dir):
        """Test scaling of a single direction.
        """
        scale = 0.7
        new_box_l = scale * self.system.box_l[dir]

        old_pos = self.system.part[:].pos
        self.system.change_volume_and_rescale_particles(new_box_l, dir=dir)
        new_pos = self.system.part[:].pos

        for i in range(3):
            if i == dir:
                max_diff = np.max(
                    np.abs(new_pos[:, i] / old_pos[:, i] - scale))
            else:
                max_diff = np.max(np.abs(new_pos[:, i] - old_pos[:, i]))

            self.assertAlmostEqual(0., max_diff, places=10)

    def test_x(self):
        self.dir_test(0)

    def test_y(self):
        self.dir_test(1)

    def test_z(self):
        self.dir_test(2)

if __name__ == "__main__":
    ut.main()
