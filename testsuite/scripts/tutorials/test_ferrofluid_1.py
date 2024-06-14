# Copyright (C) 2019-2022 The ESPResSo project
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

import unittest as ut
import importlib_wrapper
import numpy as np


tutorial, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    "@TUTORIALS_DIR@/ferrofluid/ferrofluid_part1.py",
    EQUIL_STEPS=200, EQUIL_ROUNDS=20,
    CI_DP3M_PARAMS={'cao': 3, 'r_cut': 8.34, 'mesh': [8, 8, 8], 'alpha': 0.2115, 'tune': False})


@skipIfMissingFeatures
class Tutorial(ut.TestCase):
    system = tutorial.system

    def test(self):
        self.assertEqual(
            int(np.sum(tutorial.n_clusters)), len(tutorial.cluster_sizes))
        for i in range(7):
            self.assertLess(
                tutorial.size_dist[0][i + 1],
                tutorial.size_dist[0][i])

        # check exponential decay in the tail
        xdata_sim = tutorial.xdata[4:]
        ydata_sim = tutorial.ydata[4:]
        ydata_fit = tutorial.kernel(xdata_sim, *tutorial.popt)
        np.testing.assert_allclose(ydata_sim, ydata_fit, atol=1.0, rtol=0.2)


if __name__ == "__main__":
    ut.main()
