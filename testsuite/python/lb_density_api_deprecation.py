#
# Copyright (C) 2026 The ESPResSo project
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
import unittest_decorators as utx
import warnings
import numpy as np
import espressomd
import espressomd.lb

system = espressomd.System(box_l=[8.0, 8.0, 8.0])
system.time_step = 0.01
system.cell_system.skin = 0.1


@utx.skipIfMissingFeatures("WALBERLA")
class TestLBDensityScalarDeprecation(ut.TestCase):

    def setUp(self):
        self.lbf = espressomd.lb.LBFluid(
            agrid=1.0, density=1.0, kinematic_viscosity=1.0, tau=0.01)
        system.lb = self.lbf

    def tearDown(self):
        system.lb = None

    def test_single_component_density_returns_float(self):
        rho = self.lbf[0, 0, 0].density
        self.assertAlmostEqual(float(rho), 1.0, places=10)

    def test_legacy_indexing_emits_deprecation_warning(self):
        rho = self.lbf[0, 0, 0].density
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            value = rho[0]
            self.assertAlmostEqual(float(value), 1.0, places=10)
            deprecation_warnings = [
                w for w in caught if issubclass(w.category, DeprecationWarning)
            ]
            self.assertGreater(len(deprecation_warnings), 0,
                               "expected a DeprecationWarning on density[0]")
            self.assertIn("density", str(deprecation_warnings[0].message))


class TestLBColorGradientDensityAccessor:

    lb_params = {}

    def setUp(self):
        self.lbf = espressomd.lb.LBFluid(
            agrid=1.0, density=1.0,
            kinematic_viscosity=[1.0, 1.0],
            tau=0.01, sigma=0.0, beta=0.7, **self.lb_params)
        system.lb = self.lbf
        # initialize component densities so rho_a + rho_b = 1.0
        shape = (8, 8, 8)
        half = np.full(shape, 0.5)
        self.lbf[:, :, :].density = np.stack([half, half], axis=-1)

    def tearDown(self):
        system.lb = None

    def test_density_is_length_two(self):
        rho = self.lbf[0, 0, 0].density
        self.assertEqual(len(rho), 2)
        self.assertAlmostEqual(float(np.sum(rho)), 1.0, places=10)


@utx.skipIfMissingFeatures("WALBERLA")
class TestLBColorGradientDensityAccessorCPU(TestLBColorGradientDensityAccessor,
                                            ut.TestCase):
    lb_params = {"single_precision": False, "gpu": False}


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
class TestLBColorGradientDensityAccessorGPU(TestLBColorGradientDensityAccessor,
                                            ut.TestCase):
    lb_params = {"single_precision": False, "gpu": True}


if __name__ == "__main__":
    ut.main()
