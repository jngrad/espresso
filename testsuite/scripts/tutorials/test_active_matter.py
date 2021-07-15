# Copyright (C) 2019 The ESPResSo project
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
import scipy
import setuptools

tutorial, skipIfMissingFeatures = importlib_wrapper.configure_and_import(
    "@TUTORIALS_DIR@/active_matter/active_matter.py", 
    ED_PARAMS={'n_sampling_steps': 100000,
               'time_step': 0.01,
               'box_l': 3 * [10.],
               'skin': 0.4,
               'active_velocity': 5,
               'kT': 1,
               'gamma': 1,
               'gamma_rotation': 1,
               'mass': 0.1,
               'rinertia': 3 * [1.]},
    RECT_PARAMS={'length': 100,
                 'radius': 20,
                 'funnel_inner_radius': 3,
                 'funnel_angle': np.pi / 4.0,
                 'funnel_thickness': 0.1,
                 'n_particles': 500,
                 'active_velocity': 5,
                 'steps_per_sample': 500,
                 'n_samples': 100,
                 'time_step': 0.01,
                 'wca_sigma': 0.5,
                 'wca_epsilon': 1,
                 'skin': 0.4,
                 'kT': 1.,
                 'gamma': 1.,
                 'gamma_rotation': 1},
    HYDRO_PARAMS={'box_l': 3 * [25],
                  'time_step': 0.01,
                  'run_steps': 200,
                  'skin': 1,
                  'agrid': 1,
                  'dens': 1,
                  'visc': 1,
                  'gamma': 1,
                  'mass': 5,
                  'dipole_length': 2,
                  'active_force': 0.1,
                  'mode': 'pusher'}
)


@skipIfMissingFeatures
class TestActMat(ut.TestCase):
    system = tutorial.system

    @ut.skipIf(not setuptools.version.pkg_resources.packaging.specifiers.SpecifierSet('>=1.4.0').contains(
        scipy.__version__), "Skipping test: scipy version requirement (>=1.4.0) not met")
    def test_quaternion(self):
        """ Check the quaternion function is correctly implemented
        """
        import scipy.spatial.transform as sst
        for theta in np.linspace(0, 2 * np.pi, 10):
            for phi in np.linspace(0, np.pi, 10):
                q_ref = sst.Rotation.from_euler('yz', [theta, phi]).as_quat()
                q_tut = tutorial.a2quat(theta, phi)
                q_tut = q_tut[1:] + [q_tut[0]]
                np.testing.assert_allclose(q_tut, q_ref, rtol=0., atol=1e-10)

    def test_enhanced_diffusion(self):
        """ Check that the active particle diffuses faster than the passive one
        """
        self.assertGreater(
            tutorial.msd_result[-1, 0], tutorial.msd_result[-1, 1])

    def test_rectification(self):
        """ Check that the center of mass is in the right half of the box
        """
        self.assertGreater(tutorial.com_deviations[-1], 0)

    def test_hydrodynamics(self):
        """ Check that the particle is moving up and the fluid down
        """
        self.assertGreater(
            tutorial.system.analysis.linear_momentum(
                include_lbfluid=False)[2], 0)
        self.assertLess(
            tutorial.system.analysis.linear_momentum(
                include_particles=False)[2], 0)


if __name__ == "__main__":
    ut.main() 
