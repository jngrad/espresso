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
import unittest as ut
import unittest_decorators as utx
import numpy as np 
import copy
# import tests_common

import espressomd
import espressomd.lb
import espressomd.observables

"""
Tests for the LB fluid profile observables.

"""

TIME_STEP = 0.1
AGRID = 0.7
BOX_L_X = 17.0 * AGRID
BOX_L_Y = 15.0 * AGRID
BOX_L_Z = 13.0 * AGRID
VISC = .7
DENS = 1.7
LB_PARAMS = {'agrid': AGRID,
             'dens': DENS,
             'visc': VISC,
             'tau': TIME_STEP
             }

LB_VELOCITY_PROFILE_PARAMS = {
    'n_x_bins': int(BOX_L_X / AGRID),
    'n_y_bins': int(BOX_L_Y / AGRID),
    'n_z_bins': int(BOX_L_Z / AGRID),
    'min_x': 0.0,
    'min_y': 0.0,
    'min_z': 0.0,
    'max_x': BOX_L_X,
    'max_y': BOX_L_Y,
    'max_z': BOX_L_Z,
    'sampling_delta_x': AGRID,
    'sampling_delta_y': AGRID,
    'sampling_delta_z': AGRID,
    'sampling_offset_x': 0.5 * AGRID,
    'sampling_offset_y': 0.5 * AGRID,
    'sampling_offset_z': 0.5 * AGRID,
    'allow_empty_bins': False}


class ObservableProfileLBCommon:
    lbf = None
    system = espressomd.System(box_l=[BOX_L_X, BOX_L_Y, BOX_L_Z])
    system.time_step = TIME_STEP
    system.cell_system.skin = 0.4 * AGRID

    def set_fluid_velocities(self):
        """Set an x dependent fluid velocity."""
        for x in range(int(self.system.box_l[0] / AGRID)):
            for y in range(int(self.system.box_l[1] / AGRID)):
                for z in range(int(self.system.box_l[2] / AGRID)):
                    self.lbf[x, y, z].velocity = [float(x), 0.0, 0.0]
                    np.testing.assert_allclose(np.copy(self.lbf[x, y, z].velocity), [
                                               float(x), 0.0, 0.0], atol=1E-6)

    def test_velocity_profile(self):
        self.set_fluid_velocities()
        obs = espressomd.observables.LBVelocityProfile(
            **LB_VELOCITY_PROFILE_PARAMS)
        obs_data = obs.calculate()
        for x in range(obs_data.shape[0]):
            for y in range(obs_data.shape[1]):
                for z in range(obs_data.shape[2]):
                    self.assertAlmostEqual(
                        obs_data[x, y, z, 0], float(x), places=5)
        self.assertEqual(np.prod(obs_data.shape),
                         LB_VELOCITY_PROFILE_PARAMS['n_x_bins'] *
                         LB_VELOCITY_PROFILE_PARAMS['n_y_bins'] *
                         LB_VELOCITY_PROFILE_PARAMS['n_z_bins'] * 3)

    def test_error_sampling_delta_of_0(self):
        lb_velocity_params_local = copy.copy(LB_VELOCITY_PROFILE_PARAMS)
        lb_velocity_params_local['sampling_delta_x'] = 0.0
        lb_velocity_params_local['sampling_delta_y'] = 0.0
        lb_velocity_params_local['sampling_delta_z'] = 0.0
        with self.assertRaises(RuntimeError):
            _ = espressomd.observables.LBVelocityProfile(
                **lb_velocity_params_local)

    def test_error_if_no_LB(self):
        self.system.actors.clear()
        obs = espressomd.observables.LBVelocityProfile(
            **LB_VELOCITY_PROFILE_PARAMS)
        with self.assertRaises(RuntimeError):
            obs.calculate()

    def test_error_if_empty_bin(self):
        lb_velocity_params_local = copy.copy(LB_VELOCITY_PROFILE_PARAMS)
        lb_velocity_params_local['sampling_delta_x'] = 3.0
        obs = espressomd.observables.LBVelocityProfile(
            **lb_velocity_params_local)
        with self.assertRaises(RuntimeError):
            obs.calculate()

    def test_lb_profile_interface(self):
        # test setters and getters
        params = LB_VELOCITY_PROFILE_PARAMS.copy()
        params['n_x_bins'] = 4
        params['n_y_bins'] = 6
        params['n_z_bins'] = 8
        obs = espressomd.observables.LBVelocityProfile(**params)
        # check flag
        self.assertFalse(obs.allow_empty_bins)
        obs = espressomd.observables.LBVelocityProfile(
            **{**params, 'allow_empty_bins': True})
        self.assertTrue(obs.allow_empty_bins)
        # check bins
        self.assertEqual(obs.n_x_bins, 4)
        self.assertEqual(obs.n_y_bins, 6)
        self.assertEqual(obs.n_z_bins, 8)
        obs_data = obs.calculate()
        np.testing.assert_array_equal(obs_data.shape, [4, 6, 8, 3])
        obs = espressomd.observables.LBVelocityProfile(
            **{**params, 'n_x_bins': 1, 'n_y_bins': 2, 'n_z_bins': 3})
        obs_data = obs.calculate()
        np.testing.assert_array_equal(obs_data.shape, [1, 2, 3, 3])
        # check edges lower corner
        self.assertEqual(obs.min_x, params['min_x'])
        self.assertEqual(obs.min_y, params['min_y'])
        self.assertEqual(obs.min_z, params['min_z'])
        obs = espressomd.observables.LBVelocityProfile(
            **{**params, 'min_x': 4, 'min_y': 5, 'min_z': 6})
        self.assertEqual(obs.min_x, 4)
        self.assertEqual(obs.min_y, 5)
        self.assertEqual(obs.min_z, 6)
        obs_bin_edges = obs.bin_edges()
        np.testing.assert_array_equal(obs_bin_edges[0, 0, 0], [4, 5, 6])
        # check edges upper corner
        self.assertEqual(obs.max_x, params['max_x'])
        self.assertEqual(obs.max_y, params['max_y'])
        self.assertEqual(obs.max_z, params['max_z'])
        obs = espressomd.observables.LBVelocityProfile(
            **{**params, 'max_x': 7, 'max_y': 8, 'max_z': 9})
        self.assertEqual(obs.max_x, 7)
        self.assertEqual(obs.max_y, 8)
        self.assertEqual(obs.max_z, 9)
        obs_bin_edges = obs.bin_edges()
        np.testing.assert_array_equal(obs_bin_edges[-1, -1, -1], [7, 8, 9])
        # check delta
        self.assertEqual(obs.sampling_delta_x, params['sampling_delta_x'])
        self.assertEqual(obs.sampling_delta_y, params['sampling_delta_y'])
        self.assertEqual(obs.sampling_delta_z, params['sampling_delta_z'])
        obs = espressomd.observables.LBVelocityProfile(
            **{**params, 'sampling_delta_x': 10, 'sampling_delta_y': 11, 'sampling_delta_z': 12})
        self.assertEqual(obs.sampling_delta_x, 10)
        self.assertEqual(obs.sampling_delta_y, 11)
        self.assertEqual(obs.sampling_delta_z, 12)
        # check delta
        self.assertEqual(obs.sampling_offset_x, params['sampling_offset_x'])
        self.assertEqual(obs.sampling_offset_y, params['sampling_offset_y'])
        self.assertEqual(obs.sampling_offset_z, params['sampling_offset_z'])
        obs = espressomd.observables.LBVelocityProfile(
            **{**params, 'sampling_offset_x': 13, 'sampling_offset_y': 14, 'sampling_offset_z': 15})
        self.assertEqual(obs.sampling_offset_x, 13)
        self.assertEqual(obs.sampling_offset_y, 14)
        self.assertEqual(obs.sampling_offset_z, 15)


@utx.skipIfMissingFeatures("LB_WALBERLA")
class LBWalberla(ut.TestCase, ObservableProfileLBCommon):

    """Test for the CPU implementation of the LB."""

    def setUp(self):
        self.lbf = espressomd.lb.LBFluidWalberla(**LB_PARAMS)
        self.system.actors.clear()
        self.system.actors.add(self.lbf)


if __name__ == "__main__":
    ut.main()
