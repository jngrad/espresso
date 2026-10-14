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

"""
Checkpoint/restart of a two-component color-gradient LB fluid.

Tests:
  - A simulation restarted from an LB checkpoint file continues exactly like
    the uninterrupted simulation: thermalized fluid with a droplet, a moving
    wall and a particle coupled to the fluid. The restart uses a freshly
    constructed fluid, so every piece of state that influences the dynamics
    (populations, forces, densities, phase field, velocity field, boundaries,
    RNG counter) must come from the checkpoint file.
  - Loading re-applies the boundary conditions to the ghost layers, which a
    full ghost communication would otherwise overwrite.
  - Binary files give a bit-exact restart.
  - ASCII files (16 decimal digits) give a restart within round-off.

The tests run on the CPU and on the GPU, in double precision.
"""

import unittest as ut
import unittest_decorators as utx
import numpy as np
import pathlib
import tempfile

import espressomd
import espressomd.lb
import espressomd.shapes


BOX_L = 6.0
AGRID = 0.5
TIME_STEP = 0.01
DENSITY = 1.0
VISCOSITY = [2.0, 6.0]
SIGMA = 0.005
BETA = 0.8
KT = 4e-4
LB_SEED = 7
GAMMA = 1.0
THERMOSTAT_SEED = 42
DROPLET_RADIUS = 1.5
EPSILON = 0.05
WALL_VELOCITY = [0., 1e-3, 0.]

WARMUP_STEPS = 20
CONTINUATION_STEPS = 20


class ColorGradientCheckpointTest:

    system = espressomd.System(box_l=[BOX_L] * 3)
    system.time_step = TIME_STEP
    system.cell_system.skin = 0.4

    lb_params = {}

    @classmethod
    def setUpClass(cls):
        cls.temp_dir = tempfile.TemporaryDirectory()
        cls.temp_path = pathlib.Path(cls.temp_dir.name).resolve()

    @classmethod
    def tearDownClass(cls):
        cls.temp_dir.cleanup()

    def tearDown(self):
        self.system.part.clear()
        self.system.thermostat.turn_off()
        self.system.lb = None
        self.system.time = 0.

    def make_fluid(self, walls):
        lbf = espressomd.lb.LBFluid(
            agrid=AGRID, density=DENSITY, kinematic_viscosity=VISCOSITY,
            tau=TIME_STEP, sigma=SIGMA, beta=BETA, kT=KT, seed=LB_SEED,
            **self.lb_params)
        if not walls:
            return lbf
        # a resting wall and a moving wall normal to x
        wall1 = espressomd.shapes.Wall(normal=[1., 0., 0.], dist=AGRID)
        wall2 = espressomd.shapes.Wall(normal=[-1., 0., 0.],
                                       dist=-(BOX_L - AGRID))
        lbf.add_boundary_from_shape(wall1)
        lbf.add_boundary_from_shape(wall2, WALL_VELOCITY)
        return lbf

    def init_droplet(self, lbf):
        """Droplet of component b in a matrix of component a."""
        n = lbf.shape[0]
        x = (np.arange(n) + 0.5) * AGRID
        xx, yy, zz = np.meshgrid(x, x, x, indexing='ij')
        center = BOX_L / 2.
        dist = np.sqrt((xx - center)**2 + (yy - center)**2 + (zz - center)**2)
        inside = 0.5 * (1. - np.tanh((dist - DROPLET_RADIUS) / AGRID))
        rho_b = DENSITY * (EPSILON + (1. - EPSILON) * inside)
        rho_a = DENSITY * (1. + EPSILON) - rho_b
        lbf[:, :, :].density = np.stack([rho_a, rho_b], axis=-1)
        lbf.init_two_component()

    @staticmethod
    def snapshot(lbf, p, thermostat):
        return {
            "population": np.copy(lbf[:, :, :]._population),
            "last_applied_force": np.copy(lbf[:, :, :].last_applied_force),
            "density": np.copy(lbf[:, :, :].density),
            "phasefield": np.copy(lbf[:, :, :].phasefield),
            "velocity": np.copy(lbf[:, :, :].velocity),
            "pos": np.copy(p.pos),
            "v": np.copy(p.v),
            "f": np.copy(p.f),
            "rng_state": lbf.rng_state,
            "philox_counter": thermostat.lb.philox_counter,
        }

    def check_restart(self, binary, compare):
        system = self.system
        cpt_path = self.temp_path / \
            f"lb_cg_{type(self).__name__}_{int(binary)}.cpt"

        # initial simulation
        lbf = self.make_fluid(walls=True)
        self.init_droplet(lbf)
        system.lb = lbf
        system.thermostat.set_lb(LB_fluid=lbf, gamma=GAMMA,
                                 seed=THERMOSTAT_SEED)
        # particle close to the droplet interface, so that it samples
        # velocity and density gradients of the fluid
        p = system.part.add(pos=[BOX_L / 2. + DROPLET_RADIUS, BOX_L / 2.,
                                 BOX_L / 2.], v=[0.05, 0.02, -0.01])
        system.integrator.run(WARMUP_STEPS)

        # checkpoint: LB file plus the particle and thermostat state that
        # espressomd.checkpointing would save in a real restart
        lbf.save_checkpoint(cpt_path, binary)
        saved = self.snapshot(lbf, p, system.thermostat)
        saved_time = system.time
        self.assertIsNotNone(saved["rng_state"])
        self.assertGreater(np.ptp(saved["phasefield"]), 0.5)

        # uninterrupted reference; the forces of the last step are reused
        # in both runs, as required for a restart (see the LB docs)
        system.integrator.run(CONTINUATION_STEPS, reuse_forces=True)
        ref = self.snapshot(lbf, p, system.thermostat)
        self.assertEqual(ref["rng_state"],
                         saved["rng_state"] + CONTINUATION_STEPS)
        # sanity check: the continuation actually changed the state
        self.assertGreater(
            np.max(np.abs(ref["population"] - saved["population"])), 1e-6)
        self.assertGreater(np.linalg.norm(ref["pos"] - saved["pos"]), 1e-4)

        # restart into a freshly constructed fluid without walls
        system.lb = None
        lbf_restart = self.make_fluid(walls=False)
        self.assertNotEqual(lbf_restart.rng_state, saved["rng_state"])
        lbf_restart.load_checkpoint(cpt_path, binary)
        self.assertEqual(lbf_restart.rng_state, saved["rng_state"])
        np.testing.assert_array_equal(
            np.copy(lbf_restart[:, :, :].is_boundary),
            np.copy(lbf[:, :, :].is_boundary))
        np.testing.assert_allclose(
            np.copy(lbf_restart[-1, :, :].velocity),
            np.broadcast_to(WALL_VELOCITY, (*lbf.shape[1:], 3)), atol=1e-12)
        system.lb = lbf_restart
        system.thermostat.set_lb(LB_fluid=lbf_restart, gamma=GAMMA,
                                 seed=THERMOSTAT_SEED)
        system.thermostat.lb.call_method(
            "override_philox_counter", counter=saved["philox_counter"])
        p.pos = saved["pos"]
        p.v = saved["v"]
        p.f = saved["f"]
        system.time = saved_time
        system.integrator.run(CONTINUATION_STEPS, reuse_forces=True)
        new = self.snapshot(lbf_restart, p, system.thermostat)

        self.assertEqual(new["rng_state"], ref["rng_state"])
        self.assertEqual(new["philox_counter"], ref["philox_counter"])
        for key in ("population", "last_applied_force", "density",
                    "phasefield", "velocity", "pos", "v", "f"):
            self.assertEqual(new[key].shape, ref[key].shape, msg=key)
            compare(new[key], ref[key], err_msg=f"'{key}' differs")

    def test_restart_binary(self):
        self.check_restart(True, np.testing.assert_array_equal)

    def test_restart_ascii(self):
        def compare(actual, desired, err_msg):
            np.testing.assert_allclose(actual, desired, rtol=0., atol=1e-12,
                                       err_msg=err_msg)
        self.check_restart(False, compare)


@utx.skipIfMissingFeatures(["WALBERLA"])
class ColorGradientCheckpointTestCPU(ColorGradientCheckpointTest, ut.TestCase):
    lb_params = {"single_precision": False, "gpu": False}


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
class ColorGradientCheckpointTestGPU(ColorGradientCheckpointTest, ut.TestCase):
    lb_params = {"single_precision": False, "gpu": True}


if __name__ == "__main__":
    ut.main()
