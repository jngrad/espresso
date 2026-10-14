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
Two-component color-gradient LB on the GPU compared with the CPU.

The same setup runs on the CPU and on the GPU, with the same seeds, and the
fields and particle forces are compared after a few steps. The kernels of
both architectures are generated from the same model and draw the same
Philox random numbers, so the results only differ by round-off.

Each test adds one feature to a droplet, so that the first failing test
names the part of the GPU implementation that is broken:
  - the droplet alone: stream, collide and recolour, ghost communication
  - thermal fluctuations: random numbers of the collide kernel
  - a resting and a moving wall: boundary kernel
  - an external force density: split between the components
  - the pressure tensor (node, slice, global) and the report of negative
    component densities
  - a particle with friction coupling, a particle with solvation coupling,
    particles within half a cell of the box edge and of the subdomain
    boundary, and many particles at once
  - all features together
A last test runs the GPU setup twice: the results must be bit-identical.
Both architectures run in double precision.
"""

import contextlib
import io
import re
import unittest as ut
import unittest_decorators as utx
import numpy as np

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
EPSILON = 1e-4
WALL_VELOCITY = [0., 1e-3, 0.]
EXT_FORCE_DENSITY = [0., 0., 2e-3]
N_STEPS = 50

FIELDS = ("population", "density", "phasefield", "velocity",
          "last_applied_force")


def droplet_densities(n):
    """Droplet of component b in a matrix of component a."""
    x = (np.arange(n) + 0.5) * AGRID
    xx, yy, zz = np.meshgrid(x, x, x, indexing="ij")
    center = BOX_L / 2.
    dist = np.sqrt((xx - center)**2 + (yy - center)**2 + (zz - center)**2)
    inside = 0.5 * (1. - np.tanh((dist - DROPLET_RADIUS) / AGRID))
    rho_b = DENSITY * (EPSILON + (1. - EPSILON) * inside)
    rho_a = DENSITY * (1. + EPSILON) - rho_b
    return np.stack([rho_a, rho_b], axis=-1)


def max_deviation(ref, new, min_scale=0.):
    """
    Largest deviation of each field, relative to the field's largest value,
    or to ``min_scale`` if that is larger.
    """
    out = {}
    for key, value in ref.items():
        value = np.asarray(value, dtype=float)
        scale = max(np.max(np.abs(value)), min_scale) if value.size else 0.
        diff = np.max(np.abs(np.asarray(new[key], dtype=float) - value)) \
            if value.size else 0.
        out[key] = diff / scale if scale > 0. else diff
    return out


# particle near the droplet interface, moving
INTERFACE_POS = [BOX_L / 2. + DROPLET_RADIUS, BOX_L / 2. + 0.13,
                 BOX_L / 2. - 0.07]
PARTICLE_VELOCITY = [0.05, 0.02, -0.01]


class ColorGradientGPUParity:
    """Compare a GPU fluid with a CPU fluid."""

    system = espressomd.System(box_l=[BOX_L] * 3)
    system.time_step = TIME_STEP
    system.cell_system.skin = 0.4

    ref_params = {}
    test_params = {}
    # largest deviation, relative to the largest value of each field
    rtol = 0.

    def tearDown(self):
        self.clear()

    def clear(self):
        self.system.part.clear()
        self.system.thermostat.turn_off()
        self.system.lb = None
        self.system.time = 0.

    def make_fluid(self, lb_params, kT=0., walls=False, ext_force=False):
        params = {"agrid": AGRID, "density": DENSITY,
                  "kinematic_viscosity": VISCOSITY, "tau": TIME_STEP,
                  "sigma": SIGMA, "beta": BETA}
        if kT > 0.:
            params.update(kT=kT, seed=LB_SEED)
        if ext_force:
            params["ext_force_density"] = EXT_FORCE_DENSITY
        lbf = espressomd.lb.LBFluid(**params, **lb_params)
        if walls:
            # a resting wall and a moving wall normal to x
            wall1 = espressomd.shapes.Wall(normal=[1., 0., 0.], dist=AGRID)
            wall2 = espressomd.shapes.Wall(normal=[-1., 0., 0.],
                                           dist=-(BOX_L - AGRID))
            lbf.add_boundary_from_shape(wall1)
            lbf.add_boundary_from_shape(wall2, WALL_VELOCITY)
        return lbf

    def run_setup(self, lb_params, kT=0., walls=False, ext_force=False,
                  particles=(), n_steps=N_STEPS):
        """
        Run the setup and return its state. ``particles`` holds tuples of
        position, velocity and solvation chemical potential difference.
        """
        self.clear()
        system = self.system
        lbf = self.make_fluid(lb_params, kT=kT, walls=walls,
                              ext_force=ext_force)
        system.lb = lbf
        lbf[:, :, :].density = droplet_densities(lbf.shape[0])
        lbf.init_two_component()
        if particles:
            system.thermostat.set_lb(LB_fluid=lbf, gamma=GAMMA,
                                     seed=THERMOSTAT_SEED)
            system.thermostat.lb.call_method("override_philox_counter",
                                             counter=0)
            for pos, v, delta_mu in particles:
                system.part.add(pos=pos, v=v, solvation_delta_mu=delta_mu)
        system.integrator.run(n_steps)
        state = {
            "population": np.copy(lbf[:, :, :]._population),
            "density": np.copy(lbf[:, :, :].density),
            "phasefield": np.copy(lbf[:, :, :].phasefield),
            "velocity": np.copy(lbf[:, :, :].velocity),
            "last_applied_force": np.copy(lbf[:, :, :].last_applied_force),
        }
        if particles:
            parts = system.part.all()
            state["pos"] = np.copy(parts.pos)
            state["v"] = np.copy(parts.v)
            state["f"] = np.copy(parts.f)
        return state

    def compare(self, ref, new, rtol=None, min_scale=0.):
        rtol = self.rtol if rtol is None else rtol
        self.assertEqual(sorted(ref), sorted(new))
        deviations = max_deviation(ref, new, min_scale)
        for key in ref:
            self.assertEqual(np.shape(new[key]), np.shape(ref[key]), msg=key)
            self.assertTrue(np.all(np.isfinite(new[key])), msg=key)
            self.assertLessEqual(
                deviations[key], rtol,
                msg=f"'{key}' differs between CPU and GPU by "
                    f"{deviations[key]:.3g} of its largest value")

    def check_parity(self, min_scale=0., **kwargs):
        ref = self.run_setup(self.ref_params, **kwargs)
        new = self.run_setup(self.test_params, **kwargs)
        self.compare(ref, new, min_scale=min_scale)

    def test_01_droplet(self):
        # nothing drives a flow, so the velocity is only the spurious current
        # (|u| ~ 1e-5): compare it against a scale of at least 1e-3, otherwise
        # round-off (~1e-15) exceeds the relative tolerance
        self.check_parity(min_scale=1e-3)

    def test_02_thermalized(self):
        self.check_parity(kT=KT)

    def test_03_walls(self):
        self.check_parity(walls=True)

    def test_04_external_force(self):
        self.check_parity(ext_force=True)

    def test_05_pressure_tensor(self):
        def pressure_tensors(lb_params):
            self.run_setup(lb_params, kT=KT)
            lbf = self.system.lb
            nodes = [(0, 0, 0), (3, 6, 6), (6, 6, 6), (11, 2, 9)]
            return {
                "slice": np.copy(lbf[:, :, :].pressure_tensor),
                "node": np.array([np.copy(lbf[node].pressure_tensor)
                                  for node in nodes]),
                "global": np.copy(lbf.pressure_tensor),
            }
        self.compare(pressure_tensors(self.ref_params),
                     pressure_tensors(self.test_params))

    def test_06_negative_density_report(self):
        def report(lb_params):
            self.clear()
            lbf = self.make_fluid(lb_params, kT=KT)
            self.system.lb = lbf
            n = lbf.shape[0]
            densities = np.zeros((n, n, n, 2))
            densities[..., 0] = EPSILON * DENSITY
            densities[..., 1] = DENSITY
            densities[5, 6, 7, 0] = -0.1 * DENSITY
            lbf[:, :, :].density = densities
            lbf.init_two_component()
            with contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(Exception) as context:
                    self.system.integrator.run(1)
            return str(context.exception)

        # the density value may differ in the last digits, since CPU and
        # GPU may sum the populations in a different order
        pattern = re.compile(r"\(rho = (-?[0-9.e+-]+)\)")
        message_ref = report(self.ref_params)
        message = report(self.test_params)
        self.assertIn("at node [5, 6, 7] in LB step 1", message)
        self.assertEqual(pattern.sub("(rho = X)", message),
                         pattern.sub("(rho = X)", message_ref))
        rho_ref = float(pattern.search(message_ref).group(1))
        rho = float(pattern.search(message).group(1))
        self.assertLessEqual(abs(rho - rho_ref), self.rtol * abs(rho_ref))

    def test_07_friction_particle(self):
        self.check_parity(particles=[(INTERFACE_POS, PARTICLE_VELOCITY, 0.)])

    def test_08_solvation_particle(self):
        self.check_parity(particles=[(INTERFACE_POS, [0., 0., 0.], 2.)])

    def test_09_edge_particles(self):
        # within half a cell of the box edge, and of the subdomain boundary
        # when the box is split in two along x
        edge = 0.2 * AGRID
        particles = [
            ([edge, BOX_L / 2. + DROPLET_RADIUS, BOX_L / 2.],
             PARTICLE_VELOCITY, 2.),
            ([BOX_L - edge, BOX_L / 2., BOX_L / 2. - DROPLET_RADIUS],
             PARTICLE_VELOCITY, -1.),
            ([BOX_L / 2. - edge, BOX_L / 2. - DROPLET_RADIUS, BOX_L / 2.],
             PARTICLE_VELOCITY, 1.),
            ([BOX_L / 2. + edge, BOX_L / 2., BOX_L / 2. + DROPLET_RADIUS],
             PARTICLE_VELOCITY, 0.5),
        ]
        self.check_parity(particles=particles)

    def test_10_many_particles(self):
        rng = np.random.default_rng(seed=12)
        particles = [(rng.uniform(0., BOX_L, 3), rng.uniform(-0.05, 0.05, 3),
                      rng.uniform(-2., 2.)) for _ in range(50)]
        self.check_parity(particles=particles, n_steps=10)

    def test_11_all_features(self):
        self.check_parity(
            kT=KT, walls=True, ext_force=True,
            particles=[(INTERFACE_POS, PARTICLE_VELOCITY, 2.)])

    def test_12_gpu_reproducible(self):
        kwargs = {"kT": KT, "walls": True, "ext_force": True,
                  "particles": [(INTERFACE_POS, PARTICLE_VELOCITY, 2.)]}
        first = self.run_setup(self.test_params, **kwargs)
        second = self.run_setup(self.test_params, **kwargs)
        for key in first:
            np.testing.assert_array_equal(second[key], first[key],
                                          err_msg=key)


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
class ColorGradientGPUParityTest(ColorGradientGPUParity, ut.TestCase):
    ref_params = {"single_precision": False, "gpu": False}
    test_params = {"single_precision": False, "gpu": True}
    rtol = 1e-10


if __name__ == "__main__":
    ut.main()
