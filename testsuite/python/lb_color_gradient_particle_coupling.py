#
# Copyright (C) 2024-2026 The ESPResSo project
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
Tests for particle-fluid coupling in two-component color gradient LB.

Tests:
  - Drag force on particle matches gamma * (u_fluid - v_particle)
  - Mass conservation is maintained with particle coupling active
  - Solvation force: zero at the droplet center and for delta_mu = 0,
    mainly radial at the interface, flips sign with delta_mu
  - Total momentum is conserved with solvation coupling active
  - Simulation runs without crash for multiple steps with coupling

The tests run on the CPU and on the GPU, in double precision.
"""

import unittest as ut
import unittest_decorators as utx
import numpy as np

import espressomd
import espressomd.lb


DOMAIN_SIZE = 12
AGRID = 1.0
TAU = 1.0
RHO_0 = 1.0
EPSILON = 1e-6
VISCOSITY = 1.0 / 6.0
RADIUS = 4.0
SMOOTHING_WIDTH = 2.0
GAMMA = 0.5
SIGMA = 0.005 
BETA = 0.8


def tanh_interpolation(distances, radius, smoothing_width,
                       rho_outer, rho_inner):
    """Smooth tanh interface profile between two density values."""
    return (rho_outer - rho_inner) / 2.0 * (
        np.tanh((distances - radius) / smoothing_width * 2.64665) - 1.0
    ) + rho_outer


def droplet_densities(domain_size, radius, smoothing_width, rho_0, epsilon,
                      agrid=1.0):
    """Generate initial density fields for a spherical droplet.
    All length parameters are in MD units."""
    grid_size = int(domain_size / agrid)
    x = np.arange(grid_size) + 0.5
    xx, yy, zz = np.meshgrid(x, x, x, indexing="ij")
    center = grid_size / 2.0
    distances = np.sqrt((xx - center)**2 + (yy - center)**2
                        + (zz - center)**2)
    radius = radius / agrid
    smoothing_width = smoothing_width / agrid

    rho_a = tanh_interpolation(distances, radius, smoothing_width,
                               rho_0, epsilon * rho_0)
    rho_b = tanh_interpolation(distances, radius, smoothing_width,
                               epsilon * rho_0, rho_0)
    return rho_a, rho_b

class ColorGradientParticleCouplingTest:
    """Test particle-fluid coupling for two-component LB."""

    system = espressomd.System(box_l=[DOMAIN_SIZE] * 3)
    system.time_step = TAU
    system.cell_system.skin = 0.4
    system.periodicity = [True, True, True]

    lb_params = {}

    def tearDown(self):
        self.system.lb = None
        self.system.thermostat.turn_off()
        self.system.part.clear()

    def _create_lbf(self):
        """Create a two-component LB fluid and attach it."""
        lbf = espressomd.lb.LBFluid(
            agrid=AGRID, density=RHO_0, tau=self.system.time_step,
            kinematic_viscosity=[VISCOSITY, VISCOSITY], sigma=SIGMA, beta=BETA,
            **self.lb_params)
        self.system.lb = lbf
        return lbf

    def _init_droplet(self, lbf):
        """Set densities to a spherical droplet profile and initialize the PDFs."""
        rho_a, rho_b = droplet_densities(
            DOMAIN_SIZE, RADIUS, SMOOTHING_WIDTH, RHO_0, EPSILON, agrid=AGRID)
        lbf[:, :, :].density = np.stack([rho_a, rho_b], axis=-1)
        lbf.init_two_component()

    def _test_drag_force_at_pos(self, pos, label):
        """Drag force on a particle in two-component fluid
        should be gamma * (u_fluid - v_particle) = -gamma * v_particle."""
        lbf = self._create_lbf()
        self._init_droplet(lbf)

        v_particle = np.array([0.01, -0.02, 0.015])
        p = self.system.part.add(pos=pos, v=v_particle)

        # Enable LB thermostat (kT=0, no noise)
        self.system.thermostat.set_lb(LB_fluid=lbf, gamma=GAMMA, seed=42)

        # Run one step to trigger the coupling
        self.system.integrator.run(1)

        # Expected drag force: gamma * (u_fluid - v_particle)
        # Fluid is at rest, so F = -gamma * v_particle
        expected_force = -GAMMA * v_particle

        np.testing.assert_allclose(
            np.copy(p.f), expected_force, atol=1e-10,
            err_msg=f"Drag force mismatch for particle at {label}")

    def test_drag_force_inside_droplet(self):
        """Drag force at droplet center (component b dominant)."""
        center = DOMAIN_SIZE / 2.0
        self._test_drag_force_at_pos(
            [center] * 3, "droplet center")

    def test_drag_force_outside_droplet(self):
        """Drag force outside droplet (component a dominant)."""
        center = DOMAIN_SIZE / 2.0
        self._test_drag_force_at_pos(
            [center+RADIUS+0.5, center+RADIUS+0.5, center+RADIUS+0.5], "outside droplet")

    def test_drag_force_at_interface(self):
        """Drag force at the droplet interface (mixed densities)."""
        center = DOMAIN_SIZE / 2.0
        self._test_drag_force_at_pos(
            [center + RADIUS, center, center], "interface")

    def test_mass_conservation_with_particle(self):
        """Total mass of each component should be conserved even with
        particle coupling active."""
        lbf = self._create_lbf()
        self._init_droplet(lbf)

        # Place particle inside the droplet
        center = [DOMAIN_SIZE / 2.0] * 3
        p = self.system.part.add(pos=center, v=[0.01, 0.0, 0.0])
        self.system.thermostat.set_lb(LB_fluid=lbf, gamma=GAMMA, seed=42)

        # Measure initial masses
        densities = np.copy(lbf[:, :, :].density)
        mass_a_init = np.sum(densities[:, :, :, 0])
        mass_b_init = np.sum(densities[:, :, :, 1])

        # Run with coupling
        self.system.integrator.run(20)

        # Measure final masses
        densities = np.copy(lbf[:, :, :].density)
        mass_a_final = np.sum(densities[:, :, :, 0])
        mass_b_final = np.sum(densities[:, :, :, 1])

        self.assertAlmostEqual(mass_a_init, mass_a_final, places=8,
                               msg="rho_a not conserved with particle coupling")
        self.assertAlmostEqual(mass_b_init, mass_b_final, places=8,
                               msg="rho_b not conserved with particle coupling")

    def test_solvation_force_at_center(self):
        """At the droplet center the color gradient is zero, so solvation force should be zero regardless of delta_mu."""
        lbf = self._create_lbf()
        self._init_droplet(lbf)
        # Run one step to compute the color gradient field
        self.system.integrator.run(1)

        # Particle at rest at droplet center with nonzero delta_mu
        center = [DOMAIN_SIZE / 2.0] * 3
        p = self.system.part.add(pos=center, v=[0, 0, 0],
                                 solvation_delta_mu=1.0)
        self.system.thermostat.set_lb(LB_fluid=lbf, gamma=GAMMA, seed=42)
        self.system.integrator.run(1)

        # No drag (v=0), no solvation force (grad_phi=0 at center)
        np.testing.assert_allclose(
            np.copy(p.f), [0, 0, 0], atol=1e-10,
            err_msg="Force should be zero at droplet center")

    def test_solvation_force_at_interface(self):
        """At the interface, solvation force should be nonzero and point radially (along the color gradient)."""
        lbf = self._create_lbf()
        self._init_droplet(lbf)
        # Run one step to compute the color gradient field
        self.system.integrator.run(1)

        # Particle at rest on the interface along x-axis
        center = DOMAIN_SIZE / 2.0
        pos_interface = [center + RADIUS, center, center]
        delta_mu = 2.0
        p = self.system.part.add(pos=pos_interface, v=[0, 0, 0],
                                 solvation_delta_mu=delta_mu)
        self.system.thermostat.set_lb(LB_fluid=lbf, gamma=GAMMA, seed=42)
        self.system.integrator.run(1)

        force = np.copy(p.f)
        # Solvation force should be nonzero
        self.assertGreater(np.linalg.norm(force), 1e-6,
                           "Solvation force should be nonzero at interface")
        # Force should be predominantly along x (radial direction)
        # since the color gradient points radially at this position
        self.assertGreater(abs(force[0]), abs(force[1]),
                           "Solvation force should be mainly radial (x)")
        self.assertGreater(abs(force[0]), abs(force[2]),
                           "Solvation force should be mainly radial (x)")

    def test_solvation_force_sign(self):
        """Flipping the sign of delta_mu should flip the solvation force."""
        lbf = self._create_lbf()
        self._init_droplet(lbf)

        # Run one step to compute the color gradient field
        self.system.integrator.run(1)

        center = DOMAIN_SIZE / 2.0
        pos_interface = [center + RADIUS, center, center]

        # Positive delta_mu
        p1 = self.system.part.add(pos=pos_interface, v=[0, 0, 0],
                                  solvation_delta_mu=1.0)
        self.system.thermostat.set_lb(LB_fluid=lbf, gamma=GAMMA, seed=42)
        self.system.integrator.run(1)
        force_pos = np.copy(p1.f)
        
        # Reset
        self.system.part.clear()
        self.system.thermostat.turn_off()
        self.system.lb = None

        # Negative delta_mu at same position
        lbf = self._create_lbf()
        self._init_droplet(lbf)

        # Run one step to compute the color gradient field
        self.system.integrator.run(1)

        p2 = self.system.part.add(pos=pos_interface, v=[0, 0, 0],
                                  solvation_delta_mu=-1.0)
        self.system.thermostat.set_lb(LB_fluid=lbf, gamma=GAMMA, seed=42)
        self.system.integrator.run(1)
        force_neg = np.copy(p2.f)

        # Forces should be opposite
        np.testing.assert_allclose(
            force_pos, -force_neg, atol=1e-10,
            err_msg="Flipping delta_mu should flip the solvation force")

    def test_solvation_force_zero_delta_mu(self):
        """With delta_mu=0, solvation force should be zero even at interface."""
        lbf = self._create_lbf()
        self._init_droplet(lbf)

        # Run one step to compute the color gradient field
        self.system.integrator.run(1)

        center = DOMAIN_SIZE / 2.0
        pos_interface = [center + RADIUS, center, center]
        p = self.system.part.add(pos=pos_interface, v=[0, 0, 0],
                                 solvation_delta_mu=0.0)
        self.system.thermostat.set_lb(LB_fluid=lbf, gamma=GAMMA, seed=42)
        self.system.integrator.run(1)

        np.testing.assert_allclose(
            np.copy(p.f), [0, 0, 0], atol=1e-10,
            err_msg="Force should be zero when delta_mu=0")

    def test_solvation_momentum_conservation(self):
        """Total momentum (particle + fluid) should be approximately conserved when solvation force coupling is active."""
        lbf = self._create_lbf()
        self._init_droplet(lbf)
        
        # Run one step to compute color gradient
        self.system.integrator.run(1)

        center = DOMAIN_SIZE / 2.0
        pos_interface = [center + RADIUS, center, center]
        p = self.system.part.add(pos=pos_interface, v=[0, 0, 0],
                                 solvation_delta_mu=2.0)
        self.system.thermostat.set_lb(LB_fluid=lbf, gamma=GAMMA, seed=42)

        total_momentum = []

        for i in range(1000):
            self.system.integrator.run(1)
            # Measure total momentum
            particle_momentum = np.copy(p.v) * p.mass
            densities = np.copy(lbf[:, :, :].density)
            velocities = np.copy(lbf[:, :, :].velocity)
            rho_total = densities[:, :, :, 0] + densities[:, :, :, 1]
            # density is mass/volume in MD units; multiply by cell volume to get
            # mass per cell so units match the MD particle momentum
            fluid_momentum = np.sum(
                rho_total[:, :, :, np.newaxis] * velocities * AGRID**3, axis=(0, 1, 2))
            total_momentum.append(particle_momentum + fluid_momentum)

            if i > 0:
                np.testing.assert_allclose(
                    total_momentum[-1], total_momentum[0], atol=1e-10,
                    err_msg="Momentum not conserved with solvation force coupling")
                

    def test_run_with_coupling(self):
        """Two-component LB with particle coupling should run
        without crash for multiple steps."""
        lbf = self._create_lbf()
        self._init_droplet(lbf)

        # Add particle and enable coupling
        center = [DOMAIN_SIZE / 2.0] * 3
        p = self.system.part.add(pos=center, v=[0.005, -0.005, 0.003])
        self.system.thermostat.set_lb(LB_fluid=lbf, gamma=GAMMA, seed=42)

        # Run 50 steps — should not crash
        self.system.integrator.run(50)

        # Particle should still have finite position and velocity
        self.assertTrue(np.all(np.isfinite(np.copy(p.pos))))
        self.assertTrue(np.all(np.isfinite(np.copy(p.v))))

        # Fluid velocities should be finite
        vel = np.copy(lbf[:, :, :].velocity)
        self.assertTrue(np.all(np.isfinite(vel)))


@utx.skipIfMissingFeatures(["WALBERLA"])
class ColorGradientParticleCouplingTestCPU(ColorGradientParticleCouplingTest,
                                           ut.TestCase):
    lb_params = {"single_precision": False, "gpu": False}


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
class ColorGradientParticleCouplingTestGPU(ColorGradientParticleCouplingTest,
                                           ut.TestCase):
    lb_params = {"single_precision": False, "gpu": True}


if __name__ == "__main__":
    ut.main()
