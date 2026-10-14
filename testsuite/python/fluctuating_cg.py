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
Thermal fluctuation tests for two-component color gradient LB.

Tests:
  - Total mass per component and total momentum are conserved under
    thermal fluctuations
  - Fluid temperature (global and per axis) fluctuates around the
    pre-set kT
  - Mass-weighted node velocities follow a Maxwell-Boltzmann distribution
    with variance kT

The tests run on the CPU and on the GPU, in double precision.
"""

import unittest as ut
import unittest_decorators as utx
import numpy as np

import espressomd
import espressomd.lb
import thermostats_common


BOX_L = 8.0
TIME_STEP = 0.01
AGRID = 0.5
VISCOSITY = [2.0, 6.0]
DENSITY = 0.86
SIGMA = 0.005
BETA = 0.8
TEMPERATURE = 4e-4

SMOOTHING_WIDTH = 2.0
EPSILON = 1e-12

WARMUP_STEPS = 200
STEPS_PER_FRAME = 5

# frames used for the mass/momentum/temperature checks; kept small since
# each frame performs several assertions with a fixed statistical
# confidence, so the more frames, the higher the chance of a single
# spurious failure (a multiple-comparisons effect)
N_FRAMES = 10
# additional sampling-only frames (no per-frame assertions) used to build
# up enough statistics for the velocity-distribution histogram below;
# chosen so the pooled sample has >=15 counts even in the outermost bin at
# N_BINS=11, range +-5*sqrt(kT), while keeping the added runtime to a few
# seconds
N_EXTRA_FRAMES = 116
N_BINS = 11
# flat per-bin tolerance for check_velocity_distribution, chosen
# conservatively above the ~5 sigma binomial standard error of the
# tightest (central) bin at the sample size above
ERROR_TOL = 0.004
Z = 3


def tanh_interpolation(x, y1, y2, smoothing_width):
    return (y2 - y1) / 2.0 * (np.tanh(x / smoothing_width * 2.64665) + 1.0) + y1


def plane_interface(grid_size, smoothing_width, rho_0, epsilon):
    """Generate initial density fields for a planar interface."""
    grid_size = int(grid_size)
    center = grid_size / 2.0
    # signed distance from the interface: negative on left, positive on right
    signed = np.broadcast_to(
        (np.arange(grid_size) + 0.5 - center)[:, np.newaxis, np.newaxis],
        (grid_size, grid_size, grid_size))
    rho_a = tanh_interpolation(signed, epsilon * rho_0, rho_0, smoothing_width)
    rho_b = tanh_interpolation(signed, rho_0, epsilon * rho_0, smoothing_width)
    return rho_a, rho_b


class FluctuatingColorGradientTest(thermostats_common.ThermostatsCommon):
    """Test thermal fluctuations of the two-component color gradient LB."""

    system = espressomd.System(box_l=[BOX_L] * 3)
    system.time_step = TIME_STEP
    system.cell_system.skin = 0.4

    lb_params = {}

    def setUp(self):
        self.lbf = espressomd.lb.LBFluid(
            agrid=AGRID, density=DENSITY, kinematic_viscosity=VISCOSITY,
            tau=TIME_STEP, sigma=SIGMA, beta=BETA, kT=TEMPERATURE,
            **self.lb_params)
        self.system.lb = self.lbf

        rho_b, rho_a = plane_interface(
            BOX_L / AGRID, SMOOTHING_WIDTH, DENSITY, EPSILON)
        self.lbf[:, :, :].density = np.stack([rho_a, rho_b], axis=-1)
        self.lbf.init_two_component()

        self.system.integrator.run(WARMUP_STEPS)

    def tearDown(self):
        self.system.lb = None

    def test_conservation(self):
        """Total mass per component and total momentum are conserved
        (lattice units)."""
        densities = np.copy(self.lbf[:, :, :].density)
        rho_a_init = np.sum(densities[..., 0])
        rho_b_init = np.sum(densities[..., 1])

        density = np.sum(densities, axis=-1)
        velocity = np.copy(self.lbf[:, :, :].velocity)
        momentum_init = np.sum(
            density[..., np.newaxis] * velocity, axis=(0, 1, 2))

        for _ in range(N_FRAMES):
            self.system.integrator.run(STEPS_PER_FRAME)

            # mass conservation: density per cell fluctuates, but total
            # mass per component is conserved
            densities = np.copy(self.lbf[:, :, :].density)
            rho_a = np.sum(densities[..., 0])
            rho_b = np.sum(densities[..., 1])
            np.testing.assert_allclose(rho_a_init, rho_a, atol=3e-12)
            np.testing.assert_allclose(rho_b_init, rho_b, atol=3e-12)

            # momentum conservation
            density = np.sum(densities, axis=-1)
            velocity = np.copy(self.lbf[:, :, :].velocity)
            momentum = np.sum(
                density[..., np.newaxis] * velocity, axis=(0, 1, 2))
            np.testing.assert_allclose(momentum_init, momentum, atol=3e-12)

    def _sample_frame(self):
        """Mass-weighted node velocity sqrt(rho * a^3) * v for the current
        state; should follow a Maxwell-Boltzmann distribution with unit
        mass and variance kT."""
        velocities = np.copy(self.lbf[:, :, :].velocity)
        densities = np.copy(self.lbf[:, :, :].density)
        density = np.sum(densities, axis=-1)
        return np.sqrt(density * AGRID**3)[..., np.newaxis] * velocities

    def test_fluctuations(self):
        """Fluid temperature fluctuates around the pre-set kT, and the
        mass-weighted node velocities follow a Maxwell-Boltzmann
        distribution."""
        w_samples = []
        for _ in range(N_FRAMES):
            self.system.integrator.run(STEPS_PER_FRAME)
            w = self._sample_frame()
            m_v2 = w**2

            kT_measured = np.mean(m_v2)
            kT_per_axis = np.mean(m_v2, axis=(0, 1, 2))

            # account for momentum conservation removing one degree of
            # freedom
            n_nodes = np.prod(w.shape[:-1])
            kT_measured *= n_nodes / (n_nodes - 1)
            kT_per_axis *= n_nodes / (n_nodes - 1)

            standard_error = np.sqrt(2 / (3 * (n_nodes - 1)))
            standard_error_axis = np.sqrt(2 / (n_nodes - 1))
            np.testing.assert_allclose(
                kT_measured, TEMPERATURE, rtol=Z * standard_error)
            np.testing.assert_allclose(
                kT_per_axis, [TEMPERATURE] * 3, rtol=Z * standard_error_axis)

            w_samples.append(w.reshape((-1, 3)))

        # collect further samples purely for the velocity-distribution
        # histogram below (no per-frame assertions here, so we can afford
        # many more frames without a multiple-comparisons false-failure
        # risk)
        for _ in range(N_EXTRA_FRAMES):
            self.system.integrator.run(STEPS_PER_FRAME)
            w_samples.append(self._sample_frame().reshape((-1, 3)))

        w_all = np.concatenate(w_samples, axis=0)
        v_minmax = 5.0 * np.sqrt(TEMPERATURE)
        self.check_velocity_distribution(
            w_all, v_minmax, N_BINS, ERROR_TOL, TEMPERATURE)


@utx.skipIfMissingFeatures(["WALBERLA"])
class FluctuatingColorGradientTestCPU(FluctuatingColorGradientTest,
                                      ut.TestCase):
    lb_params = {"single_precision": False, "gpu": False}


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
class FluctuatingColorGradientTestGPU(FluctuatingColorGradientTest,
                                      ut.TestCase):
    lb_params = {"single_precision": False, "gpu": True}


if __name__ == "__main__":
    ut.main()
