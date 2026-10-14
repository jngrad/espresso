#
# Copyright (C) 2024 The ESPResSo project
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
Basic tests for two-component color gradient LB.

Tests:
  - Construction with two viscosities triggers two-component mode
  - Density set/get for two-component mode (node and slice)
  - Population get/set for two-component mode (node and slice)
  - Setting velocity raises RuntimeError in two-component mode
  - Pressure tensor matches the analytical equilibrium value (node and
    slice) right after init_two_component()
  - Velocity getter returns sensible values after integration
  - init_two_component produces PDFs consistent with set densities
  - Viscosity setter works on live CG fluid
  - Bulk slice init matches per-node init
  - kT > 0 is accepted for two-component mode and thermalizes it
  - Single viscosity does not create two-component mode
  - Running integration steps without crash
  - Mass conservation (total rho_a and rho_b are conserved)
  - ext_force_density is split between components density-weighted, and
    grows total system momentum at the expected rate
  - A negative component density after streaming stops a thermalized
    integration with an error naming the node
  - A negative component density at kT = 0 gives one warning per
    integrator.run() call and the integration continues
  - A freshly constructed fluid holds zeros in every field
  - A matrix of one component with a trace of the other stays at rest
  - VTK output of phasefield, density_a and density_b matches the slice
    accessors
  - Two-layer Couette flow between a resting and a moving wall matches the
    analytical profile
  - On the GPU, the two-component model runs in double precision

The tests run on the CPU and on the GPU, in double precision.
"""

import contextlib
import io
import pathlib
import tempfile
import unittest as ut
import unittest_decorators as utx
import numpy as np

import espressomd
import espressomd.lb
import espressomd.shapes
try:
    import espressomd.io.vtk
except ImportError:
    pass


DOMAIN_SIZE = 12
AGRID = 1.0
TAU = 1.0
RHO_0 = 1.0
EPSILON = 1e-6
VISCOSITY = 1.0 / 6.0
RADIUS = 4.0
SMOOTHING_WIDTH = 2.0
# largest deviation from the analytical two-layer Couette profile, relative
# to the wall velocity; the diffuse interface gives 6.8e-3 on the CPU
COUETTE_TOLERANCE = 1e-2


def tanh_interpolation(distances, radius, smoothing_width,
                       rho_outer, rho_inner):
    """Smooth tanh interface profile between two density values."""
    return (rho_outer - rho_inner) / 2.0 * (
        np.tanh((distances - radius) / smoothing_width * 2.64665) - 1.0
    ) + rho_outer


def droplet_densities(domain_size, radius, smoothing_width, rho_0, epsilon,
                      agrid=1.0):
    """
    Generate initial density fields for a spherical droplet.

    All length parameters (domain_size, radius, smoothing_width) are in MD units.

    Returns:
        rho_a: solvent density (rho_0 outside, epsilon*rho_0 inside)
        rho_b: droplet density (epsilon*rho_0 outside, rho_0 inside)
    """
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


def couette_profile(x, wall_velocity, mu_a, mu_b, h, H):
    """
    Steady planar Couette flow with a viscosity step at ``x = h`` between a
    resting wall at ``x = 0`` and a wall moving at ``wall_velocity`` at
    ``x = H``. Without a pressure gradient, the shear stress
    ``mu * du/dx`` is the same everywhere in the gap.
    """
    stress = wall_velocity / (h / mu_a + (H - h) / mu_b)
    u_interface = stress / mu_a * h
    return np.where(x <= h, stress / mu_a * x,
                    u_interface + stress / mu_b * (x - h))


class ColorGradientLBTest:
    """Test the two-component color gradient LB method."""

    system = espressomd.System(box_l=[DOMAIN_SIZE] * 3)
    system.time_step = TAU
    system.cell_system.skin = 0.4
    system.periodicity = [True, True, True]

    lb_params = {}

    def tearDown(self):
        self.system.lb = None

    def _create_lbf(self):
        """Create a two-component LB fluid and attach it."""
        lbf = espressomd.lb.LBFluid(
            agrid=AGRID, density=RHO_0, tau=self.system.time_step,
            kinematic_viscosity=[VISCOSITY, VISCOSITY], **self.lb_params)
        self.system.lb = lbf
        return lbf

    def _init_droplet(self, lbf):
        """Set densities to a spherical droplet profile and
        initialize the PDFs."""
        rho_a, rho_b = droplet_densities(
            DOMAIN_SIZE, RADIUS, SMOOTHING_WIDTH, RHO_0, EPSILON, agrid=AGRID)
        lbf[:, :, :].density = np.stack([rho_a, rho_b], axis=-1)
        lbf.init_two_component()

    def test_construction(self):
        """Two viscosities should create a two-component LB."""
        lbf = self._create_lbf()
        visc = lbf.kinematic_viscosity
        self.assertEqual(len(visc), 2)
        self.assertAlmostEqual(visc[0], VISCOSITY, places=10)
        self.assertAlmostEqual(visc[1], VISCOSITY, places=10)

    def test_node_density_two_component(self):
        """Getting density on a node should return two values in CG mode."""
        lbf = self._create_lbf()
        dens = lbf[0, 0, 0].density
        self.assertEqual(len(dens), 2)

    def test_set_density(self):
        """Setting per-node density with two values should round-trip."""
        lbf = self._create_lbf()
        lbf[0, 0, 0].density = [0.7, 0.3]
        dens = lbf[0, 0, 0].density
        self.assertAlmostEqual(dens[0], 0.7, places=10)
        self.assertAlmostEqual(dens[1], 0.3, places=10)

    def test_slice_density_set_get(self):
        """Bulk slice density set/get should round-trip for two components."""
        lbf = self._create_lbf()
        N = int(DOMAIN_SIZE / AGRID)
        rho_a, rho_b = droplet_densities(
            DOMAIN_SIZE, RADIUS, SMOOTHING_WIDTH, RHO_0, EPSILON, agrid=AGRID)

        # Set via bulk slice
        lbf[:, :, :].density = np.stack([rho_a, rho_b], axis=-1)

        # Read back via bulk slice
        densities = np.copy(lbf[:, :, :].density)
        self.assertEqual(densities.shape, (N, N, N, 2))
        np.testing.assert_allclose(densities[:, :, :, 0], rho_a, rtol=1e-10)
        np.testing.assert_allclose(densities[:, :, :, 1], rho_b, rtol=1e-10)

        # Verify consistency with per-node access
        for x in range(0, N, N // 3):
            for y in range(0, N, N // 3):
                for z in range(0, N, N // 3):
                    d = lbf[x, y, z].density
                    self.assertAlmostEqual(d[0], rho_a[x, y, z], places=10)
                    self.assertAlmostEqual(d[1], rho_b[x, y, z], places=10)

    def test_node_population_set_get(self):
        """Node population get/set should round-trip for two components."""
        lbf = self._create_lbf()
        self._init_droplet(lbf)

        # Read populations (2 * stencil_size values)
        pop = np.copy(lbf[3, 3, 3]._population)
        stencil_size = len(pop) // 2
        self.assertEqual(len(pop), 2 * stencil_size)

        # Modify and write back
        pop_modified = pop * 0.99
        lbf[3, 3, 3]._population = pop_modified
        pop_read = np.copy(lbf[3, 3, 3]._population)
        np.testing.assert_allclose(pop_read, pop_modified, rtol=1e-10)

    def test_slice_population_set_get(self):
        """Slice population get/set should round-trip for two components."""
        lbf = self._create_lbf()
        self._init_droplet(lbf)

        N = DOMAIN_SIZE
        pop = np.copy(lbf[:, :, :]._population)
        stencil_size = pop.shape[-1] // 2
        self.assertEqual(pop.shape, (N, N, N, 2 * stencil_size))

        # Verify consistency with per-node access
        node_pop = np.copy(lbf[2, 3, 4]._population)
        np.testing.assert_allclose(pop[2, 3, 4], node_pop, rtol=1e-10)

    def test_set_velocity_raises(self):
        """Setting velocity should raise in two-component mode."""
        lbf = self._create_lbf()
        self._init_droplet(lbf)

        with self.assertRaisesRegex(RuntimeError, "not supported for two-component"):
            lbf[0, 0, 0].velocity = [0.1, 0.0, 0.0]

        with self.assertRaisesRegex(RuntimeError, "not supported for two-component"):
            lbf[:, :, :].velocity = np.zeros((DOMAIN_SIZE,) * 3 + (3,))

    def test_pressure_tensor(self):
        """
        Right after init_two_component(), the fluid is at rest (u = 0) and
        every population sits exactly at its local equilibrium, so the
        pressure tensor P = sum_i f_i c_i c_i must be purely isotropic:
        diagonal = rho_total * c_s**2, off-diagonal = 0 (color_gradient.tex).
        This holds pointwise across the whole droplet -- bulk and interface
        alike -- because both components share the same reference density
        (gamma = 1), so alpha_a = alpha_b = 1/3 everywhere and the
        density-weighted average alpha-bar collapses to 1/3 regardless of
        the local density fractions.
        """
        lbf = self._create_lbf()
        self._init_droplet(lbf)

        c_s_sq = 1.0 / 3.0
        N = int(DOMAIN_SIZE / AGRID)
        densities = np.copy(lbf[:, :, :].density)
        rho_total = densities[..., 0] + densities[..., 1]

        # Node-level: sample points spanning bulk (inside/outside the
        # droplet) and the interface.
        for x in range(0, N, N // 3):
            for y in range(0, N, N // 3):
                for z in range(0, N, N // 3):
                    P = np.copy(lbf[x, y, z].pressure_tensor)
                    self.assertEqual(P.shape, (3, 3))
                    expected_diag = rho_total[x, y, z] * c_s_sq
                    np.testing.assert_allclose(
                        np.diag(P), expected_diag, rtol=1e-9)
                    off_diag = P - np.diag(np.diag(P))
                    np.testing.assert_allclose(off_diag, 0.0, atol=1e-9)

        # Slice-level: shape, symmetry, and consistency with node access.
        P_slice = np.copy(lbf[:, :, :].pressure_tensor)
        self.assertEqual(P_slice.shape, (N, N, N, 3, 3))
        np.testing.assert_allclose(
            P_slice, np.swapaxes(P_slice, -1, -2), atol=1e-10)
        np.testing.assert_allclose(
            P_slice[3, 3, 3], np.copy(lbf[3, 3, 3].pressure_tensor),
            atol=1e-10)
        expected_diag_slice = rho_total * c_s_sq
        for i in range(3):
            np.testing.assert_allclose(
                P_slice[..., i, i], expected_diag_slice, rtol=1e-9)

    def test_velocity_getter(self):
        """Velocity getter should return finite values after integration."""
        lbf = self._create_lbf()
        self._init_droplet(lbf)
        self.system.integrator.run(5)

        # Node velocity
        vel = np.copy(lbf[3, 3, 3].velocity)
        self.assertEqual(len(vel), 3)
        self.assertTrue(np.all(np.isfinite(vel)))

        # Slice velocity
        vel_slice = np.copy(lbf[:, :, :].velocity)
        N = DOMAIN_SIZE
        self.assertEqual(vel_slice.shape, (N, N, N, 3))
        self.assertTrue(np.all(np.isfinite(vel_slice)))

        # Node and slice should be consistent
        np.testing.assert_allclose(vel_slice[3, 3, 3], vel, atol=1e-10)

    def test_init_two_component_consistency(self):
        """After init_two_component, densities should match what was set."""
        lbf = self._create_lbf()
        rho_a, rho_b = droplet_densities(
            DOMAIN_SIZE, RADIUS, SMOOTHING_WIDTH, RHO_0, EPSILON, agrid=AGRID)
        lbf[:, :, :].density = np.stack([rho_a, rho_b], axis=-1)
        lbf.init_two_component()

        # Densities should still match after init
        densities = np.copy(lbf[:, :, :].density)
        np.testing.assert_allclose(densities[:, :, :, 0], rho_a, rtol=1e-5)
        np.testing.assert_allclose(densities[:, :, :, 1], rho_b, rtol=1e-5)

    def test_viscosity_setter(self):
        """Changing viscosity on a live CG fluid should work."""
        lbf = self._create_lbf()
        self._init_droplet(lbf)

        new_visc = [VISCOSITY * 2, VISCOSITY * 3]
        lbf.kinematic_viscosity = new_visc
        visc = lbf.kinematic_viscosity
        self.assertEqual(len(visc), 2)
        self.assertAlmostEqual(visc[0], new_visc[0], places=10)
        self.assertAlmostEqual(visc[1], new_visc[1], places=10)

    def test_bulk_slice_init_matches_node_init(self):
        """Bulk slice density init should match per-node init."""
        lbf = self._create_lbf()
        N = int(DOMAIN_SIZE / AGRID)
        rho_a, rho_b = droplet_densities(
            DOMAIN_SIZE, RADIUS, SMOOTHING_WIDTH, RHO_0, EPSILON, agrid=AGRID)

        # Init via bulk slice
        lbf[:, :, :].density = np.stack([rho_a, rho_b], axis=-1)
        lbf.init_two_component()
        densities_slice = np.copy(lbf[:, :, :].density)

        # Reset and init via per-node loop (to verify both paths)
        self.system.lb = None
        lbf2 = self._create_lbf()
        for x in range(N):
            for y in range(N):
                for z in range(N):
                    lbf2[x, y, z].density = [rho_a[x, y, z],
                                             rho_b[x, y, z]]
        lbf2.init_two_component()
        densities_node = np.copy(lbf2[:, :, :].density)

        np.testing.assert_allclose(
            densities_slice, densities_node, rtol=1e-10)

    def test_thermalized_cg_accepted(self):
        """kT > 0 is accepted for two-component mode and thermalizes it.

        Only checks that the parameters arrive, that the RNG counter is
        live and that the fluid stays finite -- whether the noise satisfies
        fluctuation-dissipation is a separate question.
        """
        lbf = espressomd.lb.LBFluid(
            agrid=AGRID, density=RHO_0, tau=self.system.time_step,
            kinematic_viscosity=[VISCOSITY, VISCOSITY],
            kT=1e-4, seed=42, **self.lb_params)
        self.system.lb = lbf
        lbf[:, :, :].density = np.broadcast_to(
            [0.5 * RHO_0, 0.5 * RHO_0], (DOMAIN_SIZE,) * 3 + (2,))
        lbf.init_two_component()
        self.assertAlmostEqual(lbf.kT, 1e-4, delta=1e-10)
        self.assertEqual(lbf.seed, 42)
        # the RNG counter is exposed once thermalized, and advances per step
        self.assertEqual(lbf.rng_state, 0)
        self.system.integrator.run(2)
        self.assertEqual(lbf.rng_state, 2)
        self.assertTrue(np.all(np.isfinite(np.copy(lbf[:, :, :]._population))))

    def test_unthermalized_cg_has_no_rng_state(self):
        """kT == 0 keeps the two-component model unthermalized."""
        lbf = espressomd.lb.LBFluid(
            agrid=AGRID, density=RHO_0, tau=self.system.time_step,
            kinematic_viscosity=[VISCOSITY, VISCOSITY], **self.lb_params)
        self.system.lb = lbf
        self.assertEqual(lbf.kT, 0.0)
        self.assertIsNone(lbf.rng_state)
        with self.assertRaisesRegex(RuntimeError, "unthermalized"):
            lbf.rng_state = 5

    def test_single_viscosity_not_two_component(self):
        """Single viscosity should not create a two-component LB."""
        lbf = espressomd.lb.LBFluid(
            agrid=AGRID, density=RHO_0, tau=self.system.time_step,
            kinematic_viscosity=VISCOSITY, **self.lb_params)
        self.system.lb = lbf
        # Single-component: density should be a scalar
        dens = lbf[0, 0, 0].density
        self.assertIsInstance(dens, float)

    def test_external_force_density_weighted_split(self):
        """
        ext_force_density must be split between components in proportion to
        local rho_a/rho_b at each node (same rule as add_forces_at_pos), and
        the incremental momentum gained per step -- once the same constant
        force has been active for two consecutive steps, so the force's
        half-step velocity correction is identical on both and cancels in
        the difference -- must equal ext_force_density * volume.
        """
        ext_force_density = np.array([0.0, 0.0021, 0.0])
        lbf = espressomd.lb.LBFluid(
            agrid=AGRID, density=RHO_0, tau=self.system.time_step,
            kinematic_viscosity=[VISCOSITY, VISCOSITY],
            ext_force_density=ext_force_density, **self.lb_params)
        self.system.lb = lbf
        self._init_droplet(lbf)

        np.testing.assert_allclose(
            np.copy(lbf.ext_force_density), ext_force_density)

        # Find a genuinely mixed (interface) node to check the split rule.
        densities = np.copy(lbf[:, :, :].density)
        frac_a = densities[..., 0] / (densities[..., 0] + densities[..., 1])
        x, y, z = (int(v) for v in np.unravel_index(
            np.argmin(np.abs(frac_a - 0.5)), frac_a.shape))
        self.assertAlmostEqual(frac_a[x, y, z], 0.5, delta=0.1)

        self.system.integrator.run(1)
        laf = np.copy(lbf[x, y, z].last_applied_force)
        self.assertEqual(laf.shape, (2, 3))
        np.testing.assert_allclose(
            laf[0] + laf[1], ext_force_density, atol=1e-9)
        self.assertAlmostEqual(
            laf[0][1] / (laf[0][1] + laf[1][1]), frac_a[x, y, z], places=4)

        # Incremental momentum growth once in steady state.
        def total_momentum():
            dens = np.copy(lbf[:, :, :].density)
            rho_total = dens[..., 0] + dens[..., 1]
            velocity = np.copy(lbf[:, :, :].velocity)
            return np.sum(rho_total[..., np.newaxis] * velocity,
                          axis=(0, 1, 2)) * AGRID**3

        mom_1 = total_momentum()
        self.system.integrator.run(1)
        mom_2 = total_momentum()
        expected_delta = ext_force_density * self.system.volume()
        np.testing.assert_allclose(
            mom_2 - mom_1, expected_delta, rtol=1e-6, atol=1e-9)

    def test_sigma_unit_conversion(self):
        """
        Two setups that are identical in lattice units (agrid=1, tau=1 vs.
        agrid=0.5, tau=0.1, with all MD inputs converted accordingly) must
        evolve identically. Surface tension is energy per area [M T^-2],
        i.e. sigma_lb = sigma_md * tau^2, independent of agrid. The droplet
        is curved, so the flow it develops is driven by sigma.
        """
        n_nodes = DOMAIN_SIZE
        visc_lb = VISCOSITY
        sigma_lb = 0.005
        rho_a, rho_b = droplet_densities(
            n_nodes, RADIUS, SMOOTHING_WIDTH, RHO_0, EPSILON)
        rho_lb = np.stack([rho_a, rho_b], axis=-1)

        def restore_system():
            self.system.lb = None
            self.system.box_l = [DOMAIN_SIZE] * 3
            self.system.time_step = TAU
        self.addCleanup(restore_system)

        def run(agrid, tau):
            self.system.lb = None
            self.system.box_l = [n_nodes * agrid] * 3
            self.system.time_step = tau
            lbf = espressomd.lb.LBFluid(
                agrid=agrid, density=RHO_0 / agrid**3, tau=tau,
                kinematic_viscosity=[visc_lb * agrid**2 / tau] * 2,
                sigma=sigma_lb / tau**2, beta=0.8, **self.lb_params)
            self.system.lb = lbf
            lbf[:, :, :].density = rho_lb / agrid**3
            lbf.init_two_component()
            self.system.integrator.run(20)
            populations = np.copy(lbf[:, :, :]._population)
            velocity_lb = np.copy(lbf[:, :, :].velocity) * tau / agrid
            return populations, velocity_lb

        pop_ref, vel_ref = run(agrid=1.0, tau=1.0)
        pop_scaled, vel_scaled = run(agrid=0.5, tau=0.1)

        vel_max = np.max(np.abs(vel_ref))
        self.assertGreater(vel_max, 1e-7, "droplet flow should be nonzero")
        np.testing.assert_allclose(
            vel_scaled, vel_ref, rtol=1e-8, atol=1e-8 * vel_max)
        np.testing.assert_allclose(pop_scaled, pop_ref, rtol=1e-10)

    def test_run_steps(self):
        """Two-component LB with droplet should run without crash."""
        lbf = self._create_lbf()
        self._init_droplet(lbf)
        self.system.integrator.run(10)

    def test_mass_conservation(self):
        """Total mass of each component should be conserved."""
        lbf = self._create_lbf()
        self._init_droplet(lbf)

        # Measure initial total density per component
        densities = np.copy(lbf[:, :, :].density)
        rho_a_init = np.sum(densities[:, :, :, 0])
        rho_b_init = np.sum(densities[:, :, :, 1])

        self.assertGreater(rho_a_init, 0.0, "rho_a should be > 0 after init")
        self.assertGreater(rho_b_init, 0.0, "rho_b should be > 0 after init")

        # Run some steps
        self.system.integrator.run(20)

        # Measure final total density per component
        densities = np.copy(lbf[:, :, :].density)
        rho_a_final = np.sum(densities[:, :, :, 0])
        rho_b_final = np.sum(densities[:, :, :, 1])

        self.assertAlmostEqual(rho_a_init, rho_a_final, places=8,
                               msg="rho_a not conserved")
        self.assertAlmostEqual(rho_b_init, rho_b_final, places=8,
                               msg="rho_b not conserved")

    def _init_negative_density_node(self, lbf, node):
        """
        Uniform bulk with a trace of component a, and a negative density of
        component a at ``node``. The first stream spreads it to ``node`` and
        its 18 neighbours, with ``node`` the most negative one.
        """
        N = int(DOMAIN_SIZE / AGRID)
        densities = np.zeros((N, N, N, 2))
        densities[..., 0] = EPSILON * RHO_0
        densities[..., 1] = RHO_0
        densities[node + (0,)] = -0.1 * RHO_0
        lbf[:, :, :].density = densities
        lbf.init_two_component()

    def test_negative_density_raises_with_kT(self):
        """
        With kT > 0, a negative component density after streaming makes the
        noise amplitude sqrt(kT * rho) undefined. The integration must stop
        with an error naming the node, instead of turning the fluid into NaN.
        """
        lbf = espressomd.lb.LBFluid(
            agrid=AGRID, density=RHO_0, tau=self.system.time_step,
            kinematic_viscosity=[VISCOSITY, VISCOSITY], kT=1e-4, seed=42,
            **self.lb_params)
        self.system.lb = lbf
        node = (5, 6, 7)
        self._init_negative_density_node(lbf, node)
        with contextlib.redirect_stderr(io.StringIO()):
            with self.assertRaisesRegex(
                    Exception,
                    r"ERROR: Color-gradient LB: the populations of component "
                    r"a \(index 0\) sum to a negative density "
                    r"\(rho = -0\.0333\d*\) at node \[5, 6, 7\] in LB step 1 "
                    r"\(19 nodes affected on this MPI rank\)\. "
                    r"Negative densities are not physical\."):
                self.system.integrator.run(1)
        # the negative density stays in place for inspection
        self.assertLess(lbf[node].density[0], 0.)

    def test_negative_density_warns_without_kT(self):
        """
        With kT = 0, a negative component density is still unphysical, but
        it does not produce NaN. The integration continues with a warning,
        reported once per integrator.run() call.
        """
        lbf = self._create_lbf()
        node = (5, 6, 7)
        warning = ("WARNING: Color-gradient LB: the populations of "
                   "component a (index 0) sum to a negative density")
        for _ in range(2):
            self._init_negative_density_node(lbf, node)
            stderr = io.StringIO()
            with contextlib.redirect_stderr(stderr):
                self.system.integrator.run(5)
            self.assertEqual(stderr.getvalue().count(warning), 1)
            self.assertIn("Negative densities are not physical.",
                          stderr.getvalue())
            self.assertTrue(np.all(np.isfinite(
                np.copy(lbf[:, :, :]._population))))
    def test_fresh_fluid_is_zero(self):
        """
        Before any density is set, every field of a new fluid holds zeros.
        Device memory is not zeroed when it is allocated, so on the GPU this
        checks that every field is initialized.
        """
        lbf = self._create_lbf()
        N = DOMAIN_SIZE
        for name, shape in (("density", (N, N, N, 2)),
                            ("phasefield", (N, N, N)),
                            ("_population", (N, N, N, 38)),
                            ("velocity", (N, N, N, 3)),
                            ("last_applied_force", (N, N, N, 2, 3))):
            values = np.copy(getattr(lbf[:, :, :], name))
            self.assertEqual(values.shape, shape, msg=name)
            np.testing.assert_array_equal(values, 0., err_msg=name)

    def test_matrix_with_trace_stays_at_rest(self):
        """
        A uniform matrix of component a with a trace of component b at rest
        is an equilibrium: after a few steps every population, the velocity
        and the phase field are unchanged. (A mixture with comparable
        fractions is unstable in the color-gradient model: any perturbation
        of the phase field starts the demixing.)
        """
        lbf = self._create_lbf()
        N = DOMAIN_SIZE
        lbf[:, :, :].density = np.broadcast_to(
            [(1. - EPSILON) * RHO_0, EPSILON * RHO_0], (N, N, N, 2))
        lbf.init_two_component()
        pop_init = np.copy(lbf[:, :, :]._population)
        phasefield_init = np.copy(lbf[:, :, :].phasefield)
        self.system.integrator.run(5)
        pop = np.copy(lbf[:, :, :]._population)
        self.assertTrue(np.all(np.isfinite(pop)))
        np.testing.assert_allclose(pop, pop_init, rtol=0., atol=1e-14)
        np.testing.assert_allclose(np.copy(lbf[:, :, :].velocity), 0.,
                                   rtol=0., atol=1e-14)
        np.testing.assert_allclose(np.copy(lbf[:, :, :].phasefield),
                                   phasefield_init, rtol=0., atol=1e-14)

    @utx.skipIfMissingModules("espressomd.io.vtk")
    def test_vtk_scalar_fields(self):
        """
        VTK output of the phase field and of both component densities must
        match the slice accessors. Each field has its own output, so the
        three outputs must differ from each other.
        """
        lbf = self._create_lbf()
        self._init_droplet(lbf)
        self.system.integrator.run(2)
        observables = ["phasefield", "density_a", "density_b"]
        with tempfile.TemporaryDirectory() as tmp_directory:
            root = pathlib.Path(tmp_directory)
            vtk_obj = espressomd.lb.VTKOutput(
                identifier="cg_scalars", delta_N=0, observables=observables,
                base_folder=root)
            lbf.add_vtk_writer(vtk=vtk_obj)
            vtk_obj.write()
            filepath = root / "cg_scalars" / "simulation_step_0.vtu"
            self.assertTrue(filepath.exists())
            vtk_data = espressomd.io.vtk.VTKReader().parse(filepath)
        densities = np.copy(lbf[:, :, :].density)
        expected = {"phasefield": np.copy(lbf[:, :, :].phasefield),
                    "density_a": densities[..., 0],
                    "density_b": densities[..., 1]}
        for name in observables:
            # the VTK files are written in single precision
            np.testing.assert_allclose(
                np.reshape(vtk_data[name], expected[name].shape),
                expected[name], rtol=1e-6, atol=1e-7, err_msg=name)
        for name_1, name_2 in (("phasefield", "density_a"),
                               ("phasefield", "density_b"),
                               ("density_a", "density_b")):
            self.assertGreater(
                np.max(np.abs(vtk_data[name_1] - vtk_data[name_2])), 0.1,
                msg=f"VTK outputs '{name_1}' and '{name_2}' are equal")

    def test_two_layer_couette_flow(self):
        """
        Two fluid layers with different viscosities between a resting wall
        and a moving wall. The steady velocity profile is piecewise linear,
        with the same shear stress in both layers. The diffuse interface
        smooths the kink, hence the tolerance.
        """
        wall_velocity = 0.01
        viscosities = [1. / 6., 1. / 2.]

        def restore_system():
            self.system.lb = None
            self.system.box_l = [DOMAIN_SIZE] * 3

        self.addCleanup(restore_system)
        self.system.box_l = [4., 16., 4.]
        lbf = espressomd.lb.LBFluid(
            agrid=AGRID, density=RHO_0, tau=self.system.time_step,
            kinematic_viscosity=viscosities, **self.lb_params)
        self.system.lb = lbf
        # walls normal to y: the boundary nodes y = 0 and y = 15 put the
        # halfway bounce-back walls at y = 1 and y = 15
        lbf.add_boundary_from_shape(
            espressomd.shapes.Wall(normal=[0., 1., 0.], dist=1.))
        lbf.add_boundary_from_shape(
            espressomd.shapes.Wall(normal=[0., -1., 0.], dist=-15.),
            [wall_velocity, 0., 0.])
        # component a below y = 8, component b above
        shape = tuple(lbf.shape)
        y = np.arange(shape[1]) + 0.5
        rho_b = tanh_interpolation(y, 8., SMOOTHING_WIDTH, RHO_0,
                                   EPSILON * RHO_0)
        rho_b = np.broadcast_to(rho_b[np.newaxis, :, np.newaxis], shape)
        rho_a = RHO_0 * (1. + EPSILON) - rho_b
        lbf[:, :, :].density = np.stack([rho_a, rho_b], axis=-1)
        lbf.init_two_component()
        self.system.integrator.run(2000)

        velocity = np.copy(lbf[:, 1:15, :].velocity)
        np.testing.assert_allclose(velocity[..., 1:], 0., rtol=0.,
                                   atol=1e-10)
        profile = np.mean(velocity[..., 0], axis=(0, 2))
        x = y[1:15] - 1.
        expected = couette_profile(x, wall_velocity, RHO_0 * viscosities[0],
                                   RHO_0 * viscosities[1], h=7., H=14.)
        np.testing.assert_allclose(profile, expected, rtol=0.,
                                   atol=COUETTE_TOLERANCE * wall_velocity)


@utx.skipIfMissingFeatures(["WALBERLA"])
class ColorGradientLBTestCPU(ColorGradientLBTest, ut.TestCase):
    lb_params = {"single_precision": False, "gpu": False}


@utx.skipIfMissingGPU()
@utx.skipIfMissingFeatures(["WALBERLA", "CUDA"])
class ColorGradientLBTestGPU(ColorGradientLBTest, ut.TestCase):
    lb_params = {"single_precision": False, "gpu": True}

    def test_double_precision_by_default(self):
        """
        Single precision does not conserve mass and momentum well enough for
        the color-gradient model, so it runs in double precision on the GPU
        unless asked otherwise. The single-component model keeps single
        precision as its GPU default.
        """
        lbf = espressomd.lb.LBFluid(
            agrid=AGRID, density=RHO_0, tau=TAU,
            kinematic_viscosity=[VISCOSITY, VISCOSITY], gpu=True)
        self.assertTrue(lbf.gpu)
        self.assertFalse(lbf.single_precision)
        lbf = espressomd.lb.LBFluid(
            agrid=AGRID, density=RHO_0, tau=TAU,
            kinematic_viscosity=VISCOSITY, gpu=True)
        self.assertTrue(lbf.gpu)
        self.assertTrue(lbf.single_precision)

    def test_single_precision_rejected(self):
        """Single precision is not available for the color-gradient model
        on the GPU."""
        with self.assertRaisesRegex(
                ValueError, "single precision is not supported for the "
                "two-component LB on the GPU"):
            espressomd.lb.LBFluid(
                agrid=AGRID, density=RHO_0, tau=TAU,
                kinematic_viscosity=[VISCOSITY, VISCOSITY], gpu=True,
                single_precision=True)


if __name__ == "__main__":
    ut.main()
