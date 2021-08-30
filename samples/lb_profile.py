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
"""
Simulate the flow of a lattice-Boltzmann fluid past a cylinder,
obtain the velocity profile in polar coordinates and compare it
to the analytical solution.
"""
import numpy as np
import matplotlib.pyplot as plt

import espressomd
import espressomd.lb
import espressomd.observables
import espressomd.shapes
import espressomd.lbboundaries
import espressomd.accumulators
import espressomd.math

required_features = ["LB_WALBERLA", "LB_BOUNDARIES"]
espressomd.assert_features(required_features)

system = espressomd.System(box_l=[10.0, 10.0, 5.0])
system.time_step = 0.01
system.cell_system.skin = 0.4
n_steps = 500

lb_fluid = espressomd.lb.LBFluidWalberla(
    agrid=1.0, dens=1.0, visc=1.0, tau=0.01, ext_force_density=[0, 0, 0.15], kT=0.0)
system.actors.add(lb_fluid)
system.thermostat.set_lb(LB_fluid=lb_fluid, seed=23)
ctp = espressomd.math.CylindricalTransformationParameters(
    center=[5.0, 5.0, 0.0],
    axis=[0, 0, 1],
    orientation=[1, 0, 0])
fluid_obs = espressomd.observables.CylindricalLBVelocityProfile(
    transform_params=ctp,
    n_r_bins=100,
    max_r=4.0,
    min_z=0.0,
    max_z=10.0,
    sampling_density=0.1)
cylinder_shape = espressomd.shapes.Cylinder(
    center=[5.0, 5.0, 5.0],
    axis=[0, 0, 1],
    direction=-1,
    radius=4.0,
    length=20.0)
cylinder_boundary = espressomd.lbboundaries.LBBoundary(shape=cylinder_shape)
system.lbboundaries.add(cylinder_boundary)
system.integrator.run(n_steps)


accumulator = espressomd.accumulators.MeanVarianceCalculator(obs=fluid_obs)
system.auto_update_accumulators.add(accumulator)
system.integrator.run(n_steps)

lb_fluid_profile = accumulator.mean()


def poiseuille_flow(r, R, ext_force_density):
    return ext_force_density * 1. / 4 * (R**2.0 - r**2.0)


# Please note that due to symmetry and interpolation, a plateau is seen
# near r=0.
n_bins = len(lb_fluid_profile[:, 0, 0, 2])
r_max = 4.0
r = np.linspace(0.0, r_max, n_bins)
plt.plot(r, lb_fluid_profile[:, 0, 0, 2], label='LB profile')
plt.plot(r, poiseuille_flow(r, r_max, 0.15), label='analytical solution')
plt.legend()
plt.show()
