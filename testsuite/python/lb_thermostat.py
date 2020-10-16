# Copyright (C) 2010-2020 The ESPResSo project
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

import espressomd.lb
from tests_common import single_component_maxwell

"""
Check the lattice-Boltzmann thermostat with respect to the particle velocity
distribution.


"""

KT = 0.9 
AGRID = 0.8
VISC = 6 
DENS = 1.7
TIME_STEP = 0.008
GAMMA = 2
LB_PARAMS = {'agrid': AGRID,
             'dens': DENS,
             'visc': VISC,
             'tau': TIME_STEP,
             'kT': KT,
             'seed': 123}


class LBThermostatCommon:

    """Base class of the test that holds the test logic."""
    lbf = None
    system = espressomd.System(box_l=[AGRID * 12] * 3)
    system.time_step = TIME_STEP
    system.cell_system.skin = 0.4 * AGRID

    def test_fluid(self):
        self.prepare()
        self.system.integrator.run(100)
        vs = []
        for _ in range(100):
            for n in self.lbf.nodes():
                vs.append(n.velocity)
            self.system.integrator.run(3)

        self.assertAlmostEqual(np.var(vs), KT, delta=0.05)

    def prepare(self):
        self.system.actors.clear()
        self.system.actors.add(self.lbf)
        self.system.thermostat.set_lb(LB_fluid=self.lbf, seed=5, gamma=GAMMA)

    def test_with_particles(self):
        self.prepare()
        self.system.part.add(
            pos=np.random.random((100, 3)) * self.system.box_l)
        self.system.integrator.run(500)
        N = len(self.system.part)
        loops = 500
        v_particles = np.zeros((N * loops, 3))
        v_nodes = []
        for i in range(loops):
            self.system.integrator.run(3)
            if i % 10 == 0:
                for n in self.lbf.nodes():
                    v_nodes.append(n.velocity)
            v_particles[i * N:(i + 1) * N, :] = self.system.part[:].v
        self.assertAlmostEqual(np.var(v_nodes), KT, delta=0.01)
        np.testing.assert_allclose(np.average(v_particles), 0, atol=0.02)
        np.testing.assert_allclose(np.var(v_particles), KT, rtol=0.03)

        minmax = 3
        n_bins = 7
        for i in range(3):
            hist = np.histogram(v_particles[:, i], range=(-minmax, minmax),
                                bins=n_bins, density=False)
            data = hist[0] / float(v_particles.shape[0])
            bins = hist[1]
            expected = [single_component_maxwell(
                bins[j], bins[j + 1], KT) for j in range(n_bins)]
            np.testing.assert_allclose(data, expected, atol=0.015)


@utx.skipIfMissingFeatures("LB_WALBERLA")
class LBWalberlaThermostat(ut.TestCase, LBThermostatCommon):

    """Test for the CPU implementation of the LB."""

    def setUp(self):
        self.lbf = espressomd.lb.LBFluidWalberla(**LB_PARAMS)


if __name__ == '__main__':
    ut.main()
