#
# Copyright (C) 2013-2019 The ESPResSo project
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
import espressomd
import numpy as np
from espressomd.interactions import OifOutDirection, MembraneCollisionInteraction
import tests_common


class OifMembraneCollision(ut.TestCase):
    system = espressomd.System(box_l=[10.]*3)
    system.periodicity = 0, 0, 0
    out_dir = OifOutDirection()
    system.bonded_inter.add(out_dir)

    system.time_step = 0.01
    system.cell_system.skin = 0
    mc_a = 2.1
    mc_n = 1.3
    mc_cutoff = 2.5
    system.non_bonded_inter[1, 1].membrane_collision.set_params(
        a=mc_a, n=mc_n, cutoff=mc_cutoff)

    def create_first_group(self, system):
        """
        Prepares system to contain 4 particles. Of these, p1-p3
        form a triangle and p2,p3,p4 form a triangle sharing one edge.
        The OutDirection bond needs to be created on P4 with P1-P3 as partners.
        Returns a four instances of ParticleHandle.
        """
        base1 = np.random.random(3)
        side1_1 = np.random.random(3)
        side1_2 = np.random.random(3)
        normal1 = np.cross(side1_1, side1_2)
        normal1 /= np.linalg.norm(normal1)
        offset1 = .1 * normal1

        system.part.clear()
        triangle1 = system.part.add(pos=(
            base1,
            base1 + side1_1,
            base1 + side1_2), id=(1, 2, 3))
        extra1 = system.part.add(pos=base1 + side1_1 + side1_2 + offset1, id=4)
        return tuple(system.part[[1, 2, 3, 4]])

    def test_out_dir_handedness(self):
        """
        Checks that OutDirection is parallel to triangle normal
        and flips sign if the order of the particles forming the triangle
        is changed acyclically
        """

        system = self.system
        t1, t2, t3, extra = self.create_first_group(system)
        extra.bonds = (self.out_dir, t1.id, t2.id, t3.id)
        system.integrator.run(0)
        cross = np.cross(system.distance_vec(t1, t2),
                         system.distance_vec(t1, t3))
        normal = cross / np.linalg.norm(cross)

        np.testing.assert_allclose(
            extra.out_direction, normal)

        # Check that out_direction does not change for cyclic permutation
        extra.bonds = (self.out_dir, t2.id, t3.id, t1.id)
        system.integrator.run(0)
        np.testing.assert_allclose(
            extra.out_direction, normal)
        # Check that out_direction does not change for cyclic permutation
        extra.bonds = (self.out_dir, t3.id, t2.id, t1.id)
        system.integrator.run(0, recalc_forces=True)
        np.testing.assert_allclose(
            extra.out_direction, -normal)

    def test_membrane_collision(self):
        """
        Checks that membrane collision interaction force is calculated correctly
        """

        system = self.system
        t1_1, t1_2, t1_3, extra1 = self.create_first_group(system)
        extra1.bonds = (self.out_dir, t1_1.id, t1_2.id, t1_3.id)
        cross = np.cross(system.distance_vec(t1_1, t1_2),
                         system.distance_vec(t1_1, t1_3))
        normal1 = cross / np.linalg.norm(cross)

        # Construct shifted copy of first triangle + extra
        shift = np.array((0.1, 0.2, 0.3))
        t2_1, t2_2, t2_3, extra2 = system.part.add(
            id=(5, 6, 7, 8),
            pos=system.part[[t1_1.id, t1_2.id, t1_3.id, extra1.id]].pos + shift)
        t2_3.pos = t2_3.pos - 0.1 * shift
        extra2.bonds = (self.out_dir, t2_1.id, t2_3.id, t2_2.id)
        extra1.type = 1
        extra2.type = 1
        system.integrator.run(0, recalc_forces=True)
        force = np.copy(extra2.f)

        # Check that force does not change
        system.integrator.run(0, recalc_forces=True)
        np.testing.assert_allclose(np.copy(extra2.f), force)

        # Check that force is correct
        direction = extra1.out_direction - extra2.out_direction
        distance = np.linalg.norm(system.distance_vec(extra1, extra2))
        force_ref = direction * tests_common.membrane_collision_force(r=distance, a=self.mc_a, n=self.mc_n, dir=direction,
                                                              cutoff=self.mc_cutoff)
        np.testing.assert_allclose(force, -force_ref)

        # Check that force equals minus the counter-force
        np.testing.assert_allclose(force, -np.copy(extra1.f))


if __name__ == "__main__":
    ut.main()
