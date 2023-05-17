#
# Copyright (C) 2013-2023 The ESPResSo project
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
import numpy as np

import espressomd


class Test(ut.TestCase):
    N_PART = 200
    # Handle for espresso system
    system = espressomd.System(box_l=[10.0, 10.0, 10.0])
    partcls = system.part.add(
        id=np.arange(3, 3 + 2 * N_PART, 2),
        pos=np.random.random((N_PART, 3)) * system.box_l,
        v=np.random.random((N_PART, 3)) * 3.2 - 1,
        f=np.random.random((N_PART, 3)))

    if espressomd.has_features(["DIPOLES"]):
        partcls.dip = np.random.random((N_PART, 3)) - .3

    if espressomd.has_features(["ROTATION"]):
        direcs = np.random.random((N_PART, 3)) - 0.5
        direcs /= np.linalg.norm(direcs, axis=1)[:, None]
        partcls.director = direcs

    if espressomd.has_features("VIRTUAL_SITES"):
        p = system.part.by_id(partcls.id[8])
        p.virtual = True

    def generate_test_for_ml_observable(pprop_name):
        """Generates test cases for observables working on particle id lists.

        """

        def func(self):
            # This code is run at the execution of the generated function.
            # It will use the state of the variables in the outer function,
            # which was there, when the outer function was called
            # Randomly pick a subset of the particles
            id_list = sorted(
                np.random.choice(
                    self.system.part.all().id,
                    size=int(self.N_PART * .9),
                    replace=False))
            for id in id_list:
                self.assertTrue(self.system.part.exists(id))

            # Get data from particles
            if pprop_name == "f":
                for p_id in id_list:
                    if self.system.part.by_id(p_id).virtual:
                        id_list.remove(p_id)

            part_data = getattr(self.system.part.by_ids(id_list), pprop_name)

            # Data from observable
            ml_data = np.array(self.system.analysis.call_method(
                "particle_properties", pids=id_list,
                properties=[pprop_name])[pprop_name])

            # Check
            self.assertEqual(ml_data.shape, part_data.shape)
            np.testing.assert_allclose(ml_data, part_data, atol=1e-10)

        return func

    test_pos_folded = generate_test_for_ml_observable("pos_folded")
    test_pos_unfolded = generate_test_for_ml_observable("pos")
    test_force = generate_test_for_ml_observable("f")
    if espressomd.has_features(["ROTATION"]):
        test_director = generate_test_for_ml_observable("director")
    if espressomd.has_features(["DIPOLES"]):
        test_dipole = generate_test_for_ml_observable("dip")
    if espressomd.has_features(["ELECTROSTATICS"]):
        test_charge = generate_test_for_ml_observable("q")

    def test_write_forces(self):
        id_list = sorted(
            np.random.choice(
                self.system.part.all().id,
                size=int(self.N_PART * .9),
                replace=False))
        for id in id_list:
            self.assertTrue(self.system.part.exists(id))
        new_forces = np.random.random((len(id_list), 3))
        self.system.analysis.call_method(
            "write_ml_forces", pids=id_list, forces=new_forces)
        part_data = self.system.part.by_ids(id_list).ext_force
        np.testing.assert_allclose(part_data, new_forces, atol=1e-10)

    def test_exceptions(self):
        with self.assertRaisesRegex(RuntimeError, "Property 'unknown' is not implemented"):
            self.system.analysis.call_method(
                "particle_properties", pids=[], properties=["unknown"])


if __name__ == "__main__":
    ut.main()
